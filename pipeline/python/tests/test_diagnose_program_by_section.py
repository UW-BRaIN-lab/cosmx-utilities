#!/usr/bin/env python3
"""Tests for section derivation in diagnose_program_by_section.py (no S3 / network).

The one that matters is test_section_grouping_is_not_fooled_by_a_hot_section: two tissue
sections are mounted per CosMx slide, so grouping by slide pools two pieces. If one piece is
hot (ischaemia / slow fixation) and the letter's cells are concentrated in it, slide-level
grouping reports a confident CELL STATE that is not there.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy --with scipy \\
        python pipeline/python/tests/test_diagnose_program_by_section.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location(
    "diagnose_program_by_section", _PY_DIR / "diagnose_program_by_section.py")
dp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dp)


def test_sections_split_a_slide_at_the_fov_gap():
    """Two pieces per slide show up as two contiguous FOV runs (e.g. 1-135 and 181-205)."""
    ids = pd.Index([f"SL01A_F{f}_C{i}" for i, f in enumerate([1, 2, 3, 181, 182, 183])])
    got = dp.sections_from_cell_ids(ids, fov_gap=10)
    assert got.nunique() == 2
    assert got.iloc[0] == "SL01A:F1-3"
    assert got.iloc[-1] == "SL01A:F181-183"


def test_sections_keep_a_contiguous_run_together():
    """A slide with one continuous piece must not be chopped up by small FOV steps."""
    ids = pd.Index([f"SL02B_F{f}_C{i}" for i, f in enumerate([1, 2, 3, 4, 5])])
    assert dp.sections_from_cell_ids(ids, fov_gap=10).nunique() == 1


def test_sections_never_merge_across_slides():
    """Two slides sharing FOV numbers are still two different pieces of tissue."""
    ids = pd.Index(["SL01A_F1_C1", "SL02B_F1_C2"])
    got = dp.sections_from_cell_ids(ids, fov_gap=10)
    assert got.nunique() == 2


def _fixture_hot_section(rng):
    """One section per slide is hot; the letter is 80% concentrated in that hot section."""
    rows = []
    for s in range(10):
        slide = f"SL{s:02d}A{s:02d}B"
        for fov_lo, hot in ((1, True), (181, False)):
            n_letter = 320 if hot else 80
            for is_letter, n in ((True, n_letter), (False, 600)):
                score = rng.normal(8.0 if hot else 1.0, 0.5, size=n)
                rows.append(pd.DataFrame({
                    "score": score,
                    "t": is_letter,
                    "cell_id": [f"{slide}_F{fov_lo + (i % 100)}_C{rng.integers(1e9)}"
                                for i in range(n)],
                }))
    return pd.concat(rows, ignore_index=True)


def _verdict(df, unit):
    g = df.groupby([unit, "t"])["score"].mean().unstack()
    gap = (g[True] - g[False]).dropna()
    return float((gap > 0).mean()), float(gap.median())


def test_section_grouping_is_not_fooled_by_a_hot_section():
    """The regression this file exists for. Ground truth is HANDLING, not a cell state."""
    df = _fixture_hot_section(np.random.default_rng(3))
    df["section"] = dp.sections_from_cell_ids(pd.Index(df["cell_id"]), 10).to_numpy()
    df["slide"] = pd.Series(df["section"]).str.split(":").str[0]

    slide_pos, slide_gap = _verdict(df, "slide")
    sec_pos, sec_gap = _verdict(df, "section")

    # Grouping by slide pools the hot and cold pieces and manufactures a cell-state signal.
    assert slide_pos == 1.0, "slide-level should (wrongly) look unanimous"
    assert slide_gap > 1.0, "slide-level should (wrongly) show a large gap"
    # Grouping by section compares like with like, and the signal disappears.
    assert abs(sec_gap) < 0.2, f"section-level gap should collapse, got {sec_gap:+.3f}"
    assert 0.2 < sec_pos < 0.8, f"section-level should be a coin toss, got {sec_pos:.2f}"


def _fixture_diluted_pool(rng, n_units=20):
    """A malignant-vs-malignant deficit sitting under a large non-malignant majority.

    Ground truth: the letter carries the programme at 1.0 while the malignant comparator carries
    it at 3.0 — a real, large, cell-level deficit. But the units are 77% non-malignant cells
    scoring 0.6, so an all-other-cells pool averages to roughly the letter's own level and the
    deficit vanishes. This is the shape of the amplicon question: the reference stratifies hard
    on the programme, so the comparison group has to be chosen, not defaulted.
    """
    rows = []
    for u in range(n_units):
        for cell_type, n, mu in (("t", 60, 1.0), ("OPC-like", 60, 3.0),
                                 ("MES-like_hypoxia_MHC", 60, 1.1), ("Oligodendrocyte", 400, 0.6)):
            rows.append(pd.DataFrame({
                "score": rng.normal(mu, 0.3, size=n),
                "cell_type": cell_type,
                "slide": f"SL{u:02d}:F{u}",
            }))
    df = pd.concat(rows, ignore_index=True)
    df["is_letter"] = df["cell_type"] == "t"
    return df


def test_pooled_comparator_hides_what_compare_to_reveals():
    """The regression --compare-to exists for: dilution can erase a real deficit, even flip it."""
    df = _fixture_diluted_pool(np.random.default_rng(11))

    _, pooled = dp.section_table(df, min_cells=25)
    pooled_gap = dp.summarize(pooled)["median_gap"]

    matched = df[df["is_letter"] | (df["cell_type"] == "OPC-like")]
    _, matched_usable = dp.section_table(matched, min_cells=25)
    matched_gap = dp.summarize(matched_usable)["median_gap"]

    assert pooled_gap > 0, (
        f"the all-other-cells pool should (wrongly) put t at or above its neighbours, "
        f"got {pooled_gap:+.3f}")
    assert matched_gap < -1.5, (
        f"against OPC-like specifically the deficit must show, got {matched_gap:+.3f}")


def test_per_type_breakdown_separates_the_comparators():
    """Each comparator is its own matched comparison, so a per-type gap must track that type."""
    df = _fixture_diluted_pool(np.random.default_rng(12))
    gaps = {}
    for cell_type in ("OPC-like", "MES-like_hypoxia_MHC"):
        pair = df[df["is_letter"] | (df["cell_type"] == cell_type)]
        _, usable = dp.section_table(pair, min_cells=25)
        gaps[cell_type] = dp.summarize(usable)["median_gap"]
    # t is 1.0; OPC-like is 3.0 and MES is 1.1, so only the former is a real deficit.
    assert gaps["OPC-like"] < -1.5, gaps
    assert abs(gaps["MES-like_hypoxia_MHC"]) < 0.3, gaps


def test_letter_pct_is_a_share_of_the_pair_under_compare_to():
    """With a restricted comparator, letter_pct must not be diluted by the cells left out."""
    df = _fixture_diluted_pool(np.random.default_rng(13))
    pair = df[df["is_letter"] | (df["cell_type"] == "OPC-like")]
    _, usable = dp.section_table(pair, min_cells=25)
    assert np.allclose(usable["letter_pct"], 50.0), usable["letter_pct"].unique()


def test_resolve_comparators_expands_the_malignant_shorthand():
    present = set(dp.MALIGNANT_CORE_L4) | {"t", "Oligodendrocyte"}
    got = dp.resolve_comparators("malignant", present, "t")
    assert got == list(dp.MALIGNANT_CORE_L4), got


def test_resolve_comparators_skips_absent_types_and_the_letter_itself():
    got = dp.resolve_comparators("OPC-like,NotAType,t", {"OPC-like", "t"}, "t")
    assert got == ["OPC-like"], got


def test_resolve_comparators_defaults_to_every_other_cell():
    assert dp.resolve_comparators("", {"t", "OPC-like"}, "t") == [dp.COMPARE_ALL]


# --------------------------------------------------------------------------------------------
# Matched-null control. These two fixtures are the whole reason --null-sets exists: a depth
# difference and a real programme deficit produce the SAME negative gap, and only the null
# separates them.
# --------------------------------------------------------------------------------------------
import scipy.sparse as sp  # noqa: E402  (kept beside the tests that need it)

_PROGRAM_SIZE = 5
_SCALE = 1e4


def _simulate(rng, depth_ratio=1.0, program_knockdown=1.0,
              n_genes=400, n_cells=4000, n_units=20, letter_frac=0.35):
    """A letter and one comparator across units, with depth and programme effects separable.

    Column sums are the per-cell depth by construction, so `depth_ratio` is exactly the library
    size the letter's cells get relative to the comparator's.
    """
    gene_rate = rng.lognormal(mean=-1.0, sigma=1.5, size=n_genes)
    gene_rate /= gene_rate.sum()
    is_letter = rng.random(n_cells) < letter_frac
    depth = np.where(is_letter, depth_ratio, 1.0) * 300.0
    rate = np.outer(gene_rate, depth)

    # The programme is mid-abundance, so a matched pool actually has candidates on both sides.
    order = np.argsort(gene_rate)
    program_idx = order[int(.4 * n_genes):int(.4 * n_genes) + _PROGRAM_SIZE]
    if program_knockdown != 1.0:
        rate[np.ix_(program_idx, np.where(is_letter)[0])] *= program_knockdown

    mat = sp.csc_matrix(rng.poisson(rate).astype(np.float64))
    totals = np.asarray(mat.sum(axis=0)).ravel()
    frame = pd.DataFrame({
        "cell_type": np.where(is_letter, "t", "OPC-like"),
        "slide": [f"U{i % n_units}" for i in range(n_cells)],
        "is_letter": is_letter,
    })
    return mat, totals, program_idx, frame


def _run_null(mat, totals, program_idx, frame, n_sets=100, seed=7, min_cells=25):
    """Observed gap, the null gaps, and the empirical p — the same path main() takes."""
    comparators, pooled_only = ["OPC-like"], False
    frame = frame.copy()
    frame["score"] = dp.score_gene_set(mat, totals, np.asarray(program_idx), _SCALE)
    observed = dp.median_gaps(frame, comparators, pooled_only, min_cells)[dp.POOLED_KEY]

    mean_counts, detection = dp.gene_statistics(mat)
    pools = dp.matched_pools(mean_counts, detection, np.asarray(program_idx), 40)
    draws = dp.draw_null_sets(pools, n_sets, np.random.default_rng(seed))

    null_frame = frame.copy()
    gaps = []
    for gene_idx in draws:
        null_frame["score"] = dp.score_gene_set(mat, totals, np.asarray(gene_idx), _SCALE)
        gaps.append(dp.median_gaps(null_frame, comparators, pooled_only, min_cells)[dp.POOLED_KEY])
    gaps = np.asarray(gaps, dtype=float)
    return observed, gaps, dp.empirical_p(observed, gaps)


def test_gene_statistics_match_a_dense_computation():
    rng = np.random.default_rng(1)
    dense = rng.poisson(0.4, size=(30, 50)).astype(np.float64)
    mean_counts, detection = dp.gene_statistics(sp.csc_matrix(dense))
    assert np.allclose(mean_counts, dense.mean(axis=1))
    assert np.allclose(detection, (dense > 0).mean(axis=1))


def test_matched_pools_never_offer_a_programme_gene_as_its_own_control():
    rng = np.random.default_rng(2)
    mat, _, program_idx, _ = _simulate(rng)
    mean_counts, detection = dp.gene_statistics(mat)
    for pool in dp.matched_pools(mean_counts, detection, np.asarray(program_idx), 40):
        assert not set(pool) & set(program_idx.tolist())


def test_matched_pools_track_abundance_and_detection():
    """A control gene of the wrong abundance would not carry the dropout being controlled for."""
    rng = np.random.default_rng(3)
    mat, _, program_idx, _ = _simulate(rng)
    mean_counts, detection = dp.gene_statistics(mat)
    pools = dp.matched_pools(mean_counts, detection, np.asarray(program_idx), 40)
    for pos, pool in zip(program_idx, pools):
        # Within a factor of a few on abundance, and close on the dropout-relevant axis.
        ratio = mean_counts[pool] / max(mean_counts[pos], 1e-12)
        assert np.median(ratio) == pytest_approx_within(0.25, 4.0), np.median(ratio)
        assert abs(np.median(detection[pool]) - detection[pos]) < 0.15


def pytest_approx_within(lo, hi):
    """Tiny helper so the bound reads as a range in the assertion above."""
    class _Within:
        def __eq__(self, other):
            return lo <= other <= hi
        def __repr__(self):
            return f"within [{lo}, {hi}]"
    return _Within()


def test_draw_null_sets_are_unique_within_a_set_and_reproducible():
    rng = np.random.default_rng(4)
    mat, _, program_idx, _ = _simulate(rng)
    mean_counts, detection = dp.gene_statistics(mat)
    pools = dp.matched_pools(mean_counts, detection, np.asarray(program_idx), 40)
    first = dp.draw_null_sets(pools, 25, np.random.default_rng(11))
    again = dp.draw_null_sets(pools, 25, np.random.default_rng(11))
    assert first == again, "a reported p must be reproducible from the seed"
    for one in first:
        assert len(set(one)) == len(one)
        assert len(one) == _PROGRAM_SIZE


def test_empirical_p_is_never_zero():
    """The observed set is itself a draw, so p=0 is not an available answer."""
    assert dp.empirical_p(-9.9, np.linspace(0, 1, 100)) == 1 / 101
    assert np.isnan(dp.empirical_p(float("nan"), np.linspace(0, 1, 10)))


def test_depth_driven_deficit_is_NOT_called_a_programme_deficit():
    """The regression this control exists for.

    Ground truth: the letter has identical composition and 0.4x the depth. Every gene drops out
    more often in it, so the programme gap is negative with no programme biology at all. The
    matched null must reproduce that gap and refuse the finding.
    """
    mat, totals, program_idx, frame = _simulate(np.random.default_rng(5), depth_ratio=0.4)
    observed, gaps, p = _run_null(mat, totals, program_idx, frame)
    print(f"      observed {observed:+.3f}   null median {np.median(gaps):+.3f}   "
          f"null p05 {np.percentile(gaps, 5):+.3f}   p = {p:.3f}")

    assert observed < -0.05, f"the depth artefact should look like a deficit, got {observed:+.3f}"
    assert np.median(gaps) < -0.05, (
        f"matched genes must reproduce the artefact, null median {np.median(gaps):+.3f}")
    assert p > 0.05, f"a pure depth artefact must NOT reach significance, got p={p:.3f}"


def test_a_real_programme_deficit_survives_the_matched_null():
    """The converse: same depth, programme genes specifically knocked down. p must be small."""
    mat, totals, program_idx, frame = _simulate(
        np.random.default_rng(6), depth_ratio=1.0, program_knockdown=0.2)
    observed, gaps, p = _run_null(mat, totals, program_idx, frame)
    print(f"      observed {observed:+.3f}   null median {np.median(gaps):+.3f}   "
          f"null p05 {np.percentile(gaps, 5):+.3f}   p = {p:.3f}")

    assert observed < -0.05, f"the knockdown should show as a deficit, got {observed:+.3f}"
    assert observed < np.percentile(gaps, 5), (
        f"a real deficit must sit below the null's 5th pct "
        f"({observed:+.3f} vs {np.percentile(gaps, 5):+.3f})")
    assert p <= 0.05, f"a real programme deficit must reach significance, got p={p:.3f}"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
