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


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
