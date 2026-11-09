#!/usr/bin/env python3
"""Tests for qc_by_region.py: is ONE region of ONE donor usable?

The fixture is built so that every naive comparison gives the WRONG answer and only the
per-region one is right. The focus donor has a degraded brain (60 counts) and a normal
optic nerve (300), while optic nerve is genuinely sparser than retina cohort-wide
(300 vs 900). So:
  - against the cohort median across all regions, its good optic nerve looks degraded;
  - against its own slide median, its good optic nerve looks degraded too;
  - against OTHER DONORS' OPTIC NERVE, it is exactly normal, which is the truth.

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_qc_by_region.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location("qc_by_region", _PY_DIR / "qc_by_region.py")
qr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qr)

FOCUS = "uwa0001eyes"
PEERS = ("uwa0002eyes", "uwa0003eyes", "uwa0004eyes")
# Enough that the smallest donor-region (optic nerve, 2 FOVs) clears
# qr.MIN_CELLS_TO_RANK, as every region of every real donor does.
N_PER_FOV = 150

# region -> (fovs, cohort depth). The focus donor's brain is degraded; its eye is not.
# Retina gets twice the FOVs so that a slide-wide median sits ABOVE optic nerve even on
# the focus slide -- otherwise the degraded brain drags the slide median down to meet the
# optic nerve and the within-slide comparison accidentally gives the right answer.
REGIONS = {
    "Retina": ((1, 2, 7, 8), 900),
    "Optic nerve": ((3, 4), 300),
    "Gray matter": ((5, 6), 800),
}
FOCUS_DEPTH = {"Retina": 900, "Optic nerve": 300, "Gray matter": 60}
EXCLUDED = {(FOCUS, 5), (FOCUS, 6)}  # the degraded brain, as the SME ruled


def build_cells() -> pd.DataFrame:
    rows = []
    for slide in (FOCUS, *PEERS):
        for region, (fovs, depth) in REGIONS.items():
            counts = FOCUS_DEPTH[region] if slide == FOCUS else depth
            for fov in fovs:
                for cell in range(N_PER_FOV):
                    rows.append({"slide_id": slide, "fov": fov,
                                 "qc_gene_counts": counts,
                                 "qc_genes_detected": counts // 3,
                                 "qc_negprobe_prop": 0.01,
                                 "cell": f"{slide}_F{fov}_C{cell}"})
    return pd.DataFrame(rows).set_index("cell")


def write_annotations(path: Path) -> Path:
    rows = [{"slide_id": slide, "fov": fov, "region": region, "mixed_adjacent": 0,
             "exclude": int((slide, fov) in EXCLUDED)}
            for slide in (FOCUS, *PEERS)
            for region, (fovs, _) in REGIONS.items() for fov in fovs]
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def annotated(tmp_path: Path) -> pd.DataFrame:
    return qr.attach_regions(build_cells(), qr.load_regions(write_annotations(tmp_path / "f.csv")))


def test_region_comes_from_the_sheet_not_obs(tmp_path: Path) -> None:
    cells = annotated(tmp_path)
    assert len(cells) == len(build_cells())
    assert set(cells["region"]) == set(REGIONS)
    # Only the focus donor's brain is excluded, and it keeps its true anatomy.
    excluded = cells.loc[cells["excluded"]]
    assert set(excluded["slide_id"]) == {FOCUS}
    assert set(excluded["region"]) == {"Gray matter"}
    print("  region and exclusion come from the sheet, keyed on (slide, FOV)")


def test_a_join_that_matches_nothing_is_fatal(tmp_path: Path) -> None:
    regions = qr.load_regions(write_annotations(tmp_path / "f.csv"))
    regions["slide_id"] = "a_different_cohort"
    try:
        qr.attach_regions(build_cells(), regions)
    except SystemExit as exc:
        assert "no cell matched" in str(exc)
    else:
        raise AssertionError("a slide-id mismatch must exit, not produce an empty table")
    print("  a sheet that matches no cell exits rather than tabulating nothing")


def test_the_cohort_yardstick_excludes_excluded_cells(tmp_path: Path) -> None:
    cohort = qr.by_region(annotated(tmp_path), (50,)).set_index("region")
    # Gray matter cohort-wide must be the peers' 800, undragged by the focus donor's 60.
    assert cohort.loc["Gray matter", "counts_q50"] == 800
    assert cohort.loc["Gray matter", "n_donors"] == len(PEERS)
    print("  the yardstick is the data we intend to keep, not everything measured")


def test_a_good_region_is_not_condemned_by_its_slide(tmp_path: Path) -> None:
    """The load-bearing case: the focus donor's optic nerve is FINE and must read as fine."""
    cells = annotated(tmp_path)
    table = qr.focus(cells, FOCUS, (50, 100)).set_index("region")

    # Every naive comparison would call this region degraded.
    all_region_median = float(np.median(cells["qc_gene_counts"]))
    own_slide_median = float(np.median(
        cells.loc[cells["slide_id"] == FOCUS, "qc_gene_counts"]))
    assert 300 < all_region_median, "fixture: optic nerve must look thin cohort-wide"
    assert 300 < own_slide_median, "fixture: optic nerve must look thin within the slide"

    optic = table.loc["Optic nerve"]
    assert optic["ratio_to_others"] == 1.0, optic["ratio_to_others"]
    assert optic["other_donors"] == len(PEERS)
    assert optic["pct_under_100"] == 0.0
    print(f"  optic nerve: {optic['counts_median']:.0f} vs others' "
          f"{optic['other_counts_median']:.0f} -- normal, though it is below both the "
          f"cohort ({all_region_median:.0f}) and its own slide ({own_slide_median:.0f})")


def test_a_bad_region_still_reads_as_bad(tmp_path: Path) -> None:
    table = qr.focus(annotated(tmp_path), FOCUS, (50, 100)).set_index("region")
    gray = table.loc["Gray matter"]
    assert bool(gray["excluded"]) is True
    assert gray["ratio_to_others"] < 0.1, gray["ratio_to_others"]
    assert gray["pct_under_100"] == 100.0
    # The excluded region is the control: it must be visibly worse than the kept one.
    assert gray["ratio_to_others"] < table.loc["Optic nerve", "ratio_to_others"] / 5
    print(f"  gray matter: ratio {gray['ratio_to_others']:.2f} -- the exclusion is "
          f"separable from the regions kept on the same slide")


def test_rank_is_over_donor_medians(tmp_path: Path) -> None:
    table = qr.focus(annotated(tmp_path), FOCUS, (50,)).set_index("region")
    optic = table.loc["Optic nerve"]
    assert optic["rank_of"] == len(PEERS) + 1
    # Tied with every peer at 300, so nobody is strictly above it: rank 1 of 4.
    assert optic["rank"] == 1
    assert table.loc["Gray matter", "rank"] == len(PEERS) + 1
    print("  rank is over per-donor medians, so one donor's tail cannot set the bar")


def test_an_unknown_focus_slide_is_fatal(tmp_path: Path) -> None:
    try:
        qr.focus(annotated(tmp_path), "not_a_slide", (50,))
    except SystemExit as exc:
        assert "not in this cohort" in str(exc)
    else:
        raise AssertionError("a typo'd --focus-slide must exit, not return an empty table")
    print("  a typo'd --focus-slide exits and lists the slides there are")


def test_end_to_end_writes_the_three_tables(tmp_path: Path) -> None:
    cells = build_cells()
    n_genes = 4
    var = pd.DataFrame({"probe_type": ["gene"] * n_genes},
                       index=[f"GENE{g}" for g in range(n_genes)])
    adata = ad.AnnData(X=sp.csr_matrix(np.ones((len(cells), n_genes), dtype=np.float32)),
                       obs=cells, var=var)
    h5ad = tmp_path / "combined_qc.h5ad"
    adata.write_h5ad(h5ad)

    out = tmp_path / "out"
    argv = sys.argv  # restored, so running under pytest cannot disturb a later test
    try:
        sys.argv = ["qc_by_region.py", "--combined-h5ad", str(h5ad),
                    "--fov-annotations", str(write_annotations(tmp_path / "f2.csv")),
                    "--focus-slide", FOCUS, "--output-dir", str(out)]
        qr.main()
    finally:
        sys.argv = argv

    per_slide = pd.read_csv(out / "qc_by_region.csv")
    assert len(per_slide) == len(REGIONS) * (len(PEERS) + 1)
    assert pd.read_csv(out / "qc_region_summary.csv").shape[0] == len(REGIONS)
    assert (out / f"qc_focus_{FOCUS}.csv").exists()
    print(f"  wrote {len(per_slide)} (slide, region) rows and the focus table")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("test_qc_by_region")
        test_region_comes_from_the_sheet_not_obs(tmp)
        test_a_join_that_matches_nothing_is_fatal(tmp)
        test_the_cohort_yardstick_excludes_excluded_cells(tmp)
        test_a_good_region_is_not_condemned_by_its_slide(tmp)
        test_a_bad_region_still_reads_as_bad(tmp)
        test_rank_is_over_donor_medians(tmp)
        test_an_unknown_focus_slide_is_fatal(tmp)
        test_end_to_end_writes_the_three_tables(tmp)
    print("All tests passed.")


if __name__ == "__main__":
    main()
