#!/usr/bin/env python3
"""Tests for fov_annotations.py: the per-FOV Case/Block/Region reference and its cell-id join.

Run: python pipeline/python/tests/test_fov_annotations.py   (or under pytest)
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
import fov_annotations as fa  # noqa: E402

TWO_DONOR_SLIDE = "7495G37302G3"
EXPECTED_REGIONS = {"Tumor bulk", "Infiltrating edge", "Contralateral uninvolved"}


def _write_fixture(tmp: Path, bom: bool = False) -> tuple[Path, Path]:
    ann, xw = tmp / "ann.csv", tmp / "xw.csv"
    ann.write_text("Slide,Case,Block,Region,FOVs\n"
                   "Slide One,1,A,Tumor bulk,1\n"
                   "Slide One,2,B,Infiltrating edge,101\n"
                   "Ghost Slide,3,C,Tumor bulk,1\n",
                   encoding="utf-8-sig" if bom else "utf-8")
    xw.write_text("Canonical_slide_name,AtoMx_flatfile_folder\nSlide One,S1\n")
    return ann, xw


def test_units_are_folder_colon_fov():
    units = fa.units_from_cell_ids(pd.Index(["S1_F0007_C12", "S1_F101_C3"]))
    assert list(units) == ["S1:F7", "S1:F101"], units


def test_loader_reads_bom_and_plain_utf8():
    for bom in (False, True):
        with tempfile.TemporaryDirectory() as d:
            ann, xw = _write_fixture(Path(d), bom=bom)
            out = fa.load_fov_annotations(ann, xw)
        print(f"  bom={bom}: index={list(out.index)}")
        assert list(out.index) == ["S1:F1", "S1:F101"]
        assert list(out["Case"]) == [1, 2]


def test_slides_absent_from_the_crosswalk_are_dropped_not_kept():
    with tempfile.TemporaryDirectory() as d:
        ann, xw = _write_fixture(Path(d))
        out = fa.load_fov_annotations(ann, xw)
    assert "Ghost Slide" not in set(out["slide"]) and len(out) == 2


def test_annotate_cells_drops_unannotated_fovs():
    with tempfile.TemporaryDirectory() as d:
        ann, xw = _write_fixture(Path(d))
        out = fa.annotate_cells(pd.Index(["S1_F1_C1", "S1_F5_C1"]),
                                annotations=ann, crosswalk=xw)
    assert list(out.index) == ["S1_F1_C1"] and out.loc["S1_F1_C1", "Region"] == "Tumor bulk"


def test_committed_reference_resolves_a_real_two_donor_slide():
    """Integration check on the committed copies, using a slide that really carries two cases.

    7495G37302G3 is the slide that dominated the amplicon runs. The committed reference must
    split it into case 7495 (infiltrating edge) and case 7302 (tumour bulk) — if a future
    re-export flattens or renames it, this is what catches it.
    """
    ann = fa.load_fov_annotations()
    sub = ann[ann.index.str.startswith(f"{TWO_DONOR_SLIDE}:F")]
    by_case = sub.groupby("Case")["Region"].agg(["size", "first"])
    print(by_case)
    assert set(by_case.index) == {7495, 7302}, by_case
    assert by_case.loc[7495, "first"] == "Infiltrating edge", by_case
    assert by_case.loc[7302, "first"] == "Tumor bulk", by_case
    assert (by_case["size"] == 100).all(), by_case
    assert not ann.index.duplicated().any()


def test_committed_reference_shape_and_regions():
    ann = fa.load_fov_annotations()
    print(f"  rows={len(ann)} slides={ann['slide'].nunique()} cases={ann['Case'].nunique()}")
    assert len(ann) == 11425 and ann["slide"].nunique() == 57 and ann["Case"].nunique() == 38
    assert set(ann["Region"].unique()) == EXPECTED_REGIONS


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
