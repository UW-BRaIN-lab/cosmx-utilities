#!/usr/bin/env python3
"""Tests for fov_annotations.py: the per-FOV Case/Block/Region reference and its cell-id join.

Run: python pipeline/python/tests/test_fov_annotations.py   (or under pytest)
"""

import contextlib
import io
import sys
import tempfile
from pathlib import Path

import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
import fov_annotations as fa  # noqa: E402

TWO_DONOR_SLIDE = "7495G37302G3"
EXPECTED_REGIONS = {"Tumor bulk", "Infiltrating edge", "Contralateral uninvolved"}
# Added 2026-11-23 (flat-file export part 4): six surgical-specimen slides plus 7706-A7, whose
# FOVs 55-120 are nerve (not GBM study tissue; Case/Block/Region all "Nerve").
PART4_FOLDERS = {
    "NP17-01317 A2 NP19-00472 A3": "NP1701317A2NP1900472A3",
    "NP19-00737 A2 NP20-00960 B1": "NP1900737A2NP2000960B1",
    "NP19-00306 B1 NP19-00537 C1": "NP1900306B1NP1900537C1",
    "7706-A7 plus nerves": "7706A7plusnerves",
    "NP19.02212.B1 NP23.00712.D8": "NP1902212B1NP2300712D8",
    "NP19-00413 A2 NP20-01336 A2": "NP1900413A2NP2001336A2",
    "NP21.11263.E3 NP21.11311.A3": "NP2111263E3NP2111311A3",
}
NERVE_SLIDE, NERVE_FOVS = "7706-A7 plus nerves", set(range(55, 121))
ALL_REGIONS = EXPECTED_REGIONS | {"Surgical", "Nerve"}


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
        assert list(out["Case"]) == ["1", "2"]


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
    assert set(by_case.index) == {"7495", "7302"}, by_case
    assert by_case.loc["7495", "first"] == "Infiltrating edge", by_case
    assert by_case.loc["7302", "first"] == "Tumor bulk", by_case
    assert (by_case["size"] == 100).all(), by_case
    assert not ann.index.duplicated().any()


def test_committed_reference_shape_and_regions():
    ann = fa.load_fov_annotations()
    print(f"  rows={len(ann)} slides={ann['slide'].nunique()} cases={ann['Case'].nunique()}")
    # 38 donors plus the "Nerve" placeholder Case on 7706-A7 FOVs 55-120
    assert len(ann) == 12746 and ann["slide"].nunique() == 64 and ann["Case"].nunique() == 39
    assert set(ann["Region"].unique()) == ALL_REGIONS

def test_committed_reference_carries_the_surgical_and_nerve_rows():
    raw = pd.read_csv(fa.DEFAULT_ANNOTATIONS, dtype=str)
    print(raw.groupby(["Region"]).size())
    assert len(raw) == 12746 and raw["Slide"].nunique() == 64
    assert set(raw["Region"]) == ALL_REGIONS
    assert set(PART4_FOLDERS) <= set(raw["Slide"])
    nerve = raw[raw["Region"] == "Nerve"]
    assert set(nerve["Slide"]) == {NERVE_SLIDE}
    assert set(nerve["FOVs"].astype(int)) == NERVE_FOVS
    assert (nerve["Case"] == "Nerve").all() and (nerve["Block"] == "Nerve").all()


def test_every_surgical_case_is_also_a_postmortem_case():
    """Why stage 3a splits the PCA batch: a plain Case batch would merge these with postmortem."""
    raw = pd.read_csv(fa.DEFAULT_ANNOTATIONS, dtype=str)
    surgical = set(raw.loc[raw["Region"] == "Surgical", "Case"])
    postmortem = set(raw.loc[raw["Region"].isin(EXPECTED_REGIONS), "Case"])
    print(f"surgical cases: {sorted(surgical)}")
    assert len(surgical) == 12 and surgical <= postmortem


def test_part4_slides_resolve_to_their_flat_file_folders():
    """The surgical slides and 7706-A7 join to the folders AtoMx wrote on 2026-11-23, with no
    slide dropped for a missing folder, and the nerve FOVs stay findable for EXCLUDE_REGIONS."""
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        ann = fa.load_fov_annotations()
    out = captured.getvalue()
    print(out)
    assert "no AtoMx_flatfile_folder" not in out
    raw = pd.read_csv(fa.DEFAULT_ANNOTATIONS, dtype=str)
    for slide, folder in PART4_FOLDERS.items():
        expected = (raw["Slide"] == slide).sum()
        found = (ann["slide"] == folder).sum()
        print(f"  {slide} -> {folder}: {found}/{expected} FOVs")
        assert expected and found == expected
    nerve = ann[ann["Region"] == "Nerve"]
    assert set(nerve.index) == {f"7706A7plusnerves:F{f}" for f in NERVE_FOVS}
    assert (nerve["Case"] == "Nerve").all()


def test_text_case_reads_as_text_throughout():
    with tempfile.TemporaryDirectory() as tmp:
        ann, xw = Path(tmp) / "ann.csv", Path(tmp) / "xw.csv"
        ann.write_text("Slide,Case,Block,Region,FOVs\n"
                       "Slide One,7706,A7,Contralateral uninvolved,1\n"
                       "Slide One,Nerve,Nerve,Nerve,2\n")
        xw.write_text("Canonical_slide_name,AtoMx_flatfile_folder\nSlide One,S1\n")
        out = fa.load_fov_annotations(ann, xw)
    print(out)
    assert out["Case"].tolist() == ["7706", "Nerve"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
