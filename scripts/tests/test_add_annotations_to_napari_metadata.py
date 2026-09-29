#!/usr/bin/env python3
"""Tests for add-annotations-to-napari-metadata.py (no S3 / network).

What matters: the FOV is read from the right token of `c_<slide#>_<FOV>_<cell>`, two slides
that reuse the same FOV numbers never bleed into each other (a slide is two DIFFERENT cases and
FOV numbering repeats across slides), row order survives, and every failure mode degrades to
blank-and-reported rather than to a wrong region. A wrong Region on a slide is worse than none.

Run:  uv run --with pandas --with duckdb python scripts/tests/test_add_annotations_to_napari_metadata.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "add-annotations-to-napari-metadata.py"


def _fixture(tmp: Path, meta_by_slide: dict[str, list[str]]) -> tuple[Path, Path, Path]:
    """Two annotated slides reusing FOV numbers 1-4, with a different layout on each."""
    rows = []
    for slide, cases in (("SLIDE A", (("111", "Tumor bulk"), ("222", "Infiltrating edge"))),
                         ("SLIDE B", (("333", "Contralateral uninvolved"), ("333", "Tumor bulk")))):
        for i, (case, region) in enumerate(cases):
            for fov in (1, 2) if i == 0 else (3, 4):
                rows.append({"Slide": slide, "Case": case, "Block": "G1",
                             "Region": region, "FOVs": fov})
    ann = tmp / "ann.csv"
    # The real reference is UTF-8 with a BOM, which corrupts the first header if not handled.
    pd.DataFrame(rows).to_csv(ann, index=False, encoding="utf-8-sig")
    xw = tmp / "xw.csv"
    pd.DataFrame([{"Canonical_slide_name": "SLIDE A", "AtoMx_flatfile_folder": "SLA"},
                  {"Canonical_slide_name": "SLIDE B", "AtoMx_flatfile_folder": "SLB"}]
                 ).to_csv(xw, index=False)
    md = tmp / "md"; md.mkdir()
    for slide, ids in meta_by_slide.items():
        pd.DataFrame({"cell_ID": ids, "cell_type": "Unassigned", "hex_color": "#861933"}
                     ).to_csv(md / f"{slide}_metadata.csv", index=False)
    return ann, xw, md


def _run(tmp, ann, xw, md, *extra):
    return subprocess.run([sys.executable, str(SCRIPT), "--metadata-dir", str(md),
                           "--out-dir", str(tmp / "out"), "--annotations", str(ann),
                           "--crosswalk", str(xw), *extra], capture_output=True, text=True)


def _out(tmp, slide):
    return pd.read_csv(tmp / "out" / f"{slide}_metadata.csv", dtype=str, keep_default_na=False)


def test_region_comes_from_the_fov_token_and_slides_do_not_bleed():
    """SLA and SLB both have FOVs 1-4 but different regions; each must get its own."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ids = ["c_1_1_1", "c_1_2_5", "c_1_3_9", "c_1_4_2"]
        ann, xw, md = _fixture(tmp, {"SLA": ids, "SLB": ["c_2_1_1", "c_2_3_1"]})
        assert _run(tmp, ann, xw, md, "--columns", "Region,Case").returncode == 0
        a, b = _out(tmp, "SLA"), _out(tmp, "SLB")
        assert list(a.Region) == ["Tumor bulk", "Tumor bulk", "Infiltrating edge", "Infiltrating edge"]
        assert list(a.Case) == ["111", "111", "222", "222"]
        assert list(b.Region) == ["Contralateral uninvolved", "Tumor bulk"]   # NOT slide A's
        assert list(b.Case) == ["333", "333"]


def test_each_column_gets_a_paired_colour_and_categories_share_one():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1", "c_1_2_1", "c_1_3_1"]})
        assert _run(tmp, ann, xw, md, "--columns", "Region,Case").returncode == 0
        a = _out(tmp, "SLA")
        assert {"Region", "Region_color", "Case", "Case_color"} <= set(a.columns)
        assert a.groupby("Region").Region_color.nunique().max() == 1
        assert a.Region_color.str.fullmatch(r"#[0-9A-F]{6}").all()


def test_existing_columns_and_row_order_survive():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ids = ["c_1_4_7", "c_1_1_3", "c_1_2_1"]                    # deliberately not sorted
        ann, xw, md = _fixture(tmp, {"SLA": ids})
        assert _run(tmp, ann, xw, md).returncode == 0
        a = _out(tmp, "SLA")
        assert list(a.cell_ID) == ids
        assert list(a.cell_type) == ["Unassigned"] * 3 and list(a.hex_color) == ["#861933"] * 3


def test_a_fov_with_no_annotation_is_blank_and_reported_not_guessed():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1", "c_1_99_1"]})   # FOV 99 is not annotated
        out = _run(tmp, ann, xw, md)
        assert out.returncode == 0, out.stderr
        a = _out(tmp, "SLA")
        assert list(a.Region) == ["Tumor bulk", ""] and list(a.Region_color)[1] == ""
        assert "1" in out.stdout.split("SLA")[1].splitlines()[0]     # unannotated count shown


def test_a_slide_missing_from_the_reference_is_flagged_not_crashed():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1"], "NOSUCH": ["c_9_1_1"]})
        out = _run(tmp, ann, xw, md)
        assert out.returncode == 0, out.stderr
        assert list(_out(tmp, "NOSUCH").Region) == [""]
        assert "CHECK: none matched" in out.stdout


def test_a_malformed_cell_id_is_left_blank():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1", "not_an_id"]})
        out = _run(tmp, ann, xw, md)
        assert out.returncode == 0, out.stderr
        assert list(_out(tmp, "SLA").Region) == ["Tumor bulk", ""]
        assert "did not look like" in out.stderr


def test_it_refuses_to_overwrite_an_existing_region_column():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1"]})
        f = md / "SLA_metadata.csv"
        pd.read_csv(f).assign(Region="already here").to_csv(f, index=False)
        out = _run(tmp, ann, xw, md)
        assert out.returncode != 0 and "Refusing to overwrite" in (out.stdout + out.stderr)


def test_unknown_column_is_rejected():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw, md = _fixture(tmp, {"SLA": ["c_1_1_1"]})
        out = _run(tmp, ann, xw, md, "--columns", "Diagnosis")
        assert out.returncode != 0 and "subset of" in (out.stdout + out.stderr)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
