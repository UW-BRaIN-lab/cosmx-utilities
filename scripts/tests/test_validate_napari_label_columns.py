#!/usr/bin/env python3
"""Tests for validate-napari-label-columns.py (no S3 / network).

The important one builds a finished file with the REAL writers (labels-csv-to-typing-h5.py then
celltypes-to-napari-metadata.py --merge-into --drop-columns) and requires the validator to PASS it: a validator that
has never seen a correct file proves nothing. The rest plant one fault each (wrong label, surviving placeholder, wrong
colour, one label with two colours) and require a failure, including the "all cells mismatch" hint.

Run:  uv run --with pandas --with numpy --with h5py --with duckdb python scripts/tests/test_validate_napari_label_columns.py
"""
import gzip
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("val", SCRIPTS / "validate-napari-label-columns.py")
val = importlib.util.module_from_spec(spec); spec.loader.exec_module(val)
SLIDE = "SLX"
IDS = ["c_3_1_1", "c_3_1_2", "c_3_2_1"]


def _build_real(tmp: Path):
    """Stitched original + labels -> merged, through the real scripts. Cell SLX_F1_C2 has no label."""
    (tmp / "flat").mkdir(); (tmp / "stitched").mkdir(); (tmp / "labels").mkdir()
    with gzip.open(tmp / "flat" / f"{SLIDE}_metadata_file.csv.gz", "wt") as f:
        f.write("fov,cell_ID,cell_id\n1,1,c_3_1_1\n1,2,c_3_1_2\n2,1,c_3_2_1\n")
    pd.DataFrame({"cell_ID": IDS, "cell_type": "Unassigned", "hex_color": "#861933"}
                 ).to_csv(tmp / "stitched" / f"{SLIDE}_metadata.csv", index=False)
    pd.DataFrame({"cell_id": [f"{SLIDE}_F1_C1", f"{SLIDE}_F2_C1"], "label": ["UWA_741", "UWA_6803"]}
                 ).to_csv(tmp / "labels" / "donor.csv", index=False)
    subprocess.run([sys.executable, str(SCRIPTS / "labels-csv-to-typing-h5.py"), "--csv", str(tmp / "labels" / "donor.csv"),
                    "--label-column", "label", "--output", str(tmp / "donor.h5")], check=True, capture_output=True)
    subprocess.run([sys.executable, str(SCRIPTS / "celltypes-to-napari-metadata.py"), "--typing-h5", str(tmp / "donor.h5"),
                    "--flatfiles", str(tmp / "flat"), "--column", "Donor", "--merge-into", str(tmp / "stitched"),
                    "--drop-columns", "cell_type,hex_color", "--out-dir", str(tmp / "merged")], check=True, capture_output=True)


def _frames(tmp: Path):
    read = lambda p: pd.read_csv(p, dtype=str, keep_default_na=False)
    return (read(tmp / "merged" / f"{SLIDE}_metadata.csv"), read(tmp / "stitched" / f"{SLIDE}_metadata.csv"),
            val.load_sources(tmp / "labels", ["Donor=donor.csv"]))


def test_a_file_built_by_the_real_writers_passes():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); _build_real(tmp)
        merged, original, sources = _frames(tmp)
        print(merged)
        assert val.check_slide(SLIDE, merged, original, sources) == []
        run = subprocess.run([sys.executable, str(SCRIPTS / "validate-napari-label-columns.py"), "--merged-dir",
                              str(tmp / "merged"), "--original-dir", str(tmp / "stitched"), "--labels-dir",
                              str(tmp / "labels"), "--column", "Donor=donor.csv"], capture_output=True, text=True)
        print(run.stdout)
        assert run.returncode == 0 and "1 of 1 slides pass" in run.stdout


def test_a_wrong_label_a_surviving_placeholder_and_a_wrong_colour_each_fail():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); _build_real(tmp); merged, original, sources = _frames(tmp)
        wrong = merged.copy(); wrong.loc[0, "Donor"] = "UWA_999"
        assert any("disagree" in m for m in val.check_slide(SLIDE, wrong, original, sources))
        placeholder = merged.copy(); placeholder["hex_color"] = "#861933"
        assert any("placeholder" in m for m in val.check_slide(SLIDE, placeholder, original, sources))
        colour = merged.copy(); colour.loc[0, "Donor_color"] = "#000000"
        assert any("shared hash" in m for m in val.check_slide(SLIDE, colour, original, sources))
        two = merged.copy(); two.loc[2, "Donor"] = "UWA_741"; two.loc[2, "Donor_color"] = "#111111"
        msgs = val.check_slide(SLIDE, two, original, sources)
        assert any("more than one colour" in m for m in msgs), msgs


def test_every_cell_mismatching_prints_the_suspect_the_check_hint():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); _build_real(tmp); merged, original, sources = _frames(tmp)
        merged["Donor"] = "X"; merged["Donor_color"] = val.deterministic_color("X")
        msgs = val.check_slide(SLIDE, merged, original, sources)
        print(msgs)
        assert any("suspect the CHECK" in m for m in msgs)


def test_changed_row_order_and_an_unlabelled_cell_with_a_value_fail():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); _build_real(tmp); merged, original, sources = _frames(tmp)
        shuffled = merged.iloc[::-1].reset_index(drop=True)
        assert any("order differ" in m for m in val.check_slide(SLIDE, shuffled, original, sources))
        extra = merged.copy(); extra.loc[1, "Donor"] = "UWA_741"; extra.loc[1, "Donor_color"] = val.deterministic_color("UWA_741")
        assert any("disagree" in m for m in val.check_slide(SLIDE, extra, original, sources))   # cell with no source label


if __name__ == "__main__":
    test_a_file_built_by_the_real_writers_passes()
    test_a_wrong_label_a_surviving_placeholder_and_a_wrong_colour_each_fail()
    test_every_cell_mismatching_prints_the_suspect_the_check_hint()
    test_changed_row_order_and_an_unlabelled_cell_with_a_value_fail()
    print("OK")
