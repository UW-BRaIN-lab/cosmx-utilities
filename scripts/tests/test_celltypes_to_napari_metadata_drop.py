#!/usr/bin/env python3
"""Tests for --drop-columns in celltypes-to-napari-metadata.py (no S3 / network).

Runs the REAL script on a tiny synthetic slide: a stitched `_metadata.csv` carrying the stitcher's placeholder
`cell_type` ("Unassigned") and `hex_color`, a flat-file metadata gz to join through, and a label h5 built by the real
labels-csv-to-typing-h5.py. Pins: both placeholders are dropped, rows and order are untouched, our column arrives
with its `_color`, an untyped cell stays blank, the join key can never be dropped, and WITHOUT the flag the old
behaviour (placeholders kept) is unchanged.

Run:  uv run --with pandas --with numpy --with h5py --with duckdb python scripts/tests/test_celltypes_to_napari_metadata_drop.py
"""
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
SLIDE = "SLX"


def _setup(tmp: Path) -> dict:
    flat = tmp / "flat"; flat.mkdir()
    rows = "fov,cell_ID,cell_id\n1,1,c_3_1_1\n1,2,c_3_1_2\n2,1,c_3_2_1\n"
    with gzip.open(flat / f"{SLIDE}_metadata_file.csv.gz", "wt") as f:
        f.write(rows)
    stitched = tmp / "stitched"; stitched.mkdir()
    pd.DataFrame({"cell_ID": ["c_3_1_1", "c_3_1_2", "c_3_2_1"], "cell_type": "Unassigned",
                  "hex_color": "#861933"}).to_csv(stitched / f"{SLIDE}_metadata.csv", index=False)
    pd.DataFrame({"cell_id": [f"{SLIDE}_F1_C1", f"{SLIDE}_F2_C1"], "label": ["UWA_741", "UWA_6803"]}
                 ).to_csv(tmp / "donor.csv", index=False)           # cell F1_C2 deliberately has no label
    h5 = tmp / "donor.h5"
    done = subprocess.run([sys.executable, str(SCRIPTS / "labels-csv-to-typing-h5.py"), "--csv", str(tmp / "donor.csv"),
                           "--label-column", "label", "--output", str(h5)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return {"flat": flat, "stitched": stitched, "h5": h5}


def _run(tmp, paths, *extra):
    out = tmp / "out"
    return out, subprocess.run(
        [sys.executable, str(SCRIPTS / "celltypes-to-napari-metadata.py"), "--typing-h5", str(paths["h5"]),
         "--flatfiles", str(paths["flat"]), "--column", "Donor", "--merge-into", str(paths["stitched"]),
         "--out-dir", str(out), *extra], capture_output=True, text=True)


def test_placeholders_are_dropped_and_the_new_column_arrives():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); paths = _setup(tmp)
        out, run = _run(tmp, paths, "--drop-columns", "cell_type,hex_color")
        print(run.stdout[-300:], run.stderr[-300:])
        assert run.returncode == 0
        got = pd.read_csv(out / f"{SLIDE}_metadata.csv", dtype=str, keep_default_na=False)
        print(got)
        assert list(got.columns) == ["cell_ID", "Donor", "Donor_color"]
        assert got["cell_ID"].tolist() == ["c_3_1_1", "c_3_1_2", "c_3_2_1"]          # rows and order untouched
        assert got["Donor"].tolist() == ["UWA_741", "", "UWA_6803"]                  # unlabelled cell stays blank
        assert got.loc[1, "Donor_color"] == "" and got.loc[0, "Donor_color"].startswith("#")


def test_without_the_flag_the_placeholders_are_kept():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); paths = _setup(tmp)
        out, run = _run(tmp, paths)
        assert run.returncode == 0
        got = pd.read_csv(out / f"{SLIDE}_metadata.csv", dtype=str, keep_default_na=False)
        assert {"cell_type", "hex_color", "Donor", "Donor_color"} <= set(got.columns)


def test_the_join_key_cannot_be_dropped():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); paths = _setup(tmp)
        _, run = _run(tmp, paths, "--drop-columns", "cell_ID")
        print(run.stderr[-200:])
        assert run.returncode != 0 and "join key" in run.stderr


if __name__ == "__main__":
    test_placeholders_are_dropped_and_the_new_column_arrives()
    test_without_the_flag_the_placeholders_are_kept()
    test_the_join_key_cannot_be_dropped()
    print("OK")
