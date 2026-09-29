#!/usr/bin/env python3
"""Tests for validate-napari-anchor-metadata.py (no S3 / network).

A validator that only ever says "pass" proves nothing, so this is built around what it must REJECT.
The passing case is not hand-written: it is produced by running the REAL writer scripts end to end
(add-annotations -> celltypes gbmap -> celltypes --denovo-only -> celltypes leiden) on a tiny
synthetic slide, so the validator and the writers are tested against each other. Each fault test then
corrupts that good output in one way and requires a specific complaint.

Run:  uv run --with pandas --with numpy --with h5py --with duckdb \\
          python scripts/tests/test_validate_napari_anchor_metadata.py
"""
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
VALIDATOR = SCRIPTS / "validate-napari-anchor-metadata.py"
SLIDE = "SLA"
NAPARI_SLIDE_NUMBER = 5            # the <n> in c_<n>_<fov>_<cell>; deliberately not 1
_KEEP = tempfile.TemporaryDirectory()
BASE = Path(_KEEP.name)


def _py(script: str, *args) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)],
                          capture_output=True, text=True)


def _build_good_output() -> dict:
    """Run the real pipeline on 2 FOVs x 4 cells; four of the eight cells are 'anchor' cells."""
    cells = [(f, c) for f in (1, 2) for c in (1, 2, 3, 4)]
    (BASE / "flat").mkdir(); (BASE / "stitched").mkdir()
    with gzip.open(BASE / "flat" / f"{SLIDE}_metadata_file.csv.gz", "wt") as fh:
        fh.write("fov,cell_ID,cell_id\n" + "".join(f"{f},{c},c_{NAPARI_SLIDE_NUMBER}_{f}_{c}\n" for f, c in cells))
    ids = [f"c_{NAPARI_SLIDE_NUMBER}_{f}_{c}" for f, c in cells]
    pd.DataFrame({"cell_ID": ids, "cell_type": "Unassigned", "hex_color": "#861933"}
                 ).to_csv(BASE / "stitched" / f"{SLIDE}_metadata.csv", index=False)
    pd.DataFrame([{"Slide": "SLIDE A", "Case": 111, "Block": "G1", "Region": "Tumor bulk", "FOVs": 1},
                  {"Slide": "SLIDE A", "Case": 222, "Block": "G2", "Region": "Infiltrating edge", "FOVs": 2}]
                 ).to_csv(BASE / "ann.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([{"Canonical_slide_name": "SLIDE A", "AtoMx_flatfile_folder": SLIDE}]
                 ).to_csv(BASE / "xw.csv", index=False)

    anchor = {f"{SLIDE}_F1_C1": ("b", "3"), f"{SLIDE}_F1_C2": ("Neuron", "7"),
              f"{SLIDE}_F2_C3": ("aa", "3"), f"{SLIDE}_F2_C4": ("l", "12")}
    with h5py.File(BASE / "typing.h5", "w") as f:
        f.create_dataset("cell_id", data=np.array(list(anchor), dtype="S32"))
        f.create_dataset("cell_type", data=np.array([v[0] for v in anchor.values()], dtype="S16"))
        f.create_dataset("prob", data=np.full(len(anchor), 0.9))
    pd.DataFrame({"cell_id": list(anchor), "leiden": [v[1] for v in anchor.values()]}
                 ).to_csv(BASE / "leiden.csv", index=False)
    pd.DataFrame({"cell_id": list(anchor), "slide_id": SLIDE, "slide_cluster": f"{SLIDE}|1"}
                 ).to_csv(BASE / "anchor_cells.csv", index=False)

    steps = [
        ("add-annotations-to-napari-metadata.py", "--metadata-dir", BASE / "stitched", "--out-dir", BASE / "s1",
         "--annotations", BASE / "ann.csv", "--crosswalk", BASE / "xw.csv", "--columns", "Region,Case",
         "--drop-columns", "cell_type,hex_color"),
        ("celltypes-to-napari-metadata.py", "--typing-h5", BASE / "typing.h5", "--flatfiles", BASE / "flat",
         "--column", "gbmap_type", "--merge-into", BASE / "s1", "--out-dir", BASE / "s2"),
        ("celltypes-to-napari-metadata.py", "--typing-h5", BASE / "typing.h5", "--flatfiles", BASE / "flat",
         "--column", "gbmap_denovo", "--denovo-only", "--merge-into", BASE / "s2", "--out-dir", BASE / "s3"),
        ("labels-csv-to-typing-h5.py", "--csv", BASE / "leiden.csv", "--label-column", "leiden",
         "--label-prefix", "leiden_", "--output", BASE / "leiden.h5"),
        ("celltypes-to-napari-metadata.py", "--typing-h5", BASE / "leiden.h5", "--flatfiles", BASE / "flat",
         "--column", "leiden", "--merge-into", BASE / "s3", "--out-dir", BASE / "s4"),
    ]
    for step in steps:
        r = _py(*step)
        assert r.returncode == 0, f"{step[0]} failed:\n{r.stdout}\n{r.stderr}"
    return {"merged": BASE / "s4", "orig": BASE / "stitched"}


_GOOD = _build_good_output()


def _validate(tmp: Path, merged: Path | None = None, anchor_csv: Path | None = None):
    return subprocess.run([sys.executable, str(VALIDATOR),
                           "--merged-dir", str(merged or _GOOD["merged"]),
                           "--original-dir", str(_GOOD["orig"]), "--typing-h5", str(BASE / "typing.h5"),
                           "--leiden-csv", str(BASE / "leiden.csv"),
                           "--anchor-cells-csv", str(anchor_csv or BASE / "anchor_cells.csv"),
                           "--annotations", str(BASE / "ann.csv"), "--crosswalk", str(BASE / "xw.csv")],
                          capture_output=True, text=True)


def _corrupted(tmp: Path, edit) -> Path:
    """Copy the good merged dir and apply `edit(df) -> df` to the slide file."""
    dst = tmp / "merged"; shutil.copytree(_GOOD["merged"], dst)
    f = dst / f"{SLIDE}_metadata.csv"
    df = pd.read_csv(f, dtype=str, keep_default_na=False)
    edit(df).to_csv(f, index=False)
    return dst


def _rejects(edit, expect: str):
    with tempfile.TemporaryDirectory() as d:
        out = _validate(Path(d), merged=_corrupted(Path(d), edit))
        both = out.stdout + out.stderr
        assert out.returncode == 1, f"a corrupted file was ACCEPTED:\n{both}"
        assert expect in both, f"expected {expect!r} in:\n{both}"


def test_the_real_writers_output_passes():
    """The good case is built by the real scripts, so this pins validator <-> writers agreement."""
    with tempfile.TemporaryDirectory() as d:
        out = _validate(Path(d))
        assert out.returncode == 0, out.stdout + out.stderr
        assert "ALL 1 SLIDES PASS EVERY CHECK" in out.stdout


def test_a_wrong_gbmap_label_is_caught():
    _rejects(lambda df: df.assign(gbmap_type=df.gbmap_type.replace("Neuron", "Astrocyte")),
             "gbmap_type: 1 cells differ from a direct lookup")


def test_a_missing_leiden_prefix_triggers_the_check_bug_hint():
    """Every typed cell wrong == the typed count; the report must say to suspect the check."""
    _rejects(lambda df: df.assign(leiden=df.leiden.str.replace("leiden_", "", regex=False)),
             "signature of a bug in the CHECK")


def test_reordered_rows_are_caught():
    _rejects(lambda df: df.iloc[::-1].reset_index(drop=True), "cell_ID order differs")


def test_dropped_rows_are_caught():
    _rejects(lambda df: df.iloc[:-1], "rows 7 != original 8")


def test_a_wrong_region_is_caught():
    _rejects(lambda df: df.assign(Region=df.Region.replace("Tumor bulk", "Infiltrating edge")),
             "Region:")


def test_a_stray_label_on_an_untyped_cell_is_caught():
    def edit(df):
        df.loc[df.gbmap_type == "", "leiden"] = "leiden_99"
        return df
    _rejects(edit, "label with no gbmap_type call")


def test_a_letters_column_that_disagrees_with_the_full_column_is_caught():
    def edit(df):
        df.loc[df.gbmap_type == "Neuron", "gbmap_denovo"] = "Neuron"    # a NAMED type in the letters column
        return df
    _rejects(edit, "gbmap_denovo is not")


def test_a_colour_that_is_not_the_shared_hash_is_caught():
    def edit(df):
        df.loc[df.gbmap_type == "b", "gbmap_type_color"] = "#000000"
        return df
    _rejects(edit, "!= shared hash")


def test_a_missing_colour_column_is_caught():
    _rejects(lambda df: df.drop(columns=["leiden_color"]), "missing columns")


def test_a_typed_count_that_disagrees_with_the_anchor_is_caught():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        bad = pd.read_csv(BASE / "anchor_cells.csv")
        bad = pd.concat([bad, bad.iloc[:1].assign(cell_id="SLA_F1_C9")])          # 5 anchor cells vs 4 typed
        bad.to_csv(tmp / "anchor.csv", index=False)
        out = _validate(tmp, anchor_csv=tmp / "anchor.csv")
        assert out.returncode == 1 and "!= anchor cells on this slide" in out.stdout + out.stderr


def test_annotations_without_a_crosswalk_is_refused():
    out = subprocess.run([sys.executable, str(VALIDATOR), "--merged-dir", str(_GOOD["merged"]),
                          "--typing-h5", str(BASE / "typing.h5"), "--annotations", str(BASE / "ann.csv")],
                         capture_output=True, text=True)
    assert out.returncode != 0 and "needs --crosswalk" in out.stdout + out.stderr


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
    _KEEP.cleanup()
