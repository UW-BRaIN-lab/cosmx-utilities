#!/usr/bin/env python3
"""Tests for labels-csv-to-typing-h5.py (no S3 / network).

The one that matters is round-tripping through the SAME reader
celltypes-to-napari-metadata.py uses: if anchor_profiles.read_cell_calls cannot read what this
writes, the conversion is worthless however well-formed the file looks. The rest guard the
silent-corruption cases -- truncated ids would break the Napari join with no error, and
duplicated ids would make the join pick a row arbitrarily.

Run:  uv run --with pandas --with numpy --with h5py python scripts/tests/test_labels_csv_to_typing_h5.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "labels-csv-to-typing-h5.py"
# A long realistic pipeline index; S64 would be fine but a fixed guess is what we are avoiding.
LONG_ID = "Slide37749G17706G5_F188_C1234567"


def _csv(tmp: Path, rows) -> Path:
    p = tmp / "labels.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _run(csv: Path, out: Path, *extra) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), "--csv", str(csv),
                           "--label-column", "cluster", "--output", str(out), *extra],
                          capture_output=True, text=True)


def _read(path: Path):
    with h5py.File(path, "r") as f:
        dec = lambda d: [x.decode() if isinstance(x, bytes) else str(x) for x in d[()]]
        return dec(f["cell_id"]), dec(f["cell_type"])


def test_ids_and_labels_survive_the_round_trip():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        csv = _csv(tmp, [{"cell_id": LONG_ID, "slide_id": "x", "cluster": 12},
                         {"cell_id": "SL_F1_C2", "slide_id": "x", "cluster": 3}])
        out = tmp / "leiden.h5"
        assert _run(csv, out, "--label-prefix", "leiden_").returncode == 0
        ids, labels = _read(out)
        assert ids == [LONG_ID, "SL_F1_C2"], ids       # no truncation of the long id
        assert labels == ["leiden_12", "leiden_3"], labels


def test_it_is_readable_by_the_same_reader_the_napari_script_uses():
    """The only test that proves the output is actually usable downstream."""
    py_dir = Path(__file__).resolve().parents[2] / "pipeline" / "python"
    if not (py_dir / "anchor_profiles.py").exists():
        print("  (skipped: pipeline/python/anchor_profiles.py not present)")
        return
    sys.path.insert(0, str(py_dir))
    from anchor_profiles import read_cell_calls
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        csv = _csv(tmp, [{"cell_id": LONG_ID, "cluster": 12}, {"cell_id": "SL_F1_C2", "cluster": 3}])
        out = tmp / "leiden.h5"
        assert _run(csv, out, "--label-prefix", "leiden_").returncode == 0
        got = read_cell_calls(out)
    assert list(got["cell_id"]) == [LONG_ID, "SL_F1_C2"]
    assert list(got["cell_type"]) == ["leiden_12", "leiden_3"]


def test_duplicate_ids_are_refused():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        csv = _csv(tmp, [{"cell_id": "A_F1_C1", "cluster": 1},
                         {"cell_id": "A_F1_C1", "cluster": 2}])
        out = _run(csv, tmp / "x.h5")
    assert out.returncode != 0
    assert "duplicated cell ids" in (out.stdout + out.stderr)


def test_a_missing_column_is_named():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        csv = _csv(tmp, [{"cell_id": "A_F1_C1", "leiden": 1}])
        out = _run(csv, tmp / "x.h5")
    assert out.returncode != 0
    assert "no column 'cluster'" in (out.stdout + out.stderr)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
