"""Tests for recolor-napari-metadata.py.

Run:  uv run --with duckdb --with pandas --with pytest python -m pytest scripts/tests/test_recolor_napari_metadata.py -s
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
SCRIPT = SCRIPTS / "recolor-napari-metadata.py"
spec = importlib.util.spec_from_file_location("recolor", SCRIPT)
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

DARK = "#0A0404"
VISIBLE = "#F5C242"


def _frame():
    return pd.DataFrame({
        "cell_ID": ["1", "2", "3"],
        "Donor": ["741", "741", "6803"], "Donor_color": [VISIBLE, VISIBLE, VISIBLE],
        "InSituTree_v5": ["L4", "L4", "Microglia A"], "InSituTree_v5_color": [DARK, DARK, VISIBLE],
        "InSituTree_v5_collapsed": ["L4", "L4", "Microglia A"], "InSituTree_v5_collapsed_color": [DARK, DARK, VISIBLE],
        "cell_type": ["Unassigned"] * 3, "hex_color": [None, None, None],
    })


def test_only_dark_colours_change_in_every_colour_column():
    out, changes = rc.recolor_frame(_frame(), "s.csv")
    print(pd.DataFrame(changes).to_string())
    assert out["InSituTree_v5_color"].iloc[0] != DARK and out["InSituTree_v5_color"].iloc[2] == VISIBLE
    assert out["InSituTree_v5_color"].iloc[0] == out["InSituTree_v5_collapsed_color"].iloc[0]   # one colour per label
    assert out["Donor_color"].tolist() == [VISIBLE] * 3
    assert {c["column"] for c in changes} == {"InSituTree_v5_color", "InSituTree_v5_collapsed_color"}
    assert all(c["labels"] == "L4" and c["cells"] == 2 for c in changes)


def test_labels_and_other_columns_are_untouched():
    frame = _frame()
    out, _ = rc.recolor_frame(frame, "s.csv")
    for column in frame.columns:
        if not column.endswith("_color"):
            assert out[column].equals(frame[column]), column


def test_running_twice_changes_nothing_more():
    once, _ = rc.recolor_frame(_frame(), "s.csv")
    twice, changes = rc.recolor_frame(once, "s.csv")
    assert twice.equals(once) and changes == []


def test_cli_round_trips_files_and_lists_the_changes(tmp_path):
    in_dir, out_dir = tmp_path / "in", tmp_path / "out"
    in_dir.mkdir()
    _frame().to_csv(in_dir / "slideA_metadata.csv", index=False)
    clean = _frame().assign(InSituTree_v5_color=VISIBLE, InSituTree_v5_collapsed_color=VISIBLE)
    clean.to_csv(in_dir / "slideB_metadata.csv", index=False)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--in-dir", str(in_dir), "--out-dir", str(out_dir),
                           "--changes-csv", str(tmp_path / "changes.csv")], capture_output=True, text=True)
    print(proc.stdout, proc.stderr)
    assert proc.returncode == 0, proc.stderr
    changes = pd.read_csv(tmp_path / "changes.csv")
    assert set(changes["file"]) == {"slideA_metadata.csv"}
    assert (in_dir / "slideB_metadata.csv").read_bytes() == (out_dir / "slideB_metadata.csv").read_bytes()
    fixed = pd.read_csv(out_dir / "slideA_metadata.csv", dtype=str)
    assert fixed["InSituTree_v5_color"].iloc[0] != DARK
    assert list(fixed.columns) == list(_frame().columns)
