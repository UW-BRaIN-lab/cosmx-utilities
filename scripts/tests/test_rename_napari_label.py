"""Tests for rename-napari-label.py.

Run:  uv run --with duckdb --with pandas --with pytest python -m pytest scripts/tests/test_rename_napari_label.py -s
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
SCRIPT = SCRIPTS / "rename-napari-label.py"
spec = importlib.util.spec_from_file_location("rename_label", SCRIPT)
rl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rl)
cspec = importlib.util.spec_from_file_location("napari_colors", SCRIPTS / "napari_colors.py")
nc = importlib.util.module_from_spec(cspec)
cspec.loader.exec_module(nc)

OLD, NEW = "Excitatory_upper", "Excitatory_L2_3_mixed"
COLUMNS = ["InSituTree_v5", "InSituTree_v5_collapsed"]


def _frame():
    return pd.DataFrame({
        "cell_ID": ["1", "2", "3", "4"],
        "InSituTree_v5": [OLD, "Microglia A", OLD, ""], "InSituTree_v5_color": ["#111111", "#F5C242", "#111111", ""],
        "InSituTree_v5_collapsed": [OLD, "Microglia A", "L4 neuron", ""],
        "InSituTree_v5_collapsed_color": ["#111111", "#F5C242", "#222222", ""],
        "InSituTree_v4": [OLD, OLD, OLD, OLD], "InSituTree_v4_color": ["#111111"] * 4,
    })


def test_renames_only_the_named_columns_and_recolours_by_the_shared_rule():
    out, renamed = rl.rename_label(_frame(), COLUMNS, OLD, NEW)
    print(out.to_string(), renamed)
    assert renamed == {"InSituTree_v5": 2, "InSituTree_v5_collapsed": 1}
    assert out["InSituTree_v5"].tolist() == [NEW, "Microglia A", NEW, ""]
    assert out["InSituTree_v5_color"].tolist() == [nc.deterministic_color(NEW), "#F5C242", nc.deterministic_color(NEW), ""]
    assert out["InSituTree_v5_collapsed_color"].iloc[2] == "#222222"          # other labels keep their colour
    assert out["InSituTree_v4"].tolist() == [OLD] * 4                           # a column not named is untouched
    assert out["InSituTree_v4_color"].tolist() == ["#111111"] * 4


def test_refuses_to_merge_into_an_existing_label():
    frame = _frame()
    frame.loc[1, "InSituTree_v5"] = NEW
    with pytest.raises(ValueError, match="merge"):
        rl.rename_label(frame, COLUMNS, OLD, NEW)


def test_missing_column_is_an_error():
    with pytest.raises(KeyError):
        rl.rename_label(_frame(), ["NotAColumn"], OLD, NEW)


def test_cli_round_trip_is_byte_identical_when_the_label_is_absent(tmp_path):
    in_dir, out_dir = tmp_path / "in", tmp_path / "out"
    in_dir.mkdir()
    _frame().to_csv(in_dir / "s_metadata.csv", index=False)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--in-dir", str(in_dir), "--out-dir", str(out_dir),
                           "--column", "InSituTree_v5", "--old", "NotThere", "--new", NEW], capture_output=True, text=True)
    print(proc.stdout, proc.stderr)
    assert proc.returncode == 0, proc.stderr
    assert (in_dir / "s_metadata.csv").read_bytes() == (out_dir / "s_metadata.csv").read_bytes()
