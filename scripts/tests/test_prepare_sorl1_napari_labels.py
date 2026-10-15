#!/usr/bin/env python3
"""Tests for prepare-sorl1-napari-labels.py (no S3 / network).

Pins the label rules: variant-level SORL1 status that stops on a carrier disagreement, de novo letters turned into
their annotation text (and a missing annotation refused), the OPC confidence column that is read from the table
and never recomputed, and that only OPC-called cells appear in the OPC file.

Run:  uv run --with pandas --with numpy python scripts/tests/test_prepare_sorl1_napari_labels.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "prepare-sorl1-napari-labels.py"
spec = importlib.util.spec_from_file_location("prep", SCRIPT)
prep = importlib.util.module_from_spec(spec); spec.loader.exec_module(prep)


def _table():
    return pd.DataFrame({
        "donor": ["741", "7181", "7796", "6803"],
        "sorl1_mutation": ["AD+LATE", "AD+LATE SORL1 C467Y", "AD+LATE SORL1 R953C", "AD Pure"],
        "case_group": ["AD+LATE", "AD+LATE SORL1", "AD+LATE SORL1", "AD Pure"],
        "celltype_k21": ["o", "L4 neuron", "i", "Endothelial"],
        "celltype_tree_v2": ["Oligodendrocyte_denovo", "L4 neuron", "OPC", "OPC"],
        "opc_confident_tree_v2": [False, False, True, False],
    }, index=pd.Index(["S_F1_C1", "S_F1_C2", "S_F2_C1", "S_F2_C2"], name="cell_id"))


def _annotation(tmp: Path) -> Path:
    p = tmp / "ann.csv"
    pd.DataFrame({"denovo_label": ["o", "i"], "annotation": ["o - Oligodendrocyte", "i - Unresolved: low-information"]}
                 ).to_csv(p, index=False)
    return p


def test_sorl1_status_is_variant_level_and_checked_against_case_group():
    t = _table()
    assert prep.sorl1_status(t["sorl1_mutation"], t["case_group"]).tolist() == [
        "No SORL1 variant", "SORL1 C467Y", "SORL1 R953C", "No SORL1 variant"]
    bad = t.copy(); bad.loc["S_F1_C1", "case_group"] = "AD+LATE SORL1"      # carrier by group, not by mutation
    try:
        prep.sorl1_status(bad["sorl1_mutation"], bad["case_group"])
    except ValueError as exc:
        assert "disagree" in str(exc)
    else:
        raise AssertionError("a carrier disagreement must stop the run")


def test_k21_letters_become_their_reading_and_reference_names_stay():
    with tempfile.TemporaryDirectory() as d:
        ann = _annotation(Path(d))
        got = prep.k21_names(_table()["celltype_k21"], ann).tolist()
        print(got)
        assert got == ["o - Oligodendrocyte", "L4 neuron", "i - Unresolved: low-information", "Endothelial"]
        try:
            prep.k21_names(pd.Series(["zz"]), ann)
        except ValueError as exc:
            assert "zz" in str(exc)
        else:
            raise AssertionError("a letter with no annotation must be refused")


def test_opc_confidence_covers_only_opc_calls_and_needs_the_flag_column():
    t = _table()
    got = prep.opc_confidence(t, "tree_v2", ["OPC"])
    print(got.tolist())
    assert got.tolist()[:2] == [pd.NA, pd.NA] or got.iloc[:2].isna().all()
    assert got.iloc[2] == prep.OPC_CONFIDENT and got.iloc[3] == prep.OPC_LOW
    try:
        prep.opc_confidence(t.drop(columns="opc_confident_tree_v2"), "tree_v2", ["OPC"])
    except ValueError as exc:
        assert "opc_confident_tree_v2" in str(exc) and "OPC_FLAG" in str(exc)
    else:
        raise AssertionError("a missing flag column must say how to make it")


def test_cli_writes_five_files_and_the_opc_file_has_only_opc_cells():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d); _table().to_csv(tmp / "labels.csv.gz", compression="gzip")
        run = subprocess.run([sys.executable, str(SCRIPT), "--labels", str(tmp / "labels.csv.gz"),
                              "--k21-annotations", str(_annotation(tmp)), "--out-dir", str(tmp / "o")],
                             capture_output=True, text=True)
        print(run.stdout, run.stderr[-300:])
        assert run.returncode == 0
        assert sorted(p.name for p in (tmp / "o").glob("*.csv")) == sorted(prep.OUTPUTS.values())
        opc = pd.read_csv(tmp / "o" / "opc_confident.csv")
        assert opc["cell_id"].tolist() == ["S_F2_C1", "S_F2_C2"]
        assert pd.read_csv(tmp / "o" / "donor.csv")["label"].tolist() == ["UWA_741", "UWA_7181", "UWA_7796", "UWA_6803"]


if __name__ == "__main__":
    test_sorl1_status_is_variant_level_and_checked_against_case_group()
    test_k21_letters_become_their_reading_and_reference_names_stay()
    test_opc_confidence_covers_only_opc_calls_and_needs_the_flag_column()
    test_cli_writes_five_files_and_the_opc_file_has_only_opc_cells()
    print("OK")
