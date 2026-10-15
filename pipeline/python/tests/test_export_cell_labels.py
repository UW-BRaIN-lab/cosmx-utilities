#!/usr/bin/env python3
"""Tests for export_cell_labels.py (no S3 / network).

Plants a run that types only some cells and checks: ids are joined (not positions), the donor
annotation lands on the right cells, an untyped cell is NA with a warning, UMAP is carried, and
the CLI writes a readable gzip CSV with the expected column pairs.

Run:  uv run --with pandas --with numpy --with anndata --with scipy --with h5py python pipeline/python/tests/test_export_cell_labels.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "export_cell_labels.py"
sys.path.insert(0, str(SCRIPT.parent))
from export_cell_labels import build_table, opc_confident  # noqa: E402

CELLS = ["c0", "c1", "c2", "c3"]


def _obs():
    return pd.DataFrame({"slide_id": ["SA", "SA", "SB", "SB"], "Region": ["741", "741", "6803", "6803"],
                         "fov": ["1", "2", "1", "2"], "leiden": ["0", "1", "0", "1"],
                         "total_counts": [100.0, 200, 300, 400], "qc_genes_detected": [50.0, 60, 70, 80],
                         "qc_area": [1.0, 2, 3, 4]}, index=CELLS)


ANN = pd.DataFrame({"donor": ["741", "6803"], "case_broad": ["AD+LATE", "AD Pure"],
                    "case_group": ["AD+LATE", "AD Pure"], "sorl1_mutation": ["AD+LATE", "AD Pure"]})


def test_join_is_by_cell_id_not_position(capsys=None):
    # result lists cells in a DIFFERENT order and omits c3
    run = pd.DataFrame({"cell_type": ["x", "y", "z"], "prob": [1.0, 0.9, 0.8]},
                       index=pd.Index(["c2", "c0", "c1"], name="cell_id"))
    table = build_table(_obs(), np.arange(8.0).reshape(4, 2), ANN, {"k": run})
    print(table[["donor", "case_group", "celltype_k", "prob_k", "umap1"]])
    assert table.loc["c0", "celltype_k"] == "y" and table.loc["c2", "celltype_k"] == "x"
    assert pd.isna(table.loc["c3", "celltype_k"])
    assert table.loc["c2", "case_group"] == "AD Pure" and table.loc["c0", "case_group"] == "AD+LATE"
    assert table.loc["c1", "umap1"] == 2.0 and table.loc["c3", "umap2"] == 7.0


def test_opc_confident_needs_the_call_and_leiden_or_depth():
    table = pd.DataFrame({
        "leiden": ["0", "8", "8", "0", "2"],
        "total_counts": [100, 700, 100, 50, 900],
        "celltype_tree": ["OPC", "OPC", "OPC", "Oligodendrocyte", "OPC"]})
    flag = opc_confident(table, "tree", ["OPC"]).tolist()
    print(flag)
    # Leiden-0 shallow OPC call: yes | deep Leiden-8 OPC call: yes | shallow Leiden-8 OPC call: NO |
    # Leiden 0 but not called OPC: no | deep OPC call from Leiden 2: yes (depth alone qualifies)
    assert flag == [True, True, False, False, True]
    assert opc_confident(table, "tree", ["OPC"], min_counts=1000).tolist() == [True, False, False, False, False]


def test_opc_flag_column_and_unknown_run():
    obs = _obs()
    run = pd.DataFrame({"cell_type": ["OPC", "OPC", "x", "OPC"], "prob": 1.0}, index=pd.Index(CELLS, name="cell_id"))
    table = build_table(obs, np.zeros((4, 2)), ANN, {"tree": run}, {"tree": ["OPC"]})
    # c0 leiden 0 -> True; c1 leiden 1 and counts 200 -> False; c2 not OPC -> False; c3 leiden 1, counts 400 -> False
    assert table["opc_confident_tree"].tolist() == [True, False, False, False]
    try:
        build_table(obs, np.zeros((4, 2)), ANN, {"tree": run}, {"nope": ["OPC"]})
    except ValueError as exc:
        assert "nope" in str(exc)
    else:
        raise AssertionError("an unknown run in --opc-flag must fail loudly")


def test_cli_end_to_end():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        ad.AnnData(X=sp.csr_matrix((4, 2)), obs=_obs(), obsm={"X_umap": np.arange(8.0).reshape(4, 2)}
                   ).write_h5ad(tmp / "c.h5ad")
        with h5py.File(tmp / "r.h5", "w") as f:
            f["cell_id"] = np.array(CELLS, dtype="S")
            f["cell_type"] = np.array(["a", "b", "a", "b"], dtype="S")
            f["prob"] = np.array([1.0, 1.0, 0.5, 1.0])
        ANN.to_csv(tmp / "a.csv", index=False)
        out = subprocess.run([sys.executable, str(SCRIPT), "--clustered-h5ad", str(tmp / "c.h5ad"),
                              "--annotations", str(tmp / "a.csv"), "--output", str(tmp / "o.csv.gz"),
                              "--run", f"k21={tmp / 'r.h5'}"], capture_output=True, text=True)
        print(out.stdout, out.stderr[-300:])
        assert out.returncode == 0
        table = pd.read_csv(tmp / "o.csv.gz", index_col=0)
        assert {"celltype_k21", "prob_k21", "umap1", "case_group"} <= set(table.columns)
        assert table.loc["c2", "celltype_k21"] == "a" and len(table) == 4


if __name__ == "__main__":
    test_join_is_by_cell_id_not_position()
    test_opc_confident_needs_the_call_and_leiden_or_depth()
    test_opc_flag_column_and_unknown_run()
    test_cli_end_to_end()
    print("OK")
