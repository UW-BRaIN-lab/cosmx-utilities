#!/usr/bin/env python3
"""Tests for write_celltypes.py (no S3 / network).

Builds a 5-cell clustered AnnData and a typing result that mixes leaf labels with the tree's
node-level label, then runs the real CLI with the committed v5 rename / collapse / parent maps and
checks: cell_type is untouched, each cell_type_<name> column maps as the CSV says, the donor columns
land on the right cells, opc_confident follows the Leiden-or-depth rule, and an untyped cell is NA.

Run:  uv run --with pandas --with numpy --with anndata --with scipy --with h5py --with matplotlib python pipeline/python/tests/test_write_celltypes.py
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

PIPELINE_DIR = Path(__file__).resolve().parents[2]
SCRIPT = PIPELINE_DIR / "python" / "write_celltypes.py"
REFERENCE_DIR = PIPELINE_DIR / "reference"
RENAME_MAP = REFERENCE_DIR / "sorl1_hierarchy_v5_rename.csv"
COLLAPSE_MAP = REFERENCE_DIR / "sorl1_hierarchy_v5_collapse.csv"
PARENT_MAP = REFERENCE_DIR / "sorl1_hierarchy_v5_parent.csv"
ANNOTATIONS = REFERENCE_DIR / "fov_annotations_sorl1.csv"

CELLS = ["c0", "c1", "c2", "c3", "c4"]
OPC_LABEL = "Oligodendrocyte precursor cell"
TYPED = {  # c4 is deliberately not typed
    "c0": "Excitatory_upper",
    "c1": "Neuron_unresolved_d_denovo",
    "c2": OPC_LABEL,
    "c3": OPC_LABEL,
}


def _write_inputs(tmp: Path) -> tuple[Path, Path, str, str]:
    donors = pd.read_csv(ANNOTATIONS, dtype=str).drop_duplicates("donor")
    donor_a, donor_b = donors["donor"].iloc[0], donors["donor"].iloc[-1]
    obs = pd.DataFrame({"Region": [donor_a, donor_a, donor_b, donor_b, donor_b],
                        "leiden": ["3", "3", "0", "5", "5"],
                        "total_counts": [100.0, 100.0, 50.0, 900.0, 50.0]}, index=CELLS)
    adata = ad.AnnData(X=sp.csr_matrix(np.ones((5, 2), dtype=np.float32)), obs=obs)
    adata.obsm["X_umap"] = np.zeros((5, 2), dtype=np.float32)
    clustered = tmp / "clustered.h5ad"
    adata.write_h5ad(clustered)
    result = tmp / "result.h5"
    with h5py.File(result, "w") as f:
        f["cell_id"] = np.array(list(TYPED), dtype="S")
        f["cell_type"] = np.array(list(TYPED.values()), dtype="S")
        f["prob"] = np.array([0.9, 0.8, 0.7, 0.99])
    return clustered, result, donor_a, donor_b


def test_cli_writes_final_label_columns_donor_groups_and_opc_flag():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        clustered, result, donor_a, donor_b = _write_inputs(tmp)
        out = tmp / "typed.h5ad"
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "--clustered-h5ad", str(clustered), "--insitutype-h5", str(result),
             "--output", str(out), "--label-map", f"v5={RENAME_MAP}", "--label-map", f"v5_collapsed={COLLAPSE_MAP}",
             "--label-map", f"v5_parent={PARENT_MAP}", "--annotations", str(ANNOTATIONS),
             "--opc-labels", OPC_LABEL], capture_output=True, text=True)
        print(proc.stdout, proc.stderr)
        assert proc.returncode == 0, proc.stderr
        obs = ad.read_h5ad(out).obs
        print(obs.to_string())

        assert obs.loc["c0", "cell_type"] == "Excitatory_upper"           # the run's own label is untouched
        assert obs.loc["c0", "cell_type_v5"] == "Excitatory_L2_3_mixed"   # the final name
        assert obs.loc["c0", "cell_type_v5_collapsed"] == "Excitatory_L2_3_mixed"
        assert obs.loc["c0", "cell_type_v5_parent"] == "Excitatory_L2_3_mixed"
        assert obs.loc["c1", "cell_type_v5"] == "Neuron_unresolved_d_denovo"
        assert obs.loc["c1", "cell_type_v5_collapsed"] == "Neuron_low_signal"
        assert obs.loc["c2", "cell_type_v5_parent"] == "OPC"
        assert "Excitatory_upper" not in set(obs["cell_type_v5"].dropna())
        assert pd.isna(obs.loc["c4", "cell_type"]) and pd.isna(obs.loc["c4", "cell_type_v5"])

        donors = pd.read_csv(ANNOTATIONS, dtype=str).drop_duplicates("donor").set_index("donor")
        assert obs.loc["c0", "case_group"] == donors.loc[donor_a, "case_group"]
        assert obs.loc["c3", "sorl1_mutation"] == donors.loc[donor_b, "sorl1_mutation"]

        # c2: OPC in Leiden 0 (kept); c3: OPC, Leiden 5, 900 counts (kept); c0/c1 are not OPC; c4 untyped
        assert obs["opc_confident"].tolist() == [False, False, True, True, False]


def test_without_the_optional_flags_only_the_two_original_columns_are_added():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        clustered, result, _, _ = _write_inputs(tmp)
        out = tmp / "typed.h5ad"
        proc = subprocess.run([sys.executable, str(SCRIPT), "--clustered-h5ad", str(clustered),
                               "--insitutype-h5", str(result), "--output", str(out)], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        added = set(ad.read_h5ad(out).obs.columns) - {"Region", "leiden", "total_counts"}
        print(sorted(added))
        assert added == {"cell_type", "insitutype_prob"}


if __name__ == "__main__":
    test_cli_writes_final_label_columns_donor_groups_and_opc_flag()
    test_without_the_optional_flags_only_the_two_original_columns_are_added()
    print("OK")
