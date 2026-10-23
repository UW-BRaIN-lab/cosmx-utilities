#!/usr/bin/env python3
"""Tests for preprocess_h5ad_profiles.py on a tiny synthetic CELLxGENE-style H5AD.

Checks the per-type mean against a hand computation, duplicate-symbol summing,
Ensembl->feature_name mapping, --obs-filter, --min-cells, and that a
log-normalized X is refused when no raw slot exists.

Run:  uv run --with pandas --with numpy --with anndata --with scipy python pipeline/reference_builder/tests/test_preprocess_h5ad_profiles.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

SCRIPT = Path(__file__).resolve().parents[1] / "preprocess_h5ad_profiles.py"
GENE_SYMBOLS = ["SLC17A7", "GAD1", "MBP", "MBP"]  # MBP duplicated on purpose
ENSEMBL_IDS = ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]
CELL_TYPES = ["L2/3 IT"] * 4 + ["Sst"] * 3 + ["Rare"] * 1
DISEASE = ["normal", "normal", "dementia", "normal", "normal", "dementia", "normal", "normal"]


def make_h5ad(path: Path, counts: np.ndarray, raw_in_raw: bool = False) -> None:
    obs = pd.DataFrame({"Supertype": CELL_TYPES, "disease": DISEASE},
                       index=[f"c{i}" for i in range(len(CELL_TYPES))])
    var = pd.DataFrame({"feature_name": GENE_SYMBOLS}, index=ENSEMBL_IDS)
    adata = ad.AnnData(X=csr_matrix(counts), obs=obs, var=var)
    if raw_in_raw:
        adata.raw = ad.AnnData(X=csr_matrix(counts), var=var)
        adata.X = csr_matrix(np.log1p(counts))
    adata.write_h5ad(path)


def run(*args: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, "-I", str(SCRIPT), *args]
    print("  $", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-1500:], res.stderr[-1500:])
    return res


def main() -> None:
    rng = np.random.default_rng(0)
    counts = rng.integers(0, 20, size=(len(CELL_TYPES), len(GENE_SYMBOLS))).astype(np.float32)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        h5ad, out = tmp / "a.h5ad", tmp / "out.csv"
        make_h5ad(h5ad, counts)
        base = ["--h5ad", str(h5ad), "--cell-type-col", "Supertype", "--min-cells", "2",
                "--chunk-size", "3", "--prefix", "SEAAD_", "--output", str(out)]

        print("test_means_dups_and_min_cells")
        assert run(*base).returncode == 0
        prof = pd.read_csv(out, index_col=0)
        assert list(prof.columns) == ["SEAAD_L2/3 IT", "SEAAD_Sst"], prof.columns  # Rare dropped
        assert list(prof.index) == ["GAD1", "MBP", "SLC17A7"], prof.index
        expect_l23 = counts[:4].mean(axis=0)
        np.testing.assert_allclose(prof.loc["SLC17A7", "SEAAD_L2/3 IT"], expect_l23[0])
        np.testing.assert_allclose(prof.loc["MBP", "SEAAD_L2/3 IT"], expect_l23[2] + expect_l23[3])
        np.testing.assert_allclose(prof.loc["GAD1", "SEAAD_Sst"], counts[4:7, 1].mean())

        print("test_obs_filter")
        assert run(*base, "--obs-filter", "disease=normal").returncode == 0
        prof = pd.read_csv(out, index_col=0)
        normal = np.array(DISEASE) == "normal"
        l23 = normal & (np.array(CELL_TYPES) == "L2/3 IT")
        np.testing.assert_allclose(prof.loc["GAD1", "SEAAD_L2/3 IT"], counts[l23, 1].mean())

        print("test_bad_obs_filter_value_is_refused")
        assert run(*base, "--obs-filter", "disease=nope").returncode != 0

        print("test_raw_slot_preferred_over_lognormalized_X")
        make_h5ad(h5ad, counts, raw_in_raw=True)
        assert run(*base).returncode == 0
        prof = pd.read_csv(out, index_col=0)
        np.testing.assert_allclose(prof.loc["SLC17A7", "SEAAD_L2/3 IT"], expect_l23[0])

        print("test_lognormalized_X_without_raw_is_refused")
        make_h5ad(h5ad, np.log1p(counts))
        res = run(*base)
        assert res.returncode != 0 and "non-integer" in (res.stdout + res.stderr)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
