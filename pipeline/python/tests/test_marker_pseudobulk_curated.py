#!/usr/bin/env python3
"""Tests for the curated-gene path of marker_pseudobulk.py.

The load-bearing claim is an optimization: with a curated list, only the listed columns are
log-normalized, divided by each cell's FULL gene total. That must produce numbers identical
to normalizing the whole panel and subsetting afterwards -- otherwise the heatmap silently
changes depending on how many genes were asked for.

    uv run python pipeline/python/tests/test_marker_pseudobulk_curated.py
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
SCRIPT = _PY_DIR / "marker_pseudobulk.py"
_spec = importlib.util.spec_from_file_location("pseudobulk_core", _PY_DIR / "pseudobulk_core.py")
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)

N_CELLS, N_GENES = 400, 60
CURATED = ["GENE003", "GENE017", "GENE041", "NOT_ON_PANEL"]
MODULES = ["synaptic", "synaptic", "channel", "channel"]


def build(tmp: Path) -> tuple[Path, Path]:
    rng = np.random.default_rng(7)
    counts = sp.csr_matrix(rng.poisson(1.2, size=(N_CELLS, N_GENES)).astype(np.float32))
    var = pd.DataFrame({"probe_type": ["gene"] * N_GENES},
                       index=[f"GENE{i:03d}" for i in range(N_GENES)])
    obs = pd.DataFrame({
        "leiden": pd.Categorical(rng.choice(["0", "1", "2"], N_CELLS)),
        "Region": pd.Categorical(rng.choice(["Retina", "Optic nerve"], N_CELLS)),
    }, index=[f"c{i:04d}" for i in range(N_CELLS)])
    h5ad = tmp / "clustered.h5ad"
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(h5ad)

    genes_csv = tmp / "curated.csv"
    with genes_csv.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["gene", "module"])
        w.writerows(zip(CURATED, MODULES))
    return h5ad, genes_csv


def run(h5ad: Path, out: Path, *extra: str):
    return subprocess.run([sys.executable, str(SCRIPT), "--clustered-h5ad", str(h5ad),
                           "--output-dir", str(out), "--group-key", "leiden",
                           "--region-key", "none", "--min-group-n", "1", *extra],
                          capture_output=True, text=True)


def test_subset_normalization_matches_whole_panel():
    """The optimization is only safe if the divisor stays the full gene total."""
    rng = np.random.default_rng(3)
    counts = sp.csr_matrix(rng.poisson(1.5, size=(50, 30)).astype(np.float64))
    cols = np.array([2, 11, 29])
    whole = core.log_normalize(counts, 1000.0)[:, cols].toarray()
    totals = np.asarray(counts.sum(axis=1)).ravel()
    subset = core.log_normalize(counts[:, cols], 1000.0, totals=totals).toarray()
    assert np.allclose(whole, subset), np.abs(whole - subset).max()
    # And without `totals` it is demonstrably WRONG -- the guard earns its place.
    naive = core.log_normalize(counts[:, cols], 1000.0).toarray()
    assert not np.allclose(whole, naive)
    print("  subset-with-totals == whole-panel; subset-without-totals differs")


def test_curated_rows_in_file_order_and_missing_reported():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        h5ad, genes_csv = build(tmp)
        out = tmp / "out"
        res = run(h5ad, out, "--genes-csv", str(genes_csv), "--gene-group-column", "module")
        assert res.returncode == 0, res.stdout + res.stderr

        z = pd.read_csv(out / "marker_heatmap_zmatrix.csv", index_col=0)
        assert list(z.index) == CURATED[:3], list(z.index)

        missing = pd.read_csv(out / "genes_missing.csv")
        assert missing["gene"].tolist() == ["NOT_ON_PANEL"], missing
        assert missing["group"].tolist() == ["channel"]

        rows = pd.read_csv(out / "top_markers_per_cluster.csv")
        assert rows["cluster"].tolist() == MODULES[:3], rows
        print(f"  {len(z)} curated rows in file order, 1 off-panel gene reported")


def test_curated_values_match_the_top_n_path():
    """Same genes, two code paths: the numbers must not depend on how they were chosen."""
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        h5ad, genes_csv = build(tmp)
        curated_out, topn_out = tmp / "curated", tmp / "topn"
        assert run(h5ad, curated_out, "--genes-csv", str(genes_csv),
                   "--gene-group-column", "module").returncode == 0
        assert run(h5ad, topn_out, "--top-n", str(N_GENES)).returncode == 0

        a = pd.read_csv(curated_out / "marker_heatmap_zmatrix.csv", index_col=0)
        b = pd.read_csv(topn_out / "marker_heatmap_zmatrix.csv", index_col=0).loc[a.index, a.columns]
        assert np.allclose(a.to_numpy(), b.to_numpy()), np.abs(a.to_numpy() - b.to_numpy()).max()
        print(f"  {a.shape[0]} genes x {a.shape[1]} groups identical across both paths")


def test_no_curated_gene_on_panel_is_an_error():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        h5ad, _ = build(tmp)
        bad = tmp / "bad.csv"
        bad.write_text("gene,module\nNOPE1,x\nNOPE2,x\n")
        res = run(h5ad, tmp / "out", "--genes-csv", str(bad), "--gene-group-column", "module")
        assert res.returncode != 0
        assert "no curated gene is on the panel" in res.stdout + res.stderr
        print("  an all-off-panel list fails loudly")


if __name__ == "__main__":
    for test in (test_subset_normalization_matches_whole_panel,
                 test_curated_rows_in_file_order_and_missing_reported,
                 test_curated_values_match_the_top_n_path,
                 test_no_curated_gene_on_panel_is_an_error):
        print(test.__name__)
        test()
    print("\nall curated marker_pseudobulk tests passed")
