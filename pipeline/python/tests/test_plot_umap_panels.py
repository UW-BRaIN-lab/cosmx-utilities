#!/usr/bin/env python3
"""Tests for plot_umap_panels.py (no S3 / network).

Pins: one panel per group in natural (numeric) order, donor ids stay in numeric order, the
share table rows sum to 1 and reflect planted composition, the run works on a backed h5ad
without loading X, and a missing key fails loudly.

Run:  uv run --with pandas --with numpy --with anndata --with scipy --with matplotlib python pipeline/python/tests/test_plot_umap_panels.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "plot_umap_panels.py"
sys.path.insert(0, str(SCRIPT.parent))
from plot_umap_panels import natural_order, share_table  # noqa: E402

N_CELLS = 600


def _make(tmp: Path) -> Path:
    """3 slides x 2 clusters; slide 'S3' is planted 90% in cluster 1, the others 50/50."""
    rng = np.random.default_rng(0)
    slides = np.repeat(["S1", "S2", "S3"], N_CELLS // 3)
    cluster = np.where(rng.random(N_CELLS) < 0.5, "0", "1")
    s3 = slides == "S3"
    cluster[s3] = np.where(rng.random(s3.sum()) < 0.9, "1", "0")
    umap = rng.normal(size=(N_CELLS, 2))
    umap[0] = [5000, 5000]                       # a flung-out outlier the view must ignore
    obs = pd.DataFrame({"slide_id": slides, "leiden": cluster,
                        "Region": rng.choice(["741", "6803", "526"], N_CELLS)},
                       index=[f"c{i}" for i in range(N_CELLS)])
    adata = ad.AnnData(X=sp.csr_matrix((N_CELLS, 3)), obs=obs)
    adata.obsm["X_umap"] = umap
    path = tmp / "clustered.h5ad"
    adata.write_h5ad(path)
    return path


def _run(h5ad: Path, out: Path, *extra):
    return subprocess.run([sys.executable, str(SCRIPT), "--h5ad", str(h5ad), "--out-dir", str(out),
                           *extra], capture_output=True, text=True)


def test_natural_order_keeps_donor_ids_numeric():
    print(natural_order(["741", "6803", "526", "7316"]))
    assert natural_order(["741", "6803", "526", "7316"]) == ["526", "741", "6803", "7316"]
    assert natural_order(["S10", "S2"]) == ["S10", "S2"]      # not all ints -> lexical


def test_share_table_rows_sum_to_one_and_show_planted_composition():
    groups = pd.Series(["a"] * 10 + ["b"] * 10)
    clusters = pd.Series(["0"] * 5 + ["1"] * 5 + ["1"] * 9 + ["0"])
    shares = share_table(groups, clusters)
    print(shares)
    assert np.allclose(shares.sum(axis=1), 1.0)
    assert shares.loc["b", "1"] == 0.9


def test_end_to_end_on_backed_h5ad():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        result = _run(_make(tmp), tmp / "out")
        print(result.stdout[-500:], result.stderr[-300:])
        assert result.returncode == 0
        assert (tmp / "out" / "umap_panels_slide_id.png").stat().st_size > 5000
        shares = pd.read_csv(tmp / "out" / "slide_id_by_leiden_share.csv", index_col=0)
        print(shares)
        assert shares.index.tolist() == ["S1", "S2", "S3"]
        assert np.allclose(shares.sum(axis=1), 1.0)
        assert shares.loc["S3", "1"] > 0.8 > shares.loc["S1", "1"]
        assert "clipping the shared view" in result.stdout           # the outlier was reported


def test_donor_key_makes_numeric_panels_and_missing_key_fails():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        h5ad = _make(tmp)
        ok = _run(h5ad, tmp / "out", "--key", "Region", "--ncols", "3")
        assert ok.returncode == 0 and (tmp / "out" / "umap_panels_Region.png").exists()
        bad = _run(h5ad, tmp / "out2", "--key", "no_such_column")
        print("stderr:", bad.stderr.strip()[-160:])
        assert bad.returncode != 0 and "missing 'no_such_column'" in bad.stderr


if __name__ == "__main__":
    test_natural_order_keeps_donor_ids_numeric()
    test_share_table_rows_sum_to_one_and_show_planted_composition()
    test_end_to_end_on_backed_h5ad()
    test_donor_key_makes_numeric_panels_and_missing_key_fails()
    print("OK")
