#!/usr/bin/env python3
"""marker_pseudobulk.py with a non-tissue Region (donor ids) and --region-key none.

SORL1's Region is the donor id, which matches none of the GBM REGION_ORDER names. The stage
used to select 0 cells, write empty CSVs, report success and upload them. Pins:
  - the default region filter now FAILS LOUDLY on 0 cells instead of writing empty output;
  - --region-key none writes one z-matrix column per cluster and a marker per cluster.

Run:  uv run --with pandas --with numpy --with anndata --with scipy python pipeline/python/tests/test_marker_pseudobulk_no_region.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "marker_pseudobulk.py"
N_PER_CLUSTER = 40
GENES = [f"G{i}" for i in range(6)]
NEGS = ["NegPrb1"]


def _make(tmp: Path) -> Path:
    """3 clusters x 40 cells; cluster k over-expresses gene G(2k) and G(2k+1). Region = donor id."""
    rng = np.random.default_rng(0)
    rows, leiden, donor = [], [], []
    for k in range(3):
        for i in range(N_PER_CLUSTER):
            counts = rng.poisson(1.0, len(GENES)).astype(float)
            counts[2 * k:2 * k + 2] += 20
            rows.append(np.append(counts, 0.0))
            leiden.append(str(k))
            donor.append(str(741 + (i % 2) * 5762))   # numeric-looking donor ids
    obs = pd.DataFrame({"leiden": leiden, "Region": donor},
                       index=[f"c{i}" for i in range(len(rows))])
    var = pd.DataFrame({"probe_type": ["gene"] * len(GENES) + ["negprobe"]},
                       index=GENES + NEGS)
    path = tmp / "clustered.h5ad"
    ad.AnnData(X=sp.csr_matrix(np.vstack(rows)), obs=obs, var=var).write_h5ad(path)
    return path


def _run(h5ad: Path, out: Path, *extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--clustered-h5ad", str(h5ad), "--output-dir", str(out),
         "--top-n", "2", "--min-group-n", "5", *extra], capture_output=True, text=True)


def test_default_region_filter_fails_loudly_on_zero_cells():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        result = _run(_make(tmp), tmp / "out")
        print("returncode:", result.returncode, "| stderr:", result.stderr.strip()[-200:])
        assert result.returncode != 0
        assert "no cells selected" in result.stderr
        assert not (tmp / "out" / "marker_heatmap_zmatrix.csv").exists()


def test_region_key_none_gives_one_column_per_cluster():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        result = _run(_make(tmp), tmp / "out", "--region-key", "none")
        print(result.stdout[-400:], result.stderr[-300:])
        assert result.returncode == 0
        z = pd.read_csv(tmp / "out" / "marker_heatmap_zmatrix.csv", index_col=0)
        markers = pd.read_csv(tmp / "out" / "top_markers_per_cluster.csv")
        print("columns:", z.columns.tolist(), "| markers:", markers.values.tolist())
        assert z.columns.tolist() == ["0", "1", "2"]      # no " | " => renderer's profile mode
        assert {g for g in markers["gene"]} == set(GENES)
        # each cluster's own markers must be its highest-z column
        for cluster, gene in markers.values:
            assert z.loc[gene].idxmax() == str(cluster)


if __name__ == "__main__":
    test_default_region_filter_fails_loudly_on_zero_cells()
    test_region_key_none_gives_one_column_per_cluster()
    print("OK")
