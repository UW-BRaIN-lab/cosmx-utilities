#!/usr/bin/env python3
"""Tests for donor_depth_table.py (no S3 / network).

Plants a known structure -- carrier donor D1 is 90% in a shallow cluster, non-carrier D2 is 10%
-- and checks the shallow-cluster selection, the per-donor fractions, the joins to slide/run and
the carrier flag, and that an unannotated donor fails loudly.

Run:  uv run --with pandas --with numpy --with anndata --with scipy python pipeline/python/tests/test_donor_depth_table.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "donor_depth_table.py"
sys.path.insert(0, str(SCRIPT.parent))
from donor_depth_table import donor_table, shallow_clusters  # noqa: E402

N_PER_DONOR = 1000


def _obs() -> pd.DataFrame:
    """Clusters: '0' shallow (counts ~200), '1' deep (~2000). D1 90% shallow, D2 10%."""
    rng = np.random.default_rng(0)
    rows = []
    for donor, p_shallow in (("D1", 0.9), ("D2", 0.1)):
        shallow = rng.random(N_PER_DONOR) < p_shallow
        counts = np.where(shallow, rng.normal(200, 20, N_PER_DONOR), rng.normal(2000, 100, N_PER_DONOR))
        rows.append(pd.DataFrame({"Region": donor, "leiden": np.where(shallow, "0", "1"),
                                  "total_counts": counts, "qc_genes_detected": counts / 2,
                                  "qc_area": 10000.0}))
    return pd.concat(rows, ignore_index=True)


ANNOTATIONS = pd.DataFrame({
    "donor": ["D1", "D2"], "slide_id": ["SA", "SB"], "case_broad": ["AD+LATE", "Control"],
    "case_group": ["AD+LATE SORL1", "Control"], "sorl1_mutation": ["AD+LATE SORL1 R953C", "Control"]})
MANIFEST = pd.DataFrame({"slide_id": ["SA", "SB"], "run_date": ["20260710", "20260917"],
                         "instrument_id": ["i1", "i2"], "run_uuid": ["u1", "u2"]})


def test_shallow_cluster_selection_uses_the_median_cutoff():
    shallow = shallow_clusters(_obs(), 560.0)
    print(shallow.round(0).to_dict())
    assert list(shallow.index) == ["0"]


def test_donor_table_fractions_joins_and_carrier_flag():
    table = donor_table(_obs(), ANNOTATIONS, MANIFEST, ["0"]).set_index("donor")
    print(table[["shallow_frac", "median_total_counts", "slide_id", "run_date", "carrier"]])
    assert abs(table.loc["D1", "shallow_frac"] - 0.9) < 0.03
    assert abs(table.loc["D2", "shallow_frac"] - 0.1) < 0.03
    assert table.loc["D1", "carrier"] and not table.loc["D2", "carrier"]
    assert table.loc["D1", "run_date"] == "20260710" and table.loc["D2", "slide_id"] == "SB"
    assert table.loc["D1", "median_total_counts"] < table.loc["D2", "median_total_counts"]


def test_unannotated_donor_raises():
    try:
        donor_table(_obs(), ANNOTATIONS[ANNOTATIONS.donor == "D1"], MANIFEST, ["0"])
    except ValueError as exc:
        print("raised:", exc)
        assert "D2" in str(exc)
        return
    raise AssertionError("an unannotated donor must fail loudly")


def test_cli_end_to_end():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        obs = _obs()
        obs.index = [f"c{i}" for i in range(len(obs))]
        ad.AnnData(X=sp.csr_matrix((len(obs), 2)), obs=obs).write_h5ad(tmp / "c.h5ad")
        ANNOTATIONS.to_csv(tmp / "a.csv", index=False)
        MANIFEST.to_csv(tmp / "m.csv", index=False)
        result = subprocess.run([sys.executable, str(SCRIPT), "--h5ad", str(tmp / "c.h5ad"),
                                 "--annotations", str(tmp / "a.csv"), "--manifest", str(tmp / "m.csv"),
                                 "--output", str(tmp / "out.csv")], capture_output=True, text=True)
        print(result.stdout[-900:], result.stderr[-300:])
        assert result.returncode == 0
        out = pd.read_csv(tmp / "out.csv")
        assert out["donor"].tolist()[0] == "D1"            # sorted shallowest-first
        assert "Shallow clusters" in result.stdout


if __name__ == "__main__":
    test_shallow_cluster_selection_uses_the_median_cutoff()
    test_donor_table_fractions_joins_and_carrier_flag()
    test_unannotated_donor_raises()
    test_cli_end_to_end()
    print("OK")
