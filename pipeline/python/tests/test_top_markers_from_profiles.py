#!/usr/bin/env python3
"""Tests for top_markers_from_profiles.py (no S3 / network).

Plants: a clear marker per type, and a near-zero gene with a huge fold-change that the expression
floor must keep OFF the list (the failure of fold-enrichment with no floor).

Run:  uv run --with pandas --with numpy --with h5py python pipeline/python/tests/test_top_markers_from_profiles.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "top_markers_from_profiles.py"
sys.path.insert(0, str(SCRIPT.parent))
from top_markers_from_profiles import read_profiles, top_markers  # noqa: E402

TYPES = ["A", "B", "C"]
GENES = ["gA", "gB", "gC", "shared", "tiny", "other1"]


def _profiles() -> pd.DataFrame:
    p = pd.DataFrame(0.0, index=GENES, columns=TYPES)
    p.loc["gA", "A"] = 5.0
    p.loc["gB", "B"] = 4.0
    p.loc["gC", "C"] = 3.0
    p.loc["shared"] = 2.0                       # expressed everywhere: not a marker of anything
    p.loc["tiny", "A"] = 0.01                   # 100x any other type, but essentially absent
    p.loc["other1"] = 0.3
    return p


def test_planted_marker_is_rank_one_and_shared_gene_is_not():
    table = top_markers(_profiles(), top_n=3, min_profile=0.2)
    print(table.to_string(index=False))
    for t, g in zip(TYPES, ["gA", "gB", "gC"]):
        assert table[(table.type == t) & (table["rank"] == 1)].gene.iloc[0] == g
    top_a = table[table.type == "A"]
    assert top_a[top_a.gene == "shared"].log2fc.iloc[0] < 0.1     # ~no enrichment


def test_expression_floor_keeps_a_near_zero_gene_off_the_list():
    with_floor = top_markers(_profiles(), top_n=6, min_profile=0.2)
    no_floor = top_markers(_profiles(), top_n=6, min_profile=0.0)
    assert "tiny" not in set(with_floor[with_floor.type == "A"].gene)
    assert "tiny" in set(no_floor[no_floor.type == "A"].gene)


def test_cli_roundtrip_with_names_and_either_orientation():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        p = _profiles()
        with h5py.File(tmp / "r.h5", "w") as f:           # stored types x genes, like InSituType's file
            f["profiles"] = p.T.to_numpy()
            f["profile_genes"] = np.array(GENES, dtype="S")
            f["profile_types"] = np.array(TYPES, dtype="S")
        assert read_profiles(tmp / "r.h5").equals(p)
        pd.DataFrame({"type": ["A"], "name": ["Type A"]}).to_csv(tmp / "n.csv", index=False)
        out = subprocess.run([sys.executable, str(SCRIPT), "--profiles-h5", str(tmp / "r.h5"),
                              "--output", str(tmp / "o.csv"), "--names", str(tmp / "n.csv"),
                              "--top-n", "2"], capture_output=True, text=True)
        print(out.stdout, out.stderr[-200:])
        assert out.returncode == 0
        table = pd.read_csv(tmp / "o.csv", keep_default_na=False)
        assert table[table.type == "A"].name.iloc[0] == "Type A" and table[table.type == "B"].name.iloc[0] == ""


if __name__ == "__main__":
    test_planted_marker_is_rank_one_and_shared_gene_is_not()
    test_expression_floor_keeps_a_near_zero_gene_off_the_list()
    test_cli_roundtrip_with_names_and_either_orientation()
    print("OK")
