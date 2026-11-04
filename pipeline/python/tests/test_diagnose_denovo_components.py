#!/usr/bin/env python3
"""Tests for the de novo component decomposition.

The table exists to tell a real population apart from a low-signal pool, so the fixture
builds one of each and asserts the output separates them:

  cluster 'o'  spans Leiden 1 (deep, cohort-typical) and Leiden 2 (shallow, half the counts)
  cluster 'x'  is wholly inside Leiden 3 and is deep

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_diagnose_denovo_components.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_SCRIPT = _PY_DIR / "diagnose_denovo_components.py"
_spec = importlib.util.spec_from_file_location("diagnose_denovo_components", _SCRIPT)
ddc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ddc)

DEEP, SHALLOW = 400.0, 180.0
N_PER_GROUP = 50


# Above the script's 50-gene shared-panel floor, so the fixture exercises the real path.
N_GENES = 60
# Component '1' expresses the first third of the panel, '2' the second, '3' the last --
# so each one's nearest reference column is known in advance and the cosine is checkable.
BLOCKS = {"1": slice(0, 20), "2": slice(20, 40), "3": slice(40, 60)}
REF_TYPES = {"Fibro_like": "1", "Schwann_like": "2", "Astro_like": "3"}


def build_fixture(tmp: Path) -> tuple[Path, Path]:
    """Three Leiden components; 'o' straddles a deep one and a shallow one."""
    groups = [
        # (leiden, cell_type, total_counts, genes_detected, negprobe, region)
        ("1", "o", DEEP, 220, 0.01, "Retina"),
        ("2", "o", SHALLOW, 90, 0.06, "White matter"),
        ("3", "x", DEEP, 230, 0.01, "Gray matter"),
        ("3", "Astrocyte", DEEP, 225, 0.01, "Gray matter"),
    ]
    rows, ids, labels, expression = [], [], [], []
    for i, (leiden, ctype, tc, genes, neg, region) in enumerate(groups):
        for j in range(N_PER_GROUP):
            cid = f"slide_F1_C{i}_{j}"
            ids.append(cid)
            labels.append(ctype)
            rows.append({"leiden": leiden, "total_counts": tc,
                         "qc_genes_detected": genes, "qc_negprobe_prop": neg,
                         "qc_area": 9000.0, "Region": region})
            counts = np.zeros(N_GENES, dtype=np.float32)
            counts[BLOCKS[leiden]] = 100.0
            expression.append(counts)
    obs = pd.DataFrame(rows, index=pd.Index(ids, name="cell"))
    var = pd.DataFrame({"probe_type": ["gene"] * N_GENES},
                       index=[f"GENE{g}" for g in range(N_GENES)])
    adata = ad.AnnData(X=np.vstack(expression), obs=obs, var=var)
    h5ad = tmp / "clustered.h5ad"
    adata.write_h5ad(h5ad)

    reference = pd.DataFrame(0.0, index=var.index, columns=list(REF_TYPES))
    for ref_name, block in REF_TYPES.items():
        reference.loc[reference.index[BLOCKS[block]], ref_name] = 500.0
    reference.to_csv(tmp / "reference.csv")

    typing = tmp / "typing.h5"
    with h5py.File(typing, "w") as f:
        dt = h5py.string_dtype(encoding="utf-8")
        f.create_dataset("cell_id", data=np.array(ids, dtype=object), dtype=dt)
        f.create_dataset("cell_type", data=np.array(labels, dtype=object), dtype=dt)
        f.create_dataset("prob", data=np.full(len(ids), 0.9))
    return h5ad, typing


def test_separates_deep_from_shallow(tmp: Path) -> None:
    h5ad, typing = build_fixture(tmp)
    out = tmp / "out"
    subprocess.run(
        [sys.executable, str(_SCRIPT), "--clustered-h5ad", str(h5ad),
         "--typing-h5", str(typing), "--cluster", "o", "--cluster", "x",
         "--out-dir", str(out)],
        check=True, capture_output=True, text=True,
    )

    o = pd.read_csv(out / "o_components.csv").set_index("component")
    assert set(o.index) == {ddc.BASELINE_LABEL, "1", "2"}, o.index.tolist()

    # Composition: 'o' is split evenly, and each component is wholly 'o'.
    assert o.loc["1", "cells"] == N_PER_GROUP
    assert o.loc["1", "pct_of_cluster"] == 50.0
    assert o.loc["1", "pct_of_component"] == 100.0

    # Depth: the whole point. Leiden 2 must read as shallow against the baseline.
    assert o.loc["1", "median_total_counts"] == DEEP
    assert o.loc["2", "median_total_counts"] == SHALLOW
    assert o.loc["2", "median_total_counts"] < o.loc[ddc.BASELINE_LABEL,
                                                     "median_total_counts"]
    assert o.loc["2", "median_qc_negprobe_prop"] > o.loc["1", "median_qc_negprobe_prop"]
    assert o.loc["1", "dominant_region"] == "Retina"
    assert o.loc["1", "dominant_region_pct"] == 100.0

    # The reverse direction: 'x' fills only half of Leiden 3, which it shares with a
    # named type -- the column that says "this component is not mostly this letter".
    x = pd.read_csv(out / "x_components.csv").set_index("component")
    assert x.loc["3", "pct_of_cluster"] == 100.0
    assert x.loc["3", "pct_of_component"] == 50.0
    print("  deep/shallow separation, composition and both directions OK")


def test_tail_rollup_and_missing_cluster(tmp: Path) -> None:
    h5ad, typing = build_fixture(tmp)
    out = tmp / "out_tail"
    # 60% floor leaves neither 'o' component standing, so both roll into __other__.
    subprocess.run(
        [sys.executable, str(_SCRIPT), "--clustered-h5ad", str(h5ad),
         "--typing-h5", str(typing), "--cluster", "o",
         "--min-component-pct", "60", "--out-dir", str(out)],
        check=True, capture_output=True, text=True,
    )
    o = pd.read_csv(out / "o_components.csv").set_index("component")
    assert set(o.index) == {ddc.BASELINE_LABEL, ddc.TAIL_LABEL}, o.index.tolist()
    assert o.loc[ddc.TAIL_LABEL, "cells"] == 2 * N_PER_GROUP

    # A label that does not exist warns and skips rather than failing the run.
    r = subprocess.run(
        [sys.executable, str(_SCRIPT), "--clustered-h5ad", str(h5ad),
         "--typing-h5", str(typing), "--cluster", "zzz", "--out-dir", str(out)],
        check=True, capture_output=True, text=True,
    )
    assert "no cells labelled 'zzz'" in r.stderr, r.stderr
    assert not (out / "zzz_components.csv").exists()
    print("  tail rollup and absent-label handling OK")


def test_reference_similarity_names_the_right_branch(tmp: Path) -> None:
    """The column that decides a new leaf: each component's nearest reference type."""
    h5ad, typing = build_fixture(tmp)
    out = tmp / "out_ref"
    subprocess.run(
        [sys.executable, str(_SCRIPT), "--clustered-h5ad", str(h5ad),
         "--typing-h5", str(typing), "--cluster", "o",
         "--reference-csv", str(tmp / "reference.csv"), "--out-dir", str(out)],
        check=True, capture_output=True, text=True,
    )
    o = pd.read_csv(out / "o_components.csv").set_index("component")

    # 'o' straddles two components whose nearest reference types differ -- the signature of
    # an over-merged cluster, and the reason it should dissolve rather than gain one leaf.
    assert o.loc["1", "nearest"] == "Fibro_like", o.loc["1", "nearest"]
    assert o.loc["2", "nearest"] == "Schwann_like", o.loc["2", "nearest"]
    assert o.loc["1", "cosine"] > 0.99
    # The runner-up must be a genuinely worse match, or "nearest" means nothing.
    assert o.loc["1", "cosine_2"] < 0.5
    # The cohort baseline is not scored against the reference; it is not a component.
    assert pd.isna(o.loc[ddc.BASELINE_LABEL, "nearest"])
    print("  reference similarity picks the right branch per component")


def test_missing_column_errors(tmp: Path) -> None:
    h5ad, typing = build_fixture(tmp)
    r = subprocess.run(
        [sys.executable, str(_SCRIPT), "--clustered-h5ad", str(h5ad),
         "--typing-h5", str(typing), "--cluster", "o",
         "--partition-col", "not_a_column", "--out-dir", str(tmp / "out_err")],
        capture_output=True, text=True,
    )
    assert r.returncode != 0
    assert "missing column(s)" in r.stdout + r.stderr
    print("  missing-column guard OK")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("test_diagnose_denovo_components")
        test_separates_deep_from_shallow(tmp)
        test_tail_rollup_and_missing_cluster(tmp)
        test_reference_similarity_names_the_right_branch(tmp)
        test_missing_column_errors(tmp)
    print("All tests passed.")


if __name__ == "__main__":
    main()
