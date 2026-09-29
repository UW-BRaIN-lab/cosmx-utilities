#!/usr/bin/env python3
"""Tests for extract_anchor_leiden.py (no S3 / network).

The one that matters is test_per_slide_labels_are_refused: it plants the exact mistake that
happened -- anchor_cells.csv's `slide_cluster` ("<slide>|<n>") mistaken for the cohort-wide Leiden
-- and requires the script to refuse it. The rest pin the overlap check and that a cohort-wide
clustering is reported as spanning slides.

Run:  uv run --with pandas --with numpy --with anndata python pipeline/python/tests/test_extract_anchor_leiden.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "extract_anchor_leiden.py"
SLIDES = ["SLA", "SLB", "SLC"]


def _make(tmp: Path, leiden_of, anchor_keep=None, key="leiden"):
    """A tiny typed AnnData over 3 slides x 12 cells, and an anchor CSV over a subset."""
    ids = [f"{s}_F{f}_C{c}" for s in SLIDES for f in (1, 2) for c in range(1, 7)]
    obs = pd.DataFrame({key: [leiden_of(i) for i in ids], "cell_type": "x"}, index=pd.Index(ids))
    ad.AnnData(X=np.zeros((len(ids), 1), dtype="float32"), obs=obs).write_h5ad(tmp / "typed.h5ad")
    keep = ids if anchor_keep is None else anchor_keep
    pd.DataFrame({"cell_id": keep, "slide_id": [k.split("_F")[0] for k in keep]}
                 ).to_csv(tmp / "anchor.csv", index=False)
    return ids


def _run(tmp: Path, *extra):
    return subprocess.run([sys.executable, str(SCRIPT), "--h5ad", str(tmp / "typed.h5ad"),
                           "--cell-ids-csv", str(tmp / "anchor.csv"),
                           "--output-csv", str(tmp / "out.csv"), *extra],
                          capture_output=True, text=True)


def test_a_cohort_wide_clustering_is_extracted_and_reported_as_spanning_slides():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        # Cluster id from the cell number, so every cluster appears on every slide.
        _make(tmp, lambda i: str(int(i.rsplit("_C", 1)[1]) % 3))
        out = _run(tmp)
        assert out.returncode == 0, out.stderr
        got = pd.read_csv(tmp / "out.csv", dtype=str)
        assert list(got.columns) == ["cell_id", "leiden"] and len(got) == 36
        assert set(got.leiden) == {"0", "1", "2"}
        assert "span >= half the slides" in out.stdout
        assert "3 of 3 clusters span >= half the slides" in out.stdout


def test_per_slide_labels_are_refused():
    """The actual mistake: slide_cluster ('<slide>|<n>') is not the cohort-wide Leiden."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _make(tmp, lambda i: f"{i.split('_F')[0]}|{int(i.rsplit('_C', 1)[1]) % 3}")
        out = _run(tmp)
        assert out.returncode != 0
        assert "per-slide clustering" in (out.stdout + out.stderr)
        assert not (tmp / "out.csv").exists()


def test_only_the_anchor_cells_are_kept():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        keep = ["SLA_F1_C1", "SLB_F2_C3", "SLC_F1_C6"]
        _make(tmp, lambda i: "7", anchor_keep=keep)
        assert _run(tmp).returncode == 0
        assert set(pd.read_csv(tmp / "out.csv").cell_id) == set(keep)


def test_a_missing_leiden_column_lists_what_is_there():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _make(tmp, lambda i: "1", key="louvain")
        out = _run(tmp)
        assert out.returncode != 0
        assert "no column 'leiden'" in (out.stdout + out.stderr)
        assert "louvain" in (out.stdout + out.stderr)


def test_ids_that_do_not_overlap_fail_with_examples_from_both_sides():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _make(tmp, lambda i: "1")
        pd.DataFrame({"cell_id": ["other_F1_C1", "other_F1_C2"]}).to_csv(tmp / "anchor.csv", index=False)
        out = _run(tmp)
        assert out.returncode != 0
        both = out.stdout + out.stderr
        assert "anchor example" in both and "h5ad   example" in both


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
