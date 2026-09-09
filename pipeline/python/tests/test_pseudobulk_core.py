#!/usr/bin/env python3
"""Tests for the shared pseudobulk/marker machinery (no S3 / network).

Runnable either under pytest or directly:
    uv run --with pandas --with numpy --with scipy \\
        python pipeline/python/tests/test_pseudobulk_core.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location("pseudobulk_core", _PY_DIR / "pseudobulk_core.py")
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)


def test_log_normalize_puts_equal_compositions_on_equal_footing():
    """Two cells with the same composition but different depth must normalize identically."""
    counts = sp.csr_matrix(np.array([[1.0, 3.0], [10.0, 30.0]]))
    got = pc.log_normalize(counts, scale_factor=100.0).toarray()
    assert np.allclose(got[0], got[1])
    # log1p(100 * 1/4) on the first gene
    assert np.isclose(got[0, 0], np.log1p(25.0))


def test_log_normalize_leaves_empty_cells_at_zero():
    """A cell with no counts must not divide by zero."""
    got = pc.log_normalize(sp.csr_matrix(np.array([[0.0, 0.0], [2.0, 2.0]]))).toarray()
    assert np.all(got[0] == 0.0)


def test_onehot_and_group_means_average_within_group():
    labels = np.array(["a", "b", "a"])
    oh, cats = pc.onehot(labels)
    assert list(cats) == ["a", "b"]
    norm = sp.csr_matrix(np.array([[1.0, 0.0], [10.0, 10.0], [3.0, 4.0]]))
    means = pc.group_means(norm, oh)
    assert np.allclose(means[0], [2.0, 2.0])   # rows 0 and 2
    assert np.allclose(means[1], [10.0, 10.0])


def test_select_markers_picks_per_group_and_dedupes():
    profile = pd.DataFrame(
        {"A": [10.0, 1.0, 5.0], "B": [1.0, 10.0, 5.0]},
        index=["upA", "upB", "shared"])
    markers, gene_to_group = pc.select_markers(profile, ["A", "B"], top_n=2)
    assert markers[0] == "upA" and gene_to_group["upA"] == "A"
    assert "upB" in markers and gene_to_group["upB"] == "B"
    assert len(markers) == len(set(markers)), "a gene must not appear twice"
    # `shared` is equal in both, so whichever group claims it, it is claimed only once
    assert gene_to_group["shared"] in ("A", "B") or "shared" not in markers


def test_select_markers_differential_uses_all_columns_not_just_selected():
    """A group's markers are scored against every column, so a narrowed selection is not
    flattered by ignoring the groups it should be distinguished from."""
    profile = pd.DataFrame({"A": [10.0, 10.0], "B": [1.0, 1.0], "C": [10.0, 1.0]},
                           index=["shared_with_C", "unique_to_A"])
    markers, gene_to_group = pc.select_markers(profile, ["A"], top_n=1)
    assert markers == ["unique_to_A"], "C should suppress the gene A shares with it"


def test_zscore_rows_centres_and_scales():
    pb = pd.DataFrame({"x": [1.0, 5.0], "y": [3.0, 5.0], "z": [5.0, 5.0]}, index=["v", "flat"])
    z = pc.zscore_rows(pb)
    assert np.isclose(z.loc["v"].mean(), 0.0)
    assert np.isclose(z.loc["v", "y"], 0.0)
    assert z.loc["v", "x"] < 0 < z.loc["v", "z"]


def test_zscore_rows_keeps_constant_rows_flat_not_nan():
    """A gene identical across groups would otherwise divide by a zero sd."""
    z = pc.zscore_rows(pd.DataFrame({"x": [5.0], "y": [5.0]}, index=["flat"]))
    assert not z.isna().any().any()
    assert (z.loc["flat"] == 0.0).all()


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
