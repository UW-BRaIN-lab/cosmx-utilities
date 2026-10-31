#!/usr/bin/env python3
"""Tests for the negprobe background test (no S3 / network).

The two that matter are test_a_planted_background_gene_is_detected and
test_depth_coupling_alone_is_NOT_called_background. Depth drives log-norm expression (through
dropout) and negprobe counts alike, so an unstratified correlation is positive for EVERY gene and
would call the whole panel background. Stratifying by depth is what makes the measurement mean
anything, and the second test is the regression that proves it.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy --with scipy --with h5py \\
        python pipeline/python/tests/test_negprobe_background_test.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location(
    "negprobe_background_test", _PY_DIR / "negprobe_background_test.py")
nb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nb)

import matched_null as mn  # noqa: E402

_N_BG = 8


def _simulate(rng, n_genes=300, n_cells=6000, bg_strength=0.0, neg_tracks_depth=False):
    """Cells with a per-cell background load, optionally driving some genes and/or tracking depth.

    `bg_strength` plants genes whose counts come from the background load rather than the cell's
    depth. `neg_tracks_depth` makes the negprobe reading rise with library size, which is what
    real data does and what the stratification has to survive.
    """
    gene_rate = rng.lognormal(-1.0, 1.2, n_genes)
    gene_rate /= gene_rate.sum()
    depth = rng.lognormal(np.log(400.0), 0.6, n_cells)
    bg_load = rng.gamma(2.0, 0.5, n_cells)
    neg = bg_load * (depth / depth.mean()) if neg_tracks_depth else bg_load

    rate = np.outer(gene_rate, depth)
    bg_idx = np.arange(_N_BG)
    if bg_strength:
        rate[bg_idx, :] += bg_strength * bg_load
    mat = sp.csc_matrix(rng.poisson(rate).astype(np.float64))
    return mat, np.asarray(mat.sum(axis=0)).ravel(), neg, bg_idx


def _affinity(mat, totals, neg, n_strata):
    cells = np.arange(mat.shape[1])
    return nb.background_affinity(mat, totals, neg, cells, n_strata, 1e4, 50)[0]


def _null_p(mat, affinity, idx, n_sets=200, seed=3):
    """The main() path: matched sets, then a one-sided p for affinity ABOVE the null."""
    mean_counts, detection = mn.gene_statistics(mat)
    pools = mn.matched_pools(mean_counts, detection, np.asarray(idx), 60)
    draws = mn.draw_null_sets(pools, n_sets, np.random.default_rng(seed))
    observed = float(affinity[np.asarray(idx)].mean())
    null = np.array([affinity[np.asarray(d)].mean() for d in draws])
    return observed, null, mn.empirical_p(-observed, -null)


def test_a_planted_background_gene_is_detected():
    """Genes driven by the background load must stand out from their matched null."""
    mat, totals, neg, bg_idx = _simulate(np.random.default_rng(1), bg_strength=3.0)
    affinity = _affinity(mat, totals, neg, 10)
    observed, null, p = _null_p(mat, affinity, bg_idx)
    print(f"      planted: affinity {observed:+.4f}   null median {np.median(null):+.4f}   "
          f"p(above) {p:.3f}")
    assert observed > np.percentile(null, 95), (observed, np.percentile(null, 95))
    assert p <= 0.05, f"a planted background set must be flagged, got p={p:.3f}"


def test_depth_coupling_alone_is_NOT_called_background():
    """The regression this stratification exists for.

    Ground truth: NO gene is background-driven. But negprobe counts rise with depth and so does
    log-norm expression (shallow cells drop genes to zero), so an unstratified correlation is
    positive for the whole panel.
    """
    mat, totals, neg, _ = _simulate(np.random.default_rng(2), bg_strength=0.0,
                                    neg_tracks_depth=True)
    unstratified = np.median(_affinity(mat, totals, neg, 1))
    stratified = np.median(_affinity(mat, totals, neg, 10))
    print(f"      no background planted: unstratified median {unstratified:+.4f}   "
          f"stratified {stratified:+.4f}")
    assert unstratified > 0.05, (
        f"unstratified should (wrongly) look like panel-wide background, got {unstratified:+.4f}")
    assert abs(stratified) < 0.02, (
        f"stratifying by depth must collapse it, got {stratified:+.4f}")


def test_an_innocent_set_is_not_flagged_when_negprobes_track_depth():
    """The same fixture, now asked the question the job actually asks."""
    mat, totals, neg, _ = _simulate(np.random.default_rng(4), bg_strength=0.0,
                                    neg_tracks_depth=True)
    affinity = _affinity(mat, totals, neg, 10)
    innocent = np.arange(40, 40 + _N_BG)
    _, _, p = _null_p(mat, affinity, innocent)
    assert p > 0.05, f"a set with no background involvement must not be flagged, p={p:.3f}"


def test_background_ratio_ranks_a_dirty_population_first():
    """The ratio must compare like with like: a DEEP type can carry more neg and be cleaner."""
    totals = np.concatenate([np.full(500, 200.0), np.full(500, 2000.0)])
    neg = np.concatenate([np.full(500, 0.40), np.full(500, 0.80)])   # deep type: higher raw neg
    cell_type = np.array(["dirty_shallow"] * 500 + ["clean_deep"] * 500)
    table = nb.background_by_celltype(totals, neg, cell_type, n_genes=1000, min_cells=100)
    assert list(table["cell_type"]) == ["dirty_shallow", "clean_deep"], table
    assert table.iloc[0]["median_neg"] < table.iloc[1]["median_neg"], (
        "the fixture must be one where raw neg would give the WRONG order")


def test_background_by_celltype_drops_small_types_and_zero_depth_cells():
    totals = np.array([100.0, 100.0, 0.0, 100.0])
    neg = np.array([1.0, 1.0, 1.0, 1.0])
    cell_type = np.array(["big", "big", "big", "tiny"])
    table = nb.background_by_celltype(totals, neg, cell_type, n_genes=100, min_cells=2)
    assert list(table["cell_type"]) == ["big"]
    assert int(table.iloc[0]["n_cells"]) == 2, "the zero-depth cell must be excluded"


def test_resolve_present_drops_absent_genes_with_a_warning():
    genes = np.array(["AAA", "BBB", "CCC"])
    assert nb.resolve_present(genes, "AAA,ZZZ,CCC", "signature") == [0, 2]


def test_resolve_present_exits_when_nothing_is_on_the_panel():
    genes = np.array(["AAA"])
    try:
        nb.resolve_present(genes, "YYY,ZZZ", "signature")
    except SystemExit:
        return
    raise AssertionError("an entirely off-panel set must exit, not silently score nothing")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
