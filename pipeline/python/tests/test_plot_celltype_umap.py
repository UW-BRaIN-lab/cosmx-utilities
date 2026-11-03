#!/usr/bin/env python3
"""Tests for the cell-type UMAP's axis clipping (no S3 / network).

The stage-4 cell-type UMAP shares an embedding with the stage-3c QC plots, so it hits the
same failure: UMAP flings a few near-disconnected components hundreds of units out, the
autoscaled axes span ~60x the real manifold, and every cell lands in a handful of pixels.
plot_qc.py already solved this; the bug was that plot_celltype_umap.py never used the fix.
These tests pin the wiring and the two regimes.

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_plot_celltype_umap.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))  # plot_celltype_umap imports plot_qc as a sibling
_spec = importlib.util.spec_from_file_location(
    "plot_celltype_umap", _PY_DIR / "plot_celltype_umap.py")
pcu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pcu)


def _embedding_with_outliers(n_bulk=20_000, n_outlier=60, spread=8.0, fling=800.0):
    """A tight manifold plus a few flung points — the retina run's actual shape."""
    rng = np.random.RandomState(0)
    return np.vstack([rng.normal(0, spread, size=(n_bulk, 2)),
                      rng.uniform(-fling, fling, size=(n_outlier, 2))])


def test_umap_view_limits_is_wired_in():
    """The sibling import is the fragile part: plot_qc must resolve from plot_celltype_umap."""
    assert callable(pcu.umap_view_limits)


def test_flung_outliers_trigger_clipping():
    xy = _embedding_with_outliers()
    xlim, ylim, n_outside, ratio = pcu.umap_view_limits(xy)
    assert xlim is not None and ylim is not None, "should clip when outliers dominate the view"
    assert ratio > 10, f"expected a large overshoot, got {ratio:.1f}x"
    # The drawn view must frame the bulk, not the outliers.
    assert (xlim[1] - xlim[0]) < 200, f"clipped view still too wide: {xlim}"
    assert 0 < n_outside < len(xy) * 0.02, f"unexpected outlier count {n_outside}"


def test_outliers_are_reported_not_dropped():
    """Clipping changes the view only — every cell stays in the counts and the legend."""
    xy = _embedding_with_outliers()
    (x0, x1), (y0, y1), n_outside, _ = pcu.umap_view_limits(xy)
    actually_outside = int((
        (xy[:, 0] < x0) | (xy[:, 0] > x1) | (xy[:, 1] < y0) | (xy[:, 1] > y1)
    ).sum())
    assert n_outside == actually_outside


def test_compact_embedding_is_left_alone():
    """A manifold with no flung components keeps its autoscaled view untouched."""
    rng = np.random.RandomState(1)
    xy = rng.uniform(-20, 20, size=(20_000, 2))
    xlim, ylim, n_outside, ratio = pcu.umap_view_limits(xy)
    assert xlim is None and ylim is None, f"clipped a compact embedding (ratio {ratio:.2f})"
    assert n_outside == 0


def test_diffuse_tails_alone_can_cross_the_threshold():
    """Documents a sharp edge in AXIS_CLIP_MIN_RATIO=1.5 rather than leaving it to be
    rediscovered: a Gaussian's own tails score ~1.6, so a diffuse embedding clips even with
    no flung components, putting ~1% of cells outside the drawn view. Harmless (the count is
    printed, and real UMAPs are clumpy rather than Gaussian) but not obvious."""
    rng = np.random.RandomState(1)
    xy = rng.normal(0, 8, size=(20_000, 2))
    xlim, _, n_outside, ratio = pcu.umap_view_limits(xy)
    assert 1.5 < ratio < 2.0, f"expected a marginal ratio, got {ratio:.2f}"
    assert xlim is not None
    assert n_outside / len(xy) < 0.02


def test_napari_colours_match_the_shared_colour_rule_and_are_legible_on_black():
    sys.path.insert(0, str(_PY_DIR.parents[1] / "scripts"))
    import napari_colors as nc
    from matplotlib.colors import to_hex
    labels = ["Excitatory_L4_5_IT_like", "Microglia A", "Vascular"]
    colors = pcu.category_colors(labels, napari=True)
    print([to_hex(c) for c in colors])
    assert [to_hex(c).upper() for c in colors] == [nc.deterministic_color(label) for label in labels]
    assert all(nc.contrast_against_black(to_hex(c)) >= nc.MIN_CONTRAST_AGAINST_BLACK - 0.02 for c in colors)


def test_default_colours_are_distinct_tab20_maps():
    colors = pcu.category_colors([f"type{i}" for i in range(24)], napari=False)
    assert colors.shape == (24, 4)
    assert len({tuple(c) for c in colors}) == 24


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
