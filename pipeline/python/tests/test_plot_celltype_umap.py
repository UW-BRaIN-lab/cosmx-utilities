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


TOY_HIERARCHY = {
    "Glia": {"Astro": ["Astrocyte", "Reactive_astrocyte_denovo"], "OPC": ["OPC"]},
    "Neuron": ["Excitatory", "Inhibitory"],
    "Low_signal": ["Low_signal_denovo"],
}
TOY_COUNTS = {"Astrocyte": 300, "Reactive_astrocyte_denovo": 200, "OPC": 40,
              "Excitatory": 250, "Inhibitory": 120, "Low_signal": 400}


def test_hierarchy_order_is_depth_first_with_top_level_branch():
    order = pcu.hierarchy_order(TOY_HIERARCHY)
    print(order)
    assert [name for name, _ in order] == ["Glia", "Astro", "Astrocyte", "Reactive_astrocyte_denovo",
                                          "OPC", "OPC", "Neuron", "Excitatory", "Inhibitory",
                                          "Low_signal", "Low_signal_denovo"]
    assert dict(order)["Reactive_astrocyte_denovo"] == "Glia"


def test_single_child_node_name_maps_and_unknown_labels_fail():
    """InSituTree reports the sink by its node name (Low_signal), not its leaf; that must map.
    A label outside the tree must stop the run instead of being drawn unlabelled."""
    mapping = pcu.compartment_of(["Low_signal", "OPC", "Inhibitory"], TOY_HIERARCHY)
    assert mapping == {"Low_signal": "Low_signal", "OPC": "Glia", "Inhibitory": "Neuron"}
    try:
        pcu.compartment_of(["Astrocyte", "Microglia"], TOY_HIERARCHY)
    except KeyError as err:
        print(f"unknown label rejected: {err}")
        assert "Microglia" in str(err)
    else:
        raise AssertionError("an unknown label was accepted")


def test_compartment_colours_are_fixed_and_sink_is_grey():
    colors = pcu.compartment_colors(TOY_HIERARCHY, "Low_signal")
    assert colors == {"Glia": pcu.COMPARTMENT_COLORS[0], "Neuron": pcu.COMPARTMENT_COLORS[1],
                      "Low_signal": pcu.SINK_COLOR}
    too_many = {f"branch{i}": [f"leaf{i}"] for i in range(len(pcu.COMPARTMENT_COLORS) + 1)}
    too_many["Low_signal"] = ["Low_signal_denovo"]
    try:
        pcu.compartment_colors(too_many, "Low_signal")
    except ValueError as err:
        print(f"too many branches rejected: {err}")
    else:
        raise AssertionError("more branches than colours were accepted")


def test_k27_hierarchy_covers_every_reference_column():
    """Every profile column of the k=27 reference (the labels InSituTree can emit), plus the sink's
    node name, maps to a branch, and the branches fit the colour list."""
    import json
    import pandas as pd
    reference = _PY_DIR.parent / "reference"
    hierarchy = json.loads((reference / "insitutree_hierarchy_k27.json").read_text())
    columns = list(pd.read_csv(reference / "insitutree_profiles_k27.csv", nrows=0).columns[1:])
    mapping = pcu.compartment_of(columns + ["Low_signal"], hierarchy)
    print(pd.Series(mapping).value_counts().to_string())
    colors = pcu.compartment_colors(hierarchy, "Low_signal")
    assert len(colors) == len(hierarchy) == 8


def _write_toy_h5ad(path):
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.RandomState(2)
    labels = [label for label, n in TOY_COUNTS.items() for _ in range(n)]
    centres = {label: rng.uniform(-30, 30, size=2) for label in TOY_COUNTS}
    xy = np.vstack([rng.normal(centres[label], 2.0) for label in labels])
    obs = pd.DataFrame({"cell_type": pd.Categorical(labels)},
                       index=[f"cell{i}" for i in range(len(labels))])
    adata = ad.AnnData(X=sp.csr_matrix((len(labels), 3), dtype=np.float32), obs=obs)
    adata.obsm["X_umap"] = xy
    adata.write_h5ad(path)


def _run_main(argv):
    old_argv = sys.argv
    sys.argv = ["plot_celltype_umap.py", *argv]
    try:
        pcu.main()
    finally:
        sys.argv = old_argv


def test_hierarchy_mode_writes_compartments_panels_and_counts():
    import json
    import tempfile
    import pandas as pd
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_toy_h5ad(tmp / "typed.h5ad")
        (tmp / "hierarchy.json").write_text(json.dumps(TOY_HIERARCHY))
        _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--hierarchy", str(tmp / "hierarchy.json"),
                   "--min-panel-cells", "50", "--output", str(tmp / "umap.png")])
        assert (tmp / "umap.png").stat().st_size > 0
        assert (tmp / "umap_leaves.png").stat().st_size > 0
        counts = pd.read_csv(tmp / "umap_counts.csv")
        print(counts.to_string(index=False))
        assert counts["label"].tolist() == ["Astrocyte", "Reactive_astrocyte_denovo", "OPC",
                                            "Excitatory", "Inhibitory", "Low_signal"]
        assert dict(zip(counts["label"], counts["n_cells"])) == TOY_COUNTS
        assert counts.set_index("label").loc["OPC", "compartment"] == "Glia"


def test_flat_mode_still_writes_one_png():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_toy_h5ad(tmp / "typed.h5ad")
        _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--output", str(tmp / "flat.png")])
        assert (tmp / "flat.png").stat().st_size > 0
        assert not (tmp / "flat_counts.csv").exists()


def _write_toy_colors(path, drop=None):
    import pandas as pd
    labels = [label for label in TOY_COUNTS if label != drop]
    pd.DataFrame({"label": labels, "hex_color": ["#51F0B2", "#7449AA", "#D5DC34", "#F5C242",
                                                 "#2211D4", "#BDF613"][:len(labels)]}).to_csv(path, index=False)


def test_napari_colours_add_the_celltype_figure():
    import json
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_toy_h5ad(tmp / "typed.h5ad")
        _write_toy_colors(tmp / "colors.csv")
        (tmp / "hierarchy.json").write_text(json.dumps(TOY_HIERARCHY))
        _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--hierarchy", str(tmp / "hierarchy.json"),
                   "--label-colors-csv", str(tmp / "colors.csv"), "--min-panel-cells", "50",
                   "--output", str(tmp / "umap.png")])
        for name in ("umap.png", "umap_celltypes.png", "umap_celltypes_white.png", "umap_leaves.png",
                     "umap_counts.csv"):
            assert (tmp / name).stat().st_size > 0, name
        from PIL import Image
        corners = {name: Image.open(tmp / name).convert("RGB").getpixel((0, 0))
                   for name in ("umap_celltypes.png", "umap_celltypes_white.png")}
        print(f"top-left pixel: {corners}")
        assert corners == {"umap_celltypes.png": (0, 0, 0), "umap_celltypes_white.png": (255, 255, 255)}
        _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--label-colors-csv", str(tmp / "colors.csv"),
                   "--output", str(tmp / "flat.png")])
        assert (tmp / "flat.png").stat().st_size > 0


def test_white_legend_lists_the_types_without_branch_names():
    """The white figure's legend is the black one's minus the branch-name rows, in the same order,
    and only its swatches are outlined."""
    import pandas as pd
    table = pd.DataFrame({"compartment": ["Glia", "Glia", "Neuron", "Low_signal"],
                          "label": ["Astrocyte", "OPC", "Excitatory", "Low_signal"],
                          "n_cells": [300, 40, 250, 400]})
    colors = dict(zip(table["label"], ["#51F0B2", "#D5DC34", "#F5C242", "#BDF613"]))
    _, grouped = pcu.celltype_legend(table, colors, dark=True, branch_headers=True)
    white_handles, flat = pcu.celltype_legend(table, colors, dark=False, branch_headers=False)
    print(f"black legend: {grouped}\nwhite legend: {flat}")
    assert grouped == ["GLIA", "Astrocyte  (300)", "OPC  (40)", "NEURON", "Excitatory  (250)",
                       "LOW_SIGNAL", "Low_signal  (400)"]
    assert flat == ["Astrocyte  (300)", "OPC  (40)", "Excitatory  (250)", "Low_signal  (400)"]
    assert all(h.get_markeredgecolor() == pcu.LIGHT_SWATCH_EDGE for h in white_handles)
    dark_handles, _ = pcu.celltype_legend(table, colors, dark=True, branch_headers=False)
    assert all(h.get_markeredgewidth() == 0 for h in dark_handles)


def test_no_label_above_the_panel_minimum_skips_the_panel_figure():
    import json
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_toy_h5ad(tmp / "typed.h5ad")
        (tmp / "hierarchy.json").write_text(json.dumps(TOY_HIERARCHY))
        _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--hierarchy", str(tmp / "hierarchy.json"),
                   "--min-panel-cells", "100000", "--output", str(tmp / "umap.png")])
        assert (tmp / "umap.png").stat().st_size > 0
        assert (tmp / "umap_counts.csv").exists()
        assert not (tmp / "umap_leaves.png").exists()


def test_a_label_without_a_napari_colour_stops_the_run():
    """A cell type missing from the colour table would otherwise be drawn in some other colour
    than Napari's with no warning."""
    import json
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_toy_h5ad(tmp / "typed.h5ad")
        _write_toy_colors(tmp / "colors.csv", drop="OPC")
        (tmp / "hierarchy.json").write_text(json.dumps(TOY_HIERARCHY))
        try:
            _run_main(["--h5ad", str(tmp / "typed.h5ad"), "--hierarchy", str(tmp / "hierarchy.json"),
                       "--label-colors-csv", str(tmp / "colors.csv"), "--output", str(tmp / "umap.png")])
        except KeyError as err:
            print(f"missing colour rejected: {err}")
            assert "OPC" in str(err)
        else:
            raise AssertionError("a label without a colour was drawn")


def test_k27_colour_table_covers_every_reference_column():
    import pandas as pd
    reference = _PY_DIR.parent / "reference"
    columns = list(pd.read_csv(reference / "insitutree_profiles_k27.csv", nrows=0).columns[1:])
    colors = pcu.read_label_colors(reference / "insitutree_k27_label_colors.csv", columns + ["Low_signal"])
    print(f"{len(colors)} colours in the table")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
