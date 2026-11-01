#!/usr/bin/env python3
"""Render a UMAP coloured by cell_type with a readable, count-sorted legend.

Stage-4 figure helper. scanpy's default categorical legend is unreadable with 20–60
cell types, so this draws the UMAP directly in matplotlib: rasterised points, a distinct
colour per type (tab20 x3), and a legend placed OUTSIDE the axes, one row per type,
sorted by abundance with per-type cell counts. Optionally annotates the de-novo letters
on the fly from a denovo_annotations mapping CSV (so `a` shows as `a - MES/AC-like tumor`).

With --hierarchy (an InSituTree hierarchy JSON) the figure is organised by the tree instead,
because 60+ leaves cannot each get a distinguishable colour:
  --output                  every cell coloured by its top-level branch ("compartment"), the
                            sink branch grey underneath, each compartment named on the plot
  <output stem>_leaves.png  one panel per leaf (>= --min-panel-cells cells) in hierarchy order,
                            its cells in the compartment colour over a grey density of all cells
  <output stem>_counts.csv  cells per label with its compartment, i.e. what the figures show
Every label must belong to the hierarchy (a leaf or a node name, as InSituTree emits a
single-child node's name); an unknown label stops the run rather than being drawn unlabelled.

With --label-colors-csv (label,hex_color, e.g. pipeline/reference/insitutree_k27_label_colors.csv
from scripts/napari-label-colors.py) every cell type takes its Napari colour, so the UMAP and the
slides agree. Those colours are made legible on Napari's black background, so the figures that
use them are drawn on black: the flat UMAP, the leaf panels, and in hierarchy mode an extra
<output stem>_celltypes.png with every cell in its type's colour and the legend grouped by branch.
A label missing from the table stops the run.

Reads a stage-4c cosmx_typed.h5ad (obs has cell_type, obsm has X_umap) in backed mode, so the
expression matrix is never loaded.

Axes are clipped to a percentile view (plot_qc.umap_view_limits) when UMAP outliers would
otherwise squash the manifold into a few pixels - the same failure the stage-3c QC plots hit.
Outliers are reported on stdout, never silently dropped: they stay in the counts and the
legend, they just fall outside the drawn view.

Run in a Slurm/salloc job (e.g. via 95_celltype_umap.sh), not on a login node: the
full-cohort scatter holds millions of points.

Usage:
    uv run python pipeline/python/plot_celltype_umap.py \\
        --h5ad cosmx_typed.h5ad --output umap_cell_type_annotated.png \\
        --mapping pipeline/reference/denovo_annotations/stage4_extl3_rescale.csv \\
        --title "Extended L3 - rescale"
    uv run python pipeline/python/plot_celltype_umap.py \\
        --h5ad cosmx_typed.h5ad --output umap_cell_type_hierarchy.png \\
        --hierarchy pipeline/reference/insitutree_hierarchy_k27.json \\
        --label-colors-csv pipeline/reference/insitutree_k27_label_colors.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from plot_qc import umap_view_limits

# Fixed categorical order for the compartments (validated colour-blind-safe on adjacent pairs);
# compartments take these in hierarchy order, never cycled.
COMPARTMENT_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                      "#e87ba4", "#008300", "#6250d6", "#e34948")
SINK_LABEL = "Low_signal"
SINK_COLOR = "#cfcfcf"
SINK_PANEL_COLOR = "#555555"
# A compartment this small has no visible mass of its own, so a name at its median would sit on
# another compartment; the legend still lists it.
MIN_LABEL_FRACTION = 0.01
MIN_PANEL_CELLS = 1000
PANEL_COLUMNS = 6
PANEL_SIZE_IN = 4.0
DENSITY_BINS = 700
BACKDROP_MAX_ALPHA_FRACTION = 0.55
PANEL_POINT_SIZE = 0.4
SHUFFLE_SEED = 0
PNG_DPI = 150
# Figures in Napari colours sit on black, like the slides; the all-cell density behind each leaf
# panel then tops out at this grey level so the leaf's own colour stays the brightest thing.
DARK_BACKGROUND = "black"
DARK_INK = "white"
DARK_BACKDROP_MAX_LIGHTNESS = 0.22
LEGEND_ROWS_PER_COLUMN = 36


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True, help="cosmx_typed.h5ad (X_umap in obsm).")
    p.add_argument("--output", type=Path, required=True, help="Output PNG.")
    p.add_argument("--mapping", type=Path, default=None,
                   help="Optional denovo_annotations CSV (denovo_label,annotation) to "
                        "relabel de-novo letters before plotting.")
    p.add_argument("--hierarchy", type=Path, default=None,
                   help="InSituTree hierarchy JSON: colour by top-level branch and add the "
                        "per-leaf small multiples and counts table.")
    p.add_argument("--sink-label", default=SINK_LABEL,
                   help="Top-level branch drawn grey underneath (hierarchy mode).")
    p.add_argument("--min-panel-cells", type=int, default=MIN_PANEL_CELLS,
                   help="Labels with fewer cells get no small-multiples panel (hierarchy mode).")
    p.add_argument("--label-colors-csv", type=Path, default=None,
                   help="label,hex_color table of Napari colours; draws cell types in them, on black.")
    p.add_argument("--color-key", default="cell_type", help="obs column to colour by.")
    p.add_argument("--umap-key", default="X_umap", help="obsm key with the 2D embedding.")
    p.add_argument("--title", default=None, help="Plot title.")
    p.add_argument("--point-size", type=float, default=1.0)
    return p.parse_args()


def hierarchy_order(hierarchy: dict) -> list[tuple[str, str]]:
    """Every node and leaf name in depth-first order, each with its top-level branch."""
    order: list[tuple[str, str]] = []

    def walk(node, compartment):
        if isinstance(node, list):
            order.extend((leaf, compartment) for leaf in node)
            return
        for name, child in node.items():
            order.append((name, compartment))
            walk(child, compartment)

    for compartment, child in hierarchy.items():
        order.append((compartment, compartment))
        walk(child, compartment)
    return order


def compartment_of(labels: list[str], hierarchy: dict) -> dict[str, str]:
    """Map each observed label to its top-level branch; fail on labels outside the tree."""
    lookup = dict(hierarchy_order(hierarchy))
    unknown = sorted(set(labels) - set(lookup))
    if unknown:
        raise KeyError(f"{len(unknown)} labels are not in the hierarchy: {unknown[:10]}")
    return {label: lookup[label] for label in labels}


def compartment_colors(hierarchy: dict, sink_label: str) -> dict[str, str]:
    colored = [name for name in hierarchy if name != sink_label]
    if len(colored) > len(COMPARTMENT_COLORS):
        raise ValueError(f"{len(colored)} top-level branches but only {len(COMPARTMENT_COLORS)} "
                         f"compartment colours; group branches before plotting")
    colors = dict(zip(colored, COMPARTMENT_COLORS))
    colors[sink_label] = SINK_COLOR
    return colors


def read_label_colors(path: Path, labels) -> dict[str, str]:
    """The Napari colour of each observed label; a label the table lacks stops the run."""
    table = pd.read_csv(path, dtype=str)
    colors = dict(zip(table["label"], table["hex_color"]))
    missing = sorted(set(labels) - set(colors))
    if missing:
        raise KeyError(f"{len(missing)} labels have no colour in {path}: {missing[:10]}")
    return colors


def strip_axes(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def darken(fig, axes) -> None:
    fig.patch.set_facecolor(DARK_BACKGROUND)
    for ax in axes:
        ax.set_facecolor(DARK_BACKGROUND)
        ax.title.set_color(DARK_INK)


def apply_view(ax, xy: np.ndarray) -> None:
    xlim, ylim, n_outside, ratio = umap_view_limits(xy)
    if xlim is not None:
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        print(f"Clipped axes to the percentile view: the autoscaled view was {ratio:.0f}x "
              f"too wide, {n_outside:,} of {len(xy):,} cells fall outside it")


def plot_flat(xy: np.ndarray, s: pd.Series, args) -> None:
    """Every type its own colour (Napari's, else tab20 x3), legend sorted by abundance."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D

    codes = s.cat.codes.to_numpy()
    labels = list(s.cat.categories)
    counts = Counter(s)

    if args.label_colors_csv is None:
        # Distinct colour per category (combine three tab20 maps for up to 60).
        base = []
        for name in ("tab20", "tab20b", "tab20c"):
            base += [plt.get_cmap(name)(i) for i in range(20)]
        colors = np.array([base[i % len(base)] for i in range(len(labels))])
    else:
        napari = read_label_colors(args.label_colors_csv, labels)
        colors = np.array([to_rgba(napari[label]) for label in labels])
    point_colors = colors[codes]

    # Plot points in shuffled order so no single type is drawn on top.
    order = np.random.RandomState(SHUFFLE_SEED).permutation(len(codes))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(xy[order, 0], xy[order, 1], c=point_colors[order],
               s=args.point_size, linewidths=0, rasterized=True)
    apply_view(ax, xy)
    strip_axes(ax)
    if args.title:
        ax.set_title(args.title, fontsize=14)

    # Legend outside, one row per type, sorted by abundance (most cells first).
    by_count = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=colors[i], markeredgewidth=0) for i in by_count]
    leg = [f"{labels[i]}  ({counts[labels[i]]:,})" for i in by_count]
    ax.legend(handles, leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
              ncol=1 if len(labels) <= 36 else 2, fontsize=8, frameon=False,
              handletextpad=0.4, labelspacing=0.3,
              labelcolor=DARK_INK if args.label_colors_csv else None)
    if args.label_colors_csv is not None:
        darken(fig, [ax])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output} ({len(codes):,} cells, {len(labels)} types)")


def plot_compartments(xy, compartments: np.ndarray, colors: dict[str, str], order: list[str],
                      sink_label: str, args) -> None:
    """Cells coloured by top-level branch; the sink grey underneath, the rest shuffled on top."""
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D

    is_sink = compartments == sink_label
    typed = np.flatnonzero(~is_sink)
    typed = typed[np.random.RandomState(SHUFFLE_SEED).permutation(len(typed))]
    names = list(colors)
    rgba = np.array([to_rgba(colors[name]) for name in names])
    typed_colors = rgba[pd.Categorical(compartments[typed], categories=names).codes]

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.scatter(*xy[is_sink].T, c=SINK_COLOR, s=args.point_size, linewidths=0, rasterized=True)
    ax.scatter(xy[typed, 0], xy[typed, 1], c=typed_colors, s=args.point_size,
               linewidths=0, rasterized=True)
    apply_view(ax, xy)
    strip_axes(ax)
    for name in order:
        mask = compartments == name
        if name == sink_label or mask.sum() < MIN_LABEL_FRACTION * len(compartments):
            continue
        x, y = np.median(xy[mask], axis=0)
        ax.text(x, y, name, fontsize=11, fontweight="bold", ha="center", va="center",
                color="black", path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    if args.title:
        ax.set_title(args.title, fontsize=14)

    counts = Counter(compartments)
    present = [name for name in order if counts[name]]
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=7,
                      markerfacecolor=colors[name], markeredgewidth=0) for name in present]
    ax.legend(handles, [f"{name}  ({counts[name]:,})" for name in present],
              loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10, frameon=False,
              handletextpad=0.4, labelspacing=0.5)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output} ({len(xy):,} cells, {len(present)} compartments)")


def plot_celltypes(xy, labels: np.ndarray, table: pd.DataFrame, colors: dict[str, str],
                   sink_labels: set[str], title: str | None, path: Path, point_size: float) -> None:
    """Every cell in its type's Napari colour on black; the sink underneath, the rest shuffled on
    top; the legend lists the types under their branch, in hierarchy order."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D

    names = table["label"].tolist()
    rgba = np.array([to_rgba(colors[name]) for name in names])
    codes = pd.Categorical(labels, categories=names).codes
    is_sink = np.isin(labels, list(sink_labels))
    typed = np.flatnonzero(~is_sink)
    typed = typed[np.random.RandomState(SHUFFLE_SEED).permutation(len(typed))]

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.scatter(*xy[is_sink].T, c=rgba[codes[is_sink]], s=point_size, linewidths=0, rasterized=True)
    ax.scatter(xy[typed, 0], xy[typed, 1], c=rgba[codes[typed]], s=point_size,
               linewidths=0, rasterized=True)
    apply_view(ax, xy)
    strip_axes(ax)
    if title:
        ax.set_title(title, fontsize=14)

    handles, texts = [], []
    for compartment, group in table.groupby("compartment", sort=False):
        handles.append(Line2D([], [], linestyle=""))
        texts.append(compartment.upper())
        for row in group.itertuples():
            handles.append(Line2D([0], [0], marker="o", linestyle="", markersize=5,
                                  markerfacecolor=colors[row.label], markeredgewidth=0))
            texts.append(f"{row.label}  ({row.n_cells:,})")
    ax.legend(handles, texts, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=7, ncol=-(-len(texts) // LEGEND_ROWS_PER_COLUMN), handletextpad=0.4,
              labelspacing=0.25, labelcolor=DARK_INK)
    darken(fig, [ax])
    fig.savefig(path, dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path} ({len(xy):,} cells, {len(names)} types in Napari colours)")


def plot_leaf_panels(xy, labels: np.ndarray, panels: list[str], panel_colors: dict[str, str],
                     path: Path, dark: bool = False) -> None:
    """One panel per label: its cells over a grey density of every cell (on black when `dark`)."""
    import matplotlib.pyplot as plt

    xlim, ylim, _, _ = umap_view_limits(xy)
    xlim = xlim or (xy[:, 0].min(), xy[:, 0].max())
    ylim = ylim or (xy[:, 1].min(), xy[:, 1].max())
    density, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=DENSITY_BINS, range=[xlim, ylim])
    backdrop = np.log1p(density.T)
    if dark:
        backdrop_cmap, backdrop_vmax = "gray", backdrop.max() / DARK_BACKDROP_MAX_LIGHTNESS
    else:
        backdrop_cmap, backdrop_vmax = "Greys", backdrop.max() / BACKDROP_MAX_ALPHA_FRACTION

    rows = -(-len(panels) // PANEL_COLUMNS)
    fig, axes = plt.subplots(rows, PANEL_COLUMNS, squeeze=False,
                             figsize=(PANEL_COLUMNS * PANEL_SIZE_IN, rows * PANEL_SIZE_IN))
    for ax, label in zip(axes.flat, panels):
        ax.imshow(backdrop, origin="lower", extent=(*xlim, *ylim), aspect="auto",
                  cmap=backdrop_cmap, vmin=0, vmax=backdrop_vmax, interpolation="nearest")
        mask = labels == label
        ax.scatter(xy[mask, 0], xy[mask, 1], c=[panel_colors[label]], s=PANEL_POINT_SIZE,
                   linewidths=0, rasterized=True)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_title(f"{label}  ({int(mask.sum()):,})", fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in list(axes.flat)[len(panels):]:
        ax.axis("off")
    if dark:
        darken(fig, axes.flat)
    fig.tight_layout()
    fig.savefig(path, dpi=PNG_DPI, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Wrote {path} ({len(panels)} panels)")


def plot_hierarchy(xy: np.ndarray, s: pd.Series, args) -> None:
    hierarchy = json.loads(args.hierarchy.read_text())
    if args.sink_label not in hierarchy:
        raise KeyError(f"sink branch '{args.sink_label}' is not a top-level branch of {args.hierarchy}")
    labels = s.astype(str).to_numpy()
    label_to_compartment = compartment_of(sorted(set(labels)), hierarchy)
    colors = compartment_colors(hierarchy, args.sink_label)
    compartments = pd.Series(labels).map(label_to_compartment).to_numpy()
    plot_compartments(xy, compartments, colors, list(hierarchy), args.sink_label, args)

    rank = {name: i for i, (name, _) in enumerate(hierarchy_order(hierarchy))}
    counts = pd.Series(labels).value_counts()
    table = pd.DataFrame({"label": counts.index, "n_cells": counts.to_numpy()})
    table["compartment"] = table["label"].map(label_to_compartment)
    table = table.iloc[np.argsort(table["label"].map(rank).to_numpy(), kind="stable")]
    stem = args.output.with_suffix("")
    table[["compartment", "label", "n_cells"]].to_csv(f"{stem}_counts.csv", index=False)
    print(f"Wrote {stem}_counts.csv ({len(table)} labels)")

    panels = table.loc[table["n_cells"] >= args.min_panel_cells, "label"].tolist()
    skipped = table.loc[table["n_cells"] < args.min_panel_cells]
    print(f"Small multiples: {len(panels)} panels; skipped (< {args.min_panel_cells} cells): "
          + (", ".join(f"{r.label} ({r.n_cells})" for r in skipped.itertuples()) or "none"))
    if args.label_colors_csv is None:
        panel_colors = {label: SINK_PANEL_COLOR if label_to_compartment[label] == args.sink_label
                        else colors[label_to_compartment[label]] for label in panels}
    else:
        napari = read_label_colors(args.label_colors_csv, table["label"])
        sink_labels = {label for label, c in label_to_compartment.items() if c == args.sink_label}
        plot_celltypes(xy, labels, table, napari, sink_labels, args.title,
                       Path(f"{stem}_celltypes.png"), args.point_size)
        panel_colors = {label: napari[label] for label in panels}
    if not panels:
        print(f"No label reaches --min-panel-cells {args.min_panel_cells}; no small-multiples figure")
        return
    plot_leaf_panels(xy, labels, panels, panel_colors, Path(f"{stem}_leaves.png"),
                     dark=args.label_colors_csv is not None)


def main() -> None:
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")

    print(f"Reading {args.h5ad} (backed; obs and obsm only)")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        if args.umap_key not in adata.obsm:
            print(f"ERROR: obsm['{args.umap_key}'] missing.", file=sys.stderr); sys.exit(1)
        if args.color_key not in adata.obs:
            print(f"ERROR: obs['{args.color_key}'] missing.", file=sys.stderr); sys.exit(1)
        s = adata.obs[args.color_key].astype("category")
        xy = np.asarray(adata.obsm[args.umap_key])[:, :2]
    finally:
        adata.file.close()

    if args.mapping is not None:
        m = pd.read_csv(args.mapping, dtype=str)
        rename = {k: v for k, v in zip(m["denovo_label"].str.strip(),
                                       m["annotation"].str.strip())
                  if k in set(s.cat.categories)}
        print(f"Relabeling {len(rename)} de-novo clusters from {args.mapping.name}")
        s = s.cat.rename_categories(rename)

    if args.hierarchy is None:
        plot_flat(xy, s, args)
    else:
        plot_hierarchy(xy, s, args)


if __name__ == "__main__":
    main()
