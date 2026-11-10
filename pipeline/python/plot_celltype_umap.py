#!/usr/bin/env python3
"""Render a UMAP coloured by cell_type with a readable, count-sorted legend.

Stage-4 figure helper. scanpy's default categorical legend is unreadable with 20–60
cell types, so this draws the UMAP directly in matplotlib: rasterised points, a distinct
colour per type (tab20 x3), and a legend placed OUTSIDE the axes, one row per type,
sorted by abundance with per-type cell counts. Optionally annotates the de-novo letters
on the fly from a denovo_annotations mapping CSV (so `a` shows as `a - MES/AC-like tumor`).

Reads a stage-4c cosmx_typed.h5ad (obs has cell_type, obsm has X_umap).

Axes are clipped to a percentile view (plot_qc.umap_view_limits) when UMAP outliers would
otherwise squash the manifold into a few pixels - the same failure the stage-3c QC plots hit.
Outliers are reported on stdout, never silently dropped: they stay in the counts and the
legend, they just fall outside the drawn view. Cells with no call (NaN label) are left out of
the drawing and counted on stdout.

Only obs and obsm are read (the AnnData is opened backed), so the multi-GB expression matrix is
never loaded and this runs on a laptop.

--napari-colors gives each label the colour it has in Napari, so the UMAP and the viewer agree:
from --color-table (label,hex_color as the Napari columns carry it) where it has the label, else
re-derived with scripts/napari_colors.py (which needs duckdb). --dark draws on black, the
background those colours were made legible against.

--panels-output also writes a grid with one small panel per label (that label coloured, every other
cell faint grey), so small populations hidden under big ones in the single plot can be found.
--min-panel-cells leaves the smallest labels out of that grid.

With --hierarchy (an InSituTree hierarchy JSON) the figure is organised by the tree instead,
because 60+ leaves cannot each get a distinguishable colour:
  --output                     every cell coloured by its top-level branch ("compartment"), the
                               sink branch grey underneath, each compartment named on the plot
  <output stem>_counts.csv     cells per label with its compartment, in hierarchy order
  --panels-output              panels in hierarchy order, in the compartment colour (or the
                               Napari colour with --napari-colors)
  with --napari-colors also:
  <output stem>_celltypes.png        every cell in its Napari colour, legend grouped by branch
  <output stem>_celltypes_white.png  the same on white, branch names left out of the legend and
                                     each swatch outlined (half the Napari colours are under 3:1
                                     against white)
Every label must belong to the hierarchy (a leaf or a node name, as InSituTree emits a
single-child node's name); an unknown label stops the run rather than being drawn unlabelled.

Usage:
    uv run python pipeline/python/plot_celltype_umap.py \\
        --h5ad cosmx_typed.h5ad --output umap_cell_type_annotated.png \\
        --mapping pipeline/reference/denovo_annotations/stage4_extl3_rescale.csv \\
        --title "Extended L3 - rescale"
    uv run python pipeline/python/plot_celltype_umap.py \\
        --h5ad cosmx_typed.h5ad --output umap_cell_type_hierarchy.png \\
        --hierarchy pipeline/reference/insitutree_hierarchy_k27.json \\
        --napari-colors --color-table pipeline/reference/insitutree_k27_label_colors.csv --dark \\
        --panels-output umap_cell_type_hierarchy_leaves.png --min-panel-cells 1000
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

PANEL_COLUMNS = 6
PANEL_SIZE_INCHES = 3.0
PANEL_BACKGROUND_CELLS = 60_000     # faint context drawn in every panel; all cells would be 24x the work
PANEL_BACKGROUND_COLOR = "#303030"
PANEL_BACKGROUND_COLOR_LIGHT = "#d8d8d8"

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
SHUFFLE_SEED = 0
PNG_DPI = 150
LEGEND_ROWS_PER_COLUMN = 36
# Half the Napari colours are under 3:1 against white, so on white each legend swatch gets an outline.
LIGHT_SWATCH_EDGE = "#777777"
LIGHT_SWATCH_EDGE_WIDTH = 0.5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True, help="cosmx_typed.h5ad (X_umap in obsm).")
    p.add_argument("--output", type=Path, required=True, help="Output PNG.")
    p.add_argument("--mapping", type=Path, default=None,
                   help="Optional denovo_annotations CSV (denovo_label,annotation) to "
                        "relabel de-novo letters before plotting.")
    p.add_argument("--color-key", default="cell_type", help="obs column to colour by.")
    p.add_argument("--umap-key", default="X_umap", help="obsm key with the 2D embedding.")
    p.add_argument("--title", default=None, help="Plot title.")
    p.add_argument("--point-size", type=float, default=1.0)
    p.add_argument("--color-table", type=Path, default=None,
                   help="CSV with label,hex_color -- the colours the Napari columns actually "
                        "carry. Used with --napari-colors, and preferred over re-deriving them.")
    p.add_argument("--napari-colors", action="store_true",
                   help="Colour each label as Napari does (shared hash colour, lightened if too dark for black).")
    p.add_argument("--dark", action="store_true", help="Black background with light text.")
    p.add_argument("--panels-output", type=Path, default=None,
                   help="Also write a grid PNG with one panel per label, most cells first "
                        "(hierarchy order with --hierarchy).")
    p.add_argument("--min-panel-cells", type=int, default=0,
                   help="Labels with fewer cells get no panel in --panels-output.")
    p.add_argument("--hierarchy", type=Path, default=None,
                   help="InSituTree hierarchy JSON: colour by top-level branch, write the per-label "
                        "counts table and, with --napari-colors, the per-type figures.")
    p.add_argument("--sink-label", default=SINK_LABEL,
                   help="Top-level branch drawn grey underneath (hierarchy mode).")
    return p.parse_args()


def draw_panels(xy: np.ndarray, codes: np.ndarray, labels: list[str], colors: np.ndarray, counts: Counter,
                limits: tuple, output: Path, title: str | None, dark: bool, point_size: float,
                order: list[int] | None = None) -> int:
    """One small panel per label in `order` (default: most cells first): that label in its colour
    over faint grey context. Returns the number of panels; writes nothing when `order` is empty."""
    import matplotlib.pyplot as plt
    if order is None:
        order = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    if not order:
        print(f"No label to draw a panel for; {output.name} not written")
        return 0
    text_color = "white" if dark else "black"
    background = "black" if dark else "white"
    context_color = PANEL_BACKGROUND_COLOR if dark else PANEL_BACKGROUND_COLOR_LIGHT
    n_rows = -(-len(order) // PANEL_COLUMNS)
    fig, axes = plt.subplots(n_rows, PANEL_COLUMNS, figsize=(PANEL_COLUMNS * PANEL_SIZE_INCHES, n_rows * PANEL_SIZE_INCHES),
                             facecolor=background, squeeze=False)
    rng = np.random.RandomState(SHUFFLE_SEED)
    context = rng.choice(len(codes), size=min(PANEL_BACKGROUND_CELLS, len(codes)), replace=False)
    xlim, ylim = limits
    for ax in axes.ravel():
        ax.set_axis_off()
    for ax, label_index in zip(axes.ravel(), order):
        ax.set_facecolor(background)
        ax.scatter(xy[context, 0], xy[context, 1], c=context_color, s=point_size, linewidths=0, rasterized=True)
        in_label = np.flatnonzero(codes == label_index)
        ax.scatter(xy[in_label, 0], xy[in_label, 1], color=colors[label_index], s=point_size * 2, linewidths=0,
                   rasterized=True)
        if xlim is not None:
            ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_title(f"{labels[label_index]}\n({counts[labels[label_index]]:,})", fontsize=8, color=text_color)
        ax.set_axis_off()
    if title:
        fig.suptitle(title, fontsize=14, color=text_color)
    fig.tight_layout(rect=(0, 0, 1, 0.97 if title else 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)
    return len(order)


def load_color_table(path: Path) -> dict[str, str]:
    """{label -> hex} as published, which is not always what the rule would compute today.

    The Napari columns on S3 were written before `legible_color` lightened the near-black
    colours, so 14 of the 75 InSituTree leaves carry the raw hash there. A figure is meant to
    match the viewer someone is holding, so the published table wins over re-deriving.
    """
    table = pd.read_csv(path, dtype=str)
    missing = {"label", "hex_color"} - set(table.columns)
    if missing:
        sys.exit(f"ERROR: {path} has no {sorted(missing)} column(s)")
    return dict(zip(table["label"].str.strip(), table["hex_color"].str.strip()))


def category_colors(labels: list[str], napari: bool, color_table: Path | None = None) -> np.ndarray:
    """RGBA per label: Napari's colour for the label, or three tab20 maps combined (up to 60 distinct)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    if napari:
        published = load_color_table(color_table) if color_table else {}
        unknown = [label for label in labels if label not in published]
        if published and not unknown:
            return np.array([to_rgba(published[label]) for label in labels])
        try:
            sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
            from napari_colors import deterministic_color
            if unknown and published:
                print(f"NOTE: {len(unknown)} label(s) absent from {color_table.name}; "
                      f"hashing those: {unknown[:4]}", file=sys.stderr)
            return np.array([to_rgba(published.get(label) or deterministic_color(label))
                             for label in labels])
        except ImportError as exc:
            # duckdb supplies the hash and the rapids SIF has none. Fall back LOUDLY: silently
            # different colours between this figure and Napari is the worse failure.
            print(f"WARN: cannot compute Napari colours ({exc}); falling back to tab20. "
                  f"Colours will NOT match the Napari columns.", file=sys.stderr)
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        base += [plt.get_cmap(name)(i) for i in range(20)]
    return np.array([base[i % len(base)] for i in range(len(labels))])


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


def label_table(labels: list[str], counts: Counter, hierarchy: dict) -> pd.DataFrame:
    """compartment, label, n_cells for every observed label, in hierarchy order."""
    label_to_compartment = compartment_of(labels, hierarchy)
    rank = {name: i for i, (name, _) in reversed(list(enumerate(hierarchy_order(hierarchy))))}
    ordered = sorted(labels, key=rank.__getitem__)
    return pd.DataFrame({"compartment": [label_to_compartment[label] for label in ordered],
                         "label": ordered,
                         "n_cells": [counts[label] for label in ordered]})


def style_axes(fig, ax, title: str | None, dark: bool) -> None:
    background = "black" if dark else "white"
    fig.patch.set_facecolor(background)
    ax.set_facecolor(background)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, fontsize=14, color="white" if dark else "black")


def apply_limits(ax, limits: tuple) -> None:
    xlim, ylim = limits
    if xlim is not None:
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)


def save(fig, path: Path) -> None:
    import matplotlib.pyplot as plt
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=PNG_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_flat(xy, codes, labels, colors, counts, limits, args) -> None:
    """Every type its own colour, legend sorted by abundance (most cells first)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    point_colors = colors[codes]
    # Plot points in shuffled order so no single type is drawn on top.
    order = np.random.RandomState(SHUFFLE_SEED).permutation(len(codes))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(xy[order, 0], xy[order, 1], c=point_colors[order],
               s=args.point_size, linewidths=0, rasterized=True)
    apply_limits(ax, limits)
    style_axes(fig, ax, args.title, args.dark)

    by_count = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=colors[i], markeredgewidth=0) for i in by_count]
    leg = [f"{labels[i]}  ({counts[labels[i]]:,})" for i in by_count]
    ax.legend(handles, leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
              ncol=1 if len(labels) <= LEGEND_ROWS_PER_COLUMN else 2, fontsize=8, frameon=False,
              handletextpad=0.4, labelspacing=0.3, labelcolor="white" if args.dark else None)
    save(fig, args.output)
    print(f"Wrote {args.output} ({len(codes):,} cells, {len(labels)} types)")


def plot_compartments(xy, compartments: np.ndarray, colors: dict[str, str], order: list[str],
                      sink_label: str, limits: tuple, args) -> None:
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
    apply_limits(ax, limits)
    style_axes(fig, ax, args.title, args.dark)
    for name in order:
        mask = compartments == name
        if name == sink_label or mask.sum() < MIN_LABEL_FRACTION * len(compartments):
            continue
        x, y = np.median(xy[mask], axis=0)
        ax.text(x, y, name, fontsize=11, fontweight="bold", ha="center", va="center",
                color="black", path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    counts = Counter(compartments)
    present = [name for name in order if counts[name]]
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=7,
                      markerfacecolor=colors[name], markeredgewidth=0) for name in present]
    ax.legend(handles, [f"{name}  ({counts[name]:,})" for name in present],
              loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10, frameon=False,
              handletextpad=0.4, labelspacing=0.5, labelcolor="white" if args.dark else None)
    save(fig, args.output)
    print(f"Wrote {args.output} ({len(xy):,} cells, {len(present)} compartments)")


def celltype_legend(table: pd.DataFrame, colors: dict, dark: bool,
                    branch_headers: bool) -> tuple[list, list[str]]:
    """Legend handles and texts, one per type in table order, each branch's name above its types
    when `branch_headers`; on white the swatches are outlined."""
    from matplotlib.lines import Line2D

    edge = {"markeredgewidth": 0} if dark else {"markeredgewidth": LIGHT_SWATCH_EDGE_WIDTH,
                                                "markeredgecolor": LIGHT_SWATCH_EDGE}
    handles, texts = [], []
    for compartment, group in table.groupby("compartment", sort=False):
        if branch_headers:
            handles.append(Line2D([], [], linestyle=""))
            texts.append(compartment.upper())
        for row in group.itertuples():
            handles.append(Line2D([0], [0], marker="o", linestyle="", markersize=5,
                                  markerfacecolor=colors[row.label], **edge))
            texts.append(f"{row.label}  ({row.n_cells:,})")
    return handles, texts


def plot_celltypes(xy, point_labels: np.ndarray, table: pd.DataFrame, colors: dict,
                   sink_labels: set[str], limits: tuple, title: str | None, path: Path,
                   point_size: float, dark: bool, branch_headers: bool) -> None:
    """Every cell in its type's colour; the sink underneath, the rest shuffled on top; the legend
    lists the types in hierarchy order, under their branch name when `branch_headers`."""
    import matplotlib.pyplot as plt

    names = table["label"].tolist()
    rgba = np.array([colors[name] for name in names])
    codes = pd.Categorical(point_labels, categories=names).codes
    is_sink = np.isin(point_labels, list(sink_labels))
    typed = np.flatnonzero(~is_sink)
    typed = typed[np.random.RandomState(SHUFFLE_SEED).permutation(len(typed))]

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.scatter(*xy[is_sink].T, c=rgba[codes[is_sink]], s=point_size, linewidths=0, rasterized=True)
    ax.scatter(xy[typed, 0], xy[typed, 1], c=rgba[codes[typed]], s=point_size,
               linewidths=0, rasterized=True)
    apply_limits(ax, limits)
    style_axes(fig, ax, title, dark)
    handles, texts = celltype_legend(table, colors, dark, branch_headers)
    ax.legend(handles, texts, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=7, ncol=-(-len(texts) // LEGEND_ROWS_PER_COLUMN), handletextpad=0.4,
              labelspacing=0.25, labelcolor="white" if dark else None)
    save(fig, path)
    print(f"Wrote {path} ({len(xy):,} cells, {len(names)} types, on {'black' if dark else 'white'})")


def plot_hierarchy(xy, codes, labels, colors, counts, limits, args) -> tuple[list[int], np.ndarray]:
    """Compartment UMAP, counts table and (with --napari-colors) the per-type figures.
    Returns the panel order (hierarchy order) and the panel colour per label."""
    from matplotlib.colors import to_rgba

    hierarchy = json.loads(args.hierarchy.read_text())
    if args.sink_label not in hierarchy:
        raise KeyError(f"sink branch '{args.sink_label}' is not a top-level branch of {args.hierarchy}")
    table = label_table(labels, counts, hierarchy)
    label_to_compartment = dict(zip(table["label"], table["compartment"]))
    branch_colors = compartment_colors(hierarchy, args.sink_label)
    compartments = np.array([label_to_compartment[label] for label in labels])[codes]
    plot_compartments(xy, compartments, branch_colors, list(hierarchy), args.sink_label, limits, args)

    stem = args.output.with_suffix("")
    counts_path = Path(f"{stem}_counts.csv")
    table.to_csv(counts_path, index=False)
    print(f"Wrote {counts_path} ({len(table)} labels)")

    if args.napari_colors:
        by_label = dict(zip(labels, colors))
        sink_labels = {label for label, c in label_to_compartment.items() if c == args.sink_label}
        point_labels = np.asarray(labels, dtype=object)[codes]
        for path, dark, headers in ((Path(f"{stem}_celltypes.png"), args.dark, True),
                                    (Path(f"{stem}_celltypes_white.png"), False, False)):
            plot_celltypes(xy, point_labels, table, by_label, sink_labels, limits, args.title, path,
                           args.point_size, dark=dark, branch_headers=headers)
        panel_colors = colors
    else:
        panel_colors = np.array([to_rgba(SINK_PANEL_COLOR if label_to_compartment[label] == args.sink_label
                                         else branch_colors[label_to_compartment[label]])
                                 for label in labels])
    index = {label: i for i, label in enumerate(labels)}
    return [index[label] for label in table["label"]], panel_colors


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

    s = s.cat.remove_unused_categories()
    codes = s.cat.codes.to_numpy()
    called = codes >= 0
    if not called.all():
        print(f"Leaving out {int((~called).sum()):,} cells with no {args.color_key} call")
        xy, codes = xy[called], codes[called]
    labels = [str(label) for label in s.cat.categories]
    counts = Counter(dict(zip(labels, np.bincount(codes, minlength=len(labels)).tolist())))

    xlim, ylim, n_outside, ratio = umap_view_limits(xy)
    if xlim is not None:
        print(f"Clipped axes to the percentile view: the autoscaled view was {ratio:.0f}x "
              f"too wide, {n_outside:,} of {len(codes):,} cells fall outside it")
    limits = (xlim, ylim)

    colors = category_colors(labels, args.napari_colors, args.color_table)
    if args.hierarchy is None:
        plot_flat(xy, codes, labels, colors, counts, limits, args)
        panel_order, panel_colors = None, colors
    else:
        panel_order, panel_colors = plot_hierarchy(xy, codes, labels, colors, counts, limits, args)

    if args.panels_output is not None:
        if panel_order is None:
            panel_order = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
        kept = [i for i in panel_order if counts[labels[i]] >= args.min_panel_cells]
        skipped = [f"{labels[i]} ({counts[labels[i]]:,})" for i in panel_order if i not in kept]
        if skipped:
            print(f"No panel (< {args.min_panel_cells:,} cells): {', '.join(skipped)}")
        n_panels = draw_panels(xy, codes, labels, panel_colors, counts, limits, args.panels_output,
                               args.title, args.dark, args.point_size, order=kept)
        if n_panels:
            print(f"Wrote {args.panels_output} ({n_panels} panels)")


if __name__ == "__main__":
    main()
