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
legend, they just fall outside the drawn view.

Only obs and obsm are read (the AnnData is opened backed), so the multi-GB expression matrix is
never loaded and this runs on a laptop.

--napari-colors gives each label the colour it has in Napari (scripts/napari_colors.py, which needs
duckdb), so the UMAP and the viewer agree; --dark draws on black, the background those colours were
made legible against.

--panels-output also writes a grid with one small panel per label (that label coloured, every other cell
faint grey), so small populations hidden under big ones in the single plot can be found.

Usage:
    uv run python pipeline/python/plot_celltype_umap.py \\
        --h5ad cosmx_typed.h5ad --output umap_cell_type_annotated.png \\
        --mapping pipeline/reference/denovo_annotations/stage4_extl3_rescale.csv \\
        --title "Extended L3 - rescale"
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from plot_qc import umap_view_limits


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
    p.add_argument("--napari-colors", action="store_true",
                   help="Colour each label as Napari does (shared hash colour, lightened if too dark for black).")
    p.add_argument("--dark", action="store_true", help="Black background with light text.")
    p.add_argument("--panels-output", type=Path, default=None,
                   help="Also write a grid PNG with one panel per label, most cells first.")
    return p.parse_args()


PANEL_COLUMNS = 6
PANEL_SIZE_INCHES = 3.0
PANEL_BACKGROUND_CELLS = 60_000     # faint context drawn in every panel; all cells would be 24x the work
PANEL_BACKGROUND_COLOR = "#303030"
PANEL_BACKGROUND_COLOR_LIGHT = "#d8d8d8"


def draw_panels(xy: np.ndarray, codes: np.ndarray, labels: list[str], colors: np.ndarray, counts: Counter,
                limits: tuple, output: Path, title: str | None, dark: bool, point_size: float) -> None:
    """One small panel per label, most cells first: that label in its colour over faint grey context."""
    import matplotlib.pyplot as plt
    text_color = "white" if dark else "black"
    background = "black" if dark else "white"
    context_color = PANEL_BACKGROUND_COLOR if dark else PANEL_BACKGROUND_COLOR_LIGHT
    order = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    n_rows = -(-len(order) // PANEL_COLUMNS)
    fig, axes = plt.subplots(n_rows, PANEL_COLUMNS, figsize=(PANEL_COLUMNS * PANEL_SIZE_INCHES, n_rows * PANEL_SIZE_INCHES),
                             facecolor=background, squeeze=False)
    rng = np.random.RandomState(0)
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


def category_colors(labels: list[str], napari: bool) -> np.ndarray:
    """RGBA per label: Napari's colour for the label, or three tab20 maps combined (up to 60 distinct)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    if napari:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
        from napari_colors import deterministic_color
        return np.array([to_rgba(deterministic_color(label)) for label in labels])
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        base += [plt.get_cmap(name)(i) for i in range(20)]
    return np.array([base[i % len(base)] for i in range(len(labels))])


def main() -> None:
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    print(f"Reading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    if args.umap_key not in adata.obsm:
        print(f"ERROR: obsm['{args.umap_key}'] missing.", file=sys.stderr); sys.exit(1)
    if args.color_key not in adata.obs:
        print(f"ERROR: obs['{args.color_key}'] missing.", file=sys.stderr); sys.exit(1)

    s = adata.obs[args.color_key].astype("category")
    if args.mapping is not None:
        m = pd.read_csv(args.mapping, dtype=str)
        rename = {k: v for k, v in zip(m["denovo_label"].str.strip(),
                                       m["annotation"].str.strip())
                  if k in set(s.cat.categories)}
        print(f"Relabeling {len(rename)} de-novo clusters from {args.mapping.name}")
        s = s.cat.rename_categories(rename)

    xy = np.asarray(adata.obsm[args.umap_key])[:, :2]
    codes = s.cat.codes.to_numpy()
    labels = list(s.cat.categories)
    counts = Counter(s)

    colors = category_colors(labels, args.napari_colors)
    point_colors = colors[codes]
    text_color = "white" if args.dark else "black"

    # Plot points in shuffled order so no single type is drawn on top.
    order = np.random.RandomState(0).permutation(len(codes))
    fig, ax = plt.subplots(figsize=(12, 8), facecolor="black" if args.dark else "white")
    ax.set_facecolor("black" if args.dark else "white")
    ax.scatter(xy[order, 0], xy[order, 1], c=point_colors[order],
               s=args.point_size, linewidths=0, rasterized=True)
    xlim, ylim, n_outside, ratio = umap_view_limits(xy)
    if xlim is not None:
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        print(f"Clipped axes to the percentile view: the autoscaled view was {ratio:.0f}x "
              f"too wide, {n_outside:,} of {len(codes):,} cells fall outside it")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if args.title:
        ax.set_title(args.title, fontsize=14, color=text_color)

    # Legend outside, one row per type, sorted by abundance (most cells first).
    by_count = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=colors[i], markeredgewidth=0) for i in by_count]
    leg = [f"{labels[i]}  ({counts[labels[i]]:,})" for i in by_count]
    legend = ax.legend(handles, leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
                       ncol=1 if len(labels) <= 36 else 2, fontsize=8, frameon=False,
                       handletextpad=0.4, labelspacing=0.3)
    for text in legend.get_texts():
        text.set_color(text_color)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Wrote {args.output} ({adata.n_obs:,} cells, {len(labels)} types)")
    if args.panels_output is not None:
        draw_panels(xy, codes, labels, colors, counts, (xlim, ylim), args.panels_output, args.title,
                    args.dark, args.point_size)
        print(f"Wrote {args.panels_output} ({len(labels)} panels)")


if __name__ == "__main__":
    main()
