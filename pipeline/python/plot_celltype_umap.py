#!/usr/bin/env python3
"""Render a UMAP coloured by cell_type with a readable, count-sorted legend.

Stage-4 figure helper. scanpy's default categorical legend is unreadable with 20–60
cell types, so this draws the UMAP directly in matplotlib: rasterised points, a distinct
colour per type (tab20 x3), and a legend placed OUTSIDE the axes, one row per type,
sorted by abundance with per-type cell counts. Optionally annotates the de-novo letters
on the fly from a denovo_annotations mapping CSV (so `a` shows as `a - MES/AC-like tumor`).

Reads a stage-4c cosmx_typed.h5ad (obs has cell_type, obsm has X_umap).

Memory: loads the full typed AnnData (~3.6GB at cohort scale) — run in a Slurm/salloc
job (e.g. via 95_celltype_umap.sh), not on a login node.

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
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    print(f"Reading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
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

    # Distinct colour per category (combine three tab20 maps for up to 60).
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        base += [plt.get_cmap(name)(i) for i in range(20)]
    colors = np.array([base[i % len(base)] for i in range(len(labels))])
    point_colors = colors[codes]

    # Plot points in shuffled order so no single type is drawn on top.
    order = np.random.RandomState(0).permutation(len(codes))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(xy[order, 0], xy[order, 1], c=point_colors[order],
               s=args.point_size, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if args.title:
        ax.set_title(args.title, fontsize=14)

    # Legend outside, one row per type, sorted by abundance (most cells first).
    by_count = sorted(range(len(labels)), key=lambda i: -counts[labels[i]])
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=colors[i], markeredgewidth=0) for i in by_count]
    leg = [f"{labels[i]}  ({counts[labels[i]]:,})" for i in by_count]
    ax.legend(handles, leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
              ncol=1 if len(labels) <= 36 else 2, fontsize=8, frameon=False,
              handletextpad=0.4, labelspacing=0.3)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Wrote {args.output} ({adata.n_obs:,} cells, {len(labels)} types)")


if __name__ == "__main__":
    main()
