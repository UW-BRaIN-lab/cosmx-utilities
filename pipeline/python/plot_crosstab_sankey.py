#!/usr/bin/env python3
"""Sankey diagram from a cell-type cross-tab CSV (e.g. compare_external_typing output).

Reads a counts cross-tab (rows = left/source labeling, cols = right/destination), and
draws a two-column Sankey: left nodes -> right nodes, ribbon width = shared-cell count.
Built to visualize one cell typing vs another (e.g. Wenyu's vs our keeper) straight from
`external_vs_ours_crosstab.csv` — no need to re-touch the h5ads.

Nodes are ordered by total size; ribbons below --min-frac of all cells are accounted in
the layout but not drawn (keeps it legible). Ribbons are coloured by their source (left)
node so you can trace where each left type's cells go.

Usage (one typing vs another):
    uv run python pipeline/python/plot_crosstab_sankey.py \\
        --crosstab stage4_qc/wenyu_compare/external_vs_ours_crosstab.csv \\
        --label-left Wenyu --label-right "keeper (Ext L3 rescale)" \\
        --output stage4_qc/figures/wenyu_vs_keeper_sankey.png

Usage (unsupervised Leiden vs InSituType, off the 85d crosstab — --sort-left
natural keeps the cluster axis in 0,1,2,... order rather than by size):
    uv run python pipeline/python/plot_crosstab_sankey.py \\
        --crosstab leiden_crosstab/leiden_celltype_counts.csv \\
        --label-left "Leiden (Stage 3c)" --label-right "InSituType" \\
        --sort-left natural \\
        --output leiden_crosstab/leiden_vs_insitutype_sankey.png
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath


NATURAL_CHUNK = re.compile(r"(\d+)")


def natural_key(label: str) -> list:
    """Sort key placing '2' before '10' — Leiden clusters are numeric strings."""
    return [int(c) if c.isdigit() else c.lower()
            for c in NATURAL_CHUNK.split(str(label))]


def order_axis(totals: pd.Series, how: str) -> list:
    """Node order for one axis: by descending size, or by natural label order."""
    if how == "natural":
        return sorted(totals.index, key=natural_key)
    return list(totals.sort_values(ascending=False).index)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--crosstab", type=Path, required=True,
                   help="Counts cross-tab CSV: index = left labels, columns = right labels.")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--label-left", default="left")
    p.add_argument("--label-right", default="right")
    p.add_argument("--min-frac", type=float, default=0.0015,
                   help="Don't draw ribbons smaller than this fraction of all cells.")
    p.add_argument("--sort-left", choices=("size", "natural"), default="size",
                   help="Left node order: descending size (default), or natural label "
                        "order — use 'natural' for Leiden clusters so they read 0,1,2,...")
    p.add_argument("--sort-right", choices=("size", "natural"), default="size",
                   help="Right node order: descending size (default), or natural label order.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ct = pd.read_csv(args.crosstab, index_col=0)
    ct = ct.loc[order_axis(ct.sum(1), args.sort_left),
                order_axis(ct.sum(0), args.sort_right)]
    L, R = list(ct.index), list(ct.columns)
    total = ct.values.sum()
    thresh = args.min_frac * total

    # node colours (cycle three tab20 maps); ribbons inherit their source colour
    base = []
    for nm in ("tab20", "tab20b", "tab20c"):
        base += [plt.get_cmap(nm)(i) for i in range(20)]
    lcol = {n: base[i % len(base)] for i, n in enumerate(L)}

    gap = total * 0.012
    def layout(tots):
        y, pos = 0.0, {}
        for n, t in tots.items():
            pos[n] = (y, t); y += t + gap
        return pos, y
    Lp, yL = layout(ct.sum(1)); Rp, yR = layout(ct.sum(0))
    xL1, xR0, xm = 0.04, 0.96, 0.5
    Loff = {n: 0.0 for n in L}; Roff = {n: 0.0 for n in R}

    fig, ax = plt.subplots(figsize=(12, max(10, 0.30 * max(len(L), len(R)) + 2)))
    for s in L:                                   # ribbons, sorted big-first per source
        for d in ct.columns[ct.loc[s].argsort()[::-1]]:
            f = ct.loc[s, d]
            if f <= 0:
                continue
            y1, y2 = Lp[s][0] + Loff[s], Rp[d][0] + Roff[d]
            Loff[s] += f; Roff[d] += f
            if f < thresh:
                continue
            verts = [(xL1, y1), (xm, y1), (xm, y2), (xR0, y2),
                     (xR0, y2 + f), (xm, y2 + f), (xm, y1 + f), (xL1, y1 + f), (xL1, y1)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                     MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), facecolor=lcol[s],
                                   edgecolor="none", alpha=0.5))
    for n in L:
        ax.add_patch(Rectangle((0.0, Lp[n][0]), xL1, Lp[n][1], color=lcol[n]))
        ax.text(-0.008, Lp[n][0] + Lp[n][1] / 2, n, ha="right", va="center", fontsize=7)
    for n in R:
        ax.add_patch(Rectangle((xR0, Rp[n][0]), 1.0 - xR0, Rp[n][1], color="#888780"))
        ax.text(1.008, Rp[n][0] + Rp[n][1] / 2, n, ha="left", va="center", fontsize=7)
    ax.set_xlim(-0.30, 1.30); ax.set_ylim(0, max(yL, yR)); ax.invert_yaxis(); ax.axis("off")
    ax.text(xL1, -gap * 1.5, args.label_left, ha="center", va="bottom", fontsize=12, weight="bold")
    ax.text(xR0, -gap * 1.5, args.label_right, ha="center", va="bottom", fontsize=12, weight="bold")
    ax.set_title(f"{args.label_left} → {args.label_right}  "
                 f"({int(total):,} shared cells)", fontsize=13, pad=22)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {args.output}  ({len(L)} left x {len(R)} right nodes)")


if __name__ == "__main__":
    main()
