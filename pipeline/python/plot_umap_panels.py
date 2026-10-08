#!/usr/bin/env python3
"""One UMAP panel per group (slide, donor, ...) against a grey all-cells backdrop.

The all-in-one UMAP coloured by slide or donor cannot show batch mixing: later groups are
drawn over earlier ones and 40 donor colours are indistinguishable. Here each panel
highlights ONE group in a single colour over every other cell in grey, so a group that sits
on the same arms as everyone else (well corrected) looks different from one that forms its
own island (not corrected). Every panel shares the same, outlier-clipped view.

Also writes <key>_by_<cluster-key>_share.csv: for each group, the fraction of its cells in
each cluster -- the numbers behind the picture. Reads only obs + obsm (backed mode), so it
does not load the expression matrix.

Usage:
    python plot_umap_panels.py --h5ad cosmx_clustered.h5ad --out-dir panels --key slide_id
    python plot_umap_panels.py --h5ad cosmx_clustered.h5ad --out-dir panels --key Region --ncols 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plot_qc import umap_view_limits

BACKGROUND_COLOR = "#d4d4d4"
HIGHLIGHT_COLOR = "#c0392b"
PANEL_SIZE_IN = 3.2
DEFAULT_MAX_BACKGROUND = 150_000
DEFAULT_MAX_HIGHLIGHT = 150_000
DEFAULT_NCOLS = 5
DEFAULT_CLUSTER_KEY = "leiden"
SEED = 0


def natural_order(labels) -> list[str]:
    """Numeric order when every label is an integer (donor ids), otherwise lexical."""
    labels = [str(label) for label in labels]
    try:
        return sorted(labels, key=int)
    except ValueError:
        return sorted(labels)


def share_table(groups: pd.Series, clusters: pd.Series) -> pd.DataFrame:
    """Row-normalised crosstab: for each group, the fraction of its cells in each cluster."""
    counts = pd.crosstab(groups.astype(str), clusters.astype(str))
    counts = counts.loc[natural_order(counts.index), natural_order(counts.columns)]
    return counts.div(counts.sum(axis=1), axis=0)


def plot_panels(umap: np.ndarray, groups: pd.Series, out_path: Path, *, key: str,
                ncols: int = DEFAULT_NCOLS, max_background: int = DEFAULT_MAX_BACKGROUND,
                max_highlight: int = DEFAULT_MAX_HIGHLIGHT) -> list[str]:
    """Save the panel figure; returns the group order drawn."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(SEED)
    labels = groups.astype(str).to_numpy()
    order = natural_order(set(labels))
    n_total = len(labels)

    xlim, ylim, n_outside, ratio = umap_view_limits(umap)
    if xlim is not None:
        print(f"UMAP axis overshoot {ratio:.1f}x; clipping the shared view, "
              f"{n_outside:,} of {n_total:,} cells ({n_outside / n_total:.3%}) fall outside")

    background = umap[rng.choice(n_total, size=min(n_total, max_background), replace=False)]

    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * PANEL_SIZE_IN, nrows * PANEL_SIZE_IN),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for ax, group in zip(axes.ravel(), order):
        members = np.flatnonzero(labels == group)
        shown = rng.choice(members, size=min(len(members), max_highlight), replace=False)
        ax.scatter(background[:, 0], background[:, 1], s=0.15, c=BACKGROUND_COLOR,
                   linewidths=0, rasterized=True)
        ax.scatter(umap[shown, 0], umap[shown, 1], s=0.3, c=HIGHLIGHT_COLOR,
                   linewidths=0, rasterized=True)
        if xlim is not None:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        ax.set_title(f"{group}\nn={len(members):,} ({len(members) / n_total:.1%})", fontsize=9)
    fig.suptitle(f"UMAP by {key}: each panel highlights one group over all other cells (grey)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return order


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True,
                   help="Clustered AnnData with obsm['X_umap'] (stage 3c).")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--key", default="slide_id", help="obs column to make one panel per value of.")
    p.add_argument("--cluster-key", default=DEFAULT_CLUSTER_KEY,
                   help="obs column for the group-by-cluster share table ('' to skip).")
    p.add_argument("--ncols", type=int, default=DEFAULT_NCOLS)
    p.add_argument("--max-background", type=int, default=DEFAULT_MAX_BACKGROUND)
    p.add_argument("--max-highlight", type=int, default=DEFAULT_MAX_HIGHLIGHT)
    args = p.parse_args()

    import anndata as ad

    print(f"Reading {args.h5ad} (backed; expression matrix not loaded)")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    needed = [args.key] + ([args.cluster_key] if args.cluster_key else [])
    for key in needed:
        if key not in adata.obs:
            print(f"ERROR: obs is missing '{key}'. Available: {list(adata.obs.columns)}",
                  file=sys.stderr)
            sys.exit(1)
    if "X_umap" not in adata.obsm:
        print("ERROR: obsm['X_umap'] missing; is this a stage-3c clustered .h5ad?", file=sys.stderr)
        sys.exit(1)

    umap = np.asarray(adata.obsm["X_umap"])
    groups = adata.obs[args.key].astype(str)
    print(f"{len(groups):,} cells, {groups.nunique()} '{args.key}' groups")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_png = args.out_dir / f"umap_panels_{args.key}.png"
    plot_panels(umap, groups, out_png, key=args.key, ncols=args.ncols,
                max_background=args.max_background, max_highlight=args.max_highlight)
    print(f"Wrote {out_png}")

    if args.cluster_key:
        shares = share_table(groups, adata.obs[args.cluster_key])
        out_csv = args.out_dir / f"{args.key}_by_{args.cluster_key}_share.csv"
        shares.to_csv(out_csv, float_format="%.5f")
        print(f"Wrote {out_csv}")
    print("Done.")


if __name__ == "__main__":
    main()
