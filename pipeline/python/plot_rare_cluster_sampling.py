#!/usr/bin/env python3
"""How one small cohort-wide Leiden cluster reaches the anchor-cohort subsample (poster figure).

Three panels, one highlight colour:
  A  cohort-wide UMAP (obsm['X_umap']) with the chosen Leiden cluster highlighted
  B  one slide in tissue coordinates, the same cluster highlighted. Cluster cells that made it into the
     anchor cohort are filled; any that did not are hollow rings.
  C  for every cohort-wide Leiden cluster, the share of its cells that is in the anchor cohort. The
     dashed line is the overall share: a cluster above it was over-sampled relative to random.

The sampling strata are PER-SLIDE Leiden clusters ("<slide>|<n>" in anchor_cells.csv), not the
cohort-wide clusters drawn here, so panel C is the consequence of that stratified draw, not its
input. Everything is joined on the "<slide>_F<fov>_C<cell>" cell id, and the cohort-wide
clustering is read from the typed AnnData's obs (extract_anchor_leiden.py explains why never to
use anchor_cells.csv's own cluster column for this).

Only obs and obsm['X_umap'] are read (backed mode, no matrix). The UMAP background is thinned to
--background-cells for plotting; every cell of the highlighted cluster is always drawn.

Usage:
    uv run python pipeline/python/plot_rare_cluster_sampling.py \\
        --typed-h5ad cosmx_typed.h5ad --anchor-cells-csv anchor_cells.csv \\
        --cluster 21 --output-stem rare_cluster_21
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from plot_qc import umap_view_limits
from prep_insitucnv_input import SPATIAL_CANDIDATES

MIN_ANCHOR_OVERLAP = 0.95        # anchor cells are a subset of the typed run, so expect ~1.0
PER_SLIDE_MARK = "|"             # the slide|n separator in anchor_cells.csv's cluster column
HIGHLIGHT_COLOR = "#e0a800"
BACKGROUND_COLOR = "#d0d0d0"
OTHER_BAR_COLOR = "#8c8c8c"
UMAP_BACKGROUND_SIZE = 0.3
SLIDE_BACKGROUND_SIZE = 0.4
HIGHLIGHT_SIZE = 5.0
SHUFFLE_SEED = 0
PNG_DPI = 170
FIGURE_SIZE = (17, 6.2)
PANEL_WIDTH_RATIOS = [1, 1, 1.05]
PANEL_SPACING = 0.42


def anchor_membership(cell_ids: pd.Index, anchor_ids: pd.Series) -> np.ndarray:
    """Boolean mask over `cell_ids` of the cells that are in the anchor cohort.

    Fails loudly when the id formats differ or the files come from different runs: a silent
    mismatch would draw every cluster as unsampled.
    """
    if anchor_ids.duplicated().any():
        raise ValueError(f"{int(anchor_ids.duplicated().sum()):,} duplicated anchor cell ids")
    in_anchor = cell_ids.isin(anchor_ids)
    overlap = in_anchor.sum() / max(len(anchor_ids), 1)
    if overlap < MIN_ANCHOR_OVERLAP:
        raise ValueError(f"only {100 * overlap:.1f}% of the {len(anchor_ids):,} anchor cells are in the "
                         f"AnnData (anchor example {anchor_ids.iloc[0]!r}, AnnData example "
                         f"{cell_ids[0]!r}): different id formats or different runs")
    return np.asarray(in_anchor)


def reject_per_slide_labels(leiden: pd.Series) -> None:
    if leiden.astype(str).str.contains(PER_SLIDE_MARK, regex=False).any():
        raise ValueError(f"Leiden labels containing {PER_SLIDE_MARK!r} are per-slide strata, not the "
                         "cohort-wide clustering")


def cluster_sampling_table(leiden: pd.Series, in_anchor: np.ndarray) -> pd.DataFrame:
    """One row per cohort-wide Leiden cluster: cells, cells in the anchor cohort, and the share."""
    table = pd.DataFrame({"leiden": leiden.astype(str).to_numpy(), "in_anchor": in_anchor})
    summary = table.groupby("leiden").agg(n_cells=("in_anchor", "size"), n_sampled=("in_anchor", "sum"))
    summary["share_sampled"] = summary["n_sampled"] / summary["n_cells"]
    return summary.sort_values("n_cells", ascending=False)


def pick_slide(slide_ids: pd.Series, is_target: np.ndarray) -> str:
    """The slide holding the most cells of the target cluster."""
    return str(slide_ids[is_target].value_counts().index[0])


def spatial_columns(obs_columns: pd.Index) -> tuple[str, str]:
    for x_name, y_name in SPATIAL_CANDIDATES:
        if x_name in obs_columns and y_name in obs_columns:
            return x_name, y_name
    raise KeyError(f"no centroid columns in obs; tried {SPATIAL_CANDIDATES}")


def draw_figure(umap_xy, thin_keep, is_target, in_anchor, slide_xy, slide_target, slide_in_anchor,
                table, cluster, slide_name, path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="white")
    grid = fig.add_gridspec(1, 3, width_ratios=PANEL_WIDTH_RATIOS, wspace=PANEL_SPACING)

    ax = fig.add_subplot(grid[0])
    ax.scatter(*umap_xy[thin_keep & ~is_target].T, s=UMAP_BACKGROUND_SIZE, c=BACKGROUND_COLOR,
               linewidths=0, rasterized=True)
    ax.scatter(*umap_xy[is_target].T, s=HIGHLIGHT_SIZE, c=HIGHLIGHT_COLOR, edgecolors="black",
               linewidths=0.15, rasterized=True)
    xlim, ylim, _, _ = umap_view_limits(umap_xy)
    if xlim is not None:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
    ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_title(f"A  Cohort UMAP: Leiden {cluster} highlighted", loc="left", fontweight="bold")

    ax = fig.add_subplot(grid[1])
    ax.scatter(*slide_xy[~slide_target].T, s=SLIDE_BACKGROUND_SIZE, c=BACKGROUND_COLOR,
               linewidths=0, rasterized=True)
    kept = slide_target & slide_in_anchor
    dropped = slide_target & ~slide_in_anchor
    ax.scatter(*slide_xy[kept].T, s=HIGHLIGHT_SIZE, c=HIGHLIGHT_COLOR, edgecolors="black", linewidths=0.15)
    ax.scatter(*slide_xy[dropped].T, s=HIGHLIGHT_SIZE, facecolors="none", edgecolors=HIGHLIGHT_COLOR,
               linewidths=0.6)
    ax.set_aspect("equal")
    ax.set(xticks=[], yticks=[])
    ax.set_title(f"B  Slide {slide_name}: same cluster in tissue", loc="left", fontweight="bold")
    ax.legend(handles=[Patch(facecolor=HIGHLIGHT_COLOR, edgecolor="black", label="in anchor cohort"),
                       Patch(facecolor="none", edgecolor=HIGHLIGHT_COLOR, label="not sampled")],
              loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False, fontsize=9)

    ax = fig.add_subplot(grid[2])
    shown = table.iloc[::-1]
    colors = [HIGHLIGHT_COLOR if name == cluster else OTHER_BAR_COLOR for name in shown.index]
    ax.barh(range(len(shown)), 100 * shown["share_sampled"], color=colors, height=0.75)
    overall = 100 * table["n_sampled"].sum() / table["n_cells"].sum()
    ax.axvline(overall, color="black", linestyle="--", linewidth=1)
    ax.text(overall + 1, len(shown) - 0.5, f"overall {overall:.0f}%", fontsize=9, va="bottom")
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([f"{name}  (n={n:,})" for name, n in zip(shown.index, shown["n_cells"])],
                       fontsize=7)
    for tick, name in zip(ax.get_yticklabels(), shown.index):
        if name == cluster:
            tick.set_fontweight("bold")
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of the cluster's cells in the anchor cohort")
    ax.set_title("C  Share of each cohort Leiden cluster sampled", loc="left", fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.savefig(path, dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--typed-h5ad", type=Path, required=True)
    parser.add_argument("--anchor-cells-csv", type=Path, required=True)
    parser.add_argument("--cluster", required=True, help="cohort-wide Leiden cluster to highlight")
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--slide", default=None,
                        help="slide for panel B (default: the slide with most cells of the cluster)")
    parser.add_argument("--leiden-key", default="leiden")
    parser.add_argument("--slide-key", default="slide_id")
    parser.add_argument("--umap-key", default="X_umap")
    parser.add_argument("--background-cells", type=int, default=500_000,
                        help="UMAP background cells drawn (every cell of the cluster is always drawn)")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    try:
        for key in (args.leiden_key, args.slide_key):
            if key not in adata.obs.columns:
                sys.exit(f"ERROR: obs has no {key!r}. Present: {', '.join(map(str, adata.obs.columns))}")
        x_name, y_name = spatial_columns(adata.obs.columns)
        obs = adata.obs[[args.leiden_key, args.slide_key, x_name, y_name]].copy()
        umap_xy = np.asarray(adata.obsm[args.umap_key])[:, :2]
        cell_ids = pd.Index(adata.obs_names.astype(str))
    finally:
        adata.file.close()

    leiden = obs[args.leiden_key].astype(str).reset_index(drop=True)
    reject_per_slide_labels(leiden)
    if args.cluster not in set(leiden):
        sys.exit(f"ERROR: no cluster {args.cluster!r}; present: {sorted(set(leiden), key=int)}")

    anchor_ids = pd.read_csv(args.anchor_cells_csv, usecols=["cell_id"])["cell_id"].astype(str)
    in_anchor = anchor_membership(cell_ids, anchor_ids)
    is_target = (leiden == args.cluster).to_numpy()
    slide_ids = obs[args.slide_key].astype(str).reset_index(drop=True)

    table = cluster_sampling_table(leiden, in_anchor)
    row = table.loc[args.cluster]
    n_slides = slide_ids[is_target].nunique()
    top_slide = slide_ids[is_target].value_counts()
    print(f"Leiden {args.cluster}: {int(row['n_cells']):,} cells "
          f"({100 * row['n_cells'] / len(leiden):.2f}% of the cohort) over {n_slides} slides; "
          f"{int(row['n_sampled']):,} in the anchor cohort ({100 * row['share_sampled']:.1f}%); "
          f"overall anchor share {100 * table['n_sampled'].sum() / table['n_cells'].sum():.1f}%")
    print(f"  slide with most cells: {top_slide.index[0]} "
          f"({100 * top_slide.iloc[0] / top_slide.sum():.1f}% of the cluster)")
    print(table.to_string())

    slide_name = args.slide or pick_slide(slide_ids, is_target)
    on_slide = (slide_ids == slide_name).to_numpy()
    if not on_slide.any():
        sys.exit(f"ERROR: no cells on slide {slide_name!r}")
    slide_xy = obs.loc[on_slide, [x_name, y_name]].to_numpy(dtype=float)

    rng = np.random.default_rng(SHUFFLE_SEED)
    thin_keep = np.zeros(len(leiden), dtype=bool)
    thin_keep[rng.choice(len(leiden), size=min(args.background_cells, len(leiden)), replace=False)] = True

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(f"{args.output_stem}_sampling_by_cluster.csv")
    draw_figure(umap_xy, thin_keep, is_target, in_anchor, slide_xy, is_target[on_slide], in_anchor[on_slide],
                table, args.cluster, slide_name, Path(f"{args.output_stem}.png"))


if __name__ == "__main__":
    main()
