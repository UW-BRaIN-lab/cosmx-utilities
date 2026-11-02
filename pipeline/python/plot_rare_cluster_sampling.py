#!/usr/bin/env python3
"""How one small cohort-wide Leiden cluster reaches the anchor-cohort subsample (poster figure).

Three panels, one highlight colour:
  A  cohort-wide UMAP (obsm['X_umap']) with the chosen Leiden cluster highlighted
  B  one slide in tissue coordinates, the same cluster highlighted. Cluster cells that made it into the
     anchor cohort are filled; any that did not are hollow rings.
  C  the per-slide Leiden strata of that slide, as bars of cells sampled. The dashed line is the
     per-stratum cap (the largest stratum count in anchor_cells.csv); a stratum below it was kept
     whole. The gold part of each bar is the highlighted cohort cluster.

The sampling strata are PER-SLIDE Leiden clusters ("<slide>|<n>" in anchor_cells.csv's
slide_cluster column), not the cohort-wide clusters in panels A and B. Everything is joined on the
"<slide>_F<fov>_C<cell>" cell id, and the cohort-wide clustering is read from the typed AnnData's
obs (extract_anchor_leiden.py explains why never to use slide_cluster for that).

The per-cluster sampling table (cells and share in the anchor cohort for every cohort cluster) is
printed and written to <stem>_sampling_by_cluster.csv.

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
STRATUM_COLUMN = "slide_cluster"
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


def sampling_cap(strata: pd.Series) -> int:
    """The per-stratum cap: strata larger than it are sampled down to exactly it."""
    return int(strata.value_counts().max())


def slide_strata_table(anchor: pd.DataFrame, slide_name: str) -> pd.DataFrame:
    """Sampled cells per per-slide stratum on one slide, split into highlighted cluster and the rest.

    `anchor` has one row per anchor cell with `slide_cluster` ("<slide>|<n>") and a boolean is_target.
    """
    parts = anchor[STRATUM_COLUMN].str.rpartition(PER_SLIDE_MARK)
    on_slide = parts[0] == slide_name
    if not on_slide.any():
        raise ValueError(f"no anchor cells on slide {slide_name!r}")
    table = anchor[on_slide].groupby(parts.loc[on_slide, 2].to_numpy()).agg(
        n_sampled=("is_target", "size"), n_target=("is_target", "sum"))
    return table.sort_values("n_sampled", ascending=False)


def spatial_columns(obs_columns: pd.Index) -> tuple[str, str]:
    for x_name, y_name in SPATIAL_CANDIDATES:
        if x_name in obs_columns and y_name in obs_columns:
            return x_name, y_name
    raise KeyError(f"no centroid columns in obs; tried {SPATIAL_CANDIDATES}")


def draw_figure(umap_xy, thin_keep, is_target, in_anchor, slide_xy, slide_target, slide_in_anchor,
                strata, cap, cluster, slide_name, path):
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
    shown = strata.iloc[::-1]
    rows = range(len(shown))
    other = shown["n_sampled"] - shown["n_target"]
    ax.barh(rows, other, color=OTHER_BAR_COLOR, height=0.75)
    ax.barh(rows, shown["n_target"], left=other, color=HIGHLIGHT_COLOR, edgecolor="black",
            linewidth=0.3, height=0.75)
    ax.axvline(cap, color="black", linestyle="--", linewidth=1)
    n_whole = int((strata["n_sampled"] < cap).sum())
    ax.text(cap, len(shown) - 0.3, f"cap {cap:,} per stratum", fontsize=9, ha="right", va="bottom")
    ax.text(0.88, 0.02, f"{n_whole} of {len(strata)} strata below the cap\n(all their cells kept)",
            transform=ax.transAxes, fontsize=9, ha="right", va="bottom")
    ax.set_yticks(list(rows))
    ax.set_yticklabels([f"stratum {name}" for name in shown.index], fontsize=7)
    ax.set_xlabel(f"cells sampled on slide {slide_name}")
    ax.set_title("C  Per-slide strata sampled into the anchor cohort", loc="left", fontweight="bold")
    ax.legend(handles=[Patch(facecolor=HIGHLIGHT_COLOR, edgecolor="black", label=f"Leiden {cluster}"),
                       Patch(facecolor=OTHER_BAR_COLOR, label="other cohort clusters")],
              loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=9)
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

    anchor = pd.read_csv(args.anchor_cells_csv, usecols=["cell_id", STRATUM_COLUMN], dtype=str)
    in_anchor = anchor_membership(cell_ids, anchor["cell_id"])
    is_target = (leiden == args.cluster).to_numpy()
    anchor["is_target"] = pd.Series(is_target, index=cell_ids).reindex(anchor["cell_id"]).fillna(False).astype(bool).to_numpy()
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
    cap = sampling_cap(anchor[STRATUM_COLUMN])
    strata = slide_strata_table(anchor, slide_name)
    print(f"per-stratum cap (largest stratum count): {cap:,}; slide {slide_name}: {len(strata)} strata, "
          f"{int((strata['n_sampled'] < cap).sum())} below the cap")
    draw_figure(umap_xy, thin_keep, is_target, in_anchor, slide_xy, is_target[on_slide], in_anchor[on_slide],
                strata, cap, args.cluster, slide_name, Path(f"{args.output_stem}.png"))


if __name__ == "__main__":
    main()
