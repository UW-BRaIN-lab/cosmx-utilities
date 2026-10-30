#!/usr/bin/env python3
"""UMAP of the full cohort with the anchor cohort coloured by its InSituType call.

Every cell sits at its cohort-wide Stage-3c UMAP position (obsm['X_umap'] of the typed
AnnData). Cells NOT in the anchor cohort are drawn grey underneath; the anchor-cohort cells
are drawn on top, coloured by the label in anchor_typing.h5 (the k=27 Stage-4 fit).

Two figures are written:
  * <stem>_all_labels.png     every label its own colour
  * <stem>_denovo_letters.png de novo letters coloured, named reference types in one dark
                              neutral, so the letters being worked out stand out
  * <stem>_small_multiples.png one panel per de novo letter (plus one for the named types):
                              that group highlighted over a grey density of every cell

Only obsm['X_umap'] and the obs index are read from the typed AnnData (backed mode, no
matrix), but the 7.5M-point scatter is memory heavy: run it in a Slurm job
(pipeline/slurm/100_anchor_cohort_umap.sh), not on a login node.

Usage:
    uv run python pipeline/python/plot_anchor_cohort_umap.py \\
        --typed-h5ad cosmx_typed.h5ad --anchor-typing-h5 anchor_typing.h5 \\
        --label-colors-csv pipeline/reference/anchor_label_colors.csv \\
        --output-stem umap_anchor_cohort
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from plot_qc import umap_view_limits

OFF_COHORT_LABEL = "not in anchor cohort"
NAMED_TYPES_LABEL = "named reference types"
OFF_COHORT_COLOR = "#d9d9d9"
NAMED_TYPES_COLOR = "#6b6b6b"
OFF_COHORT_POINT_SIZE = 0.2
SHUFFLE_SEED = 0
PNG_DPI = 150
LEGEND_ROWS_PER_COLUMN = 36
COLORMAPS = ("tab20", "tab20b", "tab20c")
MIN_PANEL_CELLS = 1000
PANEL_COLUMNS = 6
PANEL_SIZE_IN = 4.0
DENSITY_BINS = 700
HIGHLIGHT_POINT_SIZE = 0.4
BACKDROP_GREYS = "Greys"
BACKDROP_MAX_ALPHA_FRACTION = 0.55


def is_denovo_letter(label: str) -> bool:
    return label.isalpha() and label.islower()


def read_anchor_labels(anchor_typing_h5: Path) -> pd.Series:
    with h5py.File(anchor_typing_h5) as handle:
        return pd.Series(handle["cell_type"][:].astype(str),
                         index=handle["cell_id"][:].astype(str), name="cell_type")


def align_to_embedding(cell_ids: pd.Index, labels: pd.Series) -> tuple[np.ndarray, pd.Series]:
    """Boolean mask of anchor-cohort cells in the embedding, plus their labels in embedding order.

    Fails loudly if an anchor-cohort cell has no UMAP position: dropping it silently would
    make the grey/coloured split lie about which cells were in the fit.
    """
    missing = labels.index.difference(cell_ids)
    if len(missing):
        raise KeyError(f"{len(missing):,} anchor-cohort cells are absent from the typed AnnData "
                       f"(first: {list(missing[:3])})")
    in_cohort = cell_ids.isin(labels.index)
    return in_cohort, labels.reindex(cell_ids[in_cohort])


def read_label_colors(label_colors_csv: Path) -> dict[str, str]:
    table = pd.read_csv(label_colors_csv, dtype=str)
    return dict(zip(table["label"], table["hex_color"]))


def label_colors(labels: list[str], denovo_only: bool,
                 fixed_colors: dict[str, str] | None = None) -> dict[str, str]:
    """One colour per label: `fixed_colors` (the Napari colours) when given, else a tab20 palette."""
    import matplotlib.pyplot as plt
    colored = [label for label in labels if is_denovo_letter(label) or not denovo_only]
    if fixed_colors is not None:
        missing = [label for label in colored if label not in fixed_colors]
        if missing:
            raise KeyError(f"no colour for labels: {missing}")
        colors = {label: fixed_colors[label] for label in colored}
    else:
        palette = [plt.get_cmap(name)(i) for name in COLORMAPS for i in range(20)]
        colors = {label: palette[i % len(palette)] for i, label in enumerate(colored)}
    colors.update({label: NAMED_TYPES_COLOR for label in labels if label not in colors})
    return colors


def draw(xy_all, in_cohort, cohort_labels, colors, title, path, denovo_only):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    counts = cohort_labels.value_counts()
    point_size = 0.5
    cohort_xy = xy_all[in_cohort]
    order = np.random.RandomState(SHUFFLE_SEED).permutation(len(cohort_xy))
    point_colors = np.array([colors[label] for label in cohort_labels], dtype=object)

    fig, ax = plt.subplots(figsize=(13, 9))
    ax.scatter(*xy_all[~in_cohort].T, c=OFF_COHORT_COLOR, s=OFF_COHORT_POINT_SIZE,
               linewidths=0, rasterized=True)
    ax.scatter(cohort_xy[order, 0], cohort_xy[order, 1],
               c=list(point_colors[order]), s=point_size, linewidths=0, rasterized=True)
    xlim, ylim, n_outside, ratio = umap_view_limits(xy_all)
    if xlim is not None:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        print(f"Clipped axes: autoscaled view was {ratio:.0f}x too wide, "
              f"{n_outside:,} of {len(xy_all):,} cells fall outside it")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=14)

    entries = [(OFF_COHORT_LABEL, OFF_COHORT_COLOR, int((~in_cohort).sum()))]
    if denovo_only:
        named_total = int(sum(n for label, n in counts.items() if not is_denovo_letter(label)))
        entries += [(label, colors[label], int(n)) for label, n in counts.items()
                    if is_denovo_letter(label)]
        entries.append((NAMED_TYPES_LABEL, NAMED_TYPES_COLOR, named_total))
    else:
        entries += [(label, colors[label], int(n)) for label, n in counts.items()]
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=color, markeredgewidth=0) for _, color, _ in entries]
    ax.legend(handles, [f"{name}  ({n:,})" for name, _, n in entries],
              loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8,
              ncol=max(1, -(-len(entries) // LEGEND_ROWS_PER_COLUMN)),
              handletextpad=0.4, labelspacing=0.3)
    fig.savefig(path, dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def select_panels(cohort_labels: pd.Series, min_cells: int) -> tuple[list[str], list[str]]:
    """De novo letters with at least `min_cells` cells (most cells first), and those skipped."""
    counts = cohort_labels[cohort_labels.map(is_denovo_letter)].value_counts()
    return list(counts[counts >= min_cells].index), list(counts[counts < min_cells].index)


def draw_small_multiples(xy_all, in_cohort, cohort_labels, colors, panels, path):
    """One panel per letter: its cells over a grey density of EVERY cell (anchor cohort or not)."""
    import matplotlib.pyplot as plt

    xlim, ylim, _, _ = umap_view_limits(xy_all)
    xlim = xlim or (xy_all[:, 0].min(), xy_all[:, 0].max())
    ylim = ylim or (xy_all[:, 1].min(), xy_all[:, 1].max())
    density, _, _ = np.histogram2d(xy_all[:, 0], xy_all[:, 1], bins=DENSITY_BINS,
                                   range=[xlim, ylim])
    backdrop = np.log1p(density.T)
    backdrop_vmax = backdrop.max() / BACKDROP_MAX_ALPHA_FRACTION

    cohort_xy = xy_all[in_cohort]
    groups = panels + [NAMED_TYPES_LABEL]
    rows = -(-len(groups) // PANEL_COLUMNS)
    fig, axes = plt.subplots(rows, PANEL_COLUMNS, squeeze=False,
                             figsize=(PANEL_COLUMNS * PANEL_SIZE_IN, rows * PANEL_SIZE_IN))
    for ax, group in zip(axes.flat, groups):
        ax.imshow(backdrop, origin="lower", extent=(*xlim, *ylim), aspect="auto",
                  cmap=BACKDROP_GREYS, vmin=0, vmax=backdrop_vmax, interpolation="nearest")
        if group == NAMED_TYPES_LABEL:
            selected = ~cohort_labels.map(is_denovo_letter).to_numpy()
            color = NAMED_TYPES_COLOR
        else:
            selected = (cohort_labels == group).to_numpy()
            color = colors[group]
        ax.scatter(cohort_xy[selected, 0], cohort_xy[selected, 1], c=[color],
                   s=HIGHLIGHT_POINT_SIZE, linewidths=0, rasterized=True)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(f"{group}  ({int(selected.sum()):,})", fontsize=13)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in list(axes.flat)[len(groups):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=PNG_DPI)
    plt.close(fig)
    print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--typed-h5ad", type=Path, required=True)
    parser.add_argument("--anchor-typing-h5", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--label-colors-csv", type=Path, default=None,
                        help="label,hex_color table (pipeline/reference/anchor_label_colors.csv) so the "
                             "UMAP matches the Napari colours; default is a generic tab20 palette.")
    parser.add_argument("--min-panel-cells", type=int, default=MIN_PANEL_CELLS,
                        help="letters with fewer anchor-cohort cells get no small-multiples panel.")
    parser.add_argument("--umap-key", default="X_umap")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    try:
        if args.umap_key not in adata.obsm:
            sys.exit(f"ERROR: obsm['{args.umap_key}'] missing from {args.typed_h5ad}")
        cell_ids = pd.Index(adata.obs_names.astype(str))
        xy_all = np.asarray(adata.obsm[args.umap_key])[:, :2]
    finally:
        adata.file.close()

    in_cohort, cohort_labels = align_to_embedding(cell_ids, read_anchor_labels(args.anchor_typing_h5))
    print(f"{len(cell_ids):,} cells in the embedding; {in_cohort.sum():,} in the anchor cohort "
          f"({cohort_labels.nunique()} labels); {(~in_cohort).sum():,} grey")
    labels = list(cohort_labels.value_counts().index)
    fixed_colors = read_label_colors(args.label_colors_csv) if args.label_colors_csv else None
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix, denovo_only in (("all_labels", False), ("denovo_letters", True)):
        draw(xy_all, in_cohort, cohort_labels, label_colors(labels, denovo_only, fixed_colors),
             f"Anchor cohort by InSituType call (k=27) - {suffix.replace('_', ' ')}",
             Path(f"{args.output_stem}_{suffix}.png"), denovo_only)

    panels, skipped = select_panels(cohort_labels, args.min_panel_cells)
    print(f"Small multiples: {len(panels)} letter panels; skipped (< {args.min_panel_cells} cells): "
          f"{', '.join(f'{l} ({(cohort_labels == l).sum()})' for l in skipped) or 'none'}")
    draw_small_multiples(xy_all, in_cohort, cohort_labels,
                         label_colors(labels, True, fixed_colors), panels,
                         Path(f"{args.output_stem}_small_multiples.png"))


if __name__ == "__main__":
    main()
