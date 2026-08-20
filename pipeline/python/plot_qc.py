#!/usr/bin/env python3
"""UMAP QC plots from a clustered AnnData, colored by discrete obs keys.

Shared by stage 3c (cluster_embedding.py imports make_qc_plots) and runnable
standalone to regenerate plots from an existing cosmx_clustered.h5ad WITHOUT
re-running the GPU pipeline — handy for recoloring/tweaking. CPU only.

    apptainer exec "$APPTAINER_RSC" python pipeline/python/plot_qc.py \\
        --h5ad cosmx_clustered.h5ad --out-dir qc_plots --color Case,Region,leiden

Each color key is cast to an unordered categorical before plotting: integer-like
IDs such as Case (donor numbers, e.g. 7134) would otherwise render on a
continuous colormap instead of as discrete groups.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Batch-correction review: Case (patient) should be well-mixed after correction;
# Region (Tumor bulk / Infiltrating edge / Contralateral uninvolved) and leiden
# should drive the structure.
DEFAULT_QC_COLOR = "Case,Region,leiden"

# UMAP routinely flings a few near-disconnected kNN components hundreds of units away.
# Autoscaled axes then span ~60x the real manifold in each direction, compressing every
# cell into a few pixels: the retina run had 98% of 849,510 cells inside a 26x21 box while
# the axes spanned 1520x1263. Clipping the view to a percentile range fixes the picture
# without touching the embedding. Outliers are reported, never silently hidden.
AXIS_CLIP_PERCENTILE = 0.5
AXIS_CLIP_MARGIN = 0.05
# Below this overshoot the autoscaled view is already fine, so leave it alone.
AXIS_CLIP_MIN_RATIO = 1.5


def umap_view_limits(coords, percentile: float = AXIS_CLIP_PERCENTILE,
                     margin: float = AXIS_CLIP_MARGIN):
    """Axis limits that frame the bulk of a UMAP, ignoring flung-out outliers.

    Returns (xlim, ylim, n_outside, ratio), or (None, None, 0, ratio) when the
    autoscaled view is already reasonable. `ratio` is the worst-axis overshoot of the
    full data range over the percentile range, i.e. how many times too wide the
    default view would be.
    """
    import numpy as np

    lims, ratio = [], 1.0
    for axis in range(2):
        c = coords[:, axis]
        lo, hi = np.percentile(c, [percentile, 100.0 - percentile])
        span = hi - lo
        if span <= 0:
            return None, None, 0, 1.0
        ratio = max(ratio, (c.max() - c.min()) / span)
        pad = span * margin
        lims.append((lo - pad, hi + pad))
    if ratio < AXIS_CLIP_MIN_RATIO:
        return None, None, 0, ratio
    (x0, x1), (y0, y1) = lims
    outside = int((
        (coords[:, 0] < x0) | (coords[:, 0] > x1)
        | (coords[:, 1] < y0) | (coords[:, 1] > y1)
    ).sum())
    return (x0, x1), (y0, y1), outside, ratio


def make_qc_plots(adata, color_keys: list[str], out_dir: Path) -> None:
    """Save one UMAP PNG per obs key, for reviewing batch correction.

    Non-critical: a missing column or plotting error is warned and skipped rather
    than failing the caller. Points are rasterized for the large cell count, and
    each key is cast to a categorical so discrete groups get distinct colors.

    Axes are clipped to a percentile view when UMAP outliers would otherwise make the
    plot unreadable; see AXIS_CLIP_PERCENTILE.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_dir
    sc.settings.set_figure_params(dpi=150, frameon=False)

    xlim = ylim = None
    if "X_umap" in adata.obsm:
        xlim, ylim, n_outside, ratio = umap_view_limits(adata.obsm["X_umap"])
        if xlim is None:
            print(f"  UMAP axis overshoot {ratio:.1f}x — autoscale is fine, not clipping")
        else:
            print(f"  UMAP axis overshoot {ratio:.1f}x — clipping view to "
                  f"p{AXIS_CLIP_PERCENTILE}-p{100 - AXIS_CLIP_PERCENTILE}; "
                  f"{n_outside:,} of {adata.n_obs:,} cells "
                  f"({n_outside / max(adata.n_obs, 1):.3%}) fall outside the frame")

    for key in color_keys:
        if key not in adata.obs:
            print(f"WARN: QC color '{key}' not in obs; skipping", file=sys.stderr)
            continue
        # Force discrete coloring: numeric IDs (e.g. Case = donor numbers) would
        # otherwise be drawn as a continuous gradient.
        adata.obs[key] = adata.obs[key].astype(str).astype("category")
        try:
            # save= is not used: scanpy writes the file before we can set the limits,
            # so take the Axes back and save it ourselves under the same filename.
            ax = sc.pl.umap(adata, color=key, show=False, size=2, legend_fontsize=6)
            if isinstance(ax, list):
                ax = ax[0]
            if xlim is not None:
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
            out_path = out_dir / f"umap_{key}.png"
            ax.figure.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(ax.figure)
            print(f"  wrote {out_path}")
        except Exception as exc:  # plotting must never sink the pipeline
            print(f"WARN: failed to plot UMAP by '{key}': {exc}", file=sys.stderr)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True,
                   help="Clustered AnnData (e.g. cosmx_clustered.h5ad) with obsm['X_umap'].")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for the PNGs.")
    p.add_argument("--color", default=DEFAULT_QC_COLOR,
                   help="Comma-separated obs keys to color UMAP plots by.")
    args = p.parse_args()

    import anndata as ad

    print(f"Reading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    if "X_umap" not in adata.obsm:
        print("ERROR: obsm['X_umap'] missing; is this a stage-3c clustered .h5ad?",
              file=sys.stderr)
        sys.exit(1)
    make_qc_plots(adata, [c.strip() for c in args.color.split(",") if c.strip()],
                  args.out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
