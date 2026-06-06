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


def make_qc_plots(adata, color_keys: list[str], out_dir: Path) -> None:
    """Save one UMAP PNG per obs key, for reviewing batch correction.

    Non-critical: a missing column or plotting error is warned and skipped rather
    than failing the caller. Points are rasterized for the large cell count, and
    each key is cast to a categorical so discrete groups get distinct colors.
    """
    import matplotlib
    matplotlib.use("Agg")
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_dir
    sc.settings.set_figure_params(dpi=150, frameon=False)
    for key in color_keys:
        if key not in adata.obs:
            print(f"WARN: QC color '{key}' not in obs; skipping", file=sys.stderr)
            continue
        # Force discrete coloring: numeric IDs (e.g. Case = donor numbers) would
        # otherwise be drawn as a continuous gradient.
        adata.obs[key] = adata.obs[key].astype(str).astype("category")
        try:
            sc.pl.umap(adata, color=key, show=False, save=f"_{key}.png",
                       size=2, legend_fontsize=6)
            print(f"  wrote {out_dir}/umap_{key}.png")
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
