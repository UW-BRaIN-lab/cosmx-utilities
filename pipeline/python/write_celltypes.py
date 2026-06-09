#!/usr/bin/env python3
"""Write InSituType cell types back into the clustered cohort AnnData.

Stage 4c of the analysis pipeline. Joins the stage-4b InSituType result onto the
stage-3c clustered AnnData by cell id, producing cosmx_typed.h5ad: the same cells and
UMAP as cosmx_clustered.h5ad, with two new obs columns

  cell_type        InSituType assignment (named GBmap types + de novo letter clusters)
  insitutype_prob  per-cell posterior confidence in that assignment

and leiden retained. Downstream marker/compare tooling then re-renders by cell type with
the existing one-flag flip (marker_pseudobulk.py --group-key cell_type). A UMAP coloured
by cell_type is emitted for review via the shared plot_qc.make_qc_plots.

Usage:
    uv run python pipeline/python/write_celltypes.py \\
        --clustered-h5ad cosmx_clustered.h5ad \\
        --insitutype-h5 insitutype_result.h5 \\
        --output cosmx_typed.h5ad \\
        --qc-plots-dir qc_plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from plot_qc import make_qc_plots


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="cosmx_clustered.h5ad from stage 3c (obs has leiden + UMAP).")
    p.add_argument("--insitutype-h5", type=Path, required=True,
                   help="insitutype_result.h5 from stage 4b.")
    p.add_argument("--output", type=Path, required=True,
                   help="Output cosmx_typed.h5ad.")
    p.add_argument("--qc-plots-dir", type=Path, default=None,
                   help="If set, write a UMAP-by-cell_type PNG here.")
    return p.parse_args()


def _decode(arr) -> np.ndarray:
    """hdf5r writes variable-length UTF-8 strings; h5py reads them as bytes or str."""
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def main() -> None:
    args = parse_args()

    print(f"Reading {args.insitutype_h5}")
    with h5py.File(args.insitutype_h5, "r") as f:
        cell_id = _decode(f["cell_id"])
        cell_type = _decode(f["cell_type"])
        prob = np.asarray(f["prob"][()], dtype=np.float64)
    typed = pd.DataFrame({"cell_type": cell_type, "insitutype_prob": prob},
                         index=pd.Index(cell_id, name="cell_id"))
    print(f"  {len(typed):,} typed cells across {typed['cell_type'].nunique()} cell types")

    print(f"Reading {args.clustered_h5ad}")
    adata = ad.read_h5ad(args.clustered_h5ad)

    aligned = typed.reindex(adata.obs.index)
    n_missing = int(aligned["cell_type"].isna().sum())
    if n_missing:
        print(f"WARN: {n_missing:,} / {adata.n_obs:,} clustered cells have no "
              f"InSituType assignment (left as NA)", file=sys.stderr)
    n_extra = int(len(typed) - aligned["cell_type"].notna().sum())
    if n_extra:
        print(f"WARN: {n_extra:,} typed cells not present in the clustered AnnData "
              f"(dropped)", file=sys.stderr)

    adata.obs["cell_type"] = aligned["cell_type"].astype("category")
    adata.obs["insitutype_prob"] = aligned["insitutype_prob"].astype(np.float32)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}")
    adata.write_h5ad(args.output, compression="gzip")

    if args.qc_plots_dir is not None:
        print(f"Rendering UMAP-by-cell_type to {args.qc_plots_dir}")
        make_qc_plots(adata, ["cell_type"], args.qc_plots_dir)

    print(f"Done. {adata.n_obs:,} cells; cell_type + insitutype_prob in obs.")


if __name__ == "__main__":
    main()
