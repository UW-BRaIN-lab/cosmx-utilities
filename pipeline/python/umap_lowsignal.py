#!/usr/bin/env python3
"""Stage 5i: UMAP of the Low_signal pool, coloured by the unsupervised (de-novo) clusters.

The Phase-3 rescue re-clustered the ~48% Low_signal pool into de-novo clusters (a..k). This
embeds those same cells (normalize -> log -> HVG -> PCA -> neighbours -> UMAP) and colours the
projection by de-novo cluster, tissue region, donor, and the CNV-malignant call, so the
structure of the flat pool is visible: discrete islands would argue for real subtypes; one
continuous manifold argues for the malignant continuum the rescue markers implied.

Runs on GPU via rapids-singlecell when available (minutes on ~1.1M cells) and falls back to
scanpy on CPU otherwise.

Inputs:
  --typed-h5ad  cosmx_typed.h5ad (raw counts; obs cell_type/Region/Case). Subset to Low_signal.
  --rescue      rescue_lowsignal.csv (cell_id, rescue_label, is_denovo) — the de-novo clusters.
  --cell-table  cell_cnv_table.csv.gz (optional; for the is_malignant_call colouring).
Writes (--output-dir): umap_denovo.png, umap_region.png, umap_donor.png, umap_cnv.png, and
  (--output-h5ad) the embedded AnnData (obsm['X_umap']).

Usage:
    python pipeline/python/umap_lowsignal.py --typed-h5ad cosmx_typed.h5ad \\
        --rescue rescue/rescue_lowsignal.csv --cell-table diagnostics/cell_cnv_table.csv.gz \\
        --output-dir umap
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# keep numba off the GPU JIT path that trips Hyak's driver (see cluster_embedding.py)
os.environ.setdefault("NUMBA_CUDA_ENABLE_PYNVJITLINK", "0")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--rescue", type=Path, required=True)
    p.add_argument("--cell-table", type=Path, default=None)
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--output-h5ad", type=Path, default=None)
    return p.parse_args()


def embed(adata, args) -> None:
    """normalize -> log1p -> HVG -> scale -> PCA -> neighbours -> UMAP, GPU if available."""
    try:
        import rmm
        from rmm.allocators.cupy import rmm_cupy_allocator
        import cupy as cp
        rmm.reinitialize(managed_memory=True, pool_allocator=True)
        cp.cuda.set_allocator(rmm_cupy_allocator)
        import rapids_singlecell as rsc
        print("Embedding on GPU (rapids-singlecell).")
        rsc.pp.normalize_total(adata, target_sum=1e4)
        rsc.pp.log1p(adata)
        rsc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
        rsc.pp.scale(adata, max_value=10)
        rsc.pp.pca(adata, n_comps=args.n_pcs)
        rsc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca")
        rsc.tl.umap(adata, random_state=args.seed)
        return adata
    except Exception as e:  # noqa: BLE001 — GPU stack absent or failed -> CPU scanpy
        print(f"GPU path unavailable ({e}); falling back to scanpy on CPU.", file=sys.stderr)
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata, n_comps=args.n_pcs)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca")
        sc.tl.umap(adata, random_state=args.seed)
        return adata


def scatter(xy, values, title, out, categorical=True, order=None) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    if categorical:
        cats = order or sorted(pd.unique(values.dropna()))
        cmap = plt.get_cmap("tab20" if len(cats) > 10 else "tab10")
        for i, c in enumerate(cats):
            m = (values == c).to_numpy()
            ax.scatter(xy[m, 0], xy[m, 1], s=1, alpha=0.4,
                       color=cmap(i % cmap.N), label=str(c), rasterized=True)
        ax.legend(markerscale=6, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  ncol=1 + len(cats) // 22)
    else:
        v = pd.to_numeric(values, errors="coerce").to_numpy()
        scv = ax.scatter(xy[:, 0], xy[:, 1], s=1, alpha=0.4, c=v, cmap="magma", rasterized=True)
        fig.colorbar(scv, ax=ax, fraction=0.03, pad=0.02)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {args.typed_h5ad}")
    adata = ad.read_h5ad(args.typed_h5ad)
    if "probe_type" in adata.var:
        adata = adata[:, (adata.var["probe_type"] == "gene").to_numpy()].copy()
    ls = (adata.obs["cell_type"].astype(str) == args.lowsignal_label).to_numpy()
    adata = adata[ls].copy()
    print(f"Low_signal pool: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    rescue = pd.read_csv(args.rescue).set_index("cell_id")
    adata.obs["denovo"] = rescue["rescue_label"].where(
        rescue["is_denovo"].astype(bool)).reindex(adata.obs_names).astype("object")
    if args.cell_table is not None and args.cell_table.exists():
        ct = pd.read_csv(args.cell_table, index_col=0)
        adata.obs["cnv_malignant"] = ct["is_malignant_call"].reindex(
            adata.obs_names).map({True: "CNV-malignant", False: "CNV-normal"}).astype("object")

    adata = embed(adata, args)
    xy = np.asarray(adata.obsm["X_umap"])

    scatter(xy, adata.obs["denovo"], "Low_signal UMAP — de-novo cluster (named-rescued = grey)",
            args.output_dir / "umap_denovo.png",
            order=sorted(x for x in adata.obs["denovo"].dropna().unique()))
    if "Region" in adata.obs:
        scatter(xy, adata.obs["Region"].astype(str), "Low_signal UMAP — tissue region",
                args.output_dir / "umap_region.png")
    if "Case" in adata.obs:
        scatter(xy, adata.obs["Case"].astype(str), "Low_signal UMAP — donor",
                args.output_dir / "umap_donor.png")
    if "cnv_malignant" in adata.obs:
        scatter(xy, adata.obs["cnv_malignant"].astype(str), "Low_signal UMAP — CNV call",
                args.output_dir / "umap_cnv.png")

    if args.output_h5ad is not None:
        keep = [c for c in ("cell_type", "Region", "Case", "denovo", "cnv_malignant")
                if c in adata.obs]
        out = ad.AnnData(X=None, obs=adata.obs[keep].copy(),
                         obsm={"X_umap": xy}, uns={"lowsignal_umap": True})
        out.write_h5ad(args.output_h5ad)
    print(f"Wrote UMAP figures to {args.output_dir}")


if __name__ == "__main__":
    main()
