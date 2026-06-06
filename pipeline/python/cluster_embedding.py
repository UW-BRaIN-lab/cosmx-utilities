#!/usr/bin/env python3
"""GPU graph clustering + UMAP on the Pearson-residual PCA embedding.

Stage 3c of the analysis pipeline. Takes the combined AnnData from stage 3a and
the batch-corrected PCA embedding from stage 3b (the R scPearsonPCA step), runs
the neighbor graph, Leiden clustering, and UMAP on the GPU via rapids-singlecell,
and writes the annotated AnnData.

This is where the L40S earns its keep: graph construction + Leiden + UMAP over ~8M
cells is minutes on a GPU versus hours on CPU. It operates only on the small
cells x npcs embedding (the dense Pearson-residual matrix never exists), so a
single L40S is sufficient — no dask-cuda cluster needed here.

Input:
  --combined-h5ad   <combined_qc>.h5ad      from stage 3a (obs + var + raw counts)
  --embedding       embedding .h5           from stage 3b (datasets: embedding,
                                            cell_id, pc, loadings, gene)
Output:
  --output          .h5ad with obs['leiden'] cluster labels, obsm['X_pearson_pca']
                    the embedding, obsm['X_umap'] the 2D projection.

Usage:
    apptainer exec --nv $APPTAINER_RSC python pipeline/python/cluster_embedding.py \\
        --combined-h5ad combined_qc.h5ad --embedding embedding.h5 \\
        --output cosmx_clustered.h5ad --resolution 1.2 --n-neighbors 15
"""

from __future__ import annotations

import os

# numba's default ctypes CUDA driver wrapper segfaults at cuCtxGetDevice on
# Hyak's 580 / CUDA-13 driver (hit inside rsc.tl.leiden). Force numba onto
# NVIDIA's official cuda-python bindings (cuda-bindings is installed), which
# work. Must be set before numba.cuda is imported — i.e., before rapids below.
os.environ.setdefault("NUMBA_CUDA_USE_NVIDIA_BINDING", "1")

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


REP_KEY = "X_pearson_pca"

# Vignette clusters at high resolution (1.2) "with the plan to condense afterwards".
DEFAULT_RESOLUTION = 1.2
DEFAULT_N_NEIGHBORS = 15

# Batch-correction review: Case (patient) should be well-mixed after correction;
# Region (Tumor bulk / Infiltrating edge / Contralateral uninvolved) and leiden
# should drive the structure.
DEFAULT_QC_COLOR = "Case,Region,leiden"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-h5ad", type=Path, required=True,
                   help="Combined AnnData from stage 3a.")
    p.add_argument("--embedding", type=Path, required=True,
                   help="Pearson-PCA embedding .h5 from stage 3b.")
    p.add_argument("--output", type=Path, required=True, help="Output .h5ad path.")
    p.add_argument("--resolution", type=float, default=DEFAULT_RESOLUTION,
                   help="Leiden resolution.")
    p.add_argument("--n-neighbors", type=int, default=DEFAULT_N_NEIGHBORS,
                   help="Neighbors for the kNN graph.")
    p.add_argument("--n-pcs", type=int, default=0,
                   help="Use the first N PCs (0 = all columns in the embedding).")
    p.add_argument("--seed", type=int, default=0, help="Random state for neighbors/UMAP.")
    p.add_argument("--qc-plots-dir", type=Path, default=None,
                   help="Directory for UMAP QC PNGs (default: alongside --output).")
    p.add_argument("--qc-color", default=DEFAULT_QC_COLOR,
                   help="Comma-separated obs keys to color UMAP QC plots by.")
    return p.parse_args()


def make_qc_plots(adata: ad.AnnData, color_keys: list[str], out_dir: Path) -> None:
    """Save one UMAP PNG per obs key, for reviewing batch correction.

    Non-critical: a missing column or plotting error is warned and skipped rather
    than failing the job. Points are rasterized for the large cell count.
    """
    import matplotlib
    matplotlib.use("Agg")
    import scanpy as sc

    out_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_dir
    sc.settings.set_figure_params(dpi=150, frameon=False)
    for key in color_keys:
        if key not in adata.obs:
            print(f"WARN: QC color '{key}' not in obs; skipping", file=sys.stderr)
            continue
        try:
            sc.pl.umap(adata, color=key, show=False, save=f"_{key}.png",
                       size=2, legend_fontsize=6)
            print(f"  wrote {out_dir}/umap_{key}.png")
        except Exception as exc:  # plotting must never sink the pipeline
            print(f"WARN: failed to plot UMAP by '{key}': {exc}", file=sys.stderr)


def load_embedding(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        emb = f["embedding"][:].astype(np.float32)
        cell_id = f["cell_id"][:]
    cell_id = np.array([c.decode() if isinstance(c, bytes) else c for c in cell_id])
    # hdf5r (which wrote this) stores R matrices with reversed dim order vs h5py,
    # so the embedding may arrive as (n_pcs, n_cells); orient to (n_cells, n_pcs).
    if emb.shape[0] != len(cell_id) and emb.shape[1] == len(cell_id):
        emb = np.ascontiguousarray(emb.T)
    return emb, cell_id


def attach_embedding(adata: ad.AnnData, emb: np.ndarray,
                     emb_cell_id: np.ndarray) -> None:
    """Place the embedding into obsm, reordered to match adata's obs index."""
    if emb.shape[0] != adata.n_obs:
        print(f"ERROR: embedding has {emb.shape[0]} cells, AnnData has "
              f"{adata.n_obs}", file=sys.stderr)
        sys.exit(1)
    order = pd_index_positions(adata.obs.index.to_numpy(), emb_cell_id)
    adata.obsm[REP_KEY] = emb[order]


def pd_index_positions(target: np.ndarray, source: np.ndarray) -> np.ndarray:
    """Positions in `source` for each id in `target`; errors if they don't match."""
    src_pos = {cid: i for i, cid in enumerate(source)}
    try:
        return np.fromiter((src_pos[c] for c in target), dtype=np.int64,
                           count=len(target))
    except KeyError as exc:
        print(f"ERROR: embedding is missing cell id {exc}; stage 3a/3b cell sets "
              f"disagree", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    args = parse_args()

    # Route GPU allocations through an RMM managed-memory pool BEFORE importing
    # rapids/cupy, so neighbors/Leiden/UMAP on millions of cells can oversubscribe
    # VRAM and spill to host RAM instead of segfaulting when a single card's
    # memory is exceeded (see containers/README "Multi-GPU memory tips"). Both
    # cupy (via the allocator) and the native cuVS/cugraph libs (via the global
    # RMM resource) then use managed memory.
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp

    rmm.reinitialize(managed_memory=True, pool_allocator=True)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    import rapids_singlecell as rsc

    print(f"Reading {args.combined_h5ad}")
    adata = ad.read_h5ad(args.combined_h5ad)

    print(f"Reading embedding {args.embedding}")
    emb, emb_cell_id = load_embedding(args.embedding)
    if args.n_pcs and args.n_pcs < emb.shape[1]:
        emb = emb[:, : args.n_pcs]
    attach_embedding(adata, emb, emb_cell_id)
    print(f"Embedding: {adata.n_obs:,} cells x {adata.obsm[REP_KEY].shape[1]} PCs")

    print(f"Neighbors (k={args.n_neighbors}) on {REP_KEY}")
    rsc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=REP_KEY,
                     random_state=args.seed)

    print(f"Leiden clustering (resolution={args.resolution})")
    rsc.tl.leiden(adata, resolution=args.resolution, key_added="leiden",
                  random_state=args.seed)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  {n_clusters} clusters")

    print("UMAP")
    rsc.tl.umap(adata, random_state=args.seed)

    # rsc leaves the UMAP embedding (obsm) and neighbor graph (obsp) as cupy
    # arrays; cupy ndarrays and cupy sparse matrices both expose .get() ->
    # host numpy/scipy. Sweep them back before writing (X was never moved to GPU).
    for store in (adata.obsm, adata.obsp):
        for key, val in list(store.items()):
            if hasattr(val, "get"):
                store[key] = val.get()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}")
    adata.write_h5ad(args.output, compression="gzip")

    qc_dir = args.qc_plots_dir or args.output.parent / "qc_plots"
    print(f"Writing UMAP QC plots to {qc_dir}")
    make_qc_plots(adata, [k.strip() for k in args.qc_color.split(",") if k.strip()],
                  qc_dir)

    print(f"Done. {adata.n_obs:,} cells, {n_clusters} Leiden clusters.")


if __name__ == "__main__":
    main()
