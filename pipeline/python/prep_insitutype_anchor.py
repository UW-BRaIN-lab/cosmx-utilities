#!/usr/bin/env python3
"""Build the InSituType ANCHOR input: a stratified subsample of the whole cohort.

USED HERE for the InSituTree REFERENCE REBUILD: a single semi-supervised InSituType
rescale over all 7.5M cells is intractable, so we fit it on this bounded, stratified
anchor (~2.3M cells, the proven-feasible Wenyu scale) and its post-rescale $profiles
become the new InSituTree reference (prep_insitutree_profiles.py). Stratification is the
point: it guarantees the 38-donor patient-private clones reach the fit, which random
subsampling would thin out. (Below describes the original two-pass framing; only the
anchor-fit "PASS 1" applies to the reference rebuild.)


Stage 4a-anchor (two-pass cell typing for large cohorts; follows the Bruker CosMx
Scratch Space "large studies" workflow). At full-cohort scale a single semi-supervised
InSituType run over every cell is intractable (the de-novo EM does not finish), so we:

  PASS 1  fit semi-supervised InSituType ONCE on a bounded, representative subsample
          (this script's output) -> a cohort-wide profile (named GBmap types + de novo
          a,b,c... clusters) via insitutype_typing.R + insitutype_profile.R
  PASS 2  assign every slide's cells to that FIXED profile in supervised mode, in
          parallel (insitutype_supervised.R, one Slurm array task per slide)

so the de-novo labels mean the same thing in every donor.

The subsample is stratified, not drawn at random: random sampling would under-represent
the small, patient-private clusters (e.g. amplicon-driven malignant programs) that are
the whole reason for the de-novo step. Capping cells per stratum guarantees every
cluster gets into the anchor, up to the cap, while bounding the total.

The strata are PER-SLIDE Leiden clusters (default --cluster-mode per-slide). We cluster
each slide independently on the stage-3b scPearsonPCA embedding already carried in
obsm[--rep-key], then namespace the labels "<slide>|<local>". Per-slide rather than one
cohort-wide partition is deliberate: a patient-private program can be a large fraction of
its one slide yet a vanishing fraction of the whole cohort (the EGFR/CDK4/MDM2 amplicon
clones are ~few % of a single donor), so a single cohort-wide Leiden would fold it into a
common type before we ever sample it. Clustering within each slide gives every locally
distinct population its own stratum, so the cap-per-stratum sampling carries it into the
anchor. (--cluster-mode column instead stratifies on an existing obs label, e.g. the
cohort-wide leiden, for back-compat.)

Output matches prep_insitutype_inputs.py exactly, so insitutype_typing.R reads it
unchanged:
  /counts/{data,indices,indptr,shape}  gene counts as CSC (genes x cells) slots
  /genes        gene-probe names (rows of /counts)
  /cell_id      cell ids (column order of /counts)
  /neg          per-cell mean negprobe count (InSituType background)
plus a sidecar anchor_cells.csv (cell_id, slide_id, cluster) recording the subsample.

Usage:
    apptainer exec --nv $APPTAINER_RSC \\
        python pipeline/python/prep_insitutype_anchor.py \\
            --clustered-h5ad cosmx_clustered.h5ad \\
            --output anchor_input.h5 \\
            --cap-per-stratum 750
"""

from __future__ import annotations

import os

# numba's default ctypes CUDA driver wrapper segfaults at cuCtxGetDevice on Hyak's
# 580 / CUDA-13 driver (hit inside rsc.tl.leiden); force numba onto NVIDIA's official
# cuda-python bindings. Must be set before numba.cuda (i.e. rapids) is imported, so set
# it at module load even though rapids is only imported in --cluster-mode per-slide.
os.environ.setdefault("NUMBA_CUDA_USE_NVIDIA_BINDING", "1")

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

REP_KEY = "X_pearson_pca"

# Match stage 3c (cluster_embedding.py): the vignette clusters at high resolution
# "with the plan to condense afterwards". Applied per-slide here.
DEFAULT_RESOLUTION = 1.2
DEFAULT_N_NEIGHBORS = 15


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="cosmx_clustered.h5ad (all probes + obs slide_id + "
                        "obsm scPearsonPCA embedding).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output anchor_input.h5 (InSituType input format).")
    p.add_argument("--cluster-mode", choices=("per-slide", "column"),
                   default="per-slide",
                   help="per-slide: Leiden-cluster each slide on obsm[--rep-key] and "
                        "stratify on those labels (default). column: stratify on an "
                        "existing obs column (--cluster-key).")
    p.add_argument("--rep-key", default=REP_KEY,
                   help=f"obsm embedding to cluster per slide (default {REP_KEY}).")
    p.add_argument("--resolution", type=float, default=DEFAULT_RESOLUTION,
                   help="Per-slide Leiden resolution (per-slide mode).")
    p.add_argument("--n-neighbors", type=int, default=DEFAULT_N_NEIGHBORS,
                   help="Neighbors for the per-slide kNN graph (per-slide mode).")
    p.add_argument("--cluster-key", default="leiden",
                   help="obs column to stratify by in --cluster-mode column "
                        "(default leiden).")
    p.add_argument("--slide-key", default="slide_id",
                   help="obs column identifying the slide (default slide_id).")
    p.add_argument("--cap-per-stratum", type=int, default=750,
                   help="Max cells sampled per (slide x cluster) stratum. Smaller "
                        "strata contribute all their cells (so rare strata are never "
                        "lost). Tune to hit a tractable anchor size (Wenyu's 2.33M "
                        "did not finish in 5.5h).")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for subsampling and per-slide clustering.")
    return p.parse_args()


def _write_strings(grp, name: str, values: np.ndarray) -> None:
    # Match concat_qc_anndata.py / prep_insitutype_inputs.py: pass an object array of
    # Python str so h5py writes a variable-length UTF-8 string the R hdf5r reader sees.
    grp.create_dataset(name, data=np.asarray(values, dtype=object),
                       dtype=h5py.string_dtype())


def per_slide_leiden(adata: ad.AnnData, slide_key: str, rep_key: str,
                     resolution: float, n_neighbors: int, seed: int) -> np.ndarray:
    """Cluster each slide independently on obsm[rep_key]; return namespaced
    "<slide>|<local>" labels as an object array in adata.obs order.

    Reuses the stage-3b scPearsonPCA embedding already in the clustered AnnData (no
    re-derivation of PCA) and runs on the GPU via rapids-singlecell, mirroring
    cluster_embedding.py (RMM managed pool + NVIDIA numba binding). The embedding lives
    in memory even when the AnnData is opened backed, so the backed X is never touched.
    """
    if rep_key not in adata.obsm:
        print(f"ERROR: obsm has no '{rep_key}'; present: {list(adata.obsm)}",
              file=sys.stderr)
        sys.exit(1)

    # Route GPU allocations through an RMM managed-memory pool BEFORE importing
    # rapids/cupy, so the per-slide graphs can spill to host RAM rather than segfault
    # (same idiom as cluster_embedding.py).
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp

    rmm.reinitialize(managed_memory=True, pool_allocator=True)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    import rapids_singlecell as rsc

    emb = np.asarray(adata.obsm[rep_key], dtype=np.float32)
    slides = adata.obs[slide_key].astype(str).to_numpy()
    labels = np.empty(adata.n_obs, dtype=object)

    for slide in pd.unique(slides):
        rows = np.flatnonzero(slides == slide)
        n = rows.size
        if n <= n_neighbors:
            # Too few cells to build a kNN graph; the whole slide is one stratum.
            labels[rows] = f"{slide}|0"
            print(f"  {slide}: {n:,} cells -> 1 cluster (below k={n_neighbors})")
            continue
        sub = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
        sub.obsm[rep_key] = emb[rows]
        rsc.pp.neighbors(sub, n_neighbors=n_neighbors, use_rep=rep_key,
                         random_state=seed)
        rsc.tl.leiden(sub, resolution=resolution, key_added="leiden",
                      random_state=seed)
        local = sub.obs["leiden"].astype(str).to_numpy()
        labels[rows] = [f"{slide}|{c}" for c in local]
        print(f"  {slide}: {n:,} cells -> {sub.obs['leiden'].nunique()} clusters")

    return labels


def stratified_indices(obs: pd.DataFrame, slide_key: str, cluster_key: str,
                       cap: int, seed: int) -> np.ndarray:
    """Row positions for a (slide x cluster)-stratified subsample, capped per stratum."""
    rng = np.random.default_rng(seed)
    pos = np.arange(len(obs))
    keep: list[np.ndarray] = []
    # Group positions by (slide, cluster); sample up to `cap` from each.
    grouper = obs.groupby([slide_key, cluster_key], observed=True)
    for _, idx in grouper.indices.items():
        if len(idx) <= cap:
            keep.append(pos[idx])
        else:
            keep.append(pos[rng.choice(idx, size=cap, replace=False)])
    out = np.sort(np.concatenate(keep))
    return out


def main() -> None:
    args = parse_args()

    print(f"Reading {args.clustered_h5ad} (backed)")
    adata = ad.read_h5ad(args.clustered_h5ad, backed="r")
    if args.slide_key not in adata.obs:
        print(f"ERROR: obs has no '{args.slide_key}'; present: "
              f"{list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing; was this built by stage 3?",
              file=sys.stderr)
        sys.exit(1)

    if args.cluster_mode == "per-slide":
        print(f"Per-slide Leiden on obsm['{args.rep_key}'] "
              f"(resolution={args.resolution}, k={args.n_neighbors})")
        adata.obs["slide_cluster"] = per_slide_leiden(
            adata, args.slide_key, args.rep_key, args.resolution,
            args.n_neighbors, args.seed)
        cluster_key = "slide_cluster"
    else:
        cluster_key = args.cluster_key
        if cluster_key not in adata.obs:
            print(f"ERROR: obs has no '{cluster_key}'; present: "
                  f"{list(adata.obs.columns)}", file=sys.stderr)
            sys.exit(1)

    idx = stratified_indices(adata.obs, args.slide_key, cluster_key,
                             args.cap_per_stratum, args.seed)
    n_strata = adata.obs.groupby([args.slide_key, cluster_key],
                                 observed=True).ngroups
    print(f"Anchor: {len(idx):,} / {adata.n_obs:,} cells "
          f"({len(idx) / adata.n_obs:.1%}) from {n_strata:,} "
          f"(slide x {cluster_key}) strata; cap {args.cap_per_stratum}")

    # Materialise only the selected rows (backed row-subset read).
    sub = adata[idx].to_memory()

    gene_mask = (sub.var["probe_type"] == "gene").to_numpy()
    neg_mask = (sub.var["probe_type"] == "negprobe").to_numpy()
    if not gene_mask.any() or not neg_mask.any():
        print("ERROR: need both gene and negprobe probe_types in var.", file=sys.stderr)
        sys.exit(1)

    # Per-stratum anchor coverage — confirm the cap left every stratum represented
    # (especially the rare, patient-private ones the de-novo step must see).
    per_stratum = sub.obs.groupby([args.slide_key, cluster_key],
                                  observed=True).size()
    print(f"Anchor cells per stratum: min {per_stratum.min():,}, "
          f"median {int(per_stratum.median()):,}, max {per_stratum.max():,}")

    X = sub.X.tocsr()
    neg = np.asarray(X[:, neg_mask].mean(axis=1)).ravel().astype(np.float64)

    # Gene counts cells x genes CSR; its (indptr, indices, data) ARE the (@p, @i, @x)
    # of a genes x cells CSC dgCMatrix (same convention as prep_insitutype_inputs.py).
    gene_counts = X[:, gene_mask].tocsr()
    gene_counts.sort_indices()
    gene_counts.eliminate_zeros()
    genes = sub.var_names[gene_mask].to_numpy()
    cell_id = sub.obs.index.to_numpy()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output} ({gene_counts.shape[1]} genes x "
          f"{gene_counts.shape[0]:,} cells, {gene_counts.nnz:,} nonzeros)")
    with h5py.File(args.output, "w") as f:
        g = f.create_group("counts")
        g.create_dataset("data", data=gene_counts.data.astype(np.int32))
        g.create_dataset("indices", data=gene_counts.indices.astype(np.int32))
        g.create_dataset("indptr", data=gene_counts.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(
            [gene_counts.shape[1], gene_counts.shape[0]], dtype=np.int64))  # genes x cells
        _write_strings(f, "genes", genes)
        _write_strings(f, "cell_id", cell_id)
        f.create_dataset("neg", data=neg)

    sidecar = args.output.with_name("anchor_cells.csv")
    sub.obs[[args.slide_key, cluster_key]].rename_axis("cell_id").to_csv(sidecar)
    print(f"Wrote {sidecar} (anchor cell -> slide/cluster record)")
    print(f"Done. anchor = {len(cell_id):,} cells x {len(genes)} genes.")


if __name__ == "__main__":
    main()
