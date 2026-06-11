#!/usr/bin/env python3
"""Build the InSituType ANCHOR input: a stratified subsample of the whole cohort.

Stage 4a-anchor (two-pass cell typing for large cohorts; follows the Bruker CosMx
Scratch Space "large studies" workflow). At full-cohort scale a single semi-supervised
InSituType run over every cell is intractable (the de-novo EM does not finish), so we:

  PASS 1  fit semi-supervised InSituType ONCE on a bounded, representative subsample
          (this script's output) -> a cohort-wide profile (named GBmap types + de novo
          a,b,c... clusters) via insitutype_typing.R + insitutype_profile.R
  PASS 2  assign every slide's cells to that FIXED profile in supervised mode, in
          parallel (insitutype_supervised.R, one Slurm array task per slide)

so the de-novo labels mean the same thing in every donor.

The subsample is stratified by (slide_id x cluster) using the Stage-3c Leiden clusters,
not drawn at random: random sampling would under-represent the small, patient-private
clusters (e.g. amplicon-driven malignant programs) that are the whole reason for the
de-novo step. Capping cells per (slide x cluster) stratum guarantees every cluster on
every slide gets into the anchor, up to the cap, while bounding the total.

Output matches prep_insitutype_inputs.py exactly, so insitutype_typing.R reads it
unchanged:
  /counts/{data,indices,indptr,shape}  gene counts as CSC (genes x cells) slots
  /genes        gene-probe names (rows of /counts)
  /cell_id      cell ids (column order of /counts)
  /neg          per-cell mean negprobe count (InSituType background)
plus a sidecar anchor_cells.csv (cell_id, slide_id, cluster) recording the subsample.

Usage:
    uv run python pipeline/python/prep_insitutype_anchor.py \\
        --clustered-h5ad cosmx_clustered.h5ad \\
        --output anchor_input.h5 \\
        --cap-per-stratum 750
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="cosmx_clustered.h5ad (all probes + obs slide_id + cluster).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output anchor_input.h5 (InSituType input format).")
    p.add_argument("--cluster-key", default="leiden",
                   help="obs column to stratify by (default leiden).")
    p.add_argument("--slide-key", default="slide_id",
                   help="obs column identifying the slide (default slide_id).")
    p.add_argument("--cap-per-stratum", type=int, default=750,
                   help="Max cells sampled per (slide x cluster) stratum. Smaller "
                        "strata contribute all their cells. Tune to hit a tractable "
                        "anchor size (Wenyu's 2.33M did not finish in 5.5h).")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for reproducible subsampling.")
    return p.parse_args()


def _write_strings(grp, name: str, values: np.ndarray) -> None:
    # Match concat_qc_anndata.py / prep_insitutype_inputs.py: pass an object array of
    # Python str so h5py writes a variable-length UTF-8 string the R hdf5r reader sees.
    grp.create_dataset(name, data=np.asarray(values, dtype=object),
                       dtype=h5py.string_dtype())


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
    for key in (args.slide_key, args.cluster_key):
        if key not in adata.obs:
            print(f"ERROR: obs has no '{key}'; present: {list(adata.obs.columns)}",
                  file=sys.stderr)
            sys.exit(1)
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing; was this built by stage 3?",
              file=sys.stderr)
        sys.exit(1)

    idx = stratified_indices(adata.obs, args.slide_key, args.cluster_key,
                             args.cap_per_stratum, args.seed)
    print(f"Anchor: {len(idx):,} / {adata.n_obs:,} cells "
          f"({len(idx) / adata.n_obs:.1%}); cap {args.cap_per_stratum} per "
          f"(slide x {args.cluster_key})")

    # Materialise only the selected rows (backed row-subset read).
    sub = adata[idx].to_memory()

    gene_mask = (sub.var["probe_type"] == "gene").to_numpy()
    neg_mask = (sub.var["probe_type"] == "negprobe").to_numpy()
    if not gene_mask.any() or not neg_mask.any():
        print("ERROR: need both gene and negprobe probe_types in var.", file=sys.stderr)
        sys.exit(1)

    # Per-cluster anchor coverage — confirm every cluster made it in (esp. the rare,
    # patient-private ones the de-novo step must see).
    cov = sub.obs[args.cluster_key].value_counts().sort_index()
    print(f"Anchor cells per {args.cluster_key} cluster (min {cov.min():,}, "
          f"max {cov.max():,}):")
    print(cov.to_string())

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
    sub.obs[[args.slide_key, args.cluster_key]].rename_axis("cell_id").to_csv(sidecar)
    print(f"Wrote {sidecar} (anchor cell -> slide/cluster record)")
    print(f"Done. anchor = {len(cell_id):,} cells x {len(genes)} genes.")


if __name__ == "__main__":
    main()
