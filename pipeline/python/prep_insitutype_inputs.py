#!/usr/bin/env python3
"""Emit the compact InSituType input from the stage-3a cohort AnnData.

Stage 4a of the analysis pipeline. InSituType (stage 4b, R) needs a gene-probe counts
matrix plus a per-cell negative-probe mean (its background estimate). Reading a full
AnnData .h5ad in R is awkward, so — exactly as stage 3a does for the PCA — this writes a
small, R-friendly HDF5 that pipeline/R/insitutype_typing.R reads with hdf5r:

  /counts/{data,indices,indptr,shape}  gene counts as CSC (genes x cells) slots
  /genes        gene-probe names (rows of /counts)
  /cell_id      cell ids (column order of /counts)
  /neg          per-cell mean negative-probe count (InSituType `neg` background)

Gene probes and negprobes are identified from var['probe_type'] (not column names), so
this is robust to CosMx export naming. The counts hold gene probes only; negprobes feed
only the /neg background vector. Reads the canonical combined_qc.h5ad, which keeps gene
+ neg + falsecode probes for precisely this.

Usage:
    uv run python pipeline/python/prep_insitutype_inputs.py \\
        --combined-h5ad combined_qc.h5ad \\
        --output insitutype_input.h5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-h5ad", type=Path, required=True,
                   help="Stage-3a combined_qc.h5ad (all probes, var['probe_type']).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output insitutype_input.h5.")
    return p.parse_args()


def _write_strings(grp, name: str, values: np.ndarray) -> None:
    # h5py cannot write numpy fixed-width unicode ('<U…') to a variable-length string
    # dtype ("No conversion path"); pass an object array of Python str. Mirrors the
    # helper in concat_qc_anndata.py so the R hdf5r reader sees the same string layout.
    grp.create_dataset(name, data=np.asarray(values, dtype=object),
                       dtype=h5py.string_dtype())


def main() -> None:
    args = parse_args()

    print(f"Reading {args.combined_h5ad}")
    adata = ad.read_h5ad(args.combined_h5ad)
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing; was this built by stage 3a?",
              file=sys.stderr)
        sys.exit(1)

    gene_mask = (adata.var["probe_type"] == "gene").to_numpy()
    neg_mask = (adata.var["probe_type"] == "negprobe").to_numpy()
    if not gene_mask.any():
        print("ERROR: no gene probes (probe_type == 'gene') in var.", file=sys.stderr)
        sys.exit(1)
    if not neg_mask.any():
        print("ERROR: no negprobes (probe_type == 'negprobe') in var; InSituType "
              "needs a per-cell negprobe mean for its background.", file=sys.stderr)
        sys.exit(1)
    print(f"{adata.n_obs:,} cells; {int(gene_mask.sum())} gene probes, "
          f"{int(neg_mask.sum())} negprobes")

    X = adata.X.tocsr()
    # Per-cell mean negprobe count = InSituType's `neg` background estimate.
    neg = np.asarray(X[:, neg_mask].mean(axis=1)).ravel().astype(np.float64)
    print(f"neg (per-cell negprobe mean): median {np.median(neg):.4f}, "
          f"{int((neg <= 0).sum())} cells at 0")

    # Gene counts as cells x genes CSR; its (indptr, indices, data) ARE the (@p, @i, @x)
    # of a genes x cells CSC dgCMatrix — no transpose needed (same as stage 3a).
    gene_counts = X[:, gene_mask].tocsr()
    gene_counts.sort_indices()
    gene_counts.eliminate_zeros()
    genes = adata.var_names[gene_mask].to_numpy()
    cell_id = adata.obs.index.to_numpy()

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

    print(f"Done. {len(cell_id):,} cells x {len(genes)} genes.")


if __name__ == "__main__":
    main()
