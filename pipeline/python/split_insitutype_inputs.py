#!/usr/bin/env python3
"""Split the QC'd cohort AnnData into per-slide InSituType/InSituTree inputs.

Full-cohort (~7.5M cell) typing shards by physical slide: InSituTree (and the flat
top-K re-score) are SUPERVISED against a fixed profile set, so every cell's label
depends only on its own counts + the fixed reference — cells are independent and can
be typed one slide at a time, then concatenated. Sharding keeps each job small and
ckpt-friendly (no single 7.5M-cell / big-memory monolith) and is embarrassingly
parallel across the ~57 slides.

This writes one <slide_id>.h5 per slide from combined_qc.h5ad — using the SAME,
already-QC'd cells the cohort clustering was built from, so typing is consistent with
the rest of the pipeline (rather than re-deriving QC from the raw per-slide Stage-1
h5ads).

Each <slide_id>.h5 matches prep_insitutype_inputs.py EXACTLY, so both
pipeline/R/insitutree_typing.R and pipeline/R/flat_posteriors.R read it unchanged:
  /counts/{data,indices,indptr,shape}  gene counts CSC (genes x cells)
  /genes /cell_id /neg

NOTE on background: the R drivers estimate InSituType's per-cell background
(estimateBackground) from the cells they are given, so a sharded run estimates it
per slide rather than cohort-wide. That is a defensible (arguably better) difference
— the slide is the natural batch — but it means a sharded result is not bit-identical
to a single monolithic run. The per-cell node assignments themselves are unchanged:
InSituTree's node gene-selection is driven by the fixed profiles, not the data.

Reads combined_qc.h5ad in backed mode and materialises one slide at a time, so memory
stays modest even at full-cohort scale.

Usage:
    uv run python pipeline/python/split_insitutype_inputs.py \\
        --combined-h5ad combined_qc.h5ad \\
        --output-dir per_slide_inputs \\
        --slide-key slide_id
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
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Directory to write one <slide_id>.h5 per slide.")
    p.add_argument("--slide-key", default="slide_id",
                   help="obs column identifying the slide (default slide_id).")
    return p.parse_args()


def _write_strings(grp, name: str, values: np.ndarray) -> None:
    # Object array of Python str -> variable-length UTF-8 (matches the other stage-4
    # writers so the R hdf5r reader sees the same string layout).
    grp.create_dataset(name, data=np.asarray(values, dtype=object),
                       dtype=h5py.string_dtype())


def write_slide_input(sub: ad.AnnData, gene_mask: np.ndarray, neg_mask: np.ndarray,
                      out_path: Path) -> int:
    """Emit one slide's InSituType input; returns the cell count."""
    X = sub.X.tocsr()
    neg = np.asarray(X[:, neg_mask].mean(axis=1)).ravel().astype(np.float64)
    gene_counts = X[:, gene_mask].tocsr()
    gene_counts.sort_indices()
    gene_counts.eliminate_zeros()
    genes = sub.var_names[gene_mask].to_numpy()
    cell_id = sub.obs.index.to_numpy()

    with h5py.File(out_path, "w") as f:
        g = f.create_group("counts")
        g.create_dataset("data", data=gene_counts.data.astype(np.int32))
        g.create_dataset("indices", data=gene_counts.indices.astype(np.int32))
        g.create_dataset("indptr", data=gene_counts.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(
            [gene_counts.shape[1], gene_counts.shape[0]], dtype=np.int64))  # genes x cells
        _write_strings(f, "genes", genes)
        _write_strings(f, "cell_id", cell_id)
        f.create_dataset("neg", data=neg)
    return len(cell_id)


def main() -> None:
    args = parse_args()

    print(f"Reading {args.combined_h5ad} (backed)")
    adata = ad.read_h5ad(args.combined_h5ad, backed="r")
    if args.slide_key not in adata.obs:
        print(f"ERROR: obs has no '{args.slide_key}'; present: "
              f"{list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing; was this built by stage 3?",
              file=sys.stderr)
        sys.exit(1)

    gene_mask = (adata.var["probe_type"] == "gene").to_numpy()
    neg_mask = (adata.var["probe_type"] == "negprobe").to_numpy()
    if not gene_mask.any() or not neg_mask.any():
        print("ERROR: need both gene and negprobe probe_types in var.", file=sys.stderr)
        sys.exit(1)

    slides = adata.obs[args.slide_key].astype(str)
    uniq = sorted(slides.unique())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"{adata.n_obs:,} cells across {len(uniq)} slides -> {args.output_dir}")

    total = 0
    for i, slide in enumerate(uniq, 1):
        mask = (slides == slide).to_numpy()
        sub = adata[mask].to_memory()
        out_path = args.output_dir / f"{slide}.h5"
        n = write_slide_input(sub, gene_mask, neg_mask, out_path)
        total += n
        print(f"  [{i:>2}/{len(uniq)}] {slide}: {n:,} cells -> {out_path.name}")

    print(f"Done. {len(uniq)} slides, {total:,} cells "
          f"({'matches' if total == adata.n_obs else 'MISMATCH vs'} "
          f"{adata.n_obs:,} cohort cells).")


if __name__ == "__main__":
    main()
