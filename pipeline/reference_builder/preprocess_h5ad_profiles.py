#!/usr/bin/env python3
"""Turn one or more scRNA-seq H5AD files into a pseudobulk profile CSV.

Generalized from the ocular/brain atlas project (HRCA v2 + Allen Brain Cell
Atlas). Large CELLxGENE H5ADs (millions of cells) can't be loaded into R, so we
read them here in backed (memory-mapped) mode, sum raw counts per cell type in
chunks, and write a small genes x cell-types CSV. `build_reference_profile.R`
then rescales and merges these CSVs with any other sources into the final
InSituType reference.

Two labeling modes (choose per invocation):
  --cell-type-col COL   label each cell by adata.obs[COL]. Use for atlases whose
                        H5AD carries per-cell annotations (e.g. HRCA v2).
  --label NAME          treat the whole file as a single cell type NAME. Use for
                        atlases split into one file per type (e.g. Allen
                        superclusters). Pass multiple --h5ad with one --label
                        each, or one --h5ad with one --label.

CELLxGENE conventions handled:
  * Gene names come from var[--gene-name-col] (default 'feature_name', HUGO
    symbols) since the index is usually Ensembl IDs and CosMx uses HUGO. Falls
    back to var_names if the column is absent.
  * Duplicate symbols (several Ensembl IDs -> one HUGO symbol) are summed.
  * Raw counts are located by probing layers['raw'], layers['counts'],
    adata.raw, then X; X is validated to be non-negative integers before use so
    a log/normalized matrix can't be mistaken for counts.

Output: genes x cell-types CSV (row index = gene symbol). Values are per-cell
mean raw counts by default (--normalize mean) — the same quantity
InSituType::getRNAprofiles approximates — so pseudobulk sources and getRNAprofiles
sources are on comparable footing before the R step's 99th-percentile rescale.

Usage (HRCA-style, per-cell labels):
    uv run python preprocess_h5ad_profiles.py \\
        --h5ad HRCA_v2_snRNA.h5ad --cell-type-col cell_type \\
        --prefix HRCA_ --output HRCA_v2_pseudobulk_profiles.csv

Usage (Allen-style, one label per file, combined into one CSV):
    uv run python preprocess_h5ad_profiles.py \\
        --h5ad Astrocyte.h5ad --label Astrocyte \\
        --h5ad Microglia.h5ad --label Microglia \\
        --prefix Allen_ --output Allen_pseudobulk_profiles.csv
"""

from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import issparse

DEFAULT_GENE_NAME_COL = "feature_name"
DEFAULT_MIN_CELLS = 15
DEFAULT_CHUNK_SIZE = 50_000
# Non-negative-integer check tolerates float storage of whole numbers.
INTEGER_ATOL = 1e-8


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--h5ad", type=Path, action="append", required=True, dest="h5ads",
        help="Path to an .h5ad. Repeatable; with --label, pair one --label per --h5ad.",
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--cell-type-col",
        help="obs column holding per-cell type labels (single --h5ad).",
    )
    mode.add_argument(
        "--label", action="append", dest="labels",
        help="Whole-file label; repeat once per --h5ad, in the same order.",
    )
    p.add_argument("--gene-name-col", default=DEFAULT_GENE_NAME_COL,
                   help=f"var column with gene symbols (default {DEFAULT_GENE_NAME_COL}; "
                        "falls back to var_names).")
    p.add_argument("--prefix", default="",
                   help="Prepended to every cell-type name (e.g. 'HRCA_') to avoid "
                        "collisions when merging datasets.")
    p.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS,
                   help=f"Drop cell types with fewer than this many cells "
                        f"(default {DEFAULT_MIN_CELLS}).")
    p.add_argument("--normalize", choices=("mean", "sum", "proportion"), default="mean",
                   help="Per-cell-type column scaling: 'mean' = per-cell mean counts "
                        "(default, matches getRNAprofiles); 'sum' = raw count totals; "
                        "'proportion' = column sums to 1. The R step rescales anyway, so "
                        "this only sets the pre-merge scale.")
    p.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
                   help=f"Cells read per chunk in backed mode (default {DEFAULT_CHUNK_SIZE}).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output genes x cell-types CSV.")
    args = p.parse_args()

    if args.labels is not None and len(args.labels) != len(args.h5ads):
        p.error(f"got {len(args.h5ads)} --h5ad but {len(args.labels)} --label; "
                "pass exactly one --label per --h5ad.")
    if args.cell_type_col is not None and len(args.h5ads) != 1:
        p.error("--cell-type-col takes a single --h5ad; for multiple files use --label.")
    return args


def gene_symbols(adata: anndata.AnnData, gene_name_col: str) -> list[str]:
    """Gene symbols from var[gene_name_col], falling back to var_names."""
    if gene_name_col in adata.var.columns:
        return adata.var[gene_name_col].astype(str).tolist()
    print(f"  note: var['{gene_name_col}'] absent; using var_names as gene symbols")
    return adata.var_names.astype(str).tolist()


def raw_counts_slot(adata: anndata.AnnData) -> str:
    """Pick the slot holding raw counts: layers['raw'|'counts'] > .raw > X."""
    if adata.layers is not None:
        if "raw" in adata.layers:
            return "raw"
        if "counts" in adata.layers:
            return "counts"
    if adata.raw is not None:
        return "adata.raw"
    return "X"


def read_chunk(adata: anndata.AnnData, sl: slice, slot: str) -> np.ndarray:
    """Dense cells x genes chunk of raw counts from the chosen slot."""
    if slot in ("raw", "counts"):
        chunk = adata.layers[slot][sl]
    elif slot == "adata.raw":
        chunk = adata.raw.X[sl]
    else:
        chunk = adata.X[sl]
    if issparse(chunk):
        chunk = chunk.toarray()
    return np.asarray(chunk, dtype=np.float64)


def validate_counts(sample: np.ndarray, slot: str) -> None:
    """Abort if the chosen slot doesn't look like non-negative integer counts."""
    nz = sample[sample != 0]
    if nz.size == 0:
        return
    if np.any(nz < 0):
        sys.exit(f"ERROR: {slot} has negative values — looks log-transformed, not raw "
                 f"counts. Point --gene-name-col/data at a raw layer.")
    if not np.allclose(nz, np.rint(nz), atol=INTEGER_ATOL):
        sys.exit(f"ERROR: {slot} has non-integer values — looks normalized, not raw "
                 f"counts. Point at a raw layer.")


def pseudobulk_one(h5ad: Path, cell_type_col: str | None, whole_file_label: str | None,
                   gene_name_col: str, chunk_size: int) -> pd.DataFrame:
    """Sum raw counts per cell type for one H5AD. Returns genes x types (raw sums)."""
    print(f"Reading {h5ad} (backed)")
    adata = anndata.read_h5ad(h5ad, backed="r")
    n_cells, n_genes = adata.shape
    print(f"  {n_cells:,} cells x {n_genes:,} genes; obs columns: "
          f"{list(adata.obs.columns)}")

    if whole_file_label is not None:
        labels = np.full(n_cells, whole_file_label, dtype=object)
    else:
        if cell_type_col not in adata.obs.columns:
            candidates = [c for c in adata.obs.columns
                          if 2 < adata.obs[c].nunique() < 200]
            sys.exit(f"ERROR: obs['{cell_type_col}'] not found. Candidate label "
                     f"columns: {candidates}")
        labels = adata.obs[cell_type_col].astype(str).to_numpy()

    slot = raw_counts_slot(adata)
    print(f"  raw counts from: {slot}")
    validate_counts(read_chunk(adata, slice(0, min(100, n_cells)), slot), slot)

    types = sorted(set(labels))
    sums = {ct: np.zeros(n_genes, dtype=np.float64) for ct in types}
    n_per_type = {ct: 0 for ct in types}
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = read_chunk(adata, slice(start, end), slot)
        chunk_labels = labels[start:end]
        for ct in types:
            mask = chunk_labels == ct
            if mask.any():
                sums[ct] += chunk[mask].sum(axis=0)
                n_per_type[ct] += int(mask.sum())
        print(f"  {end:,}/{n_cells:,} cells", end="\r", flush=True)
        del chunk
        gc.collect()
    print()

    if hasattr(adata, "file") and adata.file is not None:
        adata.file.close()

    genes = gene_symbols(adata, gene_name_col)
    profiles = pd.DataFrame(sums, index=pd.Index(genes, name="gene"))
    # Collapse duplicate HUGO symbols (several Ensembl IDs -> one symbol).
    profiles = profiles.groupby(profiles.index).sum()
    profiles.attrs["n_per_type"] = n_per_type
    return profiles


def main() -> None:
    args = parse_args()

    frames: list[pd.DataFrame] = []
    n_per_type: dict[str, int] = {}
    for i, h5ad in enumerate(args.h5ads):
        label = args.labels[i] if args.labels is not None else None
        df = pseudobulk_one(h5ad, args.cell_type_col, label,
                            args.gene_name_col, args.chunk_size)
        n_per_type.update(df.attrs.get("n_per_type", {}))
        frames.append(df)

    # Merge the per-file sums on shared genes (raw sums are additive-safe here
    # because a given cell type lives in exactly one file in --label mode; in
    # --cell-type-col mode there is a single file).
    if len(frames) == 1:
        pseudobulk = frames[0]
    else:
        shared = sorted(set.intersection(*(set(f.index) for f in frames)))
        print(f"Shared genes across {len(frames)} files: {len(shared):,}")
        pseudobulk = pd.concat([f.loc[shared] for f in frames], axis=1)

    dropped = [ct for ct, n in n_per_type.items() if n < args.min_cells]
    if dropped:
        print(f"Dropping {len(dropped)} cell types with < {args.min_cells} cells: "
              f"{dropped}")
        pseudobulk = pseudobulk.drop(columns=dropped)

    # Column scaling. getRNAprofiles yields per-cell mean counts, so 'mean' keeps
    # pseudobulk sources comparable to getRNAprofiles sources pre-rescale.
    if args.normalize == "mean":
        for ct in pseudobulk.columns:
            n = n_per_type.get(ct, 0)
            if n > 0:
                pseudobulk[ct] = pseudobulk[ct] / n
    elif args.normalize == "proportion":
        col_sums = pseudobulk.sum(axis=0)
        pseudobulk = pseudobulk.div(col_sums.replace(0, np.nan), axis=1).fillna(0.0)
    # 'sum' leaves raw totals as-is.

    if args.prefix:
        pseudobulk.columns = [f"{args.prefix}{c}" for c in pseudobulk.columns]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pseudobulk.to_csv(args.output)
    print(f"Wrote {args.output} ({pseudobulk.shape[0]:,} genes x "
          f"{pseudobulk.shape[1]} cell types)")
    for ct in pseudobulk.columns:
        raw_ct = ct[len(args.prefix):] if args.prefix else ct
        print(f"  {ct:<45} {n_per_type.get(raw_ct, 0):>10,} cells")


if __name__ == "__main__":
    main()
