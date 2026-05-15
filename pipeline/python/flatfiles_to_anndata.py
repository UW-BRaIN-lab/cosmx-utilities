#!/usr/bin/env python3
"""Convert one slide's CosMx flat files into an AnnData (.h5ad) file.

Reads the per-slide expression matrix, per-cell metadata, FOV positions, and
(optionally) a clinical-annotation sidecar (from extract_clinical_annotations.R)
plus a manifest row (from build_manifest.py), and writes a single per-slide
.h5ad with:
  - X: raw integer counts (sparse, cells x panel probes)
  - obs: per-cell metadata + clinical annotations + slide-level fields
  - var: probe name + probe_type (gene / negprobe / falsecode)
  - uns: slide_id, run_uuid, instrument_id, FOV positions, ...

Stage 1 of the analysis pipeline. Run once per slide.

Usage:
    uv run python pipeline/python/flatfiles_to_anndata.py \\
        --flatfiles-dir /path/to/<slide>_flat \\
        --clinical /path/to/<slide>_clinical.csv \\
        --manifest pipeline/manifest.csv \\
        --slide-id 7134A77439A6 \\
        --output /path/to/<slide>.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


EXPRMAT_NONGENE_COLS = ("fov", "cell_ID")

NEGPROBE_PREFIXES = ("Negative", "NegPrb", "Neg")
FALSECODE_PREFIXES = ("SystemControl", "FalseCode", "Falsecode")

MANIFEST_FIELDS = (
    "run_uuid",
    "run_date",
    "run_time",
    "instrument_id",
    "slot",
    "export_batch",
)

CLINICAL_FIELDS = ("case", "block", "region", "cell_segmentation_set")


def classify_probe(name: str) -> str:
    if name.startswith(NEGPROBE_PREFIXES):
        return "negprobe"
    if name.startswith(FALSECODE_PREFIXES):
        return "falsecode"
    return "gene"


def read_expr_mat(path: Path) -> tuple[pd.DataFrame, sp.csr_matrix, list[str]]:
    """Read a CosMx exprMat CSV(.gz). Returns (cell_index_df, counts, gene_names)."""
    df = pd.read_csv(path)
    nongene_cols = [c for c in EXPRMAT_NONGENE_COLS if c in df.columns]
    cell_index = df[nongene_cols].copy()
    cell_index["cell_ID"] = cell_index["cell_ID"].astype(str)
    gene_cols = [c for c in df.columns if c not in nongene_cols]
    counts = sp.csr_matrix(df[gene_cols].to_numpy(dtype=np.int32, copy=False))
    return cell_index, counts, gene_cols


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flatfiles-dir", type=Path, required=True,
                   help="Directory containing <slide>_exprMat_file.csv.gz etc.")
    p.add_argument("--slide-id", required=True,
                   help="Slide identifier matching the flat-file filename prefix")
    p.add_argument("--clinical", type=Path,
                   help="Per-slide clinical annotation CSV (optional)")
    p.add_argument("--manifest", type=Path,
                   help="Pipeline manifest CSV (optional)")
    p.add_argument("--output", type=Path, required=True, help="Path to write .h5ad")
    args = p.parse_args()

    slide_id: str = args.slide_id
    fdir: Path = args.flatfiles_dir
    expr_path = fdir / f"{slide_id}_exprMat_file.csv.gz"
    meta_path = fdir / f"{slide_id}_metadata_file.csv.gz"
    fov_path = fdir / f"{slide_id}_fov_positions_file.csv.gz"

    for path in (expr_path, meta_path, fov_path):
        if not path.exists():
            print(f"ERROR: missing flat file: {path}", file=sys.stderr)
            sys.exit(1)

    print(f"Reading {expr_path.name}")
    cell_index, counts, gene_names = read_expr_mat(expr_path)

    print(f"Reading {meta_path.name}")
    metadata = pd.read_csv(meta_path)
    metadata["cell_ID"] = metadata["cell_ID"].astype(str)

    # Align metadata to expression-matrix row order by (fov, cell_ID).
    cell_index["_row_idx"] = np.arange(len(cell_index))
    merged = cell_index.merge(
        metadata, on=["fov", "cell_ID"], how="left", indicator=True,
    )
    n_unmatched = int((merged["_merge"] != "both").sum())
    if n_unmatched:
        print(
            f"ERROR: {n_unmatched} cells in exprMat have no matching metadata row",
            file=sys.stderr,
        )
        sys.exit(1)
    merged = (merged.sort_values("_row_idx")
                    .drop(columns=["_row_idx", "_merge"]))

    if args.clinical and args.clinical.exists():
        print(f"Reading clinical annotations from {args.clinical}")
        clinical = pd.read_csv(args.clinical)
        clinical["cell_ID"] = clinical["cell_ID"].astype(str)
        keep_cols = ["cell_ID"] + [c for c in CLINICAL_FIELDS if c in clinical.columns]
        merged = merged.merge(clinical[keep_cols], on="cell_ID", how="left")
    elif args.clinical:
        print(f"WARN: clinical file not found: {args.clinical}", file=sys.stderr)

    merged["slide_id"] = slide_id
    uns: dict = {"slide_id": slide_id}

    if args.manifest and args.manifest.exists():
        mdf = pd.read_csv(args.manifest, dtype={"slide_id": str})
        row = mdf[mdf["slide_id"] == slide_id]
        if len(row) != 1:
            print(
                f"WARN: manifest has {len(row)} rows for slide_id={slide_id}; "
                f"skipping manifest join",
                file=sys.stderr,
            )
        else:
            for col in MANIFEST_FIELDS:
                if col in row.columns:
                    val = row.iloc[0][col]
                    merged[col] = val
                    uns[col] = val if pd.notna(val) else None
    elif args.manifest:
        print(f"WARN: manifest file not found: {args.manifest}", file=sys.stderr)

    obs = merged.set_index("cell_ID")
    obs.index.name = "cell_ID"

    var = pd.DataFrame(
        {"probe_type": [classify_probe(g) for g in gene_names]},
        index=pd.Index(gene_names, name="feature_name"),
    )

    # FOV positions kept in uns to avoid bloating obs; ndarray-of-records is the
    # simplest cross-language format.
    print(f"Reading {fov_path.name}")
    fov_df = pd.read_csv(fov_path)
    uns["fov_positions"] = fov_df.to_records(index=False)

    adata = ad.AnnData(X=counts, obs=obs, var=var, uns=uns)
    n_gene = int((var["probe_type"] == "gene").sum())
    n_neg = int((var["probe_type"] == "negprobe").sum())
    n_fc = int((var["probe_type"] == "falsecode").sum())
    print(
        f"Built AnnData: {adata.shape[0]} cells x {adata.shape[1]} probes "
        f"({n_gene} genes, {n_neg} negprobes, {n_fc} falsecodes)"
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    print("Done.")


if __name__ == "__main__":
    main()
