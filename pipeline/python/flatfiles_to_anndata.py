#!/usr/bin/env python3
"""Convert one slide's CosMx flat files into an AnnData (.h5ad) file.

Reads the per-slide expression matrix, per-cell metadata, FOV positions, and
(optionally) a manifest row (from build_manifest.py), and writes a single
per-slide .h5ad with:
  - X: raw integer counts (sparse, cells x panel probes)
  - obs: per-cell metadata, slide-prefixed FOV, and manifest fields
  - var: probe name + probe_type (gene / negprobe / falsecode)
  - uns: slide_id, run_uuid, instrument_id, FOV positions, ...

The per-FOV `cell_ID` and per-slide `fov` columns from the flat files are not
unique across slides. We follow the Bruker CosMx Scratch Space vignette and
replace them with slide-prefixed forms so per-slide AnnDatas can be safely
concatenated for cross-slide work:
  - obs index `cell` = "<slide_id>_F<fov>_C<cell_ID>"
  - obs column `FOV` = "<slide_id>_F<fov>"
The raw integer `fov` is preserved for within-slide filtering; the raw
`cell_ID` is dropped (it lives in the obs index name).

Stage 1 of the analysis pipeline. Run once per slide.

Usage:
    uv run python pipeline/python/flatfiles_to_anndata.py \\
        --flatfiles-dir /path/to/<slide>_flat \\
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


def apply_fov_annotations(obs: pd.DataFrame, table: Path, slide_id: str) -> pd.DataFrame:
    """Replace Region from an authoritative per-FOV annotation table.

    The flat-file `Region` is an AtoMx annotation that can be stale — on uwa7761eyes it
    labels 67 FOVs as Retina when only 3 are, the other 57 being ciliary body and cornea.
    Where a curated table exists it wins, and the two extra columns ride along:
      Region_mixed_adjacent  the FOV straddles a tissue boundary
      Region_excluded        the SME ruled it out of analysis

    Joins on `fov` only — cell ids differ between exports but the FOV grid does not. Every
    substitution is counted and reported; an FOV present in the data but absent from the
    table keeps its original Region and is reported rather than silently blanked, since a
    partial table must not quietly erase annotation.
    """
    ann = pd.read_csv(table, dtype={"slide_id": str})
    ann = ann[ann["slide_id"].astype(str) == slide_id]
    if ann.empty:
        print(f"WARN: no rows for slide_id={slide_id} in {table}; "
              f"leaving Region as exported", file=sys.stderr)
        return obs
    if ann["fov"].duplicated().any():
        dupes = sorted(ann.loc[ann["fov"].duplicated(), "fov"].unique())
        print(f"ERROR: {table} has duplicate fov rows for {slide_id}: {dupes}",
              file=sys.stderr)
        sys.exit(1)

    by_fov = ann.set_index(ann["fov"].astype(int))
    fov_int = obs["fov"].astype(int)
    mapped = fov_int.map(by_fov["region"])
    n_missing = int(mapped.isna().sum())
    if n_missing:
        absent = sorted(set(fov_int[mapped.isna()]))
        print(f"WARN: {n_missing:,} cells in {len(absent)} FOV(s) absent from the "
              f"annotation table, keeping exported Region: {absent[:20]}", file=sys.stderr)

    if "Region" in obs:
        changed = int((mapped.notna() & (mapped != obs["Region"].astype(str))).sum())
        print(f"FOV annotations: Region reassigned for {changed:,} of {len(obs):,} cells")
    obs["Region"] = mapped.fillna(obs["Region"] if "Region" in obs else "")
    obs["Region_mixed_adjacent"] = (
        fov_int.map(by_fov["mixed_adjacent"]).fillna(0).astype(int))
    obs["Region_excluded"] = fov_int.map(by_fov["exclude"]).fillna(0).astype(int)
    n_excl = int(obs["Region_excluded"].sum())
    if n_excl:
        print(f"FOV annotations: {n_excl:,} cells flagged Region_excluded "
              f"(drop at stage 3a with --exclude-regions)")
    print("Region composition after annotation:")
    for region, count in obs["Region"].value_counts().items():
        print(f"  {region:26s} {count:7,d}")
    return obs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flatfiles-dir", type=Path, required=True,
                   help="Directory containing <slide>_exprMat_file.csv.gz etc.")
    p.add_argument("--slide-id", required=True,
                   help="Slide identifier matching the flat-file filename prefix")
    p.add_argument("--manifest", type=Path,
                   help="Pipeline manifest CSV (optional)")
    p.add_argument("--fov-annotations", type=Path,
                   help="Authoritative per-FOV tissue table (slide_id,fov,region,...). "
                        "Overrides the flat-file Region for this slide; see "
                        "pipeline/reference/FOV_ANNOTATIONS.md.")
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

    # Slide-prefixed cell + FOV identifiers (unique across slides).
    merged["cell"] = (
        slide_id + "_F" + merged["fov"].astype(str) + "_C" + merged["cell_ID"].astype(str)
    )
    merged["FOV"] = slide_id + "_F" + merged["fov"].astype(str)
    # Keep the flat file's own `cell_id` under an unambiguous name: it is the key Napari
    # and every other AtoMx-derived artifact joins on, and it is NOT derivable from our
    # slide-prefixed index (nor the reverse). Dropping it forced downstream tools to
    # reconstruct it by re-joining the flat files on (fov, cell_ID), which is fragile and
    # needs the flat files on hand. Costs one string column per cell.
    if "cell_id" in merged:
        merged = merged.rename(columns={"cell_id": "flatfile_cell_id"})
    # The per-FOV `cell_ID` is redundant once it is in the index; `fov` (integer,
    # per-slide) is preserved for within-slide work.
    merged = merged.drop(columns=["cell_ID"], errors="ignore")

    merged["slide_id"] = slide_id
    uns: dict = {"slide_id": slide_id}

    if args.fov_annotations:
        merged = apply_fov_annotations(merged, args.fov_annotations, slide_id)

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

    obs = merged.set_index("cell")
    assert obs.index.is_unique, "obs index not unique after slide-prefix construction"

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
