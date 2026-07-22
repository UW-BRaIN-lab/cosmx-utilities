#!/usr/bin/env python3
"""Stage 5a: build per-tissue-section inputs for InSituCNV from the typed cohort.

InSituCNV (infercnvpy + spatial neighbor smoothing) needs, per cell: raw gene counts,
2D spatial coordinates, a cell-type label, and gene genomic positions. The typed cohort
AnnData (stage4_insitutree/cosmx_typed.h5ad) has the counts and labels but NOT an
obsm['spatial'] — the CosMx per-cell centroids live in obs (CenterX/Y_global_px). This
driver assembles a clean CNV input and SPLITS it into one file per tissue section.

WHY PER TISSUE SECTION: each CosMx slide co-mounts TWO donor tissues, and the global-px
centroids are per-slide (they repeat across slides). If the spatial neighbor graph were
built on the whole cohort — or even per slide — smoothing would pool counts, and
subtract the CNV reference, ACROSS patients. The physically meaningful unit is one donor
tissue on one slide, keyed here as ``tissue_section = "<slide_id>__<Case>"``. run_insitucnv.py
processes one section at a time so every neighbor graph is self-contained.

Gene positions are annotated HERE (once) from a full-genome table so every section shares
an identical gene set -> identical infercnv windows -> per-section X_cnv matrices that
concatenate cleanly in the compare step. We keep autosomes only (GBM's chr7 gain / chr10
loss / chr9p loss are all autosomal; dropping X/Y avoids a sex-mismatch baseline artifact).

Reads:
  --typed-h5ad      cosmx_typed.h5ad (obs['cell_type','Case','Region','slide_id'], X raw counts,
                    var['probe_type']).
  --gene-table      Ensembl BioMart table (Gene name, Chromosome/scaffold name, Gene start
                    (bp), Gene end (bp)); default $INSITUCNV_GENE_TABLE (baked into the SIF).
  --reference-file  the diploid reference cell_type list (only to count reference cells per
                    section, so sections that cannot be baselined are skipped).
Writes (--output-dir):
  sections/<tissue_section>.h5ad   one gene-only, position-annotated, raw-count AnnData per
                                   section (X raw, obsm['spatial'], obs cell_type/Region/...).
  sections_manifest.csv            per-section n_cells / n_reference_cells / kept.
  sections.txt                     newline-delimited KEPT section ids (drives the Slurm array).

Usage:
    python pipeline/python/prep_insitucnv_input.py \\
        --typed-h5ad cosmx_typed.h5ad \\
        --reference-file pipeline/reference/insitucnv_reference_types.txt \\
        --output-dir out --min-cells 500 --min-reference-cells 50
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Per-cell centroid columns to build obsm['spatial'] from, in priority order. Global-px
# (stitched across a slide's FOVs) is what CosMx exports; the *_local_px are per-FOV and
# would collide, so they are intentionally NOT candidates.
SPATIAL_CANDIDATES = [
    ("CenterX_global_px", "CenterY_global_px"),
    ("x_global_px", "y_global_px"),
    ("CenterX_global_mm", "CenterY_global_mm"),
    ("x_slide_mm", "y_slide_mm"),
    ("CenterX_global", "CenterY_global"),
]
AUTOSOMES = [str(i) for i in range(1, 23)]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--gene-table", type=Path,
                   default=Path(os.environ.get("INSITUCNV_GENE_TABLE", "")),
                   help="Ensembl BioMart gene-position table (default $INSITUCNV_GENE_TABLE).")
    p.add_argument("--reference-file", type=Path, required=True,
                   help="Diploid reference cell_type list (one per line; '#' comments ok).")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--donor-key", default="Case", help="obs column identifying the donor.")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--min-cells", type=int, default=500,
                   help="Skip sections with fewer cells than this.")
    p.add_argument("--min-reference-cells", type=int, default=50,
                   help="Skip sections with fewer diploid-reference cells than this "
                        "(infercnv cannot baseline them).")
    return p.parse_args()


def read_reference_types(path: Path) -> list[str]:
    out = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def pick_spatial_cols(cols) -> tuple[str, str]:
    for x, y in SPATIAL_CANDIDATES:
        if x in cols and y in cols:
            return x, y
    sys.exit(f"ERROR: no spatial centroid columns found (tried {SPATIAL_CANDIDATES}).\n"
             f"       obs columns present: {list(cols)}")


def annotate_gene_positions(var: pd.DataFrame, gene_table: Path) -> pd.DataFrame:
    """Return var with chromosome/start/end for panel genes on autosomes 1-22.

    Matches on HGNC gene symbol (the panel var index == BioMart 'Gene name'). Genes with
    no autosomal position are dropped; duplicate symbols keep the first autosomal hit.
    """
    if not gene_table or not gene_table.exists():
        sys.exit(f"ERROR: gene table not found: '{gene_table}' (set --gene-table / "
                 f"$INSITUCNV_GENE_TABLE; it is baked into insitucnv.sif).")
    bm = pd.read_csv(gene_table)
    need = ["Gene name", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"]
    missing = [c for c in need if c not in bm.columns]
    if missing:
        sys.exit(f"ERROR: gene table missing columns {missing}; has {list(bm.columns)}")

    bm = bm.rename(columns={"Gene name": "gene", "Chromosome/scaffold name": "chrom",
                            "Gene start (bp)": "start", "Gene end (bp)": "end"})
    bm["chrom"] = bm["chrom"].astype(str)
    bm = bm[bm["chrom"].isin(AUTOSOMES)]
    bm = bm.drop_duplicates(subset="gene", keep="first").set_index("gene")

    genome = bm.reindex(var.index)
    positioned = genome["chrom"].notna()
    out = var.copy()
    out["chromosome"] = ("chr" + genome["chrom"].astype("string")).to_numpy()
    out["start"] = pd.to_numeric(genome["start"], errors="coerce").to_numpy()
    out["end"] = pd.to_numeric(genome["end"], errors="coerce").to_numpy()
    out["_positioned"] = positioned.to_numpy()
    return out


def safe_id(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "-", str(text)).strip("-")


def main() -> None:
    args = parse_args()
    sections_dir = args.output_dir / "sections"
    sections_dir.mkdir(parents=True, exist_ok=True)

    reference_types = set(read_reference_types(args.reference_file))
    print(f"Reference cell types (diploid baseline): {len(reference_types)}")

    print(f"Reading {args.typed_h5ad}")
    adata = ad.read_h5ad(args.typed_h5ad)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} probes; obs cols: {list(adata.obs.columns)}")

    for key in (args.celltype_key, args.donor_key):
        if key not in adata.obs:
            sys.exit(f"ERROR: obs['{key}'] missing (needed for typing / section key).")

    # --- gene-only + genomic positions (autosomes) ---------------------------------
    if "probe_type" in adata.var:
        adata = adata[:, (adata.var["probe_type"] == "gene").to_numpy()].copy()
    var = annotate_gene_positions(adata.var, args.gene_table)
    keep_genes = var["_positioned"].to_numpy()
    print(f"Gene positions: {int(keep_genes.sum()):,} of {adata.n_vars:,} gene-probes "
          f"mapped to autosomes 1-22 ({keep_genes.mean():.0%} of the panel).")
    adata = adata[:, keep_genes].copy()
    adata.var["chromosome"] = var.loc[keep_genes, "chromosome"].to_numpy()
    adata.var["start"] = var.loc[keep_genes, "start"].to_numpy().astype(np.int64)
    adata.var["end"] = var.loc[keep_genes, "end"].to_numpy().astype(np.int64)

    # --- spatial coords + section key ----------------------------------------------
    xcol, ycol = pick_spatial_cols(adata.obs.columns)
    print(f"Building obsm['spatial'] from ('{xcol}', '{ycol}').")
    adata.obsm["spatial"] = np.column_stack([
        pd.to_numeric(adata.obs[xcol], errors="coerce").to_numpy(),
        pd.to_numeric(adata.obs[ycol], errors="coerce").to_numpy(),
    ]).astype(np.float32)

    if "slide_id" in adata.obs:
        slide = adata.obs["slide_id"].astype(str)
    else:
        slide = adata.obs.index.to_series().str.split("_F", n=1).str[0]
        print("  note: no obs['slide_id']; derived slide from the cell-id prefix.")
    donor = adata.obs[args.donor_key].astype(str)
    adata.obs["tissue_section"] = (slide.to_numpy() + "__" + donor.to_numpy())

    # keep only the obs columns downstream needs (keeps section files small)
    keep_obs = [c for c in [args.celltype_key, args.region_key, args.donor_key,
                            "slide_id", "tissue_section"] if c in adata.obs]
    adata.obs = adata.obs[keep_obs].copy()

    # --- split per section ----------------------------------------------------------
    is_ref = adata.obs[args.celltype_key].astype(str).isin(reference_types).to_numpy()
    # positional row indices per section, computed in a single pass
    pos = pd.Series(np.arange(adata.n_obs))
    groups = pos.groupby(adata.obs["tissue_section"].to_numpy(), observed=True).groups
    rows, kept = [], []
    for section in sorted(groups):
        rowpos = np.asarray(groups[section], dtype=int)
        n_cells = int(rowpos.size)
        n_ref = int(is_ref[rowpos].sum())
        ok = (n_cells >= args.min_cells) and (n_ref >= args.min_reference_cells)
        rows.append({"tissue_section": section, "n_cells": n_cells,
                     "n_reference_cells": n_ref, "kept": ok})
        if not ok:
            print(f"  SKIP {section}: n_cells={n_cells:,}, n_ref={n_ref} "
                  f"(need >= {args.min_cells} / {args.min_reference_cells})")
            continue
        fid = safe_id(section)
        adata[rowpos].copy().write_h5ad(sections_dir / f"{fid}.h5ad", compression="gzip")
        kept.append(fid)

    manifest = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    manifest.to_csv(args.output_dir / "sections_manifest.csv", index=False)
    (args.output_dir / "sections.txt").write_text("\n".join(kept) + ("\n" if kept else ""))

    n_total = len(rows)
    print(f"\nWrote {len(kept)} / {n_total} sections to {sections_dir} "
          f"({n_total - len(kept)} skipped).")
    print(f"  manifest: {args.output_dir / 'sections_manifest.csv'}")
    print(f"  array list: {args.output_dir / 'sections.txt'} "
          f"(SBATCH --array=0-{max(0, len(kept) - 1)})")


if __name__ == "__main__":
    main()
