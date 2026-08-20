#!/usr/bin/env python3
"""Restrict a genes x cell-types reference profile matrix to the CosMx panel genes.

One-time prep, run locally on a full reference export. InSituType's semi-supervised
typing (stage 4b) needs a genes x cell-types reference profile matrix on the same gene
set as the CosMx data. Atlas exports are typically ~21-28k genes and too big to commit;
this intersects them with the ~6k CosMx panel and writes a small, committable matrix.

Reference-agnostic: it only needs a genes x cell-types CSV, whatever built it (the CZI
GBmap export for GBM, pipeline/reference_builder/ for the ocular/brain atlas, ...).

Reads:
  --reference-csv  a genes x cell-types matrix of raw average expression. First
                   (unnamed) column is the gene name index; remaining columns are
                   cell types.
  the CosMx panel, as EITHER:
  --panel-h5ad     any AnnData carrying the panel in var (a per-slide stage-1 .h5ad or
                   the stage-3a combined_qc.h5ad). Gene probes are taken from
                   var['probe_type'] == 'gene' (var_names are HGNC symbols).
  --panel-genes    a plain text file, one gene-probe name per line. Use this when no
                   .h5ad exists yet (before stage 1) — the panel gene list can be read
                   straight off an exprMat header, dropping fov/cell_ID and the
                   Negative*/SystemControl/FalseCode control probes.

Writes:
  --output       genes x cell-types CSV restricted to the panel genes present in the
                 reference, ready for InSituType reference_profiles. Genes are the row
                 index; cell types are columns.

Overlap is the symbol intersection; reference rows named by Ensembl ID simply don't match
the symbol-named panel and drop out. A low overlap (<50% of the panel) is flagged loudly,
since that would point at a symbol-vs-Ensembl naming mismatch needing a gene-id map.

Usage:
    # GBM: CZI GBmap level 4, panel from a stage-3a AnnData
    uv run python pipeline/python/prep_insitutype_reference.py \\
        --reference-csv ~/keene-lab/GBM/GBmap/core_GBmap_results/\\
core_GBmap_gene_by_annotation_level_4_raw_avg_expression.csv \\
        --panel-h5ad combined_qc.h5ad \\
        --output pipeline/reference/gbmap_level4_panel.csv

    # Retina/brain: combined ocular atlas, panel from an exprMat-derived gene list
    uv run python pipeline/python/prep_insitutype_reference.py \\
        --reference-csv ~/keene-lab/retina-brain/\\
HRCA-Monavarfeshani-Allen-combined-profile.csv \\
        --panel-genes pipeline/reference/cosmx_6k_panel_genes.txt \\
        --output pipeline/reference/retina_combined_panel.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Flag a likely gene-naming mismatch rather than silently typing on a thin reference.
MIN_PANEL_OVERLAP_FRAC = 0.5

# Control / non-gene columns of a CosMx exprMat, so a --panel-genes list that still
# carries them is rejected instead of quietly skewing the overlap check. Kept in sync
# with the probe classifier in flatfiles_to_anndata.py.
CONTROL_PROBE_PREFIXES = ("Negative", "NegPrb", "Neg",
                          "SystemControl", "FalseCode", "Falsecode")
EXPRMAT_NONGENE_COLS = ("fov", "cell_ID")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference-csv", type=Path, required=True,
                   help="Reference gene-by-cell-type raw-average-expression CSV "
                        "(genes x cell types; first column is the gene index).")
    panel_src = p.add_mutually_exclusive_group(required=True)
    panel_src.add_argument("--panel-h5ad", type=Path,
                           help="AnnData with the CosMx panel in var "
                                "(per-slide or combined_qc).")
    panel_src.add_argument("--panel-genes", type=Path,
                           help="Text file of CosMx panel gene-probe names, one per "
                                "line (use before stage 1, when no .h5ad exists).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output panel-restricted reference CSV (genes x cell types).")
    return p.parse_args()


def panel_genes_from_h5ad(panel_h5ad: Path) -> pd.Index:
    """Gene-probe names from an AnnData's var (probe_type == 'gene')."""
    adata = ad.read_h5ad(panel_h5ad, backed="r")
    if "probe_type" not in adata.var:
        print(f"ERROR: var['probe_type'] missing in {panel_h5ad}; "
              f"was it built by stage 1?", file=sys.stderr)
        sys.exit(1)
    genes = adata.var_names[(adata.var["probe_type"] == "gene").to_numpy()]
    return pd.Index(genes.astype(str).str.strip(), name="gene").unique()


def panel_genes_from_txt(panel_txt: Path) -> pd.Index:
    """Gene-probe names from a one-per-line text file.

    Control probes are rejected rather than silently filtered: a list still carrying
    Negative*/SystemControl entries means the caller derived it from a raw exprMat
    header without dropping them, and every one would be a false 'absent from the
    reference' miss that inflates the overlap denominator.
    """
    names = [line.strip() for line in panel_txt.read_text().splitlines()]
    names = [n for n in names if n]
    if not names:
        print(f"ERROR: no gene names in {panel_txt}", file=sys.stderr)
        sys.exit(1)
    controls = [n for n in names if n.startswith(CONTROL_PROBE_PREFIXES)
                or n in EXPRMAT_NONGENE_COLS]
    if controls:
        print(f"ERROR: {len(controls)} control/non-gene probes in {panel_txt} "
              f"(e.g. {controls[:5]}); pass gene probes only.", file=sys.stderr)
        sys.exit(1)
    return pd.Index(names, name="gene").unique()


def panel_genes(args: argparse.Namespace) -> pd.Index:
    if args.panel_h5ad:
        return panel_genes_from_h5ad(args.panel_h5ad)
    return panel_genes_from_txt(args.panel_genes)


def main() -> None:
    args = parse_args()

    print(f"Reading reference {args.reference_csv}")
    ref = pd.read_csv(args.reference_csv, index_col=0)
    ref.index = ref.index.astype(str).str.strip()
    ref = ref[~ref.index.duplicated(keep="first")]
    print(f"  reference: {ref.shape[0]:,} genes x {ref.shape[1]} cell types")
    print(f"  cell types: {list(ref.columns)}")

    panel = panel_genes(args)
    panel_src = args.panel_h5ad or args.panel_genes
    print(f"CosMx panel: {len(panel):,} gene probes from {panel_src}")

    shared = panel.intersection(ref.index)
    frac = len(shared) / len(panel) if len(panel) else 0.0
    print(f"Overlap: {len(shared):,} / {len(panel):,} panel genes "
          f"in the reference ({frac:.1%})")
    if frac < MIN_PANEL_OVERLAP_FRAC:
        print(f"ERROR: only {frac:.1%} of the panel matched the reference — likely a "
              f"gene-naming mismatch (symbols vs Ensembl IDs). Aborting rather than "
              f"writing a thin reference.", file=sys.stderr)
        sys.exit(1)
    missing = sorted(set(panel) - set(shared))
    if missing:
        print(f"Note: {len(missing)} panel genes absent from the reference "
              f"(InSituType uses only the shared genes), e.g. {missing[:10]}")

    out = ref.loc[shared]
    out.index.name = "gene"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} cell types)")


if __name__ == "__main__":
    main()
