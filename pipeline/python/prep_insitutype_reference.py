#!/usr/bin/env python3
"""Restrict the CZI GBmap reference profiles to the CosMx panel genes.

One-time prep, run locally on the full GBmap export. InSituType's semi-supervised
typing (stage 4b) needs a genes x cell-types reference profile matrix on the same gene
set as the CosMx data. The full GBmap "gene by annotation" CSV is ~28k genes (a mix of
HGNC symbols and Ensembl IDs) x dozens of cell types and is too big to commit; this
intersects it with the ~6k CosMx panel and writes a small, committable matrix.

Reads:
  --gbmap-csv    GBmap "<...>_gene_by_annotation_level_3_raw_avg_expression.csv":
                 a genes x cell-types matrix of raw average expression. First (unnamed)
                 column is the gene name index; remaining columns are cell types.
  --panel-h5ad   any AnnData carrying the CosMx panel in var (a per-slide stage-1
                 .h5ad or the stage-3a combined_qc.h5ad). Gene probes are taken from
                 var['probe_type'] == 'gene' (var_names are HGNC symbols).

Writes:
  --output       genes x cell-types CSV restricted to the panel genes present in the
                 reference, ready for InSituType reference_profiles. Genes are the row
                 index; cell types are columns (Extended GBmap level-3 annotation).

Overlap is the symbol intersection; GBmap's Ensembl-ID rows simply don't match the
symbol-named panel and drop out. A low overlap (<50% of the panel) is flagged loudly,
since that would point at a symbol-vs-Ensembl naming mismatch needing a gene-id map.

Usage:
    uv run python pipeline/python/prep_insitutype_reference.py \\
        --gbmap-csv ~/keene-lab/GBM/GBmap/extended_GBmap_results/\\
extended_GBmap_gene_by_annotation_level_3_raw_avg_expression.csv \\
        --panel-h5ad combined_qc.h5ad \\
        --output pipeline/reference/gbmap_extended_level3_panel.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Flag a likely gene-naming mismatch rather than silently typing on a thin reference.
MIN_PANEL_OVERLAP_FRAC = 0.5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gbmap-csv", type=Path, required=True,
                   help="GBmap gene-by-annotation raw-average-expression CSV "
                        "(genes x cell types; first column is the gene index).")
    p.add_argument("--panel-h5ad", type=Path, required=True,
                   help="AnnData with the CosMx panel in var (per-slide or combined_qc).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output panel-restricted reference CSV (genes x cell types).")
    return p.parse_args()


def panel_genes(panel_h5ad: Path) -> pd.Index:
    """Gene-probe names from an AnnData's var (probe_type == 'gene')."""
    adata = ad.read_h5ad(panel_h5ad, backed="r")
    if "probe_type" not in adata.var:
        print(f"ERROR: var['probe_type'] missing in {panel_h5ad}; "
              f"was it built by stage 1?", file=sys.stderr)
        sys.exit(1)
    genes = adata.var_names[(adata.var["probe_type"] == "gene").to_numpy()]
    return pd.Index(genes.astype(str).str.strip(), name="gene").unique()


def main() -> None:
    args = parse_args()

    print(f"Reading GBmap reference {args.gbmap_csv}")
    ref = pd.read_csv(args.gbmap_csv, index_col=0)
    ref.index = ref.index.astype(str).str.strip()
    ref = ref[~ref.index.duplicated(keep="first")]
    print(f"  reference: {ref.shape[0]:,} genes x {ref.shape[1]} cell types")
    print(f"  cell types: {list(ref.columns)}")

    panel = panel_genes(args.panel_h5ad)
    print(f"CosMx panel: {len(panel):,} gene probes from {args.panel_h5ad}")

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
