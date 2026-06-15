#!/usr/bin/env python3
"""Relabel InSituType de-novo clusters with biological annotations.

InSituType's de-novo clusters come out as single letters (a, b, c, …). This rewrites the
obs cell_type column so each letter becomes a readable label, KEEPING the letter prefix
so it stays clear the cluster was de-novo discovered rather than reference-named — e.g.
`a` -> `a - MES/AC-like tumor`. Named GBmap reference types are left unchanged.

The per-run letter->annotation tables live in pipeline/reference/denovo_annotations/,
one CSV per run (columns: denovo_label, annotation, top_markers, quality). Each was
assigned from that run's de-novo marker heatmap (marker_pseudobulk.py --clusters a,…,l).

Memory: loads the full typed AnnData (~3.6GB at cohort scale) then writes a copy — run
on a node with enough RAM (e.g. the rapids-singlecell container in a Slurm/salloc job),
not on a login node.

Usage:
    uv run python pipeline/python/annotate_denovo.py \\
        --typed-h5ad cosmx_typed.h5ad \\
        --mapping pipeline/reference/denovo_annotations/stage4_extl3_rescale.csv \\
        --output cosmx_typed_annotated.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True,
                   help="cosmx_typed.h5ad from stage 4c (obs has the cell_type column).")
    p.add_argument("--mapping", type=Path, required=True,
                   help="denovo_annotations CSV with denovo_label + annotation columns.")
    p.add_argument("--output", type=Path, required=True,
                   help="Output annotated .h5ad.")
    p.add_argument("--cell-type-key", default="cell_type",
                   help="obs column holding the InSituType labels (default: cell_type).")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    m = pd.read_csv(args.mapping, dtype=str)
    for col in ("denovo_label", "annotation"):
        if col not in m.columns:
            print(f"ERROR: mapping {args.mapping} missing '{col}' column. "
                  f"Have: {list(m.columns)}", file=sys.stderr)
            sys.exit(1)
    mapping = dict(zip(m["denovo_label"].str.strip(), m["annotation"].str.strip()))

    print(f"Reading {args.typed_h5ad}")
    adata = ad.read_h5ad(args.typed_h5ad)
    key = args.cell_type_key
    if key not in adata.obs:
        print(f"ERROR: obs missing '{key}'. Have: {list(adata.obs.columns)}",
              file=sys.stderr)
        sys.exit(1)

    cats = adata.obs[key].astype("category").cat.categories
    rename = {c: mapping[c] for c in cats if c in mapping}
    absent = sorted(set(mapping) - set(cats))
    if absent:
        print(f"WARN: {len(absent)} mapping label(s) not present as '{key}' categories "
              f"(skipped): {absent}", file=sys.stderr)
    if not rename:
        print(f"ERROR: none of the mapping's de-novo labels matched a '{key}' category; "
              f"wrong run table?", file=sys.stderr)
        sys.exit(1)
    print(f"Relabeling {len(rename)} de-novo clusters; "
          f"{len(cats) - len(rename)} named/reference types left unchanged")
    for k in sorted(rename):
        print(f"  {k} -> {rename[k]}")

    adata.obs[key] = adata.obs[key].astype("category").cat.rename_categories(rename)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Done. {adata.n_obs:,} cells; '{key}' de-novo letters annotated.")


if __name__ == "__main__":
    main()
