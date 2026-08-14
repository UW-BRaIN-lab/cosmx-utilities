#!/usr/bin/env python3
"""Turn an InSituType result's /profiles into marker_heatmap.R inputs (no per-cell data needed).

The InSituType rescale $profiles matrix IS the per-cell-type pseudobulk (mean expression per
gene per type), so a de-novo triage marker heatmap can be built straight from anchor_typing.h5 --
no h5ad / Region assembly (that's marker_pseudobulk.py's job for the by-region full-cohort view).

Per de-novo cluster we pick its top specificity markers (profile / panel-wide mean, the same
fold metric as inspect_denovo_profiles.py), take the union across clusters (each gene assigned to
the cluster where it is most enriched, so rows stay unique), and z-score each gene across the
displayed cluster columns. Emits exactly what marker_heatmap.R reads:

  <out>/marker_heatmap_zmatrix.csv     genes x de-novo types, z-scored (row index = gene)
  <out>/top_markers_per_cluster.csv    columns gene, cluster (row split / left annotation)

Columns carry NO " | Region" suffix, so marker_heatmap.R renders it in its no-region profile mode.
Pass the same denovo_annotations CSV as marker_heatmap.R's label_map to relabel letters -> names.

Usage:
    uv run --no-project --with h5py --with numpy --with pandas python \\
        pipeline/python/profile_marker_inputs.py \\
            --profiles-h5 anchor_typing.h5 --output-dir denovo_markers \\
            --top-n 8 [--include-named] [--min-size 100]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# De-novo labels = InSituType cluster_name_pool = 1-2 lowercase letters (see inspect_denovo_profiles.py).
DENOVO_RE = re.compile(r"^[a-z]{1,2}$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True,
                   help="InSituType result h5 with /profiles, /profile_genes, /profile_types.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=8, help="Specificity markers per de-novo cluster.")
    p.add_argument("--include-named", action="store_true",
                   help="Also show named-type columns (context for duplicate-of-normal calls).")
    p.add_argument("--min-size", type=int, default=0,
                   help="Skip de-novo clusters with fewer than this many cells (drops junk).")
    return p.parse_args()


def _decode(arr) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]


def read_profiles(h5: Path) -> tuple[pd.DataFrame, "pd.Series | None"]:
    with h5py.File(h5, "r") as f:
        for k in ("profiles", "profile_genes", "profile_types"):
            if k not in f:
                sys.exit(f"ERROR: /{k} missing in {h5}")
        mat = f["profiles"][()]
        genes = _decode(f["profile_genes"][()])
        types = _decode(f["profile_types"][()])
        cell_type = _decode(f["cell_type"][()]) if "cell_type" in f else None
    ng, nt = len(genes), len(types)
    if mat.shape == (nt, ng):
        mat = mat.T
    elif mat.shape != (ng, nt):
        sys.exit(f"ERROR: /profiles shape {mat.shape} matches neither ({ng},{nt}) nor its transpose.")
    df = pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types).astype(float)
    sizes = pd.Series(cell_type).value_counts() if cell_type is not None else None
    return df, sizes


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    profiles, sizes = read_profiles(args.profiles_h5)
    types = list(profiles.columns)
    denovo = [t for t in types if DENOVO_RE.match(t)]
    named = [t for t in types if not DENOVO_RE.match(t)]
    if args.min_size > 0 and sizes is not None:
        denovo = [d for d in denovo if int(sizes.get(d, 0)) >= args.min_size]
    if not denovo:
        sys.exit("No de-novo clusters (after --min-size). Nothing to plot.")
    # Largest de-novo first, so the heatmap reads top-to-bottom by prevalence.
    if sizes is not None:
        denovo = sorted(denovo, key=lambda d: int(sizes.get(d, 0)), reverse=True)
    print(f"{profiles.shape[0]:,} genes x {len(types)} types: "
          f"{len(named)} named + {len(denovo)} de-novo shown {denovo}")

    panel_mean = profiles.mean(axis=1).replace(0, np.nan)
    # Per de-novo cluster: top-N genes by fold-enrichment vs the panel mean.
    fold = profiles[denovo].div(panel_mean, axis=0)
    best_cluster, best_fold, marker_rows = {}, {}, []
    for d in denovo:
        for g in fold[d].sort_values(ascending=False).head(args.top_n).index:
            marker_rows.append(g)
            # A gene topping >1 cluster is assigned to the one where it is most enriched -> unique rows.
            if g not in best_fold or fold.loc[g, d] > best_fold[g]:
                best_fold[g], best_cluster[g] = fold.loc[g, d], d
    marker_genes = list(dict.fromkeys(marker_rows))  # de-dup, preserve first-seen order

    columns = denovo + named if args.include_named else denovo
    sub = profiles.loc[marker_genes, columns]
    # z-score each gene across the displayed columns (row-wise).
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    # Order rows by their assigned cluster (in the denovo display order), then by fold within.
    order = sorted(marker_genes,
                   key=lambda g: (denovo.index(best_cluster[g]), -best_fold[g]))
    z = z.loc[order]

    zmat_path = args.output_dir / "marker_heatmap_zmatrix.csv"
    mark_path = args.output_dir / "top_markers_per_cluster.csv"
    z.to_csv(zmat_path)
    pd.DataFrame({"gene": order, "cluster": [best_cluster[g] for g in order]}).to_csv(mark_path, index=False)
    print(f"wrote {zmat_path}  ({z.shape[0]} marker genes x {z.shape[1]} columns)")
    print(f"wrote {mark_path}")
    print("render:  Rscript pipeline/R/marker_heatmap.R "
          f"{args.output_dir} {args.output_dir} <denovo_annotations.csv> "
          '"De-novo cluster markers (anchor rescale profiles, z-scored across clusters)"')


if __name__ == "__main__":
    main()
