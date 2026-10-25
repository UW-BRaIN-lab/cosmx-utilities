#!/usr/bin/env python3
"""Do an InSituType run's clusters share marker genes with a reference's cell types?

The raw-profile comparison (denovo_vs_reference_cosine.py) failed for the SEA-AD MTG reference: its
single-nucleus profiles sit at cosine/Spearman 0.2-0.4 from CosMx clusters, against 0.75-0.9 for a
CosMx-derived reference. That says the profiles are not comparable, not that the cell types differ.
Marker overlap is robust to that: it asks only whether the genes that stand out in a CosMx cluster
also stand out in a reference type, whatever the platform did to absolute levels.

For every type in the result (named and de novo):
  * cluster markers = the top genes by log2 fold-change against the other types, with an expression
    floor (top_markers_from_profiles.top_markers, the same lists the typing doc uses);
  * for each reference type, concordance = the mean percentile, across the shared genes, of those
    markers' fold-change in that reference type against the other reference types (0.5 = chance,
    1.0 = every cluster marker is among that type's most enriched genes), and
  * overlap = how many cluster markers are also in that reference type's own top markers, with a
    hypergeometric p-value against the shared-gene universe.
Named Brain_6k types in the same run are the positive control: Oligodendrocyte should land on
the reference Oligodendrocyte, L2/3 neuron on L2/3 IT, and so on. If they do not, the method is not
working on this data and the de novo rows mean nothing.

Reference fold-changes use a pseudocount set to a quarter of the reference's median non-zero value
and its markers must also exceed that median, so the answer does not depend on the reference's
scale (the CosMx side keeps top_markers_from_profiles' own 0.2 floor, which that module calibrates
to the fitted profiles).

Writes <output-dir>/marker_overlap.csv (long: cluster x reference type) and
<output-dir>/marker_concordance_matrix.csv (cluster x reference type), and prints each cluster's
three best reference types, largest clusters first.

Usage:
    marker_overlap_vs_reference.py --result-h5 insitutype_result.h5 \\
        --reference-csv seaad_mtg_panel.csv --output-dir marker_overlap/
"""

from __future__ import annotations

import argparse
import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

from anchor_profiles import is_denovo, read_cell_calls
from top_markers_from_profiles import DEFAULT_MIN_PROFILE, DEFAULT_TOP_N, read_profiles, top_markers

PSEUDOCOUNT_FRACTION = 0.25      # of the reference's median non-zero value
MIN_MARKERS_USED = 5             # a cluster with fewer shared markers is skipped, not scored
TOP_MATCHES_PRINTED = 3


def hypergeom_sf(overlap: int, universe: int, successes: int, draws: int) -> float:
    """P(X >= overlap) for X ~ Hypergeometric(universe, successes, draws)."""
    total = comb(universe, draws)
    tail = sum(comb(successes, i) * comb(universe - successes, draws - i)
               for i in range(overlap, min(successes, draws) + 1))
    return tail / total


def reference_enrichment(reference: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """genes x types log2 fold-change against the other types, and the reference's median non-zero."""
    values = reference.to_numpy(dtype=float)
    median_nonzero = float(np.median(values[values > 0]))
    pseudocount = PSEUDOCOUNT_FRACTION * median_nonzero
    n_types = reference.shape[1]
    others = (reference.sum(axis=1).to_numpy()[:, None] - values) / (n_types - 1)
    fold = np.log2((values + pseudocount) / (others + pseudocount))
    return pd.DataFrame(fold, index=reference.index, columns=reference.columns), median_nonzero


def score_clusters(markers: pd.DataFrame, reference: pd.DataFrame, top_n: int,
                   sizes: pd.Series) -> pd.DataFrame:
    """Long table: one row per (cluster, reference type)."""
    fold, median_nonzero = reference_enrichment(reference)
    percentile = fold.rank(axis=0, pct=True)
    universe = len(reference.index)
    ref_top = {}
    for t in reference.columns:
        eligible = fold[t][reference[t] >= median_nonzero]
        ref_top[t] = set(eligible.sort_values(ascending=False).head(top_n).index)

    rows = []
    for cluster, grp in markers.groupby("type", sort=False):
        genes = [g for g in grp["gene"] if g in reference.index]
        if len(genes) < MIN_MARKERS_USED:
            print(f"  skipping {cluster}: {len(genes)} markers shared with the reference "
                  f"(need {MIN_MARKERS_USED})", file=sys.stderr)
            continue
        for t in reference.columns:
            hits = len(ref_top[t].intersection(genes))
            rows.append({
                "cluster": cluster, "kind": "de novo" if is_denovo(cluster) else "named",
                "n_cells": int(sizes.get(cluster, 0)), "n_markers_used": len(genes),
                "reference_type": t,
                "concordance": round(float(percentile.loc[genes, t].mean()), 3),
                "overlap": hits, "overlap_p": hypergeom_sf(hits, universe, len(ref_top[t]), len(genes)),
            })
    table = pd.DataFrame(rows)
    table["rank"] = table.groupby("cluster")["concordance"].rank(ascending=False, method="first").astype(int)
    return table.sort_values(["n_cells", "cluster", "rank"], ascending=[False, True, True]).reset_index(drop=True)


def print_summary(table: pd.DataFrame) -> None:
    for cluster, grp in table[table["rank"] <= TOP_MATCHES_PRINTED].groupby("cluster", sort=False):
        head = grp.iloc[0]
        parts = [f"{r.reference_type} {r.concordance:.2f} ({r.overlap}/{r.n_markers_used}, p={r.overlap_p:.1e})"
                 for r in grp.itertuples()]
        print(f"{cluster:<28} {head.kind:<8} {head.n_cells:>8,}  " + "  |  ".join(parts))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--result-h5", type=Path, required=True, help="InSituType insitutype_result.h5.")
    p.add_argument("--reference-csv", type=Path, required=True, help="genes x types reference.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N, help="markers per type, both sides.")
    p.add_argument("--min-profile", type=float, default=DEFAULT_MIN_PROFILE,
                   help="expression floor for the CosMx-side markers (top_markers_from_profiles).")
    args = p.parse_args()

    profiles = read_profiles(args.result_h5)
    sizes = read_cell_calls(args.result_h5)["cell_type"].value_counts()
    reference = pd.read_csv(args.reference_csv, index_col=0)
    shared = profiles.index.intersection(reference.index)
    if len(shared) < 0.5 * len(profiles.index):
        sys.exit(f"ERROR: only {len(shared)} of {len(profiles.index)} result genes are in "
                 f"{args.reference_csv}; gene naming differs.")
    print(f"{len(shared):,} genes shared; {profiles.shape[1]} result types x "
          f"{reference.shape[1]} reference types")

    markers = top_markers(profiles, args.top_n, args.min_profile)
    table = score_clusters(markers, reference.loc[shared], args.top_n, sizes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "marker_overlap.csv", index=False)
    table.pivot(index="cluster", columns="reference_type", values="concordance").to_csv(
        args.output_dir / "marker_concordance_matrix.csv")
    print_summary(table)
    print(f"\nWrote {args.output_dir}/marker_overlap.csv and marker_concordance_matrix.csv")


if __name__ == "__main__":
    main()
