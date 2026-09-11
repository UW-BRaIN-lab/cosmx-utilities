#!/usr/bin/env python3
"""Why is a de-novo letter's cells a BAD fit to the GBmap type they were forced into?

The PI's question: given that these cells did not fit any named type well enough to be called
one, what about their expression makes them a poor match? Are they simply dimmer versions of the
named type, or are they compositionally different?

The marker heatmaps cannot answer this. Their matrix is z-scored per gene ACROSS the displayed
groups, which by construction removes amplitude — every row is forced to mean 0, sd 1. A group
that is a perfectly proportional but half-as-bright copy of its neighbour looks identical to one
that differs in composition. So this works from raw counts instead and separates the two
explanations:

  DEPTH        median counts and genes detected per cell, letter->D against D [native].
               A ratio near 0.5 means the letter's cells really do carry about half the signal.
               Note InSituType fits a per-cell scaling term, so depth alone should NOT make a
               cell fit badly — if depth is the only difference, the misfit is not explained.

  FLATNESS     per cell, the share of its counts sitting in its own top 20 genes, median over the
               group. This is what "no distinctive profile" looks like numerically: a cell whose
               expression is spread thin across the panel has nothing for a reference profile to
               latch onto, independent of how deep it is.

  COMPOSITION  Pearson r between the two groups' mean log-normalised profiles over ALL shared
               genes. Log-normalisation already divides out depth, so this isolates shape. High r
               with low depth ratio = a dimmer copy. Low r = a different thing.

  THE RESIDUAL the genes most over- and under-expressed in letter->D relative to D [native], on
               the depth-corrected scale. Under-expressed = the type's identity the cells are
               missing. Over-expressed = whatever they carry instead.

Inputs mirror denovo_vs_native_pseudobulk.py:
  --counts-h5 / --typing-h5 / --forced-csv / --letter / --destinations (or --crosstab)
Outputs:
  <out>/misfit_summary.csv        one row per letter->D pair, with the metrics above
  <out>/residual_genes.csv        per pair, the top over- and under-expressed genes

Usage:
    uv run python pipeline/python/denovo_misfit_diagnostics.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv \\
        --letter t --destinations OPC-like,AC-like,MES-like_hypoxia_MHC \\
        --output-dir t_misfit
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from anchor_profiles import read_cell_calls
from denovo_vs_native_pseudobulk import destinations_from_crosstab, read_counts
from pseudobulk_core import DEFAULT_SCALE_FACTOR, log_normalize

TOP_CONCENTRATION_GENES = 20
TOP_RESIDUAL_GENES = 12
# A pair is called "a dimmer copy" only when composition is essentially preserved AND the depth
# gap is real; "different composition" when shape itself has moved.
SAME_SHAPE_R = 0.90
DIFFERENT_SHAPE_R = 0.80
REAL_DEPTH_GAP = 0.75
# Flatness is checked first and named explicitly: a cell whose counts are spread thin across the
# panel has no peaks for any reference profile to match, which is a different failure from either
# being dim or being a different cell type — and it is the one the Low_signal sink is made of.
FLAT_SHARE_RATIO = 0.75


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True)
    p.add_argument("--typing-h5", type=Path, required=True)
    p.add_argument("--forced-csv", type=Path, required=True)
    p.add_argument("--letter", required=True)
    p.add_argument("--destinations")
    p.add_argument("--crosstab", type=Path)
    p.add_argument("--top-destinations", type=int, default=3)
    p.add_argument("--min-destination-pct", type=float, default=5.0)
    p.add_argument("--min-group-n", type=int, default=50)
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def per_cell_stats(counts_cxg: sp.csr_matrix) -> pd.DataFrame:
    """Per-cell depth and flatness from raw counts (cells x genes)."""
    totals = np.asarray(counts_cxg.sum(axis=1)).ravel()
    detected = np.diff(counts_cxg.indptr)                      # nonzero genes per cell
    # Share of a cell's counts in its own top-N genes. Done row by row on the CSR data array so
    # the full cells x genes matrix is never densified.
    top_share = np.zeros(counts_cxg.shape[0])
    data, indptr = counts_cxg.data, counts_cxg.indptr
    for i in range(counts_cxg.shape[0]):
        row = data[indptr[i]:indptr[i + 1]]
        if row.size == 0 or totals[i] <= 0:
            continue
        k = min(TOP_CONCENTRATION_GENES, row.size)
        top_share[i] = np.partition(row, -k)[-k:].sum() / totals[i]
    return pd.DataFrame({"total_counts": totals, "genes_detected": detected,
                         "top20_share": top_share})


def verdict(depth_ratio: float, comp_r: float, share_ratio: float) -> str:
    if share_ratio < FLAT_SHARE_RATIO:
        return ("FLAT - no distinctive profile, and dimmer" if depth_ratio < REAL_DEPTH_GAP
                else "FLAT - no distinctive profile")
    if comp_r >= SAME_SHAPE_R and depth_ratio < REAL_DEPTH_GAP:
        return "a dimmer copy of the type"
    if comp_r < DIFFERENT_SHAPE_R:
        return ("different composition, and dimmer" if depth_ratio < REAL_DEPTH_GAP
                else "different composition")
    if depth_ratio < REAL_DEPTH_GAP:
        return "mostly dimmer, shape partly shifted"
    return "similar depth and shape"


def main() -> None:
    args = parse_args()
    if args.destinations:
        destinations = [d.strip() for d in args.destinations.split(",") if d.strip()]
    elif args.crosstab:
        destinations = destinations_from_crosstab(
            args.crosstab, args.letter, args.top_destinations, args.min_destination_pct)
    else:
        sys.exit("ERROR: pass --destinations, or --crosstab to derive them.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]]
    forced = pd.read_csv(args.forced_csv, usecols=["cell_id", "top1_type"])
    labels = calls.merge(forced, on="cell_id", how="inner").set_index("cell_id")

    counts, genes, cell_id = read_counts(args.counts_h5)
    print(f"Counts: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")
    semisup = labels["cell_type"].reindex(cell_id)
    forced_call = labels["top1_type"].reindex(cell_id)

    summary, residuals = [], []
    for dest in destinations:
        masks = {
            f"{args.letter}->{dest}": ((semisup == args.letter) & (forced_call == dest)).to_numpy(),
            f"{dest} [native]": (semisup == dest).to_numpy(),
        }
        sizes = {k: int(v.sum()) for k, v in masks.items()}
        if min(sizes.values()) < args.min_group_n:
            print(f"  skipping {dest}: group sizes {sizes} below --min-group-n")
            continue

        stats, profiles = {}, {}
        for name, mask in masks.items():
            sub = counts[:, mask].T.tocsr()                    # cells x genes, subset first
            stats[name] = per_cell_stats(sub)
            profiles[name] = np.asarray(
                log_normalize(sub, args.scale_factor).mean(axis=0)).ravel()

        a, b = f"{args.letter}->{dest}", f"{dest} [native]"
        depth_ratio = (stats[a]["total_counts"].median()
                       / max(stats[b]["total_counts"].median(), 1))
        genes_ratio = (stats[a]["genes_detected"].median()
                       / max(stats[b]["genes_detected"].median(), 1))
        comp_r = float(np.corrcoef(profiles[a], profiles[b])[0, 1])

        summary.append({
            "pair": a, "destination": dest,
            "n_letter": sizes[a], "n_native": sizes[b],
            "counts_letter": round(stats[a]["total_counts"].median(), 1),
            "counts_native": round(stats[b]["total_counts"].median(), 1),
            "depth_ratio": round(depth_ratio, 3),
            "genes_letter": round(stats[a]["genes_detected"].median(), 1),
            "genes_native": round(stats[b]["genes_detected"].median(), 1),
            "genes_ratio": round(genes_ratio, 3),
            "top20_share_letter": round(stats[a]["top20_share"].median(), 3),
            "top20_share_native": round(stats[b]["top20_share"].median(), 3),
            "composition_r": round(comp_r, 3),
            "share_ratio": round(stats[a]["top20_share"].median()
                                 / max(stats[b]["top20_share"].median(), 1e-9), 3),
            "verdict": verdict(depth_ratio, comp_r,
                               stats[a]["top20_share"].median()
                               / max(stats[b]["top20_share"].median(), 1e-9)),
        })

        diff = pd.Series(profiles[a] - profiles[b], index=genes).sort_values()
        for gene, d in list(diff.head(TOP_RESIDUAL_GENES).items()):
            residuals.append({"pair": a, "direction": "under-expressed vs native",
                              "gene": gene, "log_diff": round(float(d), 3)})
        for gene, d in list(diff.tail(TOP_RESIDUAL_GENES)[::-1].items()):
            residuals.append({"pair": a, "direction": "over-expressed vs native",
                              "gene": gene, "log_diff": round(float(d), 3)})

    if not summary:
        sys.exit("ERROR: no destination had enough cells in both groups.")
    summ = pd.DataFrame(summary)

    print(f"\n=== {args.letter}: why these cells fit their forced type badly ===")
    print(f"{'pair':34s} {'depth':>7s} {'genes':>7s} {'top20 share':>22s} {'shape r':>8s}  verdict")
    for r in summary:
        print(f"  {r['pair']:32s} {r['depth_ratio']:>7.2f} {r['genes_ratio']:>7.2f} "
              f"{r['top20_share_letter']:>10.3f} vs {r['top20_share_native']:<8.3f} "
              f"{r['composition_r']:>8.2f}  {r['verdict']}")
    print("\n  depth / genes are letter->D medians over D [native] medians; 1.00 = same.")
    print("  top20 share = fraction of a cell's counts in its own 20 biggest genes (flatness).")
    print("  shape r = correlation of depth-corrected mean profiles; high = same shape.")
    print("  NOTE InSituType fits a per-cell scale term, so depth alone should not cause misfit —")
    print("  a pair that is ONLY dimmer has not really been explained by these numbers.")

    res = pd.DataFrame(residuals)
    print(f"\nGenes the letter's cells are MISSING relative to their forced type:")
    for pair, grp in res[res.direction.str.startswith("under")].groupby("pair", sort=False):
        print(f"  {pair:32s} " + ", ".join(grp.head(8)["gene"]))
    print(f"Genes they carry INSTEAD:")
    for pair, grp in res[res.direction.str.startswith("over")].groupby("pair", sort=False):
        print(f"  {pair:32s} " + ", ".join(grp.head(8)["gene"]))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summ.to_csv(args.output_dir / "misfit_summary.csv", index=False)
    res.to_csv(args.output_dir / "residual_genes.csv", index=False)
    print(f"\nWrote {args.output_dir}/misfit_summary.csv ({len(summ)} pairs)")
    print(f"Wrote {args.output_dir}/residual_genes.csv ({len(res)} rows)")


if __name__ == "__main__":
    main()
