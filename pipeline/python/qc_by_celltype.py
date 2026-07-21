#!/usr/bin/env python3
"""Per-cell-type QC distributions from a typed AnnData — is a given type low-quality?

Reads obs only (backed), so it's light at cohort scale. For each cell-type group it
summarizes per-cell total counts (RNA depth), detected genes, and typing probability,
and draws a boxplot of counts per type (ordered by median, the --highlight type in red)
so you can see whether a suspicious type (e.g. InSituTree's Low_signal sink) sits at the
LOW end (a QC / low-RNA artifact) or mid-range (real cells that are transcriptionally
flat on the panel). Prints a focused Low_signal-vs-rest comparison.

Reads:
  --typed-h5ad   AnnData with obs[group-key] + a per-cell counts column. Counts column is
                 auto-detected from common names unless --count-col is given.
Writes (--output-dir):
  qc_by_celltype.csv   per-type n, count quantiles, median genes, median prob
  qc_by_celltype.png   boxplot of log10(counts) per type, highlight type in red

Usage:
    python pipeline/python/qc_by_celltype.py \\
        --typed-h5ad cosmx_typed.h5ad --group-key cell_type \\
        --highlight Low_signal --output-dir out
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COUNT_CANDIDATES = ["total_counts", "totalcounts", "tc", "n_counts", "total_counts_all", "nCount"]
GENE_CANDIDATES = ["nFeature_RNA", "n_genes_by_counts", "n_genes", "nFeature", "n_genes_all"]
PROB_CANDIDATES = ["insitutype_prob", "prob", "type_prob"]


def pick(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--group-key", default="cell_type")
    p.add_argument("--count-col", default=None,
                   help="Per-cell counts column; auto-detected if omitted.")
    p.add_argument("--highlight", default="Low_signal",
                   help="Cell type to flag (red in the plot, and the vs-rest summary).")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    obs = ad.read_h5ad(args.typed_h5ad, backed="r").obs
    print(f"obs: {len(obs):,} cells; columns: {list(obs.columns)}")

    if args.group_key not in obs:
        sys.exit(f"ERROR: group-key '{args.group_key}' not in obs")
    count_col = args.count_col or pick(obs.columns, COUNT_CANDIDATES)
    if count_col is None or count_col not in obs:
        sys.exit(f"ERROR: no counts column found (tried {COUNT_CANDIDATES}); pass --count-col")
    gene_col = pick(obs.columns, GENE_CANDIDATES)
    prob_col = pick(obs.columns, PROB_CANDIDATES)
    print(f"using count-col='{count_col}', gene-col='{gene_col}', prob-col='{prob_col}'")

    df = pd.DataFrame({"grp": obs[args.group_key].astype(str).values,
                       "cnt": pd.to_numeric(obs[count_col], errors="coerce").values})
    if gene_col:
        df["genes"] = pd.to_numeric(obs[gene_col], errors="coerce").values
    if prob_col:
        df["prob"] = pd.to_numeric(obs[prob_col], errors="coerce").values
    df = df.dropna(subset=["cnt"])

    # per-type summary
    g = df.groupby("grp")
    summ = pd.DataFrame({
        "n_cells": g.size(),
        "count_median": g["cnt"].median(),
        "count_q25": g["cnt"].quantile(0.25),
        "count_q75": g["cnt"].quantile(0.75),
        "count_mean": g["cnt"].mean(),
    })
    if gene_col:
        summ["genes_median"] = g["genes"].median()
    if prob_col:
        summ["prob_median"] = g["prob"].median()
    summ = summ.sort_values("count_median")
    summ.to_csv(args.output_dir / "qc_by_celltype.csv")
    print(f"\nwrote {args.output_dir / 'qc_by_celltype.csv'}")

    # focused highlight-vs-rest comparison
    hi = df[df["grp"] == args.highlight]["cnt"]
    rest = df[df["grp"] != args.highlight]["cnt"]
    if len(hi):
        cohort_q25 = df["cnt"].quantile(0.25)
        frac_bottom = float((hi < cohort_q25).mean())
        print(f"\n=== '{args.highlight}' vs rest (counts) ===")
        print(f"  {args.highlight}: n={len(hi):,}, median={hi.median():.0f}, "
              f"q25={hi.quantile(0.25):.0f}, q75={hi.quantile(0.75):.0f}")
        print(f"  typed rest:      n={len(rest):,}, median={rest.median():.0f}, "
              f"q25={rest.quantile(0.25):.0f}, q75={rest.quantile(0.75):.0f}")
        print(f"  ratio of medians ({args.highlight}/rest): {hi.median()/rest.median():.2f}")
        print(f"  fraction of {args.highlight} below the COHORT 25th-pct count "
              f"({cohort_q25:.0f}): {frac_bottom:.1%}")
        if prob_col:
            print(f"  {args.highlight} median typing prob: "
                  f"{df[df['grp']==args.highlight]['prob'].median():.3f}")

    # boxplot of log10(counts) per type, ordered by median, highlight in red
    order = list(summ.index)
    data = [np.log10(df[df["grp"] == t]["cnt"].clip(lower=1).values) for t in order]
    fig, ax = plt.subplots(figsize=(10, max(8, 0.28 * len(order) + 2)))
    bp = ax.boxplot(data, vert=False, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black"))
    for t, box in zip(order, bp["boxes"]):
        box.set_facecolor("#d62728" if t == args.highlight else "#9ecae1")
    ax.set_yticks(range(1, len(order) + 1)); ax.set_yticklabels(order, fontsize=7)
    ax.set_xlabel(f"log10({count_col}) per cell")
    ax.axvline(np.log10(df["cnt"].median()), color="grey", ls="--", lw=1,
               label="cohort median")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Per-cell-type RNA depth ({len(df):,} cells; {args.highlight} in red)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "qc_by_celltype.png", dpi=180, bbox_inches="tight")
    print(f"wrote {args.output_dir / 'qc_by_celltype.png'}")


if __name__ == "__main__":
    main()
