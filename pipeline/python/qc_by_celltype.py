#!/usr/bin/env python3
"""Per-cell-type QC distributions from a typed AnnData — is a given type low-quality?

Reads obs only (backed), so it's light at cohort scale. For each cell-type group it
summarizes per-cell total counts (RNA depth), detected genes, and typing probability,
and draws a boxplot of counts per type (ordered by median, the --highlight type in red)
so you can see whether a suspicious type (e.g. InSituTree's Low_signal sink) sits at the
LOW end (a QC / low-RNA artifact) or mid-range (real cells that are transcriptionally
flat on the panel). Prints a focused Low_signal-vs-rest comparison (counts AND genes
detected) plus a QC-floor sensitivity table: at each candidate panel-count floor, how
much of the highlight type vs the typed rest would be dropped (a re-QC recovery curve).

Reads:
  --typed-h5ad   AnnData with obs[group-key] + a per-cell counts column (auto-detected).
Writes (--output-dir):
  qc_by_celltype.csv     per-type n, count/gene quantiles, median prob
  qc_by_celltype.png     boxplot of log10(counts) per type, highlight in red
  qc_floor_sensitivity.png  panel-count histograms (highlight vs typed) with floor lines

Usage:
    python pipeline/python/qc_by_celltype.py \\
        --typed-h5ad cosmx_typed.h5ad --group-key cell_type \\
        --highlight Low_signal --floors 50,100,150,200,300 --output-dir out
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
# panel (6k) gene counts — the metric the Stage-3a min-50 QC floor was applied to
FLOOR_CANDIDATES = ["qc_gene_counts", "nCount_RNA", "nCount", "panel_counts"]


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
    p.add_argument("--count-col", default=None, help="Per-cell counts col; auto if omitted.")
    p.add_argument("--floor-col", default=None,
                   help="Panel-count col for the QC-floor sensitivity; auto if omitted.")
    p.add_argument("--floors", default="50,100,150,200,300",
                   help="Comma-separated panel-count floors to test.")
    p.add_argument("--highlight", default="Low_signal",
                   help="Cell type to flag (red in plots, and the vs-rest summary).")
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
    floor_col = args.floor_col or pick(obs.columns, FLOOR_CANDIDATES)
    print(f"using count-col='{count_col}', gene-col='{gene_col}', "
          f"prob-col='{prob_col}', floor-col='{floor_col}'")

    df = pd.DataFrame({"grp": obs[args.group_key].astype(str).values,
                       "cnt": pd.to_numeric(obs[count_col], errors="coerce").values})
    if gene_col:
        df["genes"] = pd.to_numeric(obs[gene_col], errors="coerce").values
    if prob_col:
        df["prob"] = pd.to_numeric(obs[prob_col], errors="coerce").values
    if floor_col:
        df["floor"] = pd.to_numeric(obs[floor_col], errors="coerce").values
    df = df.dropna(subset=["cnt"])

    # per-type summary
    g = df.groupby("grp")
    summ = pd.DataFrame({"n_cells": g.size(),
                         "count_median": g["cnt"].median(),
                         "count_q25": g["cnt"].quantile(0.25),
                         "count_q75": g["cnt"].quantile(0.75),
                         "count_mean": g["cnt"].mean()})
    if gene_col:
        summ["genes_median"] = g["genes"].median()
    if prob_col:
        summ["prob_median"] = g["prob"].median()
    summ = summ.sort_values("count_median")
    summ.to_csv(args.output_dir / "qc_by_celltype.csv")
    print(f"\nwrote {args.output_dir / 'qc_by_celltype.csv'}")

    is_hi = df["grp"] == args.highlight
    hi, rest = df[is_hi], df[~is_hi]
    if len(hi):
        print(f"\n=== '{args.highlight}' vs typed rest ===")
        cq25 = df["cnt"].quantile(0.25)
        print(f"  counts  {args.highlight}: median={hi['cnt'].median():.0f} "
              f"(q25={hi['cnt'].quantile(.25):.0f}, q75={hi['cnt'].quantile(.75):.0f});  "
              f"rest: median={rest['cnt'].median():.0f}  "
              f"-> ratio {hi['cnt'].median()/rest['cnt'].median():.2f}")
        if gene_col:
            print(f"  genes   {args.highlight}: median={hi['genes'].median():.0f} "
                  f"(q25={hi['genes'].quantile(.25):.0f}, q75={hi['genes'].quantile(.75):.0f});  "
                  f"rest: median={rest['genes'].median():.0f}  "
                  f"-> ratio {hi['genes'].median()/rest['genes'].median():.2f}")
        if prob_col:
            print(f"  {args.highlight} median typing prob: {hi['prob'].median():.3f}")
        print(f"  {args.highlight} below cohort 25th-pct count ({cq25:.0f}): "
              f"{float((hi['cnt']<cq25).mean()):.1%}")

        # --- QC-floor sensitivity (on the panel-count column) --------------------
        if floor_col:
            floors = [int(x) for x in str(args.floors).split(",") if x.strip()]
            hf, rf = hi["floor"].dropna(), rest["floor"].dropna()
            print(f"\n=== QC-floor sensitivity on '{floor_col}' (cohort already >= 50) ===")
            print(f"  {'floor':>6}  {args.highlight+' dropped':>26}  {'typed dropped':>22}  enrich")
            for T in floors:
                hn, hp = int((hf < T).sum()), float((hf < T).mean())
                tn, tp = int((rf < T).sum()), float((rf < T).mean())
                enr = (hp / tp) if tp > 0 else float("nan")
                print(f"  {T:>6}  {hn:>13,} ({hp:>5.1%})  {tn:>13,} ({tp:>5.1%})  {enr:>4.1f}x")
            print(f"  (enrich = {args.highlight}%dropped / typed%dropped; >1 = floor "
                  f"preferentially removes {args.highlight})")

    # boxplot of log10(counts) per type
    order = list(summ.index)
    data = [np.log10(df[df["grp"] == t]["cnt"].clip(lower=1).values) for t in order]
    fig, ax = plt.subplots(figsize=(10, max(8, 0.28 * len(order) + 2)))
    bp = ax.boxplot(data, vert=False, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black"))
    for t, box in zip(order, bp["boxes"]):
        box.set_facecolor("#d62728" if t == args.highlight else "#9ecae1")
    ax.set_yticks(range(1, len(order) + 1)); ax.set_yticklabels(order, fontsize=7)
    ax.set_xlabel(f"log10({count_col}) per cell")
    ax.axvline(np.log10(df["cnt"].median()), color="grey", ls="--", lw=1, label="cohort median")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Per-cell-type RNA depth ({len(df):,} cells; {args.highlight} in red)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "qc_by_celltype.png", dpi=180, bbox_inches="tight")
    print(f"\nwrote {args.output_dir / 'qc_by_celltype.png'}")

    # QC-floor histogram: highlight vs typed panel-count distributions + floor lines
    if floor_col and len(hi):
        floors = [int(x) for x in str(args.floors).split(",") if x.strip()]
        fig, ax = plt.subplots(figsize=(9, 5))
        bins = np.linspace(np.log10(max(1, df["floor"].min())),
                           np.log10(max(2, df["floor"].quantile(0.995))), 60)
        ax.hist(np.log10(hi["floor"].clip(lower=1)), bins=bins, density=True,
                alpha=0.55, color="#d62728", label=args.highlight)
        ax.hist(np.log10(rest["floor"].clip(lower=1)), bins=bins, density=True,
                alpha=0.55, color="#1f77b4", label="typed rest")
        for T in floors:
            ax.axvline(np.log10(T), color="grey", ls="--", lw=0.8)
            ax.text(np.log10(T), ax.get_ylim()[1] * 0.95, str(T), fontsize=7,
                    ha="center", va="top", rotation=90)
        ax.set_xlabel(f"log10({floor_col}) per cell"); ax.set_ylabel("density")
        ax.legend(); ax.set_title(f"Panel-count distribution: {args.highlight} vs typed (floor lines)")
        fig.tight_layout()
        fig.savefig(args.output_dir / "qc_floor_sensitivity.png", dpi=180, bbox_inches="tight")
        print(f"wrote {args.output_dir / 'qc_floor_sensitivity.png'}")


if __name__ == "__main__":
    main()
