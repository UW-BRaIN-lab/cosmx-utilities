#!/usr/bin/env python3
"""Per-cluster batch-mixing QC for a clustered AnnData.

For each cluster, report how broadly it draws across the batch variable
(patient / Case). A well-batch-corrected embedding yields clusters that each
span many donors; a cluster dominated by one or two donors signals residual
patient structure (the batch correction did not fully remove the nuisance, or
the population is genuinely patient-private and worth a closer look).

Metrics per cluster:
  n_cells          cells in the cluster
  n_donors         distinct batch values present
  eff_donors       exp(Shannon entropy of the donor distribution) — the
                   "effective" number of donors; equals n_donors when cells
                   are split evenly, drops toward 1 when one donor dominates
  norm_entropy     Shannon entropy / log(total donors); 1.0 = perfectly even
                   across the whole cohort, 0 = a single donor
  top_donor_frac   fraction of the cluster from its single largest donor
  top_donor        that donor's label

Optionally also reports, per cluster, median sequencing depth and the dominant
anatomical region. Those two turn the same pass into a low-signal check: a large
cluster with the LOWEST median depth and NO dominant region is a low-signal
catch-all, not a cell type — the population that shows up on a UMAP as the hub
everything else radiates away from.

Reads obs in backed mode, so the expression matrix is never loaded — safe to
run on a login node against a multi-GB clustered .h5ad.

Usage:
    python cluster_batch_entropy.py --h5ad cosmx_clustered.h5ad
    python cluster_batch_entropy.py --h5ad cosmx_clustered.h5ad \\
        --cluster-key leiden --batch-key Case --csv out.csv
"""

from __future__ import annotations

import argparse
import sys

import anndata as ad
import numpy as np
import pandas as pd

# A cluster is flagged for review if it effectively spans too few donors or one
# donor supplies the majority of its cells.
MIN_EFF_DONORS = 3.0
MAX_TOP_DONOR_FRAC = 0.5


def per_cluster_batch_stats(
    obs: pd.DataFrame, cluster_key: str, batch_key: str
) -> tuple[pd.DataFrame, int]:
    total_donors = obs[batch_key].nunique()
    log_total = np.log(total_donors) if total_donors > 1 else 1.0

    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        counts = sub[batch_key].value_counts()
        counts = counts[counts > 0]
        proportions = counts / counts.sum()
        entropy = float(-(proportions * np.log(proportions)).sum())
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(len(sub)),
                "n_donors": int(len(counts)),
                "eff_donors": round(float(np.exp(entropy)), 2),
                "norm_entropy": round(entropy / log_total, 3),
                "top_donor_frac": round(float(proportions.iloc[0]), 3),
                "top_donor": str(counts.index[0]),
            }
        )

    stats = pd.DataFrame(rows).sort_values("eff_donors").reset_index(drop=True)
    return stats, total_donors


def per_cluster_profile(
    obs: pd.DataFrame, cluster_key: str, depth_cols: list[str], region_key: str | None
) -> pd.DataFrame:
    """Median depth and dominant region per cluster, for spotting low-signal sinks."""
    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        rec: dict = {"cluster": cluster}
        for col in depth_cols:
            rec[f"med_{col}"] = round(float(sub[col].median()), 1)
        if region_key:
            share = sub[region_key].value_counts(normalize=True)
            rec["modal_region"] = str(share.index[0])
            rec["region_share"] = round(float(share.iloc[0]), 3)
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="Clustered AnnData (.h5ad)")
    p.add_argument("--cluster-key", default="leiden", help="obs column with cluster labels")
    p.add_argument("--batch-key", default="Case", help="obs column with the batch/donor label")
    p.add_argument("--csv", default=None, help="Optional path to also write the table as CSV")
    p.add_argument("--depth-cols", default="total_counts,qc_genes_detected",
                   help="Comma-separated numeric obs columns to report a per-cluster "
                        "median for; '' to skip. Lowest-depth clusters are low-signal "
                        "candidates.")
    p.add_argument("--region-key", default="Region",
                   help="obs column with an anatomical label, reported as the dominant "
                        "region per cluster; '' to skip. A big cluster with no dominant "
                        "region is a catch-all, not a cell type.")
    args = p.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs
    for key in (args.cluster_key, args.batch_key):
        if key not in obs:
            print(f"ERROR: obs has no '{key}' column; present: {list(obs.columns)}",
                  file=sys.stderr)
            sys.exit(1)

    stats, total_donors = per_cluster_batch_stats(obs, args.cluster_key, args.batch_key)

    # Depth / region columns are optional: skipped silently when absent, so this stays
    # usable on a study whose obs lacks them.
    depth_cols = [c.strip() for c in args.depth_cols.split(",") if c.strip()]
    missing = [c for c in depth_cols if c not in obs]
    if missing:
        print(f"WARN: depth column(s) not in obs, skipped: {missing}", file=sys.stderr)
    depth_cols = [c for c in depth_cols if c in obs]
    region_key = args.region_key if args.region_key and args.region_key in obs else None
    if args.region_key and not region_key:
        print(f"WARN: region key '{args.region_key}' not in obs; skipping",
              file=sys.stderr)
    if depth_cols or region_key:
        stats = stats.merge(
            per_cluster_profile(obs, args.cluster_key, depth_cols, region_key),
            on="cluster", how="left")

    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", 200)
    print(f"{len(obs):,} cells in {len(stats)} '{args.cluster_key}' clusters; "
          f"{total_donors} total '{args.batch_key}' donors\n")
    print(stats.to_string(index=False))

    flagged = stats[(stats["eff_donors"] < MIN_EFF_DONORS)
                    | (stats["top_donor_frac"] > MAX_TOP_DONOR_FRAC)]
    print(f"\nFlagged (eff_donors < {MIN_EFF_DONORS} or top_donor_frac > "
          f"{MAX_TOP_DONOR_FRAC}) — possible residual patient structure:")
    print("  none — every cluster draws broadly across donors"
          if flagged.empty else flagged.to_string(index=False))

    # Low-signal candidates: the biggest clusters ranked by shallowest median depth.
    # A cluster that is simultaneously large, shallow and spatially diffuse is the
    # catch-all sink, and its size is the number to compare against the typing result.
    if depth_cols:
        primary = f"med_{depth_cols[0]}"
        big = stats[stats["n_cells"] >= 0.01 * len(obs)].copy()
        big["pct_cells"] = (big["n_cells"] / len(obs) * 100).round(1)
        cols = ["cluster", "n_cells", "pct_cells", primary]
        if region_key:
            cols += ["modal_region", "region_share"]
        print(f"\nLow-signal candidates — clusters >=1% of cells, shallowest first "
              f"by {primary}:")
        print(big.sort_values(primary)[cols].head(6).to_string(index=False))

    if args.csv:
        stats.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
