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


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="Clustered AnnData (.h5ad)")
    p.add_argument("--cluster-key", default="leiden", help="obs column with cluster labels")
    p.add_argument("--batch-key", default="Case", help="obs column with the batch/donor label")
    p.add_argument("--csv", default=None, help="Optional path to also write the table as CSV")
    args = p.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs
    for key in (args.cluster_key, args.batch_key):
        if key not in obs:
            print(f"ERROR: obs has no '{key}' column; present: {list(obs.columns)}",
                  file=sys.stderr)
            sys.exit(1)

    stats, total_donors = per_cluster_batch_stats(obs, args.cluster_key, args.batch_key)

    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", 160)
    print(f"{len(obs):,} cells in {len(stats)} '{args.cluster_key}' clusters; "
          f"{total_donors} total '{args.batch_key}' donors\n")
    print(stats.to_string(index=False))

    flagged = stats[(stats["eff_donors"] < MIN_EFF_DONORS)
                    | (stats["top_donor_frac"] > MAX_TOP_DONOR_FRAC)]
    print(f"\nFlagged (eff_donors < {MIN_EFF_DONORS} or top_donor_frac > "
          f"{MAX_TOP_DONOR_FRAC}) — possible residual patient structure:")
    print("  none — every cluster draws broadly across donors"
          if flagged.empty else flagged.to_string(index=False))

    if args.csv:
        stats.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
