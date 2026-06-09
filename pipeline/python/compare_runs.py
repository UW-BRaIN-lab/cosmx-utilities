#!/usr/bin/env python3
"""Compare two Stage 3 runs (e.g. different per-cell QC thresholds).

Diffs two cosmx_clustered.h5ad outputs — typically the default `stage3/` run vs a
parallel `stage3_q100/` run (see STAGE3_DIR in the slurm scripts) — into small
comparison tables: how many cells each kept (overall, per Case, per Region), how
many Leiden clusters each found, and how similar the two clusterings are on the
cells they share (adjusted Rand index + a cluster cross-tab).

Reads obs only (backed mode), so it stays light despite the multi-GB AnnDatas.

Writes to <output-dir>:
  run_comparison_summary.csv      one row per metric, columns = run A, run B, delta
  run_comparison_by_case.csv      cells per Case in each run + delta + % retained
  run_comparison_by_region.csv    cells per Region in each run + delta
  run_comparison_cluster_xtab.csv leiden(A) x leiden(B) cell counts on shared cells

Usage:
    uv run python pipeline/python/compare_runs.py \\
        --run-a combined_a/cosmx_clustered.h5ad --label-a min50 \\
        --run-b combined_b/cosmx_clustered.h5ad --label-b min100 \\
        --output-dir run_comparison
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-a", type=Path, required=True, help="Baseline clustered .h5ad.")
    p.add_argument("--run-b", type=Path, required=True, help="Comparison clustered .h5ad.")
    p.add_argument("--label-a", default="A", help="Short label for run A (e.g. min50).")
    p.add_argument("--label-b", default="B", help="Short label for run B (e.g. min100).")
    p.add_argument("--group-key", default="leiden", help="Cluster column in obs.")
    p.add_argument("--case-key", default="Case")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def load_obs(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"ERROR: {path} not found", file=sys.stderr)
        sys.exit(1)
    return ad.read_h5ad(path, backed="r").obs.copy()


def main() -> None:
    args = parse_args()
    a, b = args.label_a, args.label_b
    obs_a = load_obs(args.run_a)
    obs_b = load_obs(args.run_b)
    for key in (args.group_key, args.case_key, args.region_key):
        for label, obs in ((a, obs_a), (b, obs_b)):
            if key not in obs:
                print(f"ERROR: run {label} obs missing '{key}'", file=sys.stderr)
                sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- shared cells + clustering agreement (ARI) --------------------------------
    shared = obs_a.index.intersection(obs_b.index)
    ari = np.nan
    if len(shared) > 1:
        from sklearn.metrics import adjusted_rand_score
        ari = adjusted_rand_score(
            obs_a.loc[shared, args.group_key].astype(str),
            obs_b.loc[shared, args.group_key].astype(str),
        )

    def _median_counts(obs: pd.DataFrame) -> float:
        return float(obs["total_counts"].median()) if "total_counts" in obs else np.nan

    # --- summary table ------------------------------------------------------------
    summary = pd.DataFrame(
        {
            "metric": ["n_cells", "n_clusters", "median_total_counts",
                       "n_case", "n_region", "n_shared_cells",
                       "adjusted_rand_index"],
            a: [len(obs_a), obs_a[args.group_key].nunique(), _median_counts(obs_a),
                obs_a[args.case_key].nunique(), obs_a[args.region_key].nunique(),
                len(shared), np.nan],
            b: [len(obs_b), obs_b[args.group_key].nunique(), _median_counts(obs_b),
                obs_b[args.case_key].nunique(), obs_b[args.region_key].nunique(),
                len(shared), ari],
        }
    )
    summary["delta_b_minus_a"] = pd.to_numeric(summary[b], errors="coerce") - \
        pd.to_numeric(summary[a], errors="coerce")

    # --- per-Case and per-Region cell counts --------------------------------------
    def _counts_by(key: str) -> pd.DataFrame:
        ca = obs_a[key].astype(str).value_counts().rename(f"n_cells_{a}")
        cb = obs_b[key].astype(str).value_counts().rename(f"n_cells_{b}")
        df = pd.concat([ca, cb], axis=1).fillna(0).astype(int)
        df.index.name = key
        df[f"delta_{b}_minus_{a}"] = df[f"n_cells_{b}"] - df[f"n_cells_{a}"]
        df[f"pct_retained_{b}"] = (100 * df[f"n_cells_{b}"] /
                                   df[f"n_cells_{a}"].replace(0, np.nan)).round(1)
        return df.sort_index()

    by_case = _counts_by(args.case_key)
    by_region = _counts_by(args.region_key)

    # --- cluster cross-tab on shared cells ---------------------------------------
    xtab = pd.DataFrame()
    if len(shared) > 1:
        xtab = pd.crosstab(
            obs_a.loc[shared, args.group_key].astype(str).rename(f"{a}_cluster"),
            obs_b.loc[shared, args.group_key].astype(str).rename(f"{b}_cluster"),
        )

    summary.to_csv(args.output_dir / "run_comparison_summary.csv", index=False)
    by_case.to_csv(args.output_dir / "run_comparison_by_case.csv")
    by_region.to_csv(args.output_dir / "run_comparison_by_region.csv")
    xtab.to_csv(args.output_dir / "run_comparison_cluster_xtab.csv")

    # --- print a readable summary -------------------------------------------------
    pd.set_option("display.max_rows", None, "display.width", 200)
    print(f"\n=== Run comparison: {a} (run A) vs {b} (run B) ===")
    print(summary.to_string(index=False))
    print(f"\nCells retained per {args.region_key}:")
    print(by_region.to_string())
    print(f"\nCells per {args.case_key} (head):")
    print(by_case.head(15).to_string())
    if not np.isnan(ari):
        print(f"\nAdjusted Rand index ({a} vs {b} clusters on {len(shared):,} shared "
              f"cells): {ari:.3f}  (1=identical, ~0=random)")
    print(f"\nWrote 4 CSVs to {args.output_dir}")


if __name__ == "__main__":
    main()
