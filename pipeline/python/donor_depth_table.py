#!/usr/bin/env python3
"""Per-donor RNA depth table: is SORL1 carrier status confounded with sequencing depth?

The batch-entropy report showed that clustering is structured by depth: 57.6% of cells sit in
six shallow clusters (median counts 204-551). Carrier donors and the pilot-named slides carry
a larger share of those. This writes one row per donor with the depth metrics next to slide,
run and diagnosis, so the confound is a table Maddie can read and not just a picture.

Reads obs only (backed mode) from a stage-3c clustered AnnData, so it reflects cells that
PASSED QC: depth here is truncated at the 50-count floor and the figures understate how
shallow the lowest donors really are.

"Shallow" clusters are those whose median total_counts is below --shallow-max-median. The
default (560) selects exactly the six clusters the entropy report lists as low-signal
candidates (1, 3, 4, 5, 8, 13); the next cluster up (2) sits at 591. The clusters used and
their medians are printed, because the cutoff is a convention and the result should not
depend on a number nobody can see.

Usage:
    donor_depth_table.py --h5ad cosmx_clustered.h5ad --annotations fov_annotations_sorl1.csv \\
        --manifest manifest_sorl1_40donors.csv --output donor_depth.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_SHALLOW_MAX_MEDIAN = 560.0
CLUSTER_KEY = "leiden"
DONOR_KEY = "Region"
DEPTH_COL = "total_counts"
GENES_COL = "qc_genes_detected"
AREA_COL = "qc_area"
CARRIER_PATTERN = "SORL1"


def shallow_clusters(obs: pd.DataFrame, max_median: float) -> pd.Series:
    """Median depth per cluster, restricted to clusters under `max_median`."""
    medians = obs.groupby(CLUSTER_KEY)[DEPTH_COL].median()
    return medians[medians < max_median].sort_values()


def donor_table(obs: pd.DataFrame, annotations: pd.DataFrame, manifest: pd.DataFrame,
                shallow: list[str]) -> pd.DataFrame:
    """One row per donor: depth metrics, shallow-cluster share, and slide/run/diagnosis."""
    grouped = obs.groupby(DONOR_KEY)
    table = pd.DataFrame({
        "n_cells": grouped.size(),
        "median_total_counts": grouped[DEPTH_COL].median(),
        "median_genes_detected": grouped[GENES_COL].median(),
        "median_area": grouped[AREA_COL].median() if AREA_COL in obs else np.nan,
        "shallow_frac": obs[CLUSTER_KEY].isin(shallow).groupby(obs[DONOR_KEY]).mean(),
    })
    table.index = table.index.astype(str)
    table.index.name = "donor"

    donor_info = (annotations.astype({"donor": str})
                  .drop_duplicates("donor")
                  .set_index("donor")[["slide_id", "case_broad", "case_group", "sorl1_mutation"]])
    run_info = manifest.set_index("slide_id")[["run_date", "instrument_id", "run_uuid"]]
    table = table.join(donor_info).join(run_info, on="slide_id")
    table["carrier"] = table["case_group"].str.contains(CARRIER_PATTERN, na=False)

    unannotated = table.index[table["slide_id"].isna()].tolist()
    if unannotated:
        raise ValueError(f"donors with no annotation row: {unannotated}")
    return table.reset_index().sort_values("shallow_frac", ascending=False)


def summarize(table: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Donor-level mean/min/max shallow fraction and median depth per group."""
    return table.groupby(by).agg(
        donors=("donor", "size"),
        mean_shallow_frac=("shallow_frac", "mean"),
        min_shallow_frac=("shallow_frac", "min"),
        max_shallow_frac=("shallow_frac", "max"),
        median_counts=("median_total_counts", "median"),
    ).round(3)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True, help="Stage-3c clustered AnnData.")
    p.add_argument("--annotations", type=Path, required=True,
                   help="fov_annotations_sorl1.csv (donor, slide_id, case_*, sorl1_mutation).")
    p.add_argument("--manifest", type=Path, required=True, help="Run manifest (run_date, instrument).")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--shallow-max-median", type=float, default=DEFAULT_SHALLOW_MAX_MEDIAN)
    args = p.parse_args()

    import anndata as ad

    print(f"Reading {args.h5ad} (backed; obs only)")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    needed = [DONOR_KEY, CLUSTER_KEY, DEPTH_COL, GENES_COL]
    absent = [c for c in needed if c not in adata.obs]
    if absent:
        print(f"ERROR: obs is missing {absent}. Available: {list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)
    obs = adata.obs[[c for c in needed + [AREA_COL] if c in adata.obs]].copy()
    obs[DONOR_KEY] = obs[DONOR_KEY].astype(str)
    obs[CLUSTER_KEY] = obs[CLUSTER_KEY].astype(str)

    shallow = shallow_clusters(obs, args.shallow_max_median)
    print(f"Shallow clusters (median {DEPTH_COL} < {args.shallow_max_median:g}), "
          f"{obs[CLUSTER_KEY].isin(shallow.index).mean():.1%} of cells:")
    print(shallow.round(0).to_string())

    table = donor_table(obs, pd.read_csv(args.annotations, dtype=str),
                        pd.read_csv(args.manifest, dtype=str), list(shallow.index))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False, float_format="%.4f")
    print(f"\nWrote {args.output}: {len(table)} donors")

    pd.set_option("display.width", 200)
    print("\nShallow-cluster share by run date and carrier status (donor-level):")
    print(summarize(table, ["run_date", "carrier"]).to_string())
    print("\nBy diagnosis:")
    print(summarize(table, ["case_group"]).to_string())
    print("\nTen shallowest donors:")
    print(table.head(10)[["donor", "slide_id", "case_group", "run_date", "median_total_counts",
                          "shallow_frac"]].to_string(index=False))


if __name__ == "__main__":
    main()
