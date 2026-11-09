#!/usr/bin/env python3
"""Do the same cell type's profiles from different atlases carry REGION, or atlas signature?

The combined reference holds one astrocyte from HRCA, one from the Allen atlas and one from
Monavarfeshani; likewise microglia, endothelium, pericyte, smooth muscle, fibroblast, and the
two pigmented epithelia. Measured on the panel genes these pairs are the FARTHEST apart in the
whole reference -- Mona vs Allen pericyte sits at cosine 0.84 where two retinal ganglion
subtypes sit at 0.998. Two endothelial profiles that far apart is not plausible biology, so
one of two things is true, and they call for opposite actions:

  * they carry tissue identity -- retinal astrocyte, optic-nerve-head astrocyte and cortical
    astrocyte really are different cells. Then all three belong in the reference, and each one
    should win in its own region.
  * they carry ATLAS signature -- the dissociation, chemistry and depth of three different
    studies. Then they are three noisy copies of one profile, they compete for the same cells
    anywhere in the tissue, and InSituType's single platform-rescaling term cannot fit all
    three at once. Then keep one per class.

The data decides it. For each class this cross-tabulates the run's calls against the
neuropathologist's region, and reports ENRICHMENT -- the share of a type's cells in a region
divided by that region's share of all cells. Raw percentages cannot answer this: whichever
region holds most of the cohort is the modal region of nearly every type.

Region comes from the FOV annotation sheet joined on (slide_id, fov), NOT from obs['Region'].
A typed AnnData built before the latest annotation round carries the stale export value, and
this diagnostic exists to be run against exactly such a file.

Writes, under --output-dir:
    atlas_duplicates_summary.csv     one row per type: cells, home region, enrichment, verdict
    atlas_duplicates_<class>.csv     the full region x type table for each class

Usage:
    atlas_duplicates_by_region.py --typed-h5ad cosmx_typed.h5ad \\
        --fov-annotations fov_annotations.csv \\
        --classes atlas_duplicate_classes.csv --output-dir atlas_duplicates/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

from fov_exclusions import FOV_COL, SLIDE_COL
from qc_by_region import REGION_COL, attach_regions, load_regions

CLASS_COL = "class"
TYPE_COL = "type"
# Below this a type's regional distribution is noise; reported, never given a verdict.
MIN_CELLS_TO_JUDGE = 200
# Enrichment alone cannot carry the verdict: it is capped at 1/region_share, so in a cohort
# where gray matter holds 70% of cells NOTHING can exceed 1.43 there, and a fixed bar would
# rule out the largest region by construction. CONCENTRATION rescales it against that cap --
# (pct_of_type - region_share) / (1 - region_share) -- so 1.0 means every cell of the type is
# in that region, 0.0 means it is spread exactly like the cohort, and the bar means the same
# thing in a region of any size.
STRONG_CONCENTRATION = 0.5


def load_classes(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = [c for c in (CLASS_COL, TYPE_COL) if c not in frame.columns]
    if missing:
        sys.exit(f"ERROR: {path} has no {missing} column(s)")
    frame[CLASS_COL] = frame[CLASS_COL].astype(str).str.strip()
    frame[TYPE_COL] = frame[TYPE_COL].astype(str).str.strip()
    return frame


def region_shares(cells: pd.DataFrame) -> pd.Series:
    """Each region's share of ALL cells -- the denominator enrichment is measured against."""
    return cells[REGION_COL].value_counts(normalize=True)


def class_table(cells: pd.DataFrame, types: list[str], baseline: pd.Series,
                label_key: str) -> pd.DataFrame:
    """Region x type counts, with the row percentages and enrichment beside them."""
    present = [t for t in types if (cells[label_key] == t).any()]
    counts = pd.crosstab(cells[REGION_COL], cells[label_key])
    counts = counts.reindex(columns=present, fill_value=0)
    counts = counts.reindex(index=baseline.index, fill_value=0)
    shares = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1)
    enrichment = shares.div(baseline, axis=0)
    concentration = (shares.sub(baseline, axis=0)).div(1 - baseline, axis=0)
    out = pd.concat({"n_cells": counts, "pct_of_type": 100 * shares,
                     "enrichment": enrichment, "concentration": concentration}, axis=1)
    return out


def verdict(summary: pd.DataFrame) -> str:
    """REGIONAL, PARTIAL or MIXED for one class."""
    judged = summary[summary["n_cells"] >= MIN_CELLS_TO_JUDGE]
    if len(judged) < 2:
        return "TOO FEW CELLS"
    homes = judged["home_region"].tolist()
    if len(set(homes)) < len(homes):
        return "MIXED"
    if (judged["home_concentration"] >= STRONG_CONCENTRATION).all():
        return "REGIONAL"
    return "PARTIAL"


def cramers_v(counts: pd.DataFrame) -> float:
    """Strength of the region-by-type association, 0 (none) to 1 (deterministic)."""
    table = counts.loc[(counts.sum(axis=1) > 0), (counts.sum(axis=0) > 0)]
    if table.shape[0] < 2 or table.shape[1] < 2:
        return float("nan")
    chi2 = chi2_contingency(table.to_numpy())[0]
    n = table.to_numpy().sum()
    return float(np.sqrt(chi2 / (n * (min(table.shape) - 1))))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--typed-h5ad", type=Path, required=True,
                        help="AnnData whose obs carries the typing (read backed, obs only)")
    parser.add_argument("--fov-annotations", type=Path, required=True)
    parser.add_argument("--classes", type=Path, required=True,
                        help="class,type CSV (pipeline/reference/atlas_duplicate_classes.csv)")
    parser.add_argument("--label-key", default="cell_type",
                        help="obs column holding the InSituType call (default cell_type)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--include-excluded", action="store_true",
                        help="keep FOVs the SME ruled out (default: drop, as the refit will)")
    args = parser.parse_args(argv)

    label_key = args.label_key

    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    obs = adata.obs
    print(f"{args.typed_h5ad.name}: {len(obs):,} cells")
    if label_key not in obs.columns:
        sys.exit(f"ERROR: obs has no {label_key!r} column; it has {list(obs.columns)[:12]}")

    cells = attach_regions(obs, load_regions(args.fov_annotations), require=(label_key,))
    cells[label_key] = cells[label_key].astype(str)
    if not args.include_excluded:
        n_excluded = int(cells["excluded"].sum())
        cells = cells.loc[~cells["excluded"]]
        print(f"dropped {n_excluded:,} cells in FOVs the SME excluded")
    print(f"{len(cells):,} cells across {cells[REGION_COL].nunique()} regions, "
          f"{cells[SLIDE_COL].nunique()} slides")

    baseline = region_shares(cells)
    print("\nregion share of the cohort (the enrichment denominator):")
    for region, share in baseline.items():
        print(f"  {region:<24} {100 * share:5.1f}%")

    classes = load_classes(args.classes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for class_name, group in classes.groupby(CLASS_COL, sort=True):
        types = group[TYPE_COL].tolist()
        absent = [t for t in types if not (cells[label_key] == t).any()]
        table = class_table(cells, types, baseline, label_key)
        table.to_csv(args.output_dir / f"atlas_duplicates_{class_name}.csv")

        summary = []
        for type_name in types:
            if type_name in absent:
                summary.append({"class": class_name, TYPE_COL: type_name, "n_cells": 0,
                                "home_region": "", "home_concentration": np.nan,
                                "home_enrichment": np.nan, "home_pct_of_type": np.nan,
                                "modal_region": ""})
                continue
            column = table["concentration"][type_name]
            home = column.idxmax()
            summary.append({
                "class": class_name, TYPE_COL: type_name,
                "n_cells": int(table["n_cells"][type_name].sum()),
                "home_region": home,
                "home_concentration": float(column.loc[home]),
                "home_enrichment": float(table["enrichment"][type_name].loc[home]),
                "home_pct_of_type": float(table["pct_of_type"][type_name].loc[home]),
                # Where the RAW percentages would have put it. Differs from home_region
                # whenever one region dominates the cohort, which is the whole point.
                "modal_region": table["pct_of_type"][type_name].idxmax(),
            })
        summary = pd.DataFrame(summary)
        call = verdict(summary)
        summary["verdict"] = call
        summary["cramers_v"] = cramers_v(table["n_cells"])
        rows.append(summary)

        print(f"\n================ {class_name}  --  {call} "
              f"(Cramer's V {summary['cramers_v'].iloc[0]:.3f}) ================")
        if absent:
            print(f"  not called at all: {', '.join(absent)}")
        print("  concentration by region (1.0 = every cell of the type is here, 0.0 = "
              "spread like\n  the cohort, negative = under-represented there):")
        with pd.option_context("display.width", 170, "display.max_columns", None):
            print(table["concentration"].round(2).to_string())
        print(summary[[TYPE_COL, "n_cells", "home_region", "home_concentration",
                       "home_enrichment", "home_pct_of_type", "modal_region"]]
              .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(args.output_dir / "atlas_duplicates_summary.csv", index=False)
    print(f"\nWrote {1 + len(rows)} tables to {args.output_dir}")

    print("\nREAD:")
    print("  REGIONAL -- each atlas copy owns its own tissue. Keep them all; the split is "
          "biology.")
    print("  MIXED    -- two copies compete for the same tissue. Keep ONE per class, the one "
          "from\n              the atlas whose tissue this cohort actually sampled.")
    print("  PARTIAL  -- the copies lean different ways without owning a tissue. Decide with "
          "the\n              calibration table: a copy with no pinned cells is not earning "
          "its column.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
