#!/usr/bin/env python3
"""Per-(donor, region) QC from a cohort AnnData — is ONE region of ONE donor usable?

The slide-level QC table cannot answer this. A slide carrying both eye and brain averages
them, so a donor whose brain is degraded and whose eye is fine looks unremarkable overall:
uwa7761eyes sits mid-pack at 164 median counts while its brain is the reason 100 of its
FOVs are excluded. The decision in front of us is per region, so the table has to be too.

Reads obs only (backed), so it is light at cohort scale.

The comparison that matters is a region against ITSELF across donors, never against the
cohort as a whole. Optic nerve is sparser than retina everywhere, so a donor's optic nerve
measured against an all-region median would look degraded in every donor, and measured
against its own slide's median would look degraded wherever the eye is good. Both are the
same mistake. `--focus-slide` therefore reports each of that donor's regions against the
other donors' cells in THAT region, and gives its rank among them.

Region comes from the FOV annotation sheet, joined on (slide_id, fov), not from obs —
obs['Region'] carries whatever the upstream export said, which is what the SME's rulings
correct. Rows the sheet marks exclude=1 are reported too, under their true region with
`excluded` set: they are the control. A donor whose excluded brain is visibly worse than
its kept eye is evidence the exclusion was drawn in the right place; one where both look
alike means the exclusion is either too wide or too narrow.

Reads:
  --combined-h5ad   stage-3 combined_qc.h5ad (obs: slide_id, fov, qc_gene_counts,
                    qc_genes_detected, qc_negprobe_prop).
  --fov-annotations fov_annotations.csv (slide_id, fov, region, exclude).
Writes (--output-dir):
  qc_by_region.csv        every (slide, region): n, count/gene quantiles, floor failures
  qc_region_summary.csv   per region: the cohort distribution each donor is judged against
  qc_focus_<slide>.csv    --focus-slide only: that donor's regions vs the other donors',
                          carrying both depth (counts) and complexity (genes per 100 counts)
                          because a region can be normal on one and not the other

Usage:
    python pipeline/python/qc_by_region.py \\
        --combined-h5ad combined_qc.h5ad \\
        --fov-annotations pipeline/reference/fov_annotations.csv \\
        --focus-slide uwa7761eyes --output-dir out
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from fov_exclusions import FOV_EXCLUDE_COL, FOV_FOV_COL, FOV_SLIDE_COL, FOV_COL, SLIDE_COL

REGION_COL = "region"
COUNT_COL = "qc_gene_counts"
GENES_COL = "qc_genes_detected"
NEGPROP_COL = "qc_negprobe_prop"
# Genes detected per 100 panel counts -- the measure that separates "fewer molecules"
# from "less distinguishable molecules". See describe().
COMPLEXITY_COL = "genes_per_100_counts"
# The Stage-3a floor this cohort was QC'd at, plus the candidates a re-QC would consider.
DEFAULT_FLOORS = (50, 100, 150, 200)
QUANTILES = (0.25, 0.5, 0.75)
# A donor region this small cannot be judged against anything; reported, never ranked.
MIN_CELLS_TO_RANK = 200


def load_regions(annotations_csv: Path) -> pd.DataFrame:
    """The (slide, FOV) -> (region, excluded) sheet, validated."""
    sheet = pd.read_csv(annotations_csv, dtype={FOV_SLIDE_COL: str})
    required = (FOV_SLIDE_COL, FOV_FOV_COL, REGION_COL, FOV_EXCLUDE_COL)
    missing = [c for c in required if c not in sheet.columns]
    if missing:
        sys.exit(f"ERROR: {annotations_csv} has no {missing} column(s)")
    return pd.DataFrame({
        FOV_SLIDE_COL: sheet[FOV_SLIDE_COL].astype(str).str.strip(),
        FOV_FOV_COL: sheet[FOV_FOV_COL].astype(int),
        REGION_COL: sheet[REGION_COL].astype(str).str.strip(),
        "excluded": sheet[FOV_EXCLUDE_COL].astype(int) == 1,
    })


def attach_regions(obs: pd.DataFrame, regions: pd.DataFrame,
                   require: tuple[str, ...] = (COUNT_COL,)) -> pd.DataFrame:
    """obs with region/excluded joined on (slide, FOV). Exits if the join matches nothing.

    A join that matches nothing means the slide ids or FOV numbering disagree, and every
    table downstream would then be empty or silently partial -- the same failure mode
    fov_exclusions.py guards, for the same reason.

    `require` is what the CALLER needs beyond the join keys: QC needs the counts column, a
    caller counting labels needs its label column instead. The join itself needs neither.
    """
    missing = [c for c in (SLIDE_COL, FOV_COL, *require) if c not in obs.columns]
    if missing:
        sys.exit(f"ERROR: obs is missing {missing}; is this a stage-3 combined_qc.h5ad?")

    keyed = pd.DataFrame({
        FOV_SLIDE_COL: obs[SLIDE_COL].astype(str).str.strip(),
        FOV_FOV_COL: obs[FOV_COL].astype(int),
    }, index=obs.index)
    joined = keyed.merge(regions, on=[FOV_SLIDE_COL, FOV_FOV_COL], how="left")
    joined.index = obs.index

    matched = joined[REGION_COL].notna()
    if not matched.any():
        sys.exit(
            f"ERROR: no cell matched any (slide, FOV) in the annotation sheet.\n"
            f"  obs slides:   {sorted(keyed[FOV_SLIDE_COL].unique())[:6]}\n"
            f"  sheet slides: {sorted(regions[FOV_SLIDE_COL].unique())[:6]}")
    if not matched.all():
        unmatched = keyed.loc[~matched, FOV_SLIDE_COL].value_counts()
        print(f"WARN: {int((~matched).sum()):,} cells have no annotated FOV and are dropped "
              f"from this table; worst slides: {dict(unmatched.head(3))}", file=sys.stderr)

    out = obs.loc[matched].copy()
    out[REGION_COL] = joined.loc[matched, REGION_COL].to_numpy()
    # Cast, do not inherit: the left merge above leaves NaN in `excluded` wherever the sheet
    # has no row, which makes the whole column object dtype. `~` on an object column is
    # bitwise negation -- -1, -2 -- and .loc then reads those as row LABELS and raises.
    # Matched rows are never NaN, so this is safe as well as necessary.
    out["excluded"] = joined.loc[matched, "excluded"].to_numpy(dtype=bool)
    return out


def describe(group: pd.DataFrame, floors: tuple[int, ...]) -> dict:
    """The QC summary for one set of cells."""
    counts = group[COUNT_COL].to_numpy()
    row = {"n_cells": len(group)}
    for q in QUANTILES:
        row[f"counts_q{int(q * 100)}"] = float(np.quantile(counts, q))
    if GENES_COL in group:
        genes = group[GENES_COL].to_numpy()
        row["genes_median"] = float(np.median(genes))
        # Transcriptional complexity at fixed depth. Degradation shows up HERE before it
        # shows up in counts: a cell can keep its depth and still lose the specificity that
        # separates a mature oligodendrocyte from an OPC, which is what a brain full of
        # implausible OPCs looks like from the typing side. The ratio saturates, so it falls
        # as depth rises -- compare it within ONE region across donors, never across regions.
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(counts > 0, genes / np.maximum(counts, 1) * 100, np.nan)
        row[COMPLEXITY_COL] = float(np.nanmedian(ratio))
    if NEGPROP_COL in group:
        row["negprobe_prop_median"] = float(np.median(group[NEGPROP_COL].to_numpy()))
    for floor in floors:
        row[f"pct_under_{floor}"] = float((counts < floor).mean() * 100)
    return row


def by_slide_region(cells: pd.DataFrame, floors: tuple[int, ...]) -> pd.DataFrame:
    rows = []
    for (slide, region, excluded), group in cells.groupby(
            [SLIDE_COL, REGION_COL, "excluded"], observed=True, sort=True):
        rows.append({SLIDE_COL: slide, REGION_COL: region, "excluded": bool(excluded),
                     **describe(group, floors)})
    return pd.DataFrame(rows)


def by_region(cells: pd.DataFrame, floors: tuple[int, ...]) -> pd.DataFrame:
    """The cohort distribution per region — what an individual donor is judged against.

    Excluded cells are left out: the yardstick has to be the data we intend to keep.
    """
    kept = cells.loc[~cells["excluded"]]
    rows = []
    for region, group in kept.groupby(REGION_COL, observed=True, sort=True):
        rows.append({REGION_COL: region, "n_donors": group[SLIDE_COL].nunique(),
                     **describe(group, floors)})
    return pd.DataFrame(rows)


def focus(cells: pd.DataFrame, slide: str, floors: tuple[int, ...]) -> pd.DataFrame:
    """One donor's regions against the OTHER donors' cells in the same region.

    Holding the region fixed is the whole point: optic nerve is sparser than retina in
    every donor, so any comparison that crosses regions measures anatomy, not quality.
    """
    if slide not in set(cells[SLIDE_COL].astype(str)):
        sys.exit(f"ERROR: --focus-slide '{slide}' is not in this cohort: "
                 f"{sorted(cells[SLIDE_COL].astype(str).unique())}")
    mine = cells[SLIDE_COL].astype(str) == slide

    rows = []
    for (region, excluded), group in cells.loc[mine].groupby(
            [REGION_COL, "excluded"], observed=True, sort=True):
        others = cells.loc[~mine & (cells[REGION_COL] == region) & ~cells["excluded"]]
        mine_median = float(np.median(group[COUNT_COL].to_numpy()))
        row = {REGION_COL: region, "excluded": bool(excluded), "n_cells": len(group),
               "counts_median": mine_median,
               "other_donors": others[SLIDE_COL].nunique(),
               "other_counts_median": float(np.median(others[COUNT_COL].to_numpy()))
                                      if len(others) else float("nan")}
        row["ratio_to_others"] = (row["counts_median"] / row["other_counts_median"]
                                  if row["other_counts_median"] else float("nan"))

        # Depth and complexity can disagree, and which one is off changes the diagnosis:
        # normal counts with low complexity is degradation, low counts with normal
        # complexity is thin tissue or undersegmentation.
        row[COMPLEXITY_COL] = describe(group, ())[COMPLEXITY_COL]
        row[f"other_{COMPLEXITY_COL}"] = (describe(others, ())[COMPLEXITY_COL]
                                          if len(others) else float("nan"))

        # Rank among the per-donor medians, so one donor's extreme tail cannot set the bar.
        per_donor = others.groupby(SLIDE_COL, observed=True)[COUNT_COL].median()
        per_donor = per_donor[
            others.groupby(SLIDE_COL, observed=True)[COUNT_COL].size() >= MIN_CELLS_TO_RANK]
        if len(group) >= MIN_CELLS_TO_RANK and len(per_donor):
            # 1 = the deepest donor in this region.
            row["rank_of"] = len(per_donor) + 1
            row["rank"] = int((per_donor > mine_median).sum()) + 1
        else:
            row["rank"], row["rank_of"] = None, None
        for floor in floors:
            row[f"pct_under_{floor}"] = float((group[COUNT_COL].to_numpy() < floor).mean() * 100)
        rows.append(row)
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-h5ad", type=Path, required=True,
                   help="Stage-3 combined_qc.h5ad; obs is read backed.")
    p.add_argument("--fov-annotations", type=Path, required=True,
                   help="fov_annotations.csv — the region source, joined on (slide_id, fov).")
    p.add_argument("--focus-slide", default=None,
                   help="Report this donor's regions against the other donors' same region.")
    p.add_argument("--floors", default=",".join(str(f) for f in DEFAULT_FLOORS),
                   help="Comma-separated panel-count floors for the failure columns.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    floors = tuple(int(f) for f in args.floors.split(",") if f.strip())

    adata = ad.read_h5ad(args.combined_h5ad, backed="r")
    obs = adata.obs
    print(f"Read obs for {len(obs):,} cells from {args.combined_h5ad}")

    cells = attach_regions(obs, load_regions(args.fov_annotations))
    print(f"Annotated {len(cells):,} cells across {cells[SLIDE_COL].nunique()} slides "
          f"and {cells[REGION_COL].nunique()} regions "
          f"({int(cells['excluded'].sum()):,} in excluded FOVs)")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_slide = by_slide_region(cells, floors)
    per_slide.to_csv(args.output_dir / "qc_by_region.csv", index=False)
    cohort = by_region(cells, floors)
    cohort.to_csv(args.output_dir / "qc_region_summary.csv", index=False)

    print("\nCohort, by region (excluded FOVs left out):")
    print(cohort.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    if args.focus_slide:
        table = focus(cells, args.focus_slide, floors)
        table.to_csv(args.output_dir / f"qc_focus_{args.focus_slide}.csv", index=False)
        print(f"\n{args.focus_slide}, each region against the OTHER donors' same region:")
        print(table.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
        print("\nrank 1 = the deepest donor in that region. A kept region ranking mid-pack "
              "is usable; one ranking last on the same measure its excluded regions rank "
              "last on is the same problem, not a different one.")

    print(f"\nWrote {args.output_dir}/qc_by_region.csv, qc_region_summary.csv"
          + (f", qc_focus_{args.focus_slide}.csv" if args.focus_slide else ""))


if __name__ == "__main__":
    main()
