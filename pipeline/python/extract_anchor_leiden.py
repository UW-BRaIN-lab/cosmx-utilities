#!/usr/bin/env python3
"""The COHORT-WIDE Leiden cluster of every anchor-cohort cell, as a small per-cell CSV.

WHY THIS EXISTS. The anchor cohort's own sidecar, anchor_cells.csv, has a column `slide_cluster`
whose values look like "7104A297104A23|16". That is a PER-SLIDE Leiden used only to stratify the
subsample (13-34 clusters on each of 57 slides, 1,075 in all): cluster 16 on one slide has no
relationship to cluster 16 on another, and none to the cohort-wide Stage-3c clusters the team
talks about (12, 19, 21, ...). Colouring a slide in Napari by it would give every slide its own
private, uncomparable labels -- plausible-looking and wrong.

The cohort-wide clustering lives in the 7.5M-cell typed AnnData's obs, keyed by the same
"<slide>_F<fov>_C<cell>" index. This reads only that obs (backed mode, no matrix) and keeps the
anchor cells.

Two guards, both aimed at that exact mistake:
  * labels containing "|" are refused as per-slide;
  * the report states how many slides each cluster spans. A cohort-wide clustering has clusters
    spanning most slides; a per-slide one has each cluster on exactly one.

Usage:
    uv run python pipeline/python/extract_anchor_leiden.py \\
        --h5ad cosmx_typed.h5ad --cell-ids-csv anchor_cells.csv --output-csv anchor_leiden.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

MIN_OVERLAP = 0.95           # anchor cells are a subset of the typed run, so expect ~1.0
PER_SLIDE_MARK = "|"         # the slide|n separator in anchor_cells.csv's slide_cluster
CELL_ID_RE = re.compile(r"^(?P<slide>.+)_F\d+_C\d+$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True,
                   help="AnnData whose obs carries the cohort-wide Leiden (cosmx_typed.h5ad).")
    p.add_argument("--cell-ids-csv", type=Path, required=True,
                   help="CSV listing the anchor cells (anchor_cells.csv).")
    p.add_argument("--id-column", default="cell_id", help="Cell-id column in --cell-ids-csv.")
    p.add_argument("--leiden-key", default="leiden", help="obs column holding the Leiden cluster.")
    p.add_argument("--output-csv", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad, backed="r")            # obs only; never touches the matrix
    if args.leiden_key not in adata.obs.columns:
        sys.exit(f"ERROR: obs has no column {args.leiden_key!r}. Present: "
                 f"{', '.join(map(str, adata.obs.columns))}")
    leiden = pd.Series(adata.obs[args.leiden_key].to_numpy(), index=adata.obs.index.to_numpy(),
                       name="leiden")
    print(f"{args.h5ad.name}: {len(leiden):,} cells; obs column {args.leiden_key!r}")

    anchor = pd.read_csv(args.cell_ids_csv, usecols=[args.id_column])[args.id_column].astype(str)
    if anchor.duplicated().any():
        sys.exit(f"ERROR: {int(anchor.duplicated().sum()):,} duplicated ids in {args.cell_ids_csv}.")
    print(f"anchor cells listed: {len(anchor):,}")

    found = anchor[anchor.isin(leiden.index)]
    frac = len(found) / max(len(anchor), 1)
    print(f"  found in the AnnData: {len(found):,} ({100 * frac:.2f}%)")
    if frac < MIN_OVERLAP:
        sys.exit(f"ERROR: only {100 * frac:.1f}% of anchor cells are in the AnnData. The id "
                 f"formats differ or these are different runs.\n  anchor example: {anchor.iloc[0]!r}"
                 f"\n  h5ad   example: {leiden.index[0]!r}")

    out = pd.DataFrame({"cell_id": found.to_numpy(), "leiden": leiden.reindex(found).to_numpy()})
    n_missing = int(out["leiden"].isna().sum())
    if n_missing:
        print(f"  NOTE: {n_missing:,} anchor cells have no Leiden value and are dropped.")
        out = out.dropna(subset=["leiden"])
    out["leiden"] = out["leiden"].astype(str)

    # THE GUARD. This is the label-shape check that would have caught slide_cluster.
    per_slide = out["leiden"].str.contains(PER_SLIDE_MARK, regex=False)
    if per_slide.any():
        sys.exit(f"ERROR: {int(per_slide.sum()):,} labels contain {PER_SLIDE_MARK!r} (e.g. "
                 f"{out.loc[per_slide, 'leiden'].iloc[0]!r}); that is a per-slide clustering, "
                 f"not the cohort-wide Leiden this script is for.")

    slides = out["cell_id"].str.extract(CELL_ID_RE)["slide"]
    n_slides = slides.nunique()
    span = out.assign(slide=slides).groupby("leiden")["slide"].nunique().sort_values(ascending=False)
    sizes = out["leiden"].value_counts()
    print(f"\n{out['leiden'].nunique()} clusters across {n_slides} slides. Slides spanned per cluster "
          f"(cohort-wide => most clusters span most slides):")
    print(f"  median {int(span.median())} of {n_slides}; "
          f"{int((span >= 0.5 * n_slides).sum())} of {len(span)} clusters span >= half the slides")
    print(f"  largest: " + ", ".join(f"{k} ({v:,})" for k, v in sizes.head(5).items()))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)
    print(f"\nWrote {args.output_csv} ({len(out):,} cells)")


if __name__ == "__main__":
    main()
