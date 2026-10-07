#!/usr/bin/env python3
"""Fold Maddie's four per-slide annotation CSVs into one pipeline FOV table.

Each SORL1 pilot slide carries four donors in contiguous FOV blocks, and the AtoMx
export left `UWA` / `Case_broad` / `SORL1_mutation` blank on three of the four slides.
This rebuilds those assignments from the curated per-slide CSVs into the schema
`flatfiles_to_anndata.py --fov-annotations` reads, so stage 1 can stamp donor and
diagnosis onto every cell.

`region` is set to the donor id: one MTG section per donor, so donor *is* the tissue
block, and it is the level at which RNA yield varies (11-fold across this cohort).
The donor / case columns ride alongside for the same reason the retina table carries
`mixed_adjacent` and `exclude`.

Usage:
    build_fov_annotations.py --annotations-dir DIR --output CSV
"""

import argparse
import pathlib
import sys

import pandas as pd

SLIDE_PREFIX = "20260708_UWA_"
SOURCE_COLUMNS = {
    "Flow Cells": "slide_id",
    "FOVs": "fov",
    "UWA": "donor",
    "Case_broad": "case_broad",
    "Case_specific_SORL1": "case_group",
    "SORL1_mutation": "sorl1_mutation",
}
# Written for the existing --fov-annotations reader; this study has no tissue-boundary
# or excluded FOVs, so both flags are constant.
MIXED_ADJACENT = 0
EXCLUDE = 0

OUTPUT_COLUMNS = ["slide_id", "fov", "region", "mixed_adjacent", "exclude",
                  "donor", "case_broad", "case_group", "sorl1_mutation"]


def report_conflicts(table):
    """Flag FOVs whose case_group disagrees with sorl1_mutation.

    A SORL1 carrier's `case_group` should say so; where it does not, the row is
    reported rather than rewritten -- which of the two columns is authoritative is
    a curation decision, not one to make silently here.
    """
    carries_variant = table["sorl1_mutation"].str.contains("SORL1", na=False)
    group_says_variant = table["case_group"].str.contains("SORL1", na=False)
    conflicts = table[carries_variant != group_says_variant]
    if conflicts.empty:
        return
    print(f"WARN: {len(conflicts)} FOV(s) where case_group and sorl1_mutation disagree "
          f"-- written as given, please confirm:", file=sys.stderr)
    for _, row in conflicts.iterrows():
        print(f"  {row.slide_id} FOV {row.fov} donor {row.donor}: "
              f"case_group={row.case_group!r} sorl1_mutation={row.sorl1_mutation!r}",
              file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.annotations_dir.glob("*_annotations.csv"))
    if not paths:
        parser.error(f"no *_annotations.csv under {args.annotations_dir}")

    frames = []
    for path in paths:
        frame = pd.read_csv(path)[list(SOURCE_COLUMNS)].rename(columns=SOURCE_COLUMNS)
        slides = frame["slide_id"].unique()
        if len(slides) != 1:
            raise ValueError(f"{path.name}: expected one slide, found {list(slides)}")
        if not slides[0].startswith(SLIDE_PREFIX):
            raise ValueError(f"{path.name}: unexpected slide id {slides[0]!r}")
        frames.append(frame)

    table = pd.concat(frames, ignore_index=True)
    duplicated = table.duplicated(subset=["slide_id", "fov"])
    if duplicated.any():
        raise ValueError(f"duplicate slide_id/fov rows: "
                         f"{table[duplicated][['slide_id', 'fov']].to_dict('records')}")

    table["donor"] = table["donor"].astype(int)
    table["region"] = table["donor"].astype(str)
    table["mixed_adjacent"] = MIXED_ADJACENT
    table["exclude"] = EXCLUDE
    report_conflicts(table)

    table = table[OUTPUT_COLUMNS].sort_values(["slide_id", "fov"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    print(f"wrote {args.output}: {len(table)} FOVs, {table.slide_id.nunique()} slides, "
          f"{table.donor.nunique()} donors")
    counts = table.groupby("case_group")["donor"].nunique().sort_index()
    print("donors per case group:")
    for group, n in counts.items():
        print(f"  {group:16s} {n}")


if __name__ == "__main__":
    main()
