#!/usr/bin/env python3
"""Turn Maddie's combined SORL1 annotation sheet into the pipeline FOV table.

Each slide carries four donors in contiguous FOV blocks, and the AtoMx export leaves
`UWA` / `Case_*` / `SORL1_mutation` blank on 9 of the 10 slides. This reads the curated
sheet (one row per FOV) and writes the schema `flatfiles_to_anndata.py --fov-annotations`
reads, so stage 1 can stamp the donor onto every cell.

`region` is set to the donor id: one MTG section per donor, so donor *is* the tissue
block, and it is the level at which RNA yield varies. The donor / case columns ride
alongside for the same reason the retina table carries `mixed_adjacent` and `exclude`.

The sheet names the new slides with spaces ("6562 A7 576 A6 ..."); the export directories
drop them, so `Flow Cells` is matched to the manifest after removing whitespace. The
output is checked against the manifest: every slide present, every FOV 1..N covered
exactly once -- a slide the sheet misses would silently keep a blank `Region`.

Usage:
    build_fov_annotations_sorl1.py --annotations-csv CSV --manifest MANIFEST --output CSV
"""

import argparse
import pathlib
import sys

import pandas as pd

SOURCE_COLUMNS = {
    "Flow Cells": "slide_id",
    "FOVs": "fov",
    "UWA": "donor",
    "Case_broad": "case_broad",
    "Case_specific": "case_group",
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

    Reported rather than rewritten: which column is authoritative is a curation call.
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


def check_against_manifest(table, manifest_slides, fovs_per_slide):
    """Raise unless the table covers exactly the manifest's slides, FOVs 1..N once each."""
    unknown = set(table["slide_id"]) - set(manifest_slides)
    missing = set(manifest_slides) - set(table["slide_id"])
    if unknown or missing:
        raise ValueError(f"slides not in manifest: {sorted(unknown)}; "
                         f"manifest slides absent from sheet: {sorted(missing)}")
    for slide_id, group in table.groupby("slide_id"):
        expected = set(range(1, fovs_per_slide[slide_id] + 1))
        found = set(group["fov"])
        if found != expected:
            raise ValueError(f"{slide_id}: FOVs missing {sorted(expected - found)}, "
                             f"unexpected {sorted(found - expected)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations-csv", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--fov-counts", type=pathlib.Path, required=True,
                        help="CSV with slide_id,n_fovs (FOVs per slide in the flat files)")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.annotations_csv, dtype=str, encoding="utf-8-sig")
    table = table[list(SOURCE_COLUMNS)].rename(columns=SOURCE_COLUMNS)
    table["slide_id"] = table["slide_id"].str.replace(r"\s+", "", regex=True)
    table["fov"] = table["fov"].astype(int)
    table["donor"] = table["donor"].astype(int)

    duplicated = table.duplicated(subset=["slide_id", "fov"])
    if duplicated.any():
        raise ValueError(f"duplicate slide_id/fov rows: "
                         f"{table[duplicated][['slide_id', 'fov']].to_dict('records')}")

    manifest_slides = pd.read_csv(args.manifest)["slide_id"].tolist()
    counts = pd.read_csv(args.fov_counts).set_index("slide_id")["n_fovs"].to_dict()
    check_against_manifest(table, manifest_slides, counts)

    table["region"] = table["donor"].astype(str)
    table["mixed_adjacent"] = MIXED_ADJACENT
    table["exclude"] = EXCLUDE
    report_conflicts(table)

    table = table[OUTPUT_COLUMNS].sort_values(["slide_id", "fov"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    print(f"wrote {args.output}: {len(table)} FOVs, {table.slide_id.nunique()} slides, "
          f"{table.donor.nunique()} donors")
    per_group = table.groupby("case_group")["donor"].nunique().sort_index()
    print("donors per case group:")
    for group, n in per_group.items():
        print(f"  {group:16s} {n}")


if __name__ == "__main__":
    main()
