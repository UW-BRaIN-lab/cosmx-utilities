#!/usr/bin/env python3
"""Check whether one cohort's slides were segmented the same way.

Three slides scored at pixel level showed that AtoMx changed its compartment
output mid-cohort: slides segmented 2026-04-16 carry a `cytoplasm` compartment,
slides segmented 2026-05-07 and later do not, from an identical config. If that
update changed only compartment LABELLING it is cosmetic; if it moved cell
BOUNDARIES then cell area shifts with the segmentation date, and every cross-slide
comparison in the cohort carries a batch effect keyed to when segmentation ran
rather than to biology.

Cell geometry is already in the flat-file metadata, so the question is answerable
without touching the images. Per slide this reports the Area, AspectRatio and
Width/Height distributions plus cells-per-FOV, alongside whatever `version` and
`cellSegmentationSetId` the export carries, and then compares the distributions
between groups.

`version` is undocumented in the AtoMx flat-file readme, so it is reported rather
than trusted: if it varies across slides it is offered as the grouping, and
--group-by accepts an external mapping (e.g. segmentation dates read off the AtoMx
UI) when it does not.

Usage:
    uv run python pipeline/python/segmentation_homogeneity.py \
        --flatfiles-dir staged_flatfiles --manifest pipeline/manifest.csv \
        --output segmentation_by_slide.csv --figure segmentation_area.png
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AREA_COLUMN = "Area"
FOV_COLUMN = "fov"
GEOMETRY_COLUMNS = ("Area", "AspectRatio", "Width", "Height")
PROVENANCE_COLUMNS = ("version", "cellSegmentationSetId", "Run_name")
METADATA_SUFFIX_PATTERN = re.compile(r"^(?P<slide>.+)_metadata_file\.csv(?:\.gz)?$")

# A slide whose geometry median sits this far from the cohort median is called out.
# 10% of a linear dimension is roughly 20% of an area, well beyond segmentation noise.
OUTLIER_FRACTION = 0.10
MIN_CELLS = 100


class MissingColumnError(Exception):
    """A structural column is absent from a slide's metadata file."""


def metadata_paths(flatfiles_dir: Path | None,
                   explicit: list[Path]) -> list[tuple[str, Path]]:
    paths = list(explicit)
    if flatfiles_dir:
        paths.extend(sorted(flatfiles_dir.rglob("*_metadata_file.csv*")))
    resolved, seen = [], set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        match = METADATA_SUFFIX_PATTERN.match(path.name)
        if match:
            resolved.append((match.group("slide"), path))
        else:
            print(f"WARN: skipping {path.name}", file=sys.stderr)
    return resolved


def single_value(frame: pd.DataFrame, column: str) -> str:
    """The column's one value, or a joined list when a slide is not internally uniform."""
    if column not in frame.columns:
        return ""
    values = sorted({str(v).strip() for v in frame[column].dropna().unique()})
    if len(values) > 1:
        print(f"  WARN: {column} is not uniform within the slide: {values[:4]}",
              file=sys.stderr)
    return "|".join(values)


def slide_summary(slide_id: str, path: Path) -> dict:
    """Geometry distributions and provenance for one slide."""
    header = pd.read_csv(path, nrows=0).columns
    wanted = [c for c in (FOV_COLUMN, *GEOMETRY_COLUMNS, *PROVENANCE_COLUMNS)
              if c in header]
    if AREA_COLUMN not in wanted or FOV_COLUMN not in wanted:
        raise MissingColumnError(f"{path.name} lacks {AREA_COLUMN} or {FOV_COLUMN}")

    frame = pd.read_csv(path, usecols=wanted)
    frame = frame.loc[frame[AREA_COLUMN] > 0]
    if len(frame) < MIN_CELLS:
        print(f"WARN: {slide_id} has {len(frame)} usable cells, skipping",
              file=sys.stderr)
        return {}

    n_fovs = int(frame[FOV_COLUMN].nunique())
    row = {
        "slide_id": slide_id,
        "n_cells": int(len(frame)),
        "n_fovs": n_fovs,
        "cells_per_fov": round(len(frame) / n_fovs, 1) if n_fovs else float("nan"),
    }
    for column in GEOMETRY_COLUMNS:
        if column not in frame.columns:
            continue
        values = frame[column].astype(float)
        low, median, high = (float(np.percentile(values, q)) for q in (5, 50, 95))
        row[f"{column}_p05"] = round(low, 2)
        row[f"{column}_p50"] = round(median, 2)
        row[f"{column}_p95"] = round(high, 2)
        # Spread relative to the median, so slides of different cell sizes compare.
        row[f"{column}_spread"] = round((high - low) / median, 3) if median else float("nan")
    for column in PROVENANCE_COLUMNS:
        row[column] = single_value(frame, column)
    return row


def choose_grouping(frame: pd.DataFrame, requested: str | None) -> str | None:
    """The column that actually separates slides, preferring an explicit request."""
    if requested:
        if requested not in frame.columns:
            raise MissingColumnError(
                f"--group-by {requested!r} is not a column; have {list(frame.columns)}")
        return requested
    for column in ("version", "cellSegmentationSetId"):
        if column in frame.columns and 1 < frame[column].nunique() < len(frame):
            return column
    return None


def report_groups(frame: pd.DataFrame, grouping: str) -> pd.DataFrame:
    columns = [c for c in frame.columns if c.endswith(("_p50", "_spread"))]
    summary = frame.groupby(grouping)[columns].median()
    summary.insert(0, "n_slides", frame.groupby(grouping).size())
    return summary


def report_outliers(frame: pd.DataFrame) -> pd.DataFrame:
    """Slides whose geometry sits far from the cohort median."""
    flagged = []
    for column in (f"{c}_p50" for c in GEOMETRY_COLUMNS):
        if column not in frame.columns:
            continue
        cohort = frame[column].median()
        if not cohort:
            continue
        deviation = (frame[column] - cohort).abs() / cohort
        for slide, value, dev in zip(frame["slide_id"], frame[column], deviation):
            if dev > OUTLIER_FRACTION:
                flagged.append({"slide_id": slide, "metric": column,
                                "value": value, "cohort_median": round(cohort, 2),
                                "deviation": round(float(dev), 3)})
    return pd.DataFrame(flagged).sort_values("deviation", ascending=False) \
        if flagged else pd.DataFrame()


def write_figure(frame: pd.DataFrame, grouping: str | None, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = frame.sort_values("Area_p50").reset_index(drop=True)
    groups = (list(ordered[grouping].fillna("")) if grouping
              else [""] * len(ordered))
    palette = {g: c for g, c in zip(sorted(set(groups)),
                                    ["#4c72b0", "#dd8452", "#55a868", "#c44e52",
                                     "#8172b3", "#937860"])}

    height = max(4.0, 0.16 * len(ordered))
    fig, axis = plt.subplots(figsize=(9, height))
    positions = np.arange(len(ordered))
    axis.hlines(positions, ordered["Area_p05"], ordered["Area_p95"],
                color="0.75", linewidth=1.2, zorder=1)
    for group in palette:
        mask = [g == group for g in groups]
        axis.scatter(ordered.loc[mask, "Area_p50"], positions[mask],
                     s=26, color=palette[group], zorder=2,
                     label=str(group) if grouping else None)
    axis.axvline(ordered["Area_p50"].median(), color="0.4",
                 linestyle="--", linewidth=0.9)
    axis.set_yticks(positions)
    axis.set_yticklabels(ordered["slide_id"], fontsize=6)
    axis.set_xlabel("Cell area (px) — dot is the median, bar spans p05–p95")
    axis.set_title("Per-slide cell area" + (f", coloured by {grouping}" if grouping else ""))
    if grouping:
        axis.legend(fontsize=7, title=grouping, title_fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--flatfiles-dir", type=Path)
    parser.add_argument("--metadata", type=Path, action="append", default=[])
    parser.add_argument("--manifest", type=Path,
                        help="Joins run_date and instrument_id onto the per-slide table")
    parser.add_argument("--group-by",
                        help="Column to compare across, e.g. 'version' or a column "
                             "brought in with --slide-groups")
    parser.add_argument("--slide-groups", type=Path,
                        help="CSV with slide_id plus any grouping columns — use it to "
                             "supply segmentation dates read off the AtoMx UI, which are "
                             "not present in the flat files")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()

    if not args.flatfiles_dir and not args.metadata:
        parser.error("give --flatfiles-dir or at least one --metadata")

    targets = metadata_paths(args.flatfiles_dir, args.metadata)
    if not targets:
        print("ERROR: no metadata files found", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(targets)} metadata file(s)")

    rows = []
    for slide_id, path in targets:
        print(f"Reading {path.name}")
        summary = slide_summary(slide_id, path)
        if summary:
            rows.append(summary)
    if not rows:
        print("ERROR: no slide produced a summary", file=sys.stderr)
        sys.exit(1)

    frame = pd.DataFrame(rows)
    for extra, key in ((args.manifest, ("run_date", "instrument_id")),
                       (args.slide_groups, None)):
        if extra and extra.exists():
            other = pd.read_csv(extra, dtype={"slide_id": str})
            columns = (["slide_id", *[c for c in key if c in other.columns]] if key
                       else list(other.columns))
            frame = frame.merge(other[columns], on="slide_id", how="left")

    frame.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

    print(f"\nCohort geometry across {len(frame)} slides:")
    for column in (f"{c}_p50" for c in GEOMETRY_COLUMNS):
        if column in frame:
            values = frame[column]
            print(f"  {column:18s} median {values.median():9.2f}   "
                  f"range {values.min():.2f}–{values.max():.2f}   "
                  f"fold {values.max() / values.min():.2f}x"
                  if values.min() else f"  {column}: undefined")

    for column in PROVENANCE_COLUMNS:
        if column in frame:
            distinct = frame[column].nunique()
            note = "" if distinct <= 1 else "   <-- NOT uniform across the cohort"
            print(f"  {column:18s} {distinct} distinct value(s){note}")

    grouping = choose_grouping(frame, args.group_by)
    if grouping:
        print(f"\nGrouped by {grouping}:")
        print(report_groups(frame, grouping).round(3).to_string())
    else:
        print("\nNo column separates the slides into groups. Supply segmentation "
              "dates with --slide-groups + --group-by to test the AtoMx-version "
              "split directly.")

    outliers = report_outliers(frame)
    if len(outliers):
        print(f"\nSlides deviating >{OUTLIER_FRACTION:.0%} from the cohort median:")
        print(outliers.head(20).to_string(index=False))
    else:
        print(f"\nNo slide deviates more than {OUTLIER_FRACTION:.0%} from the cohort "
              f"median on any geometry metric.")

    if args.figure:
        write_figure(frame, grouping, args.figure)


if __name__ == "__main__":
    main()
