#!/usr/bin/env python3
"""Generate _metadata.csv for a single slide, matching the segmentation version
used for stitching.

When a slide has been resegmented, multiple AtoMx runs may exist under the same
experiment, each with its own flatFiles metadata tied to a specific segmentation.
This script finds the correct metadata source by matching the
cellSegmentationSetId UUID and produces a _metadata.csv whose cell IDs align
with the CellLabels TIFFs selected by the stitcher.

The output always begins with a cell_ID column, followed by one or more
categorical annotation columns. Each annotation column is paired with a
`<name>_color` column holding a deterministic per-value hex color, so the same
category gets the same color across every slide (and across annotation columns
that share label names, e.g. two cell-typing runs). The napari-cosmx reader
prefers `<name>_color` for a column and falls back to a generic `hex_color`
column for legacy single-annotation files.

Usage (single cell-type column, legacy — writes cell_type + hex_color):
    uv run python scripts/generate-slide-metadata.py \
        --bucket my-bucket \
        --experiment-prefix CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26 \
        --slide-name UWA7522G2G5Glioblastoma \
        --seg-id 12d18c13-3b25-4cbf-be1a-24d6c24703d5 \
        --output /tmp/slide/output/_metadata.csv

Usage (multiple annotation columns — writes <name> + <name>_color per column):
    uv run python scripts/generate-slide-metadata.py \
        --bucket my-bucket \
        --experiment-prefix CosMx-retina/... \
        --slide-name UWA7575eyes \
        --column celltype_norefit='RNA_RNA_Cell.Typing.InSituType.No.Refit_1_clusters' \
        --column Region=Region \
        --output /tmp/slide/output/_metadata.csv

    When --seg-id is omitted, all cells from the first metadata file found are
    used (backwards-compatible with single-segmentation slides).
"""

import argparse
import csv
import gzip
import os
import re
import sys
import tempfile
from dataclasses import dataclass

import boto3
import duckdb
from botocore.exceptions import ClientError

CELL_ID_COLUMN = "cell_ID"
SEG_ID_COLUMN = "cellSegmentationSetId"
LEGACY_OUTPUT_NAME = "cell_type"
LEGACY_COLOR_COLUMN = "hex_color"


@dataclass(frozen=True)
class ColumnSpec:
    """One output annotation column: `out` holds the source column `src`'s value,
    and `color_col` holds its deterministic per-value hex color."""
    out: str
    src: str
    color_col: str


def s3_ls_prefixes(s3, bucket: str, prefix: str) -> list[str]:
    """List immediate subdirectory prefixes under an S3 path."""
    response = s3.list_objects_v2(
        Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/",
    )
    return [
        p["Prefix"].rstrip("/").rsplit("/", 1)[-1]
        for p in response.get("CommonPrefixes", [])
    ]


def s3_file_exists(s3, bucket: str, key: str) -> bool:
    """Check if an S3 object exists."""
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError:
        return False


def s3_download(s3, bucket: str, key: str, local_path: str) -> bool:
    """Download an S3 object to a local path. Returns True on success."""
    try:
        s3.download_file(bucket, key, local_path)
        return True
    except ClientError:
        return False


def deterministic_color(value: str) -> str:
    """Generate a deterministic hex color from a string using DuckDB's hash()."""
    result = duckdb.sql(
        f"SELECT printf('#%06X', abs(hash($1)) % 16777216)",
        params=[value],
    ).fetchone()
    return result[0]


def find_metadata_file(
    s3, bucket: str, experiment_prefix: str, slide_name: str
) -> list[tuple[str, str]]:
    """Find all metadata files for a slide across all AtoMx runs.

    Returns list of (atomx_run_name, s3_key) tuples.
    """
    atomx_runs = s3_ls_prefixes(s3, bucket, experiment_prefix)
    results = []
    for run in atomx_runs:
        key = f"{experiment_prefix}/{run}/flatFiles/{slide_name}/{slide_name}_metadata_file.csv.gz"
        if s3_file_exists(s3, bucket, key):
            results.append((run, key))
    return results


CELL_TYPE_COLUMN_PATTERN = re.compile(r"RNA_RNA_Cell\.Typing\.InSituType\..*_clusters$")


def _detect_cell_type_column(headers: list[str]) -> str | None:
    """Auto-detect the cell type column from CSV headers."""
    matches = [h for h in headers if CELL_TYPE_COLUMN_PATTERN.match(h)]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        print(f"  Multiple InSituType columns found: {matches}")
        print(f"  Using first match: {matches[0]}")
        return matches[0]
    return None


def read_rows(
    input_path: str,
    seg_id_set: set[str] | None,
    specs: list[ColumnSpec],
) -> tuple[list[tuple[str, list[str]]], set[str]]:
    """Read a gzipped metadata CSV, optionally filtering by segmentation ID.

    Returns (rows, seg_ids_seen), where each row is
    (cell_id, [value for each spec, in spec order]). A source column absent from
    this file yields an empty string for that spec (with a one-time warning).
    """
    rows: list[tuple[str, list[str]]] = []
    seg_ids_seen: set[str] = set()

    with gzip.open(input_path, "rt") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        missing = [spec for spec in specs if spec.src not in headers]
        for spec in missing:
            print(f"  WARNING: source column '{spec.src}' not found; "
                  f"'{spec.out}' will be empty")

        for row in reader:
            row_seg_id = row.get(SEG_ID_COLUMN, "").strip()
            seg_ids_seen.add(row_seg_id)
            if seg_id_set is not None and row_seg_id not in seg_id_set:
                continue
            cell_id = row.get("cell_id", "")
            values = [row.get(spec.src, "") for spec in specs]
            rows.append((cell_id, values))

    return rows, seg_ids_seen


def write_output(
    output_path: str,
    specs: list[ColumnSpec],
    rows: list[tuple[str, list[str]]],
) -> dict:
    """Write cell_ID + (value, color) columns for each spec.

    Colors are deterministic per value, built over the full set of collected
    rows so they are stable across slides. Returns a stats dict.
    """
    # One color map per column, over that column's unique non-empty values.
    color_maps: list[dict[str, str]] = []
    for i, _spec in enumerate(specs):
        values = sorted({r[1][i] for r in rows if r[1][i]})
        color_maps.append({v: deterministic_color(v) for v in values})

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    header = [CELL_ID_COLUMN]
    for spec in specs:
        header += [spec.out, spec.color_col]
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for cell_id, values in rows:
            out_row = [cell_id]
            for i, spec in enumerate(specs):
                value = values[i]
                out_row += [value, color_maps[i].get(value, "")]
            writer.writerow(out_row)

    return {
        "total_written": len(rows),
        "value_counts": {spec.out: len(color_maps[i]) for i, spec in enumerate(specs)},
    }


def build_specs(args) -> list[ColumnSpec]:
    """Build the ordered list of output columns from CLI args.

    --column OUT=SRC (repeatable) defines the new multi-column format, with each
    column paired to a `<OUT>_color`. When no --column is given, fall back to the
    legacy single cell_type column paired to `hex_color` (source column from
    --cell-type-column, else auto-detected from the file).
    """
    if args.column:
        specs = []
        for item in args.column:
            if "=" not in item:
                print(f"ERROR: --column must be OUTNAME=SOURCE_HEADER, got '{item}'",
                      file=sys.stderr)
                sys.exit(2)
            out, src = item.split("=", 1)
            out, src = out.strip(), src.strip()
            specs.append(ColumnSpec(out=out, src=src, color_col=f"{out}_color"))
        return specs
    # Legacy single-column mode; src resolved per-file if not given (None here).
    return [ColumnSpec(out=LEGACY_OUTPUT_NAME, src=args.cell_type_column,
                       color_col=LEGACY_COLOR_COLUMN)]


def resolve_legacy_src(specs: list[ColumnSpec], input_path: str) -> list[ColumnSpec]:
    """For legacy mode with no explicit --cell-type-column, auto-detect the source
    column from the file header. No-op when a source is already set."""
    if len(specs) == 1 and specs[0].src is None:
        with gzip.open(input_path, "rt") as f:
            reader = csv.DictReader(f)
            detected = _detect_cell_type_column(reader.fieldnames or [])
        if detected:
            print(f"  Auto-detected cell type column: {detected}")
        else:
            print("  WARNING: No InSituType column found, cell_type will be empty")
            detected = ""
        return [ColumnSpec(out=specs[0].out, src=detected, color_col=specs[0].color_col)]
    return specs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate _metadata.csv for a slide, matching the segmentation version.",
    )
    parser.add_argument(
        "--bucket", required=True,
        help="S3 bucket name",
    )
    parser.add_argument(
        "--experiment-prefix", required=True,
        help="S3 prefix for the experiment (e.g. CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26)",
    )
    parser.add_argument(
        "--slide-name", required=True,
        help="Slide name (e.g. UWA7522G2G5Glioblastoma)",
    )
    parser.add_argument(
        "--seg-id", default=None,
        help="cellSegmentationSetId UUID to filter by (extracted from CellLabels subdir name). "
             "When omitted, all cells are included.",
    )
    parser.add_argument(
        "--cell-type-column",
        default=None,
        help="Legacy single-column mode: source column for cell type annotations. "
             "Auto-detected from InSituType columns if omitted. Ignored when --column is used.",
    )
    parser.add_argument(
        "--column",
        action="append",
        default=None,
        metavar="OUTNAME=SOURCE_HEADER",
        help="Add an output annotation column named OUTNAME sourced from SOURCE_HEADER "
             "(repeatable). Each gets a deterministic <OUTNAME>_color column.",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output path for _metadata.csv",
    )
    args = parser.parse_args()

    specs = build_specs(args)
    print(f"  Slide:     {args.slide_name}")
    print(f"  Seg ID:    {args.seg_id or '(all)'}")
    print(f"  Columns:   {', '.join(f'{s.out}<-{s.src}' for s in specs)}")

    s3 = boto3.client("s3")

    # Find all metadata files for this slide
    sources = find_metadata_file(
        s3, args.bucket, args.experiment_prefix, args.slide_name,
    )

    if not sources:
        print(f"ERROR: No metadata files found for slide {args.slide_name}", file=sys.stderr)
        print(f"  Searched: s3://{args.bucket}/{args.experiment_prefix}/*/flatFiles/{args.slide_name}/", file=sys.stderr)
        sys.exit(1)

    print(f"  Found {len(sources)} metadata source(s):")
    for run_name, key in sources:
        print(f"    - {run_name}")

    seg_id_set = (
        set(s.strip() for s in args.seg_id.split(",")) if args.seg_id else None
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        local_gz = os.path.join(tmpdir, "metadata.csv.gz")

        if seg_id_set is None:
            # No filtering — use the first source.
            run_name, s3_key = sources[0]
            print(f"  Downloading from: {run_name}")
            if not s3_download(s3, args.bucket, s3_key, local_gz):
                print(f"ERROR: Failed to download s3://{args.bucket}/{s3_key}", file=sys.stderr)
                sys.exit(1)
            specs = resolve_legacy_src(specs, local_gz)
            rows, _ = read_rows(local_gz, None, specs)
            stats = write_output(args.output, specs, rows)
            print(f"  Generated: {stats['total_written']} cells, "
                  f"value counts {stats['value_counts']}")
            return

        # Multiple seg IDs may span multiple AtoMx runs (e.g. two-step
        # resegmentation). Collect matching rows from ALL sources, then write a
        # single merged _metadata.csv with colors built over the full set.
        seg_ids_requested = set(seg_id_set)
        seg_ids_found: set[str] = set()
        all_rows: list[tuple[str, list[str]]] = []
        all_seg_ids_seen: set[str] = set()
        specs_resolved = False

        for run_name, s3_key in sources:
            print(f"  Trying: {run_name} ...")
            if not s3_download(s3, args.bucket, s3_key, local_gz):
                print(f"    Download failed, skipping")
                continue
            if not specs_resolved:
                specs = resolve_legacy_src(specs, local_gz)
                specs_resolved = True

            rows, seg_ids_seen = read_rows(local_gz, seg_id_set, specs)
            all_seg_ids_seen.update(seg_ids_seen)

            if rows:
                all_rows += rows
                matched_ids = seg_ids_requested & seg_ids_seen
                seg_ids_found.update(matched_ids)
                print(f"    Matched {len(rows)} cells (seg IDs: {matched_ids})")
                if seg_ids_found >= seg_ids_requested:
                    break
            else:
                print(f"    No matching cells (seg IDs in file: {seg_ids_seen})")

        if not all_rows:
            print(f"ERROR: No metadata source contains cells for "
                  f"segmentation IDs {seg_ids_requested}", file=sys.stderr)
            print(f"  Seg IDs found across all sources: {all_seg_ids_seen}",
                  file=sys.stderr)
            sys.exit(1)

        stats = write_output(args.output, specs, all_rows)
        print(f"  Merged metadata: {stats['total_written']} cells, "
              f"value counts {stats['value_counts']} "
              f"from {len(seg_ids_found)} segmentation(s)")


if __name__ == "__main__":
    main()
