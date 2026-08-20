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

A source column may list fallbacks with `|`, taking the first header that the
chosen metadata file actually has. This covers studies that renamed a column
between AtoMx runs:

        --column 'Case Specific=Case_specific_SORL1|Case_specific'

When an AtoMx export is re-run (e.g. cell typing was redone because the first
pass used the wrong Refit setting), the newer export lands *nested* under the
original run's flatFiles:

    <run>/flatFiles/<rerun>/flatFiles/<slide>/<slide>_metadata_file.csv.gz

Both nested and top-level exports are discovered, and the one that actually
contains the requested source columns is preferred — the seg ID alone cannot
tell them apart, since a re-export carries the same segmentation.

Annotation columns that AtoMx left blank on one study can be filled from
another study of the same physical slide with --fill-from, joining on FOV:

        --fill-from CosMx-Maddie/<other AtoMx run>

Only empty values are filled; existing values are never overwritten.
"""

import argparse
import csv
import gzip
import os
import re
import sys
import tempfile
import zlib
from collections import defaultdict
from dataclasses import dataclass

import boto3
import duckdb
from botocore.exceptions import ClientError

CELL_ID_COLUMN = "cell_ID"
SEG_ID_COLUMN = "cellSegmentationSetId"
FOV_COLUMN = "fov"
LEGACY_OUTPUT_NAME = "cell_type"
LEGACY_COLOR_COLUMN = "hex_color"
SOURCE_FALLBACK_SEPARATOR = "|"
# Column naming AtoMx annotation sheets use for the FOV number, most specific first.
ANNOTATION_FOV_COLUMNS = ("FOVs", "FOV", "fov")
# Column naming the slide an annotation sheet describes, used to catch a wrong file.
ANNOTATION_SLIDE_COLUMNS = ("Flow Cells", "Flow Cell", "slide", "Run_Tissue_name")
# 64 KB of a gzip stream decompresses to far more than one CSV header row.
HEADER_PROBE_BYTES = 1 << 16


@dataclass(frozen=True)
class ColumnSpec:
    """One output annotation column: `out` holds the source column `src`'s value,
    and `color_col` holds its deterministic per-value hex color.

    `src` may name several candidate headers separated by `|`; the first one a
    given metadata file actually has wins, so one command can serve studies that
    renamed a column between AtoMx runs.
    """
    out: str
    src: str
    color_col: str

    @property
    def candidates(self) -> list[str]:
        if not self.src:
            return []
        return [c.strip() for c in self.src.split(SOURCE_FALLBACK_SEPARATOR) if c.strip()]

    def resolve(self, headers) -> str:
        """The first candidate header present in `headers`, else ""."""
        available = set(headers)
        for candidate in self.candidates:
            if candidate in available:
                return candidate
        return ""


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


def read_gz_header(s3, bucket: str, key: str) -> list[str]:
    """Column headers of a gzipped CSV in S3, without downloading the whole object.

    Metadata files run to tens of MB, and source selection only needs the first
    row, so this ranged-GETs the head of the object and decompresses just that.
    Returns [] if the object cannot be read.
    """
    try:
        body = s3.get_object(
            Bucket=bucket, Key=key, Range=f"bytes=0-{HEADER_PROBE_BYTES - 1}",
        )["Body"].read()
    except ClientError:
        return []
    # wbits 16+MAX_WBITS selects gzip framing; a truncated tail is expected here.
    text = zlib.decompressobj(zlib.MAX_WBITS | 16).decompress(body).decode(
        "utf-8", errors="replace")
    first_line, _, _ = text.partition("\n")
    if not first_line:
        return []
    return next(csv.reader([first_line.rstrip("\r")]), [])


def find_metadata_file(
    s3, bucket: str, experiment_prefix: str, slide_name: str
) -> list[tuple[str, str]]:
    """Find all metadata files for a slide across all AtoMx runs.

    Looks both directly under a run's flatFiles and one level deeper, since
    re-running an AtoMx export nests the new flatFiles under the original run:
    ``<run>/flatFiles/<rerun>/flatFiles/<slide>/``. Nested exports are listed
    before their parent so a re-export — which is what carries corrected or
    additional analysis columns — is preferred on ties.

    Returns list of (source_label, s3_key) tuples.
    """
    def metadata_key(prefix: str) -> str:
        return f"{prefix}/flatFiles/{slide_name}/{slide_name}_metadata_file.csv.gz"

    results = []
    for run in s3_ls_prefixes(s3, bucket, experiment_prefix):
        run_prefix = f"{experiment_prefix}/{run}"
        for nested in s3_ls_prefixes(s3, bucket, f"{run_prefix}/flatFiles"):
            if nested == slide_name:
                continue  # the run's own per-slide directory, handled below
            nested_key = metadata_key(f"{run_prefix}/flatFiles/{nested}")
            if s3_file_exists(s3, bucket, nested_key):
                results.append((f"{run}/{nested}", nested_key))
        key = metadata_key(run_prefix)
        if s3_file_exists(s3, bucket, key):
            results.append((run, key))
    return results


def rank_sources(
    s3, bucket: str, sources: list[tuple[str, str]], specs: list[ColumnSpec],
) -> list[tuple[str, str]]:
    """Order sources so ones carrying every requested column come first.

    A re-export and its parent share a segmentation ID, so seg-ID matching alone
    cannot pick between them; what distinguishes them is which analysis columns
    they contain. Sources are otherwise left in discovery order.
    """
    wanted = [spec for spec in specs if spec.candidates]
    if not wanted or len(sources) < 2:
        return sources

    complete, partial = [], []
    for label, key in sources:
        headers = read_gz_header(s3, bucket, key)
        missing = [spec.out for spec in wanted if not spec.resolve(headers)]
        if missing:
            partial.append((label, key))
            print(f"    {label}: missing {', '.join(missing)}")
        else:
            complete.append((label, key))
            print(f"    {label}: has all requested columns")
    return complete + partial


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
) -> tuple[list[tuple[str, str, list[str]]], set[str]]:
    """Read a gzipped metadata CSV, optionally filtering by segmentation ID.

    Returns (rows, seg_ids_seen), where each row is
    (cell_id, fov, [value for each spec, in spec order]). The FOV travels with
    each row so annotations can later be filled in from another study of the
    same slide. A source column absent from this file yields an empty string
    for that spec (with a one-time warning).
    """
    rows: list[tuple[str, str, list[str]]] = []
    seg_ids_seen: set[str] = set()

    with gzip.open(input_path, "rt") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        resolved = [spec.resolve(headers) for spec in specs]
        for spec, source_column in zip(specs, resolved):
            if not source_column:
                print(f"  WARNING: source column '{spec.src}' not found; "
                      f"'{spec.out}' will be empty")
            elif source_column != spec.src:
                print(f"  '{spec.out}' <- {source_column}")

        for row in reader:
            row_seg_id = row.get(SEG_ID_COLUMN, "").strip()
            seg_ids_seen.add(row_seg_id)
            if seg_id_set is not None and row_seg_id not in seg_id_set:
                continue
            cell_id = row.get("cell_id", "")
            fov = row.get(FOV_COLUMN, "").strip()
            values = [row.get(col, "") if col else "" for col in resolved]
            rows.append((cell_id, fov, values))

    return rows, seg_ids_seen


def fov_annotation_map(
    rows: list[tuple[str, str, list[str]]], specs: list[ColumnSpec],
) -> tuple[dict[int, dict[str, str]], list[str]]:
    """Collapse per-cell rows into {spec index: {fov: value}}, FOV-constant only.

    An FOV join is only meaningful for columns that describe the *tissue* rather
    than the cell: AtoMx case annotations are per-case, a case occupies a
    contiguous block of FOVs, so every cell in an FOV carries the same value.
    Per-cell columns (cell typing, QC flags) vary within an FOV and belong to
    the segmentation that produced them — transferring those across studies
    would invent data.

    So a column qualifies only if every FOV has at most one distinct non-empty
    value. Returns (fillable, skipped) where skipped names the columns that
    varied within an FOV.
    """
    values_by_column: list[dict[str, set[str]]] = [defaultdict(set) for _ in specs]
    for _cell_id, fov, values in rows:
        if not fov:
            continue
        for index, value in enumerate(values):
            if value:
                values_by_column[index][fov].add(value)

    fillable: dict[int, dict[str, str]] = {}
    skipped: list[str] = []
    for index, spec in enumerate(specs):
        per_fov = values_by_column[index]
        if not per_fov:
            continue
        if any(len(vals) > 1 for vals in per_fov.values()):
            skipped.append(spec.out)
            continue
        fillable[index] = {fov: next(iter(vals)) for fov, vals in per_fov.items()}
    return fillable, skipped


def fill_missing_annotations(
    rows: list[tuple[str, str, list[str]]],
    donor_rows: list[tuple[str, str, list[str]]],
    specs: list[ColumnSpec],
) -> tuple[dict[str, int], list[str]]:
    """Fill empty annotation values from another study of the same slide, by FOV.

    Two AtoMx studies over one physical slide segment the same FOVs, so cell IDs
    differ but the FOV grid does not. Only FOV-constant columns are transferred
    (see fov_annotation_map), and only into empty values — anything already
    present is left alone.

    Returns ({output column: cells filled}, [columns skipped as per-cell]).
    """
    donor_by_column, skipped = fov_annotation_map(donor_rows, specs)
    filled: dict[str, int] = defaultdict(int)
    for _cell_id, fov, values in rows:
        if not fov:
            continue
        for index, donor_by_fov in donor_by_column.items():
            if values[index]:
                continue
            donor_value = donor_by_fov.get(fov)
            if donor_value:
                values[index] = donor_value
                filled[specs[index].out] += 1
    return dict(filled), skipped


def write_output(
    output_path: str,
    specs: list[ColumnSpec],
    rows: list[tuple[str, str, list[str]]],
) -> dict:
    """Write cell_ID + (value, color) columns for each spec.

    Colors are deterministic per value, built over the full set of collected
    rows so they are stable across slides. Returns a stats dict.
    """
    # One color map per column, over that column's unique non-empty values.
    color_maps: list[dict[str, str]] = []
    for i, _spec in enumerate(specs):
        values = sorted({r[2][i] for r in rows if r[2][i]})
        color_maps.append({v: deterministic_color(v) for v in values})

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    header = [CELL_ID_COLUMN]
    for spec in specs:
        header += [spec.out, spec.color_col]
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for cell_id, _fov, values in rows:
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


def read_annotation_csv(
    path: str, slide_name: str, specs: list[ColumnSpec],
) -> list[tuple[str, str, list[str]]]:
    """Read a per-FOV annotation sheet into the same row shape as metadata.

    These sheets carry one row per FOV rather than per cell — the case-level
    annotations AtoMx should have held. Returning the (cell_id, fov, values)
    shape lets them flow through the same FOV join and the same FOV-constant
    guard as a donor study, with a blank cell_id since there are no cells here.

    Raises ValueError if the sheet names a different slide, which is the likely
    failure when per-slide files are wired up by prefix.
    """
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        fov_column = next((c for c in ANNOTATION_FOV_COLUMNS if c in headers), "")
        if not fov_column:
            raise ValueError(
                f"{path}: no FOV column; expected one of {ANNOTATION_FOV_COLUMNS}, "
                f"got {headers}")
        slide_column = next((c for c in ANNOTATION_SLIDE_COLUMNS if c in headers), "")
        resolved = [spec.resolve(headers) for spec in specs]

        rows: list[tuple[str, str, list[str]]] = []
        slides_seen: set[str] = set()
        for row in reader:
            if slide_column:
                slides_seen.add((row.get(slide_column) or "").strip())
            fov = (row.get(fov_column) or "").strip()
            if not fov:
                continue
            rows.append(("", fov, [(row.get(col) or "") if col else "" for col in resolved]))

    off_slide = {s for s in slides_seen if s and s != slide_name}
    if off_slide:
        raise ValueError(
            f"{path}: describes slide(s) {sorted(off_slide)}, not {slide_name}")

    filled = [spec.out for spec, col in zip(specs, resolved) if col]
    print(f"  Annotation sheet: {len(rows)} FOVs, provides {', '.join(filled) or '(nothing)'}")
    return rows


def load_annotation_rows(
    s3, bucket: str, source: str, slide_name: str,
    specs: list[ColumnSpec], tmpdir: str,
) -> list[tuple[str, str, list[str]]] | None:
    """Load an annotation sheet from a local path or an s3:// URI."""
    local = source
    if source.startswith("s3://"):
        rest = source[5:]
        src_bucket, _, key = rest.partition("/")
        local = os.path.join(tmpdir, "annotations.csv")
        if not s3_download(s3, src_bucket, key, local):
            print(f"  WARNING: could not download annotation sheet {source}")
            return None
    elif not os.path.exists(local):
        print(f"  WARNING: annotation sheet not found: {local}")
        return None
    return read_annotation_csv(local, slide_name, specs)


def load_donor_rows(
    s3, bucket: str, donor_prefix: str, slide_name: str,
    specs: list[ColumnSpec], tmpdir: str,
) -> list[tuple[str, str, list[str]]] | None:
    """Read the same slide's metadata from another study, for annotation filling.

    Segmentation is deliberately not filtered here: the donor is only consulted
    for per-FOV annotation values, which are identical across a slide's cells
    regardless of which segmentation produced them.
    """
    donor_sources = find_metadata_file(s3, bucket, donor_prefix, slide_name)
    if not donor_sources:
        print(f"  WARNING: --fill-from found no metadata for {slide_name} "
              f"under {donor_prefix}")
        return None

    donor_sources = rank_sources(s3, bucket, donor_sources, specs)
    donor_label, donor_key = donor_sources[0]
    donor_gz = os.path.join(tmpdir, "donor_metadata.csv.gz")
    print(f"  Filling blanks from: {donor_label}")
    if not s3_download(s3, bucket, donor_key, donor_gz):
        print(f"  WARNING: --fill-from download failed for s3://{bucket}/{donor_key}")
        return None

    donor_rows, _ = read_rows(donor_gz, None, specs)
    return donor_rows


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
             "(repeatable). Each gets a deterministic <OUTNAME>_color column. "
             "SOURCE_HEADER may list fallbacks separated by '|' — the first header "
             "present in the chosen metadata file wins.",
    )
    parser.add_argument(
        "--annotations-csv",
        default=None,
        metavar="PATH_OR_S3URI",
        help="Per-FOV annotation sheet (one row per FOV) to fill blank annotation "
             "values from. Use when AtoMx never captured the annotations. "
             "Applied before --fill-from.",
    )
    parser.add_argument(
        "--fill-from",
        default=None,
        metavar="EXPERIMENT_PREFIX",
        help="Fill annotation values AtoMx left blank from another study of the same "
             "physical slide, joining on FOV. Only empty values are filled.",
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

    # A re-export shares its parent's segmentation ID, so prefer whichever
    # source actually carries the requested columns.
    sources = rank_sources(s3, args.bucket, sources, specs)

    seg_id_set = (
        set(s.strip() for s in args.seg_id.split(",")) if args.seg_id else None
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        local_gz = os.path.join(tmpdir, "metadata.csv.gz")

        def _report_fill(source, filled, skipped):
            if skipped:
                print(f"  Not filled (varies within an FOV, so per-cell not "
                      f"per-case): {', '.join(skipped)}")
            if filled:
                print(f"  Filled from {source}: "
                      + ", ".join(f"{col} ({n} cells)" for col, n in filled.items()))
            else:
                print(f"  Nothing to fill from {source}")

        def apply_fill(rows):
            """Fill blank annotations, sheet first then donor study.

            An explicit --annotations-csv is the more authoritative source, so it
            goes first; --fill-from then covers anything it left blank. Both only
            ever write into empty values.
            """
            if args.annotations_csv:
                sheet_rows = load_annotation_rows(
                    s3, args.bucket, args.annotations_csv, args.slide_name,
                    specs, tmpdir,
                )
                if sheet_rows:
                    filled, skipped = fill_missing_annotations(rows, sheet_rows, specs)
                    _report_fill(args.annotations_csv, filled, skipped)
            if not args.fill_from:
                return
            donor_rows = load_donor_rows(
                s3, args.bucket, args.fill_from, args.slide_name, specs, tmpdir,
            )
            if donor_rows is None:
                return
            filled, skipped = fill_missing_annotations(rows, donor_rows, specs)
            _report_fill(args.fill_from, filled, skipped)

        if seg_id_set is None:
            # No filtering — use the first source.
            run_name, s3_key = sources[0]
            print(f"  Downloading from: {run_name}")
            if not s3_download(s3, args.bucket, s3_key, local_gz):
                print(f"ERROR: Failed to download s3://{args.bucket}/{s3_key}", file=sys.stderr)
                sys.exit(1)
            specs = resolve_legacy_src(specs, local_gz)
            rows, _ = read_rows(local_gz, None, specs)
            apply_fill(rows)
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

        apply_fill(all_rows)
        stats = write_output(args.output, specs, all_rows)
        print(f"  Merged metadata: {stats['total_written']} cells, "
              f"value counts {stats['value_counts']} "
              f"from {len(seg_ids_found)} segmentation(s)")


if __name__ == "__main__":
    main()
