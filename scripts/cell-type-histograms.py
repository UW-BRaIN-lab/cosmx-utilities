#!/usr/bin/env python3
"""Generate annotation abundance histograms for each CosMx slide from stitched
`_metadata.csv` files in S3.

Reads the multi-column `_metadata.csv` files produced by
generate-slide-metadata.py (each cell has one or more annotation columns, e.g.
`Cell Type`, each paired with a `<column>_color`) and writes one horizontal bar
chart per slide per requested column.

Usage (multi-column metadata — quote names containing spaces):
    uv run python scripts/cell-type-histograms.py \
        --bucket keene-cosmx-data \
        --prefix napari-stitched/CosMx-Maddie/20260720_MV_SORL1_Pilot_1_12_08_2026_11_57_08_503 \
        --column "Cell Type" \
        --output-dir histograms

Usage (legacy single-column metadata — defaults to cell_type + hex_color):
    uv run python scripts/cell-type-histograms.py \
        --bucket keene-cosmx-data \
        --prefix napari-stitched/CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26/Resegmentationcosmxretinabrain22626_01_04_2026_15_10_18_504 \
        --output-dir histograms
"""

import argparse
import csv
import hashlib
import io
import os
import re
from collections import Counter

import boto3
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

METADATA_FILENAME = "_metadata.csv"
LEGACY_COLUMN = "cell_type"
LEGACY_COLOR_COLUMN = "hex_color"
COLOR_COLUMN_SUFFIX = "_color"
UNSLUGGABLE_PATTERN = re.compile(r"[^0-9A-Za-z]+")


class MissingColumnError(Exception):
    """A requested annotation column is absent from a slide's metadata file."""


def deterministic_color(value: str) -> str:
    """Fallback color when a <column>_color column is absent: stable per value."""
    digest = hashlib.md5(value.encode()).hexdigest()
    return f"#{digest[:6]}"


def slugify(value: str) -> str:
    """Filename-safe form of a column name (e.g. 'Cell Type' -> 'Cell_Type')."""
    return UNSLUGGABLE_PATTERN.sub("_", value).strip("_") or "column"


def annotation_columns(headers: list[str]) -> list[str]:
    """Columns that look like annotations, i.e. those paired with a <name>_color."""
    return [h for h in headers if f"{h}{COLOR_COLUMN_SUFFIX}" in headers]


def list_slides(s3, bucket: str, prefix: str) -> list[str]:
    """List slide subdirectories under the given S3 prefix."""
    paginator = s3.get_paginator("list_objects_v2")
    slides = []
    for page in paginator.paginate(
        Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/",
    ):
        slides.extend(
            p["Prefix"].rstrip("/").rsplit("/", 1)[-1]
            for p in page.get("CommonPrefixes", [])
        )
    return slides


def load_metadata(
    s3, bucket: str, prefix: str, slide: str, columns: list[str],
) -> dict[str, tuple[Counter, dict]]:
    """Count values for each requested column in one pass over a slide's metadata.

    Returns {column: (value counts, {value: hex color})}. Colors come from the
    paired `<column>_color`, falling back to a generic `hex_color` column and
    then to a deterministic per-value color.
    """
    key = f"{prefix.rstrip('/')}/{slide}/{METADATA_FILENAME}"
    obj = s3.get_object(Bucket=bucket, Key=key)
    body = obj["Body"].read().decode("utf-8")

    reader = csv.DictReader(io.StringIO(body))
    headers = reader.fieldnames or []
    missing = [column for column in columns if column not in headers]
    if missing:
        available = annotation_columns(headers) or headers
        raise MissingColumnError(
            f"column(s) {missing} not in {key}; available: {available}"
        )

    results: dict[str, tuple[Counter, dict]] = {
        column: (Counter(), {}) for column in columns
    }
    for row in reader:
        for column in columns:
            counts, color_map = results[column]
            value = (row.get(column) or "").strip()
            if not value:
                continue
            counts[value] += 1
            if value not in color_map:
                color_map[value] = (
                    (row.get(f"{column}{COLOR_COLUMN_SUFFIX}") or "").strip()
                    or (row.get(LEGACY_COLOR_COLUMN) or "").strip()
                    or deterministic_color(value)
                )
    return results


def plot_histogram(
    slide: str, column: str, counts: Counter, color_map: dict, output_path: str,
) -> None:
    """Create and save a horizontal bar chart of per-value cell counts."""
    sorted_values = counts.most_common()
    labels = [value for value, _ in sorted_values]
    values = [n for _, n in sorted_values]
    colors = [color_map.get(label, "#4C72B0") for label in labels]

    fig_height = max(4, len(labels) * 0.35)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    ax.barh(labels, values, color=colors, edgecolor="white", linewidth=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Cell Count")
    ax.set_title(f"{column} Distribution — {slide}")
    ax.tick_params(axis="y", labelsize=8)

    for i, v in enumerate(values):
        ax.text(v + max(values) * 0.01, i, str(v), va="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate annotation abundance histograms from stitched "
            "_metadata.csv files in S3."
        ),
    )
    parser.add_argument("--bucket", required=True, help="S3 bucket name")
    parser.add_argument("--prefix", required=True, help="S3 prefix containing slide directories")
    parser.add_argument(
        "--column",
        action="append",
        dest="columns",
        metavar="NAME",
        help=(
            "Annotation column to histogram, repeatable. Quote names containing "
            f"spaces (e.g. --column 'Cell Type'). Default: {LEGACY_COLUMN}."
        ),
    )
    parser.add_argument("--output-dir", default="histograms", help="Local directory for output PNGs")
    args = parser.parse_args()

    columns = args.columns or [LEGACY_COLUMN]

    s3 = boto3.client("s3")
    os.makedirs(args.output_dir, exist_ok=True)

    slides = list_slides(s3, args.bucket, args.prefix)
    print(f"Found {len(slides)} slides; columns: {columns}")

    for slide in sorted(slides):
        print(f"  {slide} ...", end=" ", flush=True)
        try:
            results = load_metadata(s3, args.bucket, args.prefix, slide, columns)
        except MissingColumnError as e:
            print(f"SKIPPED ({e})")
            continue
        except Exception as e:
            print(f"SKIPPED ({type(e).__name__}: {e})")
            continue

        summaries = []
        for column, (counts, color_map) in results.items():
            output_path = os.path.join(
                args.output_dir, f"{slide}_{slugify(column)}.png"
            )
            plot_histogram(slide, column, counts, color_map, output_path)
            summaries.append(
                f"{column}: {sum(counts.values())} cells, {len(counts)} values"
            )
        print("; ".join(summaries))

    print(f"\nHistograms saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
