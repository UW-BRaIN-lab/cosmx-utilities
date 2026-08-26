#!/usr/bin/env python3
"""Cross-tabulate AtoMx's Leiden clustering against its InSituType calls, from
CosMx flat files, and emit counts cross-tabs ready for plot_crosstab_sankey.py.

AtoMx exports both an unsupervised Leiden clustering (built in expression space
by its Neighbor-network module) and a supervised InSituType call for every cell,
side by side in `*_metadata_file.csv.gz`. Crossing them answers "does the
reference-free clustering agree with the reference-based typing, and where does
it disagree?" without any pipeline run -- the stitched Napari `_metadata.csv`
carries the typing but not the Leiden, so this reads the flat files directly.

Two AtoMx quirks this handles, both of which silently pick the wrong run:

  - A study can hold a *nested* re-export (`<study>/flatFiles/<rerun>/flatFiles/...`)
    that supersedes the outer flat files. The deepest, newest file per slide wins.
  - A re-export can carry several typing runs (InSituType.1_1, .2_1, ...). The
    highest run index wins -- header order is not a recency signal.

Writes, per study:
  <slide>_leiden_x_insitutype.csv   per-slide counts contingency (Leiden x type)
  pooled_leiden_x_insitutype.csv    all slides pooled
  leiden_purity.csv                 per-cluster dominant type + its share

Render any of those with:
    uv run python pipeline/python/plot_crosstab_sankey.py \\
        --crosstab <dir>/pooled_leiden_x_insitutype.csv \\
        --label-left "Leiden (AtoMx)" --label-right "InSituType (AtoMx)" \\
        --sort-left natural --output <dir>/pooled_sankey.png

Usage:
    uv run python scripts/atomx-leiden-typing-crosstab.py \\
        --bucket keene-cosmx-data \\
        --prefix CosMx-Maddie/20260813_MV_SORL1_pilot_3D_resegmentation_18_08_2026_9_51_12_114 \\
        --output-dir crosstabs/maddie-3D
"""

import argparse
import gzip
import io
import os
import re

import boto3
import pandas as pd

METADATA_SUFFIX = "_metadata_file.csv.gz"
FLATFILES_SEGMENT = "flatFiles"
TYPING_PATTERN = re.compile(
    r"^RNA_RNA_Cell\.Typing\.InSituType\.(\d+)_(\d+)_clusters$"
)
LEIDEN_PATTERN = re.compile(r"Leiden\.Clustering")
UNASSIGNED_LABEL = "Unassigned"


class ColumnDetectionError(Exception):
    """A slide's flat file lacks a Leiden or InSituType column."""


def newest_metadata_per_slide(s3, bucket: str, prefix: str) -> dict[str, str]:
    """Map slide name -> key of its authoritative metadata file.

    A nested re-export supersedes the outer flat files, so among the candidates
    for one slide the deepest key wins, and among equal depths the newest.
    """
    paginator = s3.get_paginator("list_objects_v2")
    candidates: dict[str, list[tuple[int, object, str]]] = {}
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix.rstrip("/") + "/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not key.endswith(METADATA_SUFFIX):
                continue
            slide = key.rsplit("/", 1)[0].rsplit("/", 1)[-1]
            depth = key.count(f"/{FLATFILES_SEGMENT}/")
            candidates.setdefault(slide, []).append(
                (depth, obj["LastModified"], key)
            )
    return {
        slide: max(entries)[2] for slide, entries in candidates.items()
    }


def detect_columns(
    headers: list[str],
    typing_version: str | None = None,
    leiden_column: str | None = None,
) -> tuple[str, str]:
    """Return (leiden column, InSituType column) for a flat file.

    With no --typing-version, the run with the highest (major, minor) index wins;
    AtoMx orders headers arbitrarily, so position must not decide which run is
    current. Pass a version (e.g. "2_1") to pin one run across a whole study.
    """
    leiden = ([leiden_column] if leiden_column
              else [h for h in headers if LEIDEN_PATTERN.search(h)])
    typing = [(m, h) for m, h in
              ((TYPING_PATTERN.match(h), h) for h in headers) if m]
    if not leiden or leiden[0] not in headers:
        raise ColumnDetectionError(
            f"Leiden column {leiden[0]!r} absent" if leiden_column
            else "no Leiden.Clustering column")
    if not typing:
        raise ColumnDetectionError("no InSituType *_clusters column")

    if typing_version:
        wanted = f"RNA_RNA_Cell.Typing.InSituType.{typing_version}_clusters"
        if wanted not in headers:
            found = sorted(f"{m.group(1)}_{m.group(2)}" for m, _ in typing)
            raise ColumnDetectionError(
                f"InSituType version {typing_version!r} absent; present: {found}")
        return leiden[0], wanted

    newest = max(typing, key=lambda mh: (int(mh[0].group(1)), int(mh[0].group(2))))
    return leiden[0], newest[1]


def read_slide(
    s3, bucket: str, key: str,
    typing_version: str | None = None,
    leiden_column: str | None = None,
) -> tuple[pd.DataFrame, str, str]:
    """Read just the Leiden and chosen-typing columns from one slide's flat file."""
    body = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
    with gzip.open(io.BytesIO(body), "rt") as f:
        headers = f.readline().rstrip("\n").split(",")
    leiden_col, typing_col = detect_columns(headers, typing_version, leiden_column)
    with gzip.open(io.BytesIO(body), "rt") as f:
        df = pd.read_csv(f, usecols=[leiden_col, typing_col], low_memory=False)
    df = df.rename(columns={leiden_col: "leiden", typing_col: "cell_type"})
    df["leiden"] = df["leiden"].astype(str)
    df["cell_type"] = df["cell_type"].fillna(UNASSIGNED_LABEL)
    return df, leiden_col, typing_col


def leiden_purity(counts: pd.DataFrame) -> pd.DataFrame:
    """Per Leiden cluster: size, dominant InSituType, and that type's share.

    A cluster that is 90% one type says the two labelings agree there; a cluster
    split across many types is where the reference-based call is doing work the
    unsupervised clustering does not support (or vice versa).
    """
    rows = []
    for cluster, row in counts.iterrows():
        total = row.sum()
        dominant = row.idxmax()
        rows.append({
            "leiden": cluster,
            "n_cells": int(total),
            "dominant_type": dominant,
            "dominant_frac": round(float(row.max()) / total, 4) if total else 0.0,
            "n_types_present": int((row > 0).sum()),
        })
    return pd.DataFrame(rows).sort_values(
        "n_cells", ascending=False, ignore_index=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-tabulate AtoMx Leiden vs InSituType from CosMx flat files.",
    )
    parser.add_argument("--bucket", required=True, help="S3 bucket name")
    parser.add_argument("--prefix", required=True, help="S3 study prefix (holds flatFiles/)")
    parser.add_argument("--output-dir", required=True, help="Local directory for the CSVs")
    parser.add_argument(
        "--typing-version", metavar="MAJOR_MINOR",
        help="Pin one InSituType run across the study, e.g. 2_1 for the 3D "
             "resegmentation. Default: the highest run index present per slide.")
    parser.add_argument(
        "--leiden-column",
        help="Exact Leiden header, if a study has more than one clustering.")
    args = parser.parse_args()

    s3 = boto3.client("s3")
    os.makedirs(args.output_dir, exist_ok=True)

    keys = newest_metadata_per_slide(s3, args.bucket, args.prefix)
    if not keys:
        raise SystemExit(f"no {METADATA_SUFFIX} found under s3://{args.bucket}/{args.prefix}")
    print(f"Found {len(keys)} slides")

    pooled = []
    for slide in sorted(keys):
        print(f"  {slide} ...", end=" ", flush=True)
        try:
            df, leiden_col, typing_col = read_slide(
                s3, args.bucket, keys[slide], args.typing_version, args.leiden_column)
        except (ColumnDetectionError, ValueError) as e:
            print(f"SKIPPED ({e})")
            continue

        counts = pd.crosstab(df["leiden"], df["cell_type"])
        counts.to_csv(os.path.join(args.output_dir, f"{slide}_leiden_x_insitutype.csv"))
        pooled.append(df.assign(slide=slide))
        print(
            f"{len(df)} cells, {counts.shape[0]} Leiden x {counts.shape[1]} types "
            f"(typing: {typing_col.split('.InSituType.')[-1]})"
        )

    if not pooled:
        raise SystemExit("no slide yielded a crosstab")

    all_cells = pd.concat(pooled, ignore_index=True)
    pooled_counts = pd.crosstab(all_cells["leiden"], all_cells["cell_type"])
    pooled_counts.to_csv(os.path.join(args.output_dir, "pooled_leiden_x_insitutype.csv"))

    purity = leiden_purity(pooled_counts)
    purity.to_csv(os.path.join(args.output_dir, "leiden_purity.csv"), index=False)

    weighted = (purity["dominant_frac"] * purity["n_cells"]).sum() / purity["n_cells"].sum()
    print(f"\nPooled: {len(all_cells)} cells, "
          f"{pooled_counts.shape[0]} Leiden clusters x {pooled_counts.shape[1]} types")
    print(f"Weighted mean cluster purity (dominant type's share): {weighted:.1%}")
    print(f"\nWrote crosstabs to {args.output_dir}/")


if __name__ == "__main__":
    main()
