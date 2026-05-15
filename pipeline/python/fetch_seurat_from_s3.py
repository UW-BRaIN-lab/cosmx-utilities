#!/usr/bin/env python3
"""Download one slide's Seurat .RDS from AWS S3 to a local path.

Used by the one-time clinical-annotation extraction job. The .RDS files live
on AWS S3 (never copied to Kopah); this script locates the matching file by
slide_id and downloads it to a path the R script can read.

AtoMx exports name files like `seuratObject_7134.A7.7439.A6.RDS` — stripping
dots from the suffix after `seuratObject_` yields the slide_id.

Usage:
    uv run python pipeline/python/fetch_seurat_from_s3.py \\
        --slide-id 7134A77439A6 \\
        --batch-prefix CosMx-GBM/Wenyu-cohort/Wenyupart142926.../ \\
        --output /path/to/local/7134A77439A6.RDS
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import boto3
from dotenv import load_dotenv

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"


def make_source_client():
    access_key = os.environ.get("AWS_SOURCE_ACCESS_KEY_ID")
    secret_key = os.environ.get("AWS_SOURCE_SECRET_ACCESS_KEY")
    if access_key and secret_key:
        return boto3.client(
            "s3",
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
        )
    profile = os.environ.get("AWS_SOURCE_PROFILE")
    if profile:
        return boto3.Session(profile_name=profile).client("s3")
    return boto3.client("s3")


def find_seurat_key(s3, bucket: str, batch_prefix: str, slide_id: str) -> str | None:
    if not batch_prefix.endswith("/"):
        batch_prefix += "/"
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=batch_prefix, Delimiter="/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            filename = key.rsplit("/", 1)[-1]
            if not filename.startswith("seuratObject_") or not filename.endswith(".RDS"):
                continue
            stripped = filename[len("seuratObject_"):-len(".RDS")].replace(".", "")
            if stripped == slide_id:
                return key
    return None


def main() -> None:
    if ENV_PATH.exists():
        load_dotenv(ENV_PATH)

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--slide-id", required=True)
    p.add_argument(
        "--batch-prefix", required=True,
        help="S3 prefix for the export batch (the directory containing flatFiles/, DecodedFiles/, seuratObject_*.RDS)",
    )
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--source-bucket", default=os.environ.get("SOURCE_S3_BUCKET"),
        help="AWS S3 bucket holding the source data (env: SOURCE_S3_BUCKET)",
    )
    args = p.parse_args()

    if not args.source_bucket:
        print("ERROR: --source-bucket required (or set SOURCE_S3_BUCKET)", file=sys.stderr)
        sys.exit(1)

    s3 = make_source_client()
    key = find_seurat_key(s3, args.source_bucket, args.batch_prefix, args.slide_id)
    if not key:
        print(
            f"ERROR: no seuratObject_*.RDS at s3://{args.source_bucket}/{args.batch_prefix} "
            f"matching slide_id={args.slide_id}",
            file=sys.stderr,
        )
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading s3://{args.source_bucket}/{key} -> {args.output}")
    s3.download_file(args.source_bucket, key, str(args.output))
    print(f"Done ({args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
