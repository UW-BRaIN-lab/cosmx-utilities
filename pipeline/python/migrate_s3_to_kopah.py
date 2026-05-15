#!/usr/bin/env python3
"""Mirror CosMx flat files from AWS S3 to UW Kopah.

For each slide in the manifest, copies:
  - <slide>_exprMat_file.csv.gz
  - <slide>_metadata_file.csv.gz
  - <slide>_fov_positions_file.csv.gz
  - <slide>-polygons.csv.gz
  - optionally <slide>_tx_file.csv.gz (1-2 GB per slide; skipped unless --include-transcripts)

The Seurat .RDS files are deliberately NOT migrated — the clinical annotations
we originally thought were Seurat-only (Case / Block / Region) turned out to be
present in the flat-file metadata too, so the .RDS files are unused by this
pipeline. Flat files + AnnData are the working format throughout.

Two boto3 clients in one process:
  - source: AWS S3 — credentials from AWS_SOURCE_ACCESS_KEY_ID/_SECRET (or
            AWS_SOURCE_PROFILE, or boto3's default credential chain).
  - dest:   Kopah  — credentials from KOPAH_ACCESS_KEY_ID/_SECRET.
Explicit source creds prevent confusion when AWS_* env vars are temporarily
set to Kopah keys (e.g. by the `kopah_brain` shell helper on klone).

The script is idempotent: if the destination object already exists and matches
the source size, it is skipped. Re-run safely after interruptions.

Destination layout on Kopah mirrors the source S3 path under the Kopah bucket,
so the manifest column `flat_files_prefix` refers to the Kopah location directly.

Usage:
    uv run python pipeline/python/migrate_s3_to_kopah.py
    uv run python pipeline/python/migrate_s3_to_kopah.py --include-transcripts
    uv run python pipeline/python/migrate_s3_to_kopah.py --only-slides 7134A77439A6
    uv run python pipeline/python/migrate_s3_to_kopah.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterable

from botocore.exceptions import ClientError
from dotenv import load_dotenv

from _clients import make_kopah_client, make_source_client

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"

# Files we always copy per slide (small to medium-sized, all needed downstream).
ESSENTIAL_FLAT_SUFFIXES = (
    "_exprMat_file.csv.gz",
    "_metadata_file.csv.gz",
    "_fov_positions_file.csv.gz",
    "-polygons.csv.gz",
)
# Heavy file, ~1-2 GB per slide; only copy on demand.
TRANSCRIPTS_SUFFIX = "_tx_file.csv.gz"


def _env(key: str) -> str:
    value = os.environ.get(key)
    if not value:
        print(f"ERROR: missing required env var {key}", file=sys.stderr)
        sys.exit(1)
    return value


def head_size(s3, bucket: str, key: str) -> int | None:
    try:
        resp = s3.head_object(Bucket=bucket, Key=key)
        return int(resp["ContentLength"])
    except ClientError as e:
        code = e.response.get("Error", {}).get("Code")
        if code in ("404", "NoSuchKey", "NotFound"):
            return None
        raise


def copy_object(
    src_s3, src_bucket: str, src_key: str,
    dst_s3, dst_bucket: str, dst_key: str,
    *, dry_run: bool,
) -> tuple[str, int]:
    """Returns ('copied'|'skipped'|'missing-src', size_bytes)."""
    src_size = head_size(src_s3, src_bucket, src_key)
    if src_size is None:
        return ("missing-src", 0)

    dst_size = head_size(dst_s3, dst_bucket, dst_key)
    if dst_size is not None and dst_size == src_size:
        return ("skipped", src_size)

    if dry_run:
        return ("copied", src_size)

    with tempfile.NamedTemporaryFile(delete=True, prefix="kopah_migrate_") as tmp:
        src_s3.download_fileobj(src_bucket, src_key, tmp)
        tmp.flush()
        tmp.seek(0)
        dst_s3.upload_fileobj(tmp, dst_bucket, dst_key)

    return ("copied", src_size)


def fmt_bytes(n: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    f = float(n)
    for u in units:
        if f < 1024 or u == units[-1]:
            return f"{f:.1f} {u}"
        f /= 1024
    return f"{n} B"


def migrate(
    rows: Iterable[dict], *,
    src_s3, src_bucket: str,
    dst_s3, dst_bucket: str,
    include_transcripts: bool, dry_run: bool,
) -> None:
    totals = {"copied": 0, "skipped": 0, "missing-src": 0}
    total_bytes = 0
    t0 = time.time()

    for row in rows:
        slide_id = row["slide_id"]
        flat_prefix = row["flat_files_prefix"]
        print(f"\n=== {slide_id} (batch={row['export_batch']}) ===")

        suffixes = list(ESSENTIAL_FLAT_SUFFIXES)
        if include_transcripts:
            suffixes.append(TRANSCRIPTS_SUFFIX)

        for suf in suffixes:
            src_key = f"{flat_prefix}{slide_id}{suf}"
            dst_key = src_key
            status, size = copy_object(
                src_s3, src_bucket, src_key,
                dst_s3, dst_bucket, dst_key,
                dry_run=dry_run,
            )
            totals[status] += 1
            if status == "copied":
                total_bytes += size
            print(f"  [{status:>11}] {fmt_bytes(size):>10}  {src_key}")

    elapsed = time.time() - t0
    print(
        f"\nDone. {totals['copied']} copied, {totals['skipped']} skipped, "
        f"{totals['missing-src']} missing at source. "
        f"Transferred {fmt_bytes(total_bytes)} in {elapsed:.0f}s."
    )


def main() -> None:
    if ENV_PATH.exists():
        load_dotenv(ENV_PATH)

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--manifest", type=Path,
        default=Path(__file__).resolve().parent.parent / "manifest.csv",
        help="Path to manifest.csv produced by build_manifest.py",
    )
    p.add_argument("--include-transcripts", action="store_true",
                   help="Also copy the heavy _tx_file.csv.gz per slide")
    p.add_argument("--only-slides", nargs="+",
                   help="Migrate only these slide_ids (default: all in manifest)")
    p.add_argument("--dry-run", action="store_true",
                   help="Report what would be copied without transferring bytes")
    args = p.parse_args()

    if not args.manifest.exists():
        print(f"ERROR: manifest not found at {args.manifest}", file=sys.stderr)
        sys.exit(1)

    src_bucket = _env("SOURCE_S3_BUCKET")
    dst_bucket = _env("KOPAH_BUCKET")

    src_s3 = make_source_client()
    dst_s3 = make_kopah_client(
        endpoint_url=_env("KOPAH_ENDPOINT_URL"),
        access_key=_env("KOPAH_ACCESS_KEY_ID"),
        secret_key=_env("KOPAH_SECRET_ACCESS_KEY"),
    )

    with args.manifest.open() as f:
        rows = list(csv.DictReader(f))
    if args.only_slides:
        keep = set(args.only_slides)
        rows = [r for r in rows if r["slide_id"] in keep]
        if not rows:
            print(f"ERROR: no manifest rows matched --only-slides", file=sys.stderr)
            sys.exit(1)

    print(
        f"Migrating flat files for {len(rows)} slide(s) "
        f"from s3://{src_bucket}/ to s3://{dst_bucket}/"
    )
    if args.dry_run:
        print("(DRY RUN — no objects will be transferred)")

    migrate(
        rows,
        src_s3=src_s3, src_bucket=src_bucket,
        dst_s3=dst_s3, dst_bucket=dst_bucket,
        include_transcripts=args.include_transcripts,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
