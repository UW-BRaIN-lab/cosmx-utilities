#!/usr/bin/env python3
"""Build a per-slide manifest from the CosMx S3 directory structure.

Crawls SOURCE_S3_BUCKET/SOURCE_S3_PREFIX, extracts slide_id, run_uuid, run_date,
instrument_id, slot, and export_batch for each slide, and writes a CSV manifest.

The manifest is the source of truth for downstream pipeline stages: which slides
exist, when/where they were run, and what batch/instrument context to join into
AnnData obs (so we can revisit cohort granularity later if needed).

Instrument and run metadata come from RunSummary filenames, which look like:
    Run_<uuid>_<YYYYMMDD>_<HHMMSS>_S<slot>_<instrument>_<suffix>
    Run_430bb75a-4832-4787-ac6a-06eddf7c82e3_20250320_234427_S2_2309H0195_ExptConfig.txt

Usage:
    uv run python pipeline/python/build_manifest.py
    uv run python pipeline/python/build_manifest.py --output pipeline/manifest.csv
    uv run python pipeline/python/build_manifest.py \\
        --source-bucket my-bucket --source-prefix CosMx-GBM/Wenyu-cohort
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Iterator

from dotenv import load_dotenv

from _clients import make_source_client

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"

RUNSUMMARY_FILE_RE = re.compile(
    r"^Run_(?P<run_uuid>[0-9a-fA-F-]+)"
    r"_(?P<run_date>\d{8})_(?P<run_time>\d{6})"
    r"_S(?P<slot>\d+)"
    r"_(?P<instrument_id>[A-Za-z0-9]+)"
    r"_"
)


@dataclass
class SlideManifestRow:
    export_batch: str
    slide_id: str
    run_uuid: str
    run_date: str
    run_time: str
    instrument_id: str
    slot: int
    flat_files_prefix: str
    decoded_prefix: str


def _list_subprefixes(s3, bucket: str, prefix: str) -> list[str]:
    if not prefix.endswith("/"):
        prefix += "/"
    paginator = s3.get_paginator("list_objects_v2")
    out: list[str] = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
        for cp in page.get("CommonPrefixes", []):
            out.append(cp["Prefix"])
    return out


def _first_runsummary_key(s3, bucket: str, decoded_scan_prefix: str) -> str | None:
    resp = s3.list_objects_v2(
        Bucket=bucket,
        Prefix=decoded_scan_prefix + "RunSummary/Run_",
        MaxKeys=1,
    )
    contents = resp.get("Contents") or []
    return contents[0]["Key"] if contents else None


def discover_slides(s3, bucket: str, prefix: str) -> Iterator[SlideManifestRow]:
    """Yield one SlideManifestRow per slide found under prefix.

    Expected layout:
        <prefix>/<export_batch>/DecodedFiles/<slide_id>/<YYYYMMDD_HHMMSS_S{slot}>/RunSummary/Run_...
        <prefix>/<export_batch>/flatFiles/<slide_id>/
    """
    for batch_prefix in _list_subprefixes(s3, bucket, prefix):
        export_batch = batch_prefix.rstrip("/").rsplit("/", 1)[-1]
        decoded_root = batch_prefix + "DecodedFiles/"
        flat_root = batch_prefix + "flatFiles/"

        for slide_prefix in _list_subprefixes(s3, bucket, decoded_root):
            slide_id = slide_prefix.rstrip("/").rsplit("/", 1)[-1]

            scan_dirs = [
                d for d in _list_subprefixes(s3, bucket, slide_prefix)
                if not d.rstrip("/").endswith("/Logs")
            ]
            if not scan_dirs:
                print(f"WARN: no scan dir for slide {slide_id}", file=sys.stderr)
                continue
            if len(scan_dirs) > 1:
                print(
                    f"WARN: slide {slide_id} has {len(scan_dirs)} scan dirs; "
                    f"using {scan_dirs[0]}",
                    file=sys.stderr,
                )
            scan_prefix = scan_dirs[0]

            run_key = _first_runsummary_key(s3, bucket, scan_prefix)
            if not run_key:
                print(
                    f"WARN: no RunSummary file for slide {slide_id} at {scan_prefix}",
                    file=sys.stderr,
                )
                continue
            m = RUNSUMMARY_FILE_RE.match(run_key.rsplit("/", 1)[-1])
            if not m:
                print(
                    f"WARN: unparseable RunSummary filename: {run_key}",
                    file=sys.stderr,
                )
                continue

            yield SlideManifestRow(
                export_batch=export_batch,
                slide_id=slide_id,
                run_uuid=m.group("run_uuid"),
                run_date=m.group("run_date"),
                run_time=m.group("run_time"),
                instrument_id=m.group("instrument_id"),
                slot=int(m.group("slot")),
                flat_files_prefix=f"{flat_root}{slide_id}/",
                decoded_prefix=scan_prefix,
            )


def main() -> None:
    if ENV_PATH.exists():
        load_dotenv(ENV_PATH)

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--source-bucket",
        default=os.environ.get("SOURCE_S3_BUCKET"),
        help="S3 bucket containing the CosMx data (env: SOURCE_S3_BUCKET)",
    )
    p.add_argument(
        "--source-prefix",
        default=os.environ.get("SOURCE_S3_PREFIX"),
        help="Top-level prefix within the bucket (env: SOURCE_S3_PREFIX)",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "manifest.csv",
        help="Path to write the manifest CSV",
    )
    args = p.parse_args()

    if not args.source_bucket or not args.source_prefix:
        print(
            "ERROR: --source-bucket and --source-prefix are required "
            "(set SOURCE_S3_BUCKET / SOURCE_S3_PREFIX in .env, or pass on CLI).",
            file=sys.stderr,
        )
        sys.exit(1)

    s3 = make_source_client()
    rows = list(discover_slides(s3, args.source_bucket, args.source_prefix))
    if not rows:
        print("ERROR: no slides found.", file=sys.stderr)
        sys.exit(1)

    field_names = [f.name for f in fields(SlideManifestRow)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    by_instrument: dict[str, int] = {}
    for r in rows:
        by_instrument[r.instrument_id] = by_instrument.get(r.instrument_id, 0) + 1
    print(f"Wrote {len(rows)} slides to {args.output}")
    print("Slides per instrument:")
    for k, v in sorted(by_instrument.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
