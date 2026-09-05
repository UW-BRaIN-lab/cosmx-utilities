#!/usr/bin/env python3
"""Recover each slide's AtoMx segmentation version and parameters from S3.

The flat files do not record which AtoMx build segmented a slide -- `version` is
a single value across the whole cohort -- so a version-driven batch effect cannot
be tested from them. The SegmentationManifest_Parameters JSON in each slide's
CellStatsDir does record it, two ways:

  * any creation/date field it carries, and
  * the PRESENCE of `Run3DSegmentation` / `assignRNAtoNearest`, which AtoMx gained
    between 2026-04-16 and 2026-05-07. Those two keys appear in profiles from the
    later build and not the earlier one, so presence alone fingerprints the build
    even when no date is stored.

Emits a slide_groups.csv to feed back into segmentation_homogeneity.py:

    uv run python pipeline/python/segmentation_versions.py \
        --manifest pipeline/manifest.csv --output slide_groups.csv
    uv run python pipeline/python/segmentation_homogeneity.py \
        --flatfiles-dir staged --slide-groups slide_groups.csv \
        --group-by atomx_build --output by_slide.csv

Reads from the AWS SOURCE bucket (CellStatsDir was never migrated to Kopah).
Small enough to run on a login node; no Slurm needed.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _clients import make_source_client  # noqa: E402

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"

# Keys AtoMx gained in the newer build. Presence, not value, is the fingerprint.
FINGERPRINT_KEYS = ("Run3DSegmentation", "assignRNAtoNearest")
BUILD_WITH = "post-2026-05"
BUILD_WITHOUT = "pre-2026-05"

# Parameters worth carrying alongside, to prove configs really were identical.
PARAMETER_KEYS = (
    "NuclearDiameterUm", "CellDiameterUm", "CellDilationUm", "MinCellSizeUm",
    "MaxCellSizeUm", "NucleiModel", "CytoplasmModel", "ForegroundThreshold",
)
DATE_KEY_PATTERN = re.compile(r"(date|created|timestamp)", re.I)
PARAMETERS_FILE_PATTERN = re.compile(r"SegmentationManifest_Parameters_.*\.json$")
SEGMENTATION_DIR_PATTERN = re.compile(r"(Segmentation_[^/]+)/")


class MissingParametersError(Exception):
    """No SegmentationManifest_Parameters JSON under a slide's CellStatsDir."""


def walk_scalars(payload, found: dict) -> dict:
    """First scalar value seen for each key, at any nesting depth."""
    if isinstance(payload, dict):
        for key, value in payload.items():
            if not isinstance(value, (dict, list)) and key not in found:
                found[key] = value
            walk_scalars(value, found)
    elif isinstance(payload, list):
        for item in payload:
            walk_scalars(item, found)
    return found


def extract_facts(payload) -> dict:
    """Version fingerprint plus key parameters from one parameters payload."""
    scalars = walk_scalars(payload, {})
    lowered = {k.lower(): k for k in scalars}

    present = [k for k in FINGERPRINT_KEYS if k.lower() in lowered]
    facts: dict = {
        "atomx_build": BUILD_WITH if present else BUILD_WITHOUT,
        "fingerprint_keys": "|".join(present),
    }
    for key in PARAMETER_KEYS:
        actual = lowered.get(key.lower())
        facts[key] = scalars.get(actual, "") if actual else ""
    for key, original in lowered.items():
        if DATE_KEY_PATTERN.search(key):
            facts.setdefault("date_field", f"{original}={scalars[original]}")
    facts.setdefault("date_field", "")
    return facts


def cellstats_prefix(decoded_prefix: str) -> str:
    """CellStatsDir sits INSIDE the slide/scan directory decoded_prefix names."""
    return f"{decoded_prefix.rstrip('/')}/CellStatsDir"


def find_parameters_keys(client, bucket: str, prefix: str) -> list[str]:
    """Every SegmentationManifest_Parameters JSON under a CellStatsDir."""
    keys, token = [], None
    while True:
        kwargs = {"Bucket": bucket, "Prefix": f"{prefix}/"}
        if token:
            kwargs["ContinuationToken"] = token
        page = client.list_objects_v2(**kwargs)
        for item in page.get("Contents", []):
            if PARAMETERS_FILE_PATTERN.search(item["Key"]):
                keys.append(item["Key"])
        if not page.get("IsTruncated"):
            return keys
        token = page["NextContinuationToken"]


def read_manifest(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str)
    # csv.DictWriter writes CRLF, so the LAST column carries a trailing \r.
    for column in frame.columns:
        frame[column] = frame[column].str.strip().str.strip("\r")
    return frame


def main() -> None:
    if ENV_PATH.exists():
        load_dotenv(ENV_PATH)

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--source-bucket", default=None,
                        help="Defaults to SOURCE_S3_BUCKET from pipeline/.env")
    parser.add_argument("--seg-uuid-column", default="cellSegmentationSetId",
                        help="Ignored unless --flatfiles-dir is given")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import os
    bucket = args.source_bucket or os.environ.get("SOURCE_S3_BUCKET")
    if not bucket:
        print("ERROR: set --source-bucket or SOURCE_S3_BUCKET", file=sys.stderr)
        sys.exit(1)

    manifest = read_manifest(args.manifest)
    if "decoded_prefix" not in manifest.columns:
        print("ERROR: manifest has no decoded_prefix column", file=sys.stderr)
        sys.exit(1)

    client = make_source_client()
    rows = []
    for _, entry in manifest.iterrows():
        slide_id = entry["slide_id"]
        prefix = cellstats_prefix(entry["decoded_prefix"])
        print(f"{slide_id} ...", flush=True)
        try:
            keys = find_parameters_keys(client, bucket, prefix)
        except Exception as exc:                      # noqa: BLE001 - report and continue
            print(f"  WARN: listing failed: {exc}", file=sys.stderr)
            continue
        if not keys:
            print(f"  WARN: no parameters JSON under {prefix}", file=sys.stderr)
            continue

        # A slide can carry several segmentations; record each, newest suffix last.
        for key in sorted(keys):
            body = client.get_object(Bucket=bucket, Key=key)["Body"].read()
            try:
                payload = json.loads(body)
            except ValueError as exc:
                print(f"  WARN: {key} is not JSON: {exc}", file=sys.stderr)
                continue
            match = SEGMENTATION_DIR_PATTERN.search(key)
            row = {
                "slide_id": slide_id,
                "segmentation_dir": match.group(1) if match else "",
            }
            row.update(extract_facts(payload))
            rows.append(row)
            print(f"  {row['segmentation_dir']}: {row['atomx_build']}"
                  f"{' ' + row['date_field'] if row['date_field'] else ''}")

    if not rows:
        print("ERROR: no segmentation parameters recovered", file=sys.stderr)
        sys.exit(1)

    frame = pd.DataFrame(rows)
    frame.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}  ({len(frame)} segmentation(s), "
          f"{frame['slide_id'].nunique()} slides)")

    print("\nBuilds across the cohort:")
    print(frame.groupby("atomx_build")["slide_id"].nunique()
          .rename("slides").to_string())
    for key in PARAMETER_KEYS:
        if key in frame and frame[key].astype(str).nunique() > 1:
            print(f"  {key} is NOT uniform: "
                  f"{sorted(frame[key].astype(str).unique())[:6]}")


if __name__ == "__main__":
    main()
