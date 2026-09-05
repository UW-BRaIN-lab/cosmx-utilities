#!/usr/bin/env python3
"""Recover each slide's AtoMx segmentation version and parameters from S3.

The flat files do not record which AtoMx build segmented a slide -- `version` is
a single value across the whole cohort -- so a version-driven batch effect cannot
be tested from them. The SegmentationManifest_Parameters JSON in each slide's
CellStatsDir does record it, two ways:

  * a `Datecreated` field, in UTC, which is what this uses.

A key-presence fingerprint was tried first -- `Run3DSegmentation` and
`assignRNAtoNearest` appear in the AtoMx UI for later-build profiles and not
earlier ones -- and it does NOT work: the JSON stores neither key, so all 57 slides
came back as the earlier build including two known to be later. The UI renders
those fields from its own current defaults rather than from the saved record. The
stored parameters are absent too; the date is the only usable signal.

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

import numpy as np
import pandas as pd
from dotenv import load_dotenv

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _clients import make_source_client  # noqa: E402

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"

# Kept only as a diagnostic: these are absent from every JSON we have seen, so the
# build is derived from the stored date instead. See the module docstring.
FINGERPRINT_KEYS = ("Run3DSegmentation", "assignRNAtoNearest")
# The compartment change lands between 2026-04-16 (cytoplasm present) and
# 2026-05-07 (absent), so any cutoff in that gap separates the two behaviours.
DEFAULT_BUILD_CUTOFF = "2026-05-01"
DATE_VALUE_PATTERN = re.compile(r"(\d{4}-\d{2}-\d{2})")

# Parameters worth carrying alongside, to prove configs really were identical.
PARAMETER_KEYS = (
    "NuclearDiameterUm", "CellDiameterUm", "CellDilationUm", "MinCellSizeUm",
    "MaxCellSizeUm", "NucleiModel", "CytoplasmModel", "ForegroundThreshold",
)
DATE_KEY_PATTERN = re.compile(r"(date|created|timestamp)", re.I)
PARAMETERS_FILE_PATTERN = re.compile(r"SegmentationManifest_Parameters_.*\.json$")
SEG_UUID_COLUMN = "cellSegmentationSetId"
# Enough of the gzipped metadata to reach its header and first data row.
METADATA_PEEK_BYTES = 131072
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

    facts: dict = {"fingerprint_keys": "|".join(
        k for k in FINGERPRINT_KEYS if k.lower() in lowered)}
    for key in PARAMETER_KEYS:
        actual = lowered.get(key.lower())
        facts[key] = scalars.get(actual, "") if actual else ""
    for key, original in lowered.items():
        if DATE_KEY_PATTERN.search(key):
            facts.setdefault("date_field", f"{original}={scalars[original]}")
    facts.setdefault("date_field", "")
    match = DATE_VALUE_PATTERN.search(facts["date_field"])
    facts["segmentation_date"] = match.group(1) if match else ""
    return facts


def active_segmentation_uuid(client, bucket: str, flat_prefix: str,
                             slide_id: str) -> str:
    """The cellSegmentationSetId the flat files were built from.

    A slide carries the ORIGINAL segmentation from acquisition alongside any later
    resegmentation, so dates alone mix 2024 originals with 2026 resegmentations.
    Only the segmentation named here produced the data everything downstream uses.

    Range-requests the head of the gzipped metadata rather than the whole file --
    the header and one data row are all that is needed.
    """
    import zlib

    key = f"{flat_prefix.rstrip('/')}/{slide_id}_metadata_file.csv.gz"
    try:
        raw = client.get_object(Bucket=bucket, Key=key,
                                Range=f"bytes=0-{METADATA_PEEK_BYTES - 1}")["Body"].read()
    except Exception as exc:                          # noqa: BLE001
        print(f"  WARN: could not read {key}: {exc}", file=sys.stderr)
        return ""
    # A truncated gzip stream raises at the end; whatever decompressed is enough.
    decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
    try:
        text = decompressor.decompress(raw).decode("utf-8", "ignore")
    except zlib.error as exc:
        print(f"  WARN: could not decompress {key}: {exc}", file=sys.stderr)
        return ""
    lines = text.splitlines()
    if len(lines) < 2:
        print(f"  WARN: {key} decompressed to {len(lines)} line(s)", file=sys.stderr)
        return ""
    header = [h.strip().strip('"') for h in lines[0].split(",")]
    if SEG_UUID_COLUMN not in header:
        print(f"  WARN: {key} has no {SEG_UUID_COLUMN} column "
              f"(saw {len(header)} columns, e.g. {header[:4]})", file=sys.stderr)
        return ""
    values = [v.strip().strip('"') for v in lines[1].split(",")]
    index = header.index(SEG_UUID_COLUMN)
    return values[index] if index < len(values) else ""


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
    parser.add_argument("--all-segmentations", action="store_true",
                        help="Keep every Segmentation_* dir. By default only the one "
                             "the flat files were built from is kept, so the output "
                             "joins one-to-one and 2024 originals do not contaminate "
                             "the 2026 resegmentation dates.")
    parser.add_argument("--build-cutoff", default=DEFAULT_BUILD_CUTOFF,
                        help="Segmentations on or after this date are labelled the "
                             f"later build (default {DEFAULT_BUILD_CUTOFF})")
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
        active_uuid = ""
        if not args.all_segmentations and "flat_files_prefix" in manifest.columns:
            active_uuid = active_segmentation_uuid(
                client, bucket, entry["flat_files_prefix"], slide_id)
        try:
            keys = find_parameters_keys(client, bucket, prefix)
        except Exception as exc:                      # noqa: BLE001 - report and continue
            print(f"  WARN: listing failed: {exc}", file=sys.stderr)
            continue
        if not keys:
            print(f"  WARN: no parameters JSON under {prefix}", file=sys.stderr)
            continue

        # A slide carries its acquisition-time segmentation plus any resegmentation.
        if active_uuid:
            matching = [k for k in keys if active_uuid in k]
            if matching:
                keys = matching
            else:
                print(f"  WARN: no segmentation matches {active_uuid}, keeping all",
                      file=sys.stderr)
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
                "is_active": bool(active_uuid and active_uuid in key),
            }
            row.update(extract_facts(payload))
            rows.append(row)
            print(f"  {row['segmentation_dir']}: "
                  f"{row['segmentation_date'] or 'no date'}")

    if not rows:
        print("ERROR: no segmentation parameters recovered", file=sys.stderr)
        sys.exit(1)

    frame = pd.DataFrame(rows)
    dated = frame["segmentation_date"].astype(str)
    frame["atomx_build"] = np.where(
        dated == "", "unknown",
        np.where(dated >= args.build_cutoff,
                 f"on/after {args.build_cutoff}", f"before {args.build_cutoff}"))
    frame.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}  ({len(frame)} segmentation(s), "
          f"{frame['slide_id'].nunique()} slides)")

    undated = int((frame["segmentation_date"] == "").sum())
    if undated:
        print(f"\nWARN: {undated} segmentation(s) carry no date", file=sys.stderr)
    print(f"\nSegmentation dates: {dated[dated != ''].min()} to "
          f"{dated[dated != ''].max()}")
    resolved = int(frame["is_active"].sum()) if "is_active" in frame else 0
    slides = frame["slide_id"].nunique()
    print(f"\nActive segmentation resolved for {resolved} of {slides} slides")
    if resolved < slides and not args.all_segmentations:
        print(f"  {slides - resolved} slide(s) fell back to keeping every "
              f"segmentation -- see the WARN lines above for why. Redirect stderr "
              f"to a file (2> versions.err) to read them.", file=sys.stderr)

    counts = frame["slide_id"].value_counts()
    extra = counts[counts > 1]
    if len(extra):
        print(f"\nWARN: {len(extra)} slide(s) contribute more than one row "
              f"({int(counts.sum() - slides)} extra rows); a join on slide_id will "
              f"duplicate their geometry.", file=sys.stderr)
        print(f"  most segmentations on one slide: {int(counts.max())} "
              f"({counts.idxmax()})", file=sys.stderr)

    print("\nBuilds across the cohort (segmentations, not slides):")
    print(frame["atomx_build"].value_counts().rename("n").to_string())
    for key in PARAMETER_KEYS:
        if key in frame and frame[key].astype(str).nunique() > 1:
            print(f"  {key} is NOT uniform: "
                  f"{sorted(frame[key].astype(str).unique())[:6]}")


if __name__ == "__main__":
    main()
