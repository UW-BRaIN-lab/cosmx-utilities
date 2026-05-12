#!/usr/bin/env python3
"""Process a single CosMx slide: detect segmentation, download, stitch, and upload.

Replaces process-slide-wrapper.sh with a Python implementation that uses DuckDB
to query segmentation manifests directly in S3 (zero pre-downloads) and does
targeted S3 syncs to download only the files stitch-images actually needs.

Usage:
    uv run python scripts/process-slide.py my-bucket CosMx-GBM/.../SlideName/ScanId
    uv run python scripts/process-slide.py --whatif my-bucket CosMx-GBM/.../SlideName/ScanId
    uv run python scripts/process-slide.py --segmentation-version Segmentation_uuid_003 my-bucket CosMx-GBM/.../SlideName/ScanId
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import boto3
from botocore.exceptions import ClientError


# ── Path parsing ────────────────────────────────────────────────────────────


@dataclass
class SlideContext:
    """All path components derived from the S3 bucket + slide base path."""

    bucket: str
    slide_base_path: str  # Study/[Experiment/]AtoMxRun/DecodedFiles/SlideName/ScanId
    work_dir: Path = Path("/tmp/slide")

    def __post_init__(self) -> None:
        parts = self.slide_base_path.rstrip("/").split("/")
        # Two supported layouts:
        #   5-level: Study/AtoMxRun/DecodedFiles/SlideName/ScanId            (typical)
        #   6-level: Study/Experiment/AtoMxRun/DecodedFiles/SlideName/ScanId (when
        #            grouping multiple resegmentations of the same scans under a
        #            shared experiment wrapper for comparison)
        self.slide_name = parts[-2]
        self.scan_id = parts[-1]
        self.atomx_run = parts[-4]
        self.study = parts[0]
        if len(parts) >= 6:
            self.experiment = parts[-5]
            self.experiment_prefix = f"{self.study}/{self.experiment}"
            self.output_prefix = (
                f"napari-stitched/{self.study}/{self.experiment}"
                f"/{self.atomx_run}/{self.slide_name}"
            )
        else:
            self.experiment = ""
            self.experiment_prefix = self.study
            self.output_prefix = (
                f"napari-stitched/{self.study}"
                f"/{self.atomx_run}/{self.slide_name}"
            )

    def s3(self, suffix: str = "") -> str:
        return f"s3://{self.bucket}/{self.slide_base_path}/{suffix}"


# ── Benchmarking ────────────────────────────────────────────────────────────


@dataclass
class Benchmark:
    """Tracks step timings and writes a benchmark CSV on success or failure."""

    ctx: SlideContext
    whatif: bool
    start_time: str = ""
    status: str = "failed"
    failed_step: str = ""
    _uploaded: bool = False
    _timings: dict[str, tuple[float | None, float | None]] = field(
        default_factory=lambda: {
            "download": (None, None),
            "stitch": (None, None),
            "read_targets": (None, None),
            "metadata": (None, None),
            "upload": (None, None),
        }
    )

    def start(self, step: str) -> float:
        self.failed_step = step
        t = time.time()
        self._timings[step] = (t, None)
        return t

    def end(self, step: str) -> float:
        t = time.time()
        start, _ = self._timings[step]
        self._timings[step] = (start, t)
        return t

    def _duration(self, step: str) -> str:
        start, end = self._timings.get(step, (None, None))
        if start is not None and end is not None:
            return str(int(end - start))
        return ""

    def _duration_seconds(self, step: str) -> int | None:
        start, end = self._timings.get(step, (None, None))
        if start is not None and end is not None:
            return int(end - start)
        return None

    def _instance_metadata(self) -> tuple[str, str]:
        vcpus = str(multiprocessing.cpu_count())
        mem_mb = ""
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        mem_mb = str(int(line.split()[1]) // 1024)
                        break
        except FileNotFoundError:
            pass
        return vcpus, mem_mb

    def _fov_stats(self) -> tuple[str, str]:
        cell_stats = self.ctx.work_dir / "CellStatsDir"
        if not cell_stats.is_dir():
            return "", ""
        fov_count = sum(1 for d in cell_stats.iterdir() if d.is_dir() and d.name.startswith("FOV"))
        total_bytes = sum(f.stat().st_size for f in cell_stats.rglob("*") if f.is_file())
        return str(fov_count), str(total_bytes)

    def write_and_upload(self) -> None:
        if self._uploaded:
            return
        self._uploaded = True

        csv_path = self.ctx.work_dir / "benchmark.csv"
        vcpus, mem_mb = self._instance_metadata()
        fov_count, fov_bytes = self._fov_stats()

        download_start, _ = self._timings.get("download", (None, None))
        _, upload_end = self._timings.get("upload", (None, None))
        total = ""
        if download_start is not None and upload_end is not None:
            total = str(int(upload_end - download_start))

        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "slide_name", "study", "experiment", "atomx_run",
                "vcpus", "memory_mb", "fov_count", "fov_bytes", "start_time",
                "duration_download_s", "duration_stitch_s", "duration_read_targets_s",
                "duration_metadata_s", "duration_upload_s", "duration_total_s",
                "status", "failed_step",
            ])
            writer.writerow([
                self.ctx.slide_name, self.ctx.study, self.ctx.experiment, self.ctx.atomx_run,
                vcpus, mem_mb, fov_count, fov_bytes, self.start_time,
                self._duration("download"), self._duration("stitch"),
                self._duration("read_targets"), self._duration("metadata"),
                self._duration("upload"), total,
                self.status, self.failed_step,
            ])

        bench_name = f"benchmark_{vcpus}cpu_{mem_mb}mb.csv"
        if self.whatif:
            log(f"  Benchmark (whatif, local):")
            log(csv_path.read_text())
        else:
            _get_s3().upload_file(
                str(csv_path), self.ctx.bucket,
                f"{self.ctx.output_prefix}/{bench_name}",
            )
            log(f"  Benchmark: s3://{self.ctx.bucket}/{self.ctx.output_prefix}/{bench_name}")


# ── Helpers ─────────────────────────────────────────────────────────────────


def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(msg: str) -> None:
    print(msg, flush=True)


_s3 = None


def _get_s3():
    global _s3
    if _s3 is None:
        _s3 = boto3.client("s3")
    return _s3


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, **kwargs)


def s3_sync(
    src: str,
    dst: str,
    *,
    exclude: str | None = None,
    includes: list[str] | None = None,
    dryrun: bool = False,
) -> None:
    cmd = ["aws", "s3", "sync", src, dst]
    if exclude:
        cmd += ["--exclude", exclude]
    for inc in (includes or []):
        cmd += ["--include", inc]
    if dryrun:
        cmd.append("--dryrun")
    run(cmd)


def s3_ls_prefixes(bucket: str, prefix: str) -> list[str]:
    """List immediate subdirectory prefixes under an S3 path."""
    response = _get_s3().list_objects_v2(
        Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/",
    )
    return [
        p["Prefix"].rstrip("/").rsplit("/", 1)[-1]
        for p in response.get("CommonPrefixes", [])
    ]


# ── DuckDB S3 credentials ───────────────────────────────────────────────────


def _duckdb_create_s3_secret(con: "duckdb.DuckDBPyConnection") -> None:
    """Configure DuckDB S3 access using boto3-resolved credentials.

    Resolves credentials through boto3's credential chain (env vars, SSO,
    profiles, IAM role). This ensures DuckDB can access S3 even when using
    SSO sessions that DuckDB's credential_chain provider can't auto-detect.
    """
    session = boto3.Session()
    credentials = session.get_credentials()
    if credentials is None:
        con.sql("INSTALL httpfs; LOAD httpfs;")
        return

    creds = credentials.get_frozen_credentials()
    region = (
        session.region_name
        or os.environ.get("AWS_DEFAULT_REGION")
        or os.environ.get("AWS_REGION", "us-west-2")
    )

    secret_sql = f"""
        CREATE OR REPLACE SECRET secret (
            TYPE s3,
            PROVIDER config,
            KEY_ID '{creds.access_key}',
            SECRET '{creds.secret_key}',
            REGION '{region}'
    """
    if creds.token:
        secret_sql += f",\n        SESSION_TOKEN '{creds.token}'"
    secret_sql += "\n    );"

    con.sql("INSTALL httpfs; LOAD httpfs;")
    con.sql(secret_sql)


# ── Step 0: Detect segmentation version ────────────────────────────────────


def _read_seg_id_from_flatfiles(ctx: SlideContext) -> str:
    """Read the cellSegmentationSetId from this AtoMx run's flatFiles metadata.

    Downloads the first row of the run's metadata CSV to extract the
    segmentation UUID that was used for this specific AtoMx run's analysis.
    Returns the UUID string, or "" if not found.
    """
    import csv
    import gzip
    import tempfile

    s3_key = (
        f"{ctx.experiment_prefix}/{ctx.atomx_run}"
        f"/flatFiles/{ctx.slide_name}/{ctx.slide_name}_metadata_file.csv.gz"
    )

    with tempfile.NamedTemporaryFile(suffix=".csv.gz", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        try:
            _get_s3().download_file(ctx.bucket, s3_key, tmp_path)
        except ClientError:
            log(f"  WARNING: Could not download s3://{ctx.bucket}/{s3_key}")
            return ""

        with gzip.open(tmp_path, "rt") as f:
            reader = csv.DictReader(f)
            row = next(reader, None)
            if row:
                return row.get("cellSegmentationSetId", "").strip()
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return ""


def _find_seg_subdir_by_uuid(ctx: SlideContext, seg_id: str) -> str:
    """Find the Segmentation_* subdirectory in CellStatsDir matching a UUID.

    Queries S3 for Segmentation_* prefixes and returns the one whose name
    contains the given UUID. Returns "" if no match.
    """
    seg_dirs = s3_ls_prefixes(
        ctx.bucket,
        f"{ctx.slide_base_path}/CellStatsDir",
    )
    for dirname in seg_dirs:
        if seg_id in dirname:
            return dirname
    return ""


def _detect_highest_version(ctx: SlideContext) -> tuple[str, str]:
    """Fallback: query S3 manifest JSONs via DuckDB to find the highest version.

    Returns (celllabels_subdir, seg_uuid) or ("", "") if none found.
    """
    import duckdb

    con = duckdb.connect()
    _duckdb_create_s3_secret(con)

    s3_glob = (
        f"s3://{ctx.bucket}/{ctx.slide_base_path}"
        f"/CellStatsDir/Segmentation_*"
        f"/SegmentationManifest_Parameters_*.json"
    )
    try:
        result = con.sql(f"""
            SELECT
                regexp_extract(filename, '(Segmentation_[^/]+)', 1) AS seg_dir
            FROM read_json_auto(
                '{s3_glob}',
                filename=true,
                union_by_name=true
            )
            ORDER BY CAST(Version AS INTEGER) DESC
            LIMIT 1
        """).fetchone()
        if result and result[0]:
            subdir = result[0]
            uuid = re.sub(r"^Segmentation_", "", subdir)
            uuid = re.sub(r"_\d+$", "", uuid)
            return subdir, uuid
    except Exception as e:
        print(f"WARNING: DuckDB manifest query failed: {e}", file=sys.stderr)

    return "", ""


def detect_segmentation(ctx: SlideContext) -> tuple[str, str]:
    """Detect the correct segmentation version for this AtoMx run.

    Strategy: read cellSegmentationSetId from the AtoMx run's own flatFiles
    metadata, then find the matching Segmentation_* subdirectory in
    CellStatsDir. This ensures each run uses its own segmentation version
    when multiple resegmentations exist.

    Falls back to picking the highest version from manifest JSONs if the
    flatFiles lookup fails.

    Returns (celllabels_subdir, seg_uuid) or ("", "") if none found.
    """
    # Primary: match via flatFiles cellSegmentationSetId
    seg_id = _read_seg_id_from_flatfiles(ctx)
    if seg_id:
        log(f"  flatFiles cellSegmentationSetId: {seg_id}")
        subdir = _find_seg_subdir_by_uuid(ctx, seg_id)
        if subdir:
            return subdir, seg_id
        log(f"  WARNING: No CellStatsDir subdir matches UUID {seg_id}, "
            "falling back to highest version")

    # Fallback: highest version from manifest JSONs
    return _detect_highest_version(ctx)


def _seg_version(dirname: str) -> int:
    """Extract the numeric version suffix from a Segmentation_* directory name."""
    m = re.search(r"_(\d+)$", dirname)
    return int(m.group(1)) if m else 0


def _seg_uuid(dirname: str) -> str:
    """Extract the UUID from a Segmentation_* directory name."""
    uuid = re.sub(r"^Segmentation_", "", dirname)
    return re.sub(r"_\d+$", "", uuid)


def _list_seg_fovs(ctx: SlideContext, seg_dir: str) -> set[int]:
    """List which FOV numbers a Segmentation_* subdir contains (via S3)."""
    prefix = f"{ctx.slide_base_path}/CellStatsDir/{seg_dir}/"
    fov_dirs = s3_ls_prefixes(ctx.bucket, prefix)
    fovs = set()
    for d in fov_dirs:
        m = re.match(r"FOV(\d+)$", d)
        if m:
            fovs.add(int(m.group(1)))
    return fovs


def detect_all_segmentations(ctx: SlideContext) -> list[tuple[str, str]]:
    """Detect the segmentation versions needed to cover all FOVs for this run.

    Strategy:
    1. Identify this run's primary segmentation via flatFiles metadata
       (cellSegmentationSetId). This is the segmentation the run's analysis
       was performed against.
    2. Check if it covers all FOVs (using the original _001 as reference).
    3. If not, fill gaps by working backward from the primary version through
       older versions — never newer, since newer versions didn't exist when
       this run's analysis was done.

    Returns list of (celllabels_subdir, seg_uuid) tuples, ordered by
    version descending (primary first). Returns [] if no Segmentation_*
    subdirs exist.
    """
    cellstats_prefix = f"{ctx.slide_base_path}/CellStatsDir"
    seg_dirs = s3_ls_prefixes(ctx.bucket, cellstats_prefix)
    seg_dirs = [d for d in seg_dirs if d.startswith("Segmentation_")]
    if not seg_dirs:
        return []

    # Find this run's primary segmentation
    primary_seg_id = _read_seg_id_from_flatfiles(ctx)
    if not primary_seg_id:
        log("  WARNING: Could not read cellSegmentationSetId from flatFiles")
        return _detect_highest_version_as_list(ctx, seg_dirs)

    primary_dir = _find_seg_subdir_by_uuid(ctx, primary_seg_id)
    if not primary_dir:
        log(f"  WARNING: No CellStatsDir subdir matches UUID {primary_seg_id}")
        return _detect_highest_version_as_list(ctx, seg_dirs)

    primary_version = _seg_version(primary_dir)
    log(f"  Primary segmentation: {primary_dir} (UUID: {primary_seg_id})")

    # Reference FOVs from the original segmentation (lowest version = all FOVs)
    seg_dirs_sorted = sorted(seg_dirs, key=_seg_version)
    original_dir = seg_dirs_sorted[0]
    all_fovs = _list_seg_fovs(ctx, original_dir)
    log(f"  Total FOVs in slide: {len(all_fovs)} (from {original_dir})")

    # Check primary coverage
    primary_fovs = _list_seg_fovs(ctx, primary_dir)
    if primary_fovs >= all_fovs:
        log(f"  {primary_dir}: covers all {len(all_fovs)} FOVs")
        return [(primary_dir, primary_seg_id)]

    # Primary doesn't cover all FOVs — fill gaps from older versions
    selected = [(primary_dir, primary_seg_id)]
    covered_fovs = set(primary_fovs)
    log(f"  {primary_dir}: {len(primary_fovs)} FOVs "
        f"({len(covered_fovs)}/{len(all_fovs)} covered)")

    # Only consider versions older than the primary, newest first
    older_dirs = [d for d in seg_dirs if _seg_version(d) < primary_version]
    older_dirs.sort(key=_seg_version, reverse=True)

    for dirname in older_dirs:
        if covered_fovs >= all_fovs:
            break
        fovs = _list_seg_fovs(ctx, dirname)
        new_fovs = fovs - covered_fovs
        if new_fovs:
            selected.append((dirname, _seg_uuid(dirname)))
            covered_fovs.update(new_fovs)
            log(f"  {dirname}: {len(new_fovs)} new FOVs "
                f"({len(covered_fovs)}/{len(all_fovs)} covered)")
        else:
            log(f"  {dirname}: skipped (all {len(fovs)} FOVs already covered)")

    if covered_fovs < all_fovs:
        missing = all_fovs - covered_fovs
        log(f"  WARNING: {len(missing)} FOVs still uncovered: {sorted(missing)[:20]}...")

    return selected


def _detect_highest_version_as_list(
    ctx: SlideContext, seg_dirs: list[str],
) -> list[tuple[str, str]]:
    """Fallback: return the highest-version segmentation as a single-element list."""
    result = _detect_highest_version(ctx)
    if result[0]:
        return [result]
    return []


# ── Step 1: Download slide data ────────────────────────────────────────────


def download_slide(ctx: SlideContext, celllabels_subdirs: str, dryrun: bool = False) -> str:
    """Download only the files needed for processing. Returns the AnalysisResults subdir name."""
    # CellStatsDir: CellLabels TIFs + Morphology2D images only
    s3_sync(
        ctx.s3("CellStatsDir/"), str(ctx.work_dir / "CellStatsDir/"),
        exclude="*",
        includes=[
            "FOV*/CellLabels_F*.tif", "FOV*/CellLabels_F*.TIF",
            "Morphology2D/*.TIF", "Morphology2D/*.tif",
        ],
        dryrun=dryrun,
    )

    # Selected segmentation version(s)' CellLabels
    if celllabels_subdirs:
        includes = []
        for subdir in celllabels_subdirs.split(","):
            subdir = subdir.strip()
            includes += [
                f"{subdir}/FOV*/CellLabels_F*.tif",
                f"{subdir}/FOV*/CellLabels_F*.TIF",
            ]
        s3_sync(
            ctx.s3("CellStatsDir/"), str(ctx.work_dir / "CellStatsDir/"),
            exclude="*",
            includes=includes,
            dryrun=dryrun,
        )

    # RunSummary: only FOV location CSVs
    s3_sync(
        ctx.s3("RunSummary/"), str(ctx.work_dir / "RunSummary/"),
        exclude="*",
        includes=["*.fovs.csv", "*FOV_Locations*"],
        dryrun=dryrun,
    )

    # AnalysisResults: full sync (needed for read-targets)
    subdirs = s3_ls_prefixes(ctx.bucket, f"{ctx.slide_base_path}/AnalysisResults")
    if len(subdirs) != 1:
        raise RuntimeError(
            f"Expected exactly 1 AnalysisResults subdirectory, found {len(subdirs)}: {subdirs}"
        )
    analysis_subdir = subdirs[0]
    s3_sync(
        ctx.s3(f"AnalysisResults/{analysis_subdir}/"),
        str(ctx.work_dir / "AnalysisResults" / analysis_subdir / ""),
        dryrun=dryrun,
    )

    return analysis_subdir


# ── Step 2: Stitch images ──────────────────────────────────────────────────


def stitch_images(
    ctx: SlideContext,
    celllabels_subdirs: str,
    channel_names: list[str] | None = None,
) -> None:
    cmd = [
        "uv", "run", "stitch-images",
        "-i", str(ctx.work_dir / "CellStatsDir"),
        "-f", str(ctx.work_dir / "RunSummary"),
        "-o", str(ctx.work_dir / "output"),
    ]
    if celllabels_subdirs:
        cmd += ["--celllabels-subdir", celllabels_subdirs]
    for pair in channel_names or []:
        cmd += ["--channel-name", pair]
    run(cmd)


# ── Step 3: Read targets ───────────────────────────────────────────────────


def read_targets(ctx: SlideContext, analysis_subdir: str) -> None:
    run([
        "uv", "run", "read-targets",
        str(ctx.work_dir / "AnalysisResults" / analysis_subdir),
        "-o", str(ctx.work_dir / "output"),
        "--filename", "targets.hdf5",
    ])


# ── Step 4: Generate metadata CSV ──────────────────────────────────────────


def generate_metadata(ctx: SlideContext, seg_uuids: list[str]) -> None:
    script_dir = Path(__file__).resolve().parent
    metadata_script = script_dir / "generate-slide-metadata.py"
    if not metadata_script.exists():
        metadata_script = Path("/app/scripts/generate-slide-metadata.py")

    cmd = [
        "uv", "run", "python", str(metadata_script),
        "--bucket", ctx.bucket,
        "--experiment-prefix", ctx.experiment_prefix,
        "--slide-name", ctx.slide_name,
        "--output", str(ctx.work_dir / "output" / "_metadata.csv"),
    ]
    if seg_uuids:
        cmd += ["--seg-id", ",".join(seg_uuids)]
    run(cmd)


# ── Step 5: Upload results ─────────────────────────────────────────────────


def upload_results(ctx: SlideContext) -> None:
    s3_sync(
        str(ctx.work_dir / "output/"),
        f"s3://{ctx.bucket}/{ctx.output_prefix}/",
    )


def _write_success_marker(ctx: SlideContext) -> None:
    """Write a zero-byte _SUCCESS marker to the S3 output prefix.

    The bucket policy on keene-cosmx-data denies s3:DeleteObject to everything
    except the DELETER role, so we can't toggle a separate _FAILED marker —
    PutObject is the only mutation available. Absence of _SUCCESS means the
    slide isn't done (failed, in-flight, or never run); the ECS task history
    is the authoritative source for the "why" of any failure.
    """
    marker_path = ctx.work_dir / "_SUCCESS"
    marker_path.touch()
    _get_s3().upload_file(str(marker_path), ctx.bucket, f"{ctx.output_prefix}/_SUCCESS")
    log(f"  Status marker: s3://{ctx.bucket}/{ctx.output_prefix}/_SUCCESS")


# ── Main ────────────────────────────────────────────────────────────────────


def process_slide(
    ctx: SlideContext,
    *,
    whatif: bool = False,
    seg_version_override: str = "",
    channel_names: list[str] | None = None,
) -> None:
    bench = Benchmark(ctx=ctx, whatif=whatif)

    log("=== Process Slide ===")
    if whatif:
        log("  Mode:       whatif (dry run)")
    log(f"  Study:      {ctx.study}")
    log(f"  Experiment: {ctx.experiment}")
    log(f"  AtoMx run:  {ctx.atomx_run}")
    log(f"  Slide:      {ctx.slide_name}")
    log(f"  Output:     s3://{ctx.bucket}/{ctx.output_prefix}/")
    log("")

    bench.start_time = now_iso()

    # Clean up any previous run
    if ctx.work_dir.exists():
        shutil.rmtree(ctx.work_dir)
    ctx.work_dir.mkdir(parents=True)

    try:
        if seg_version_override:
            celllabels_subdirs = seg_version_override
            seg_uuids = [_seg_uuid(s.strip()) for s in celllabels_subdirs.split(",")]
            log(f"  Using segmentation override: {celllabels_subdirs}")
        else:
            bench.failed_step = "detect_segmentation"
            log(f"[{now_iso()}] Detecting segmentation versions ...")
            segmentations = detect_all_segmentations(ctx)
            if segmentations:
                celllabels_subdirs = ",".join(s for s, _ in segmentations)
                seg_uuids = [uuid for _, uuid in segmentations]
            else:
                celllabels_subdirs = ""
                seg_uuids = []
                log("  No Segmentation_* subdirs found, using base CellStatsDir")

        bench.start("download")
        log(f"[{now_iso()}] Downloading slide data ...")
        analysis_subdir = download_slide(ctx, celllabels_subdirs, dryrun=whatif)
        bench.end("download")
        dur = bench._duration_seconds("download")
        log(f"[{now_iso()}] Download complete ({dur}s)")

        if whatif:
            # In whatif mode, we've shown what would be synced — skip processing steps
            for step in ("stitch", "read_targets", "metadata", "upload"):
                bench.start(step)
                log(f"[{now_iso()}] [whatif] {step} ... skipped")
                bench.end(step)
        else:
            bench.start("stitch")
            log(f"[{now_iso()}] Stitching images ...")
            stitch_images(ctx, celllabels_subdirs, channel_names=channel_names)
            bench.end("stitch")
            log(f"[{now_iso()}] Stitch complete ({bench._duration_seconds('stitch')}s)")

            bench.start("read_targets")
            log(f"[{now_iso()}] Reading targets ...")
            read_targets(ctx, analysis_subdir)
            bench.end("read_targets")
            log(f"[{now_iso()}] Read targets complete ({bench._duration_seconds('read_targets')}s)")

            bench.start("metadata")
            log(f"[{now_iso()}] Generating metadata CSV ...")
            generate_metadata(ctx, seg_uuids)
            bench.end("metadata")
            log(f"[{now_iso()}] Metadata complete ({bench._duration_seconds('metadata')}s)")

            bench.start("upload")
            log(f"[{now_iso()}] Uploading results ...")
            upload_results(ctx)
            bench.end("upload")
            log(f"[{now_iso()}] Upload complete ({bench._duration_seconds('upload')}s)")

        bench.status = "success"
        bench.failed_step = ""

    except Exception:
        log(f"\n=== FAILED at step: {bench.failed_step or 'unknown'} ===")
        bench.write_and_upload()
        raise

    bench.write_and_upload()
    if not whatif:
        _write_success_marker(ctx)

    download_start, _ = bench._timings["download"]
    _, upload_end = bench._timings["upload"]
    if download_start and upload_end:
        log(f"\n=== Done (total {int(upload_end - download_start)}s) ===")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Process a single CosMx slide: detect segmentation, download, stitch, and upload.",
    )
    parser.add_argument("s3_bucket", help="S3 bucket name")
    parser.add_argument(
        "slide_base_path",
        help="Path within bucket to slide (Study/Experiment/AtoMxRun/DecodedFiles/SlideName/ScanId)",
    )
    parser.add_argument("--whatif", action="store_true", help="Dry run: detect segmentation and show what would be synced")
    parser.add_argument(
        "--segmentation-version",
        default="",
        help="Override segmentation version subdirectory (e.g. Segmentation_uuid_003)",
    )
    parser.add_argument(
        "--channel-name",
        action="append",
        default=[],
        help="Override the kit's channel-to-marker mapping. Repeatable. "
             "Format: CH=MARKER (CH in B/G/Y/R/U). "
             "Example: --channel-name B=AT8 --channel-name G=6E10",
    )
    args = parser.parse_args()

    ctx = SlideContext(bucket=args.s3_bucket, slide_base_path=args.slide_base_path)
    process_slide(
        ctx,
        whatif=args.whatif,
        seg_version_override=args.segmentation_version,
        channel_names=args.channel_name,
    )


if __name__ == "__main__":
    main()
