#!/usr/bin/env python3
"""Discover and process all CosMx slides under an S3 experiment directory.

Discovers slides by walking the S3 hierarchy, then launches one Fargate task
per slide to process them in parallel.

Usage:
    uv run python scripts/process-slides.py s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/
    uv run python scripts/process-slides.py s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/ --skip
    uv run python scripts/process-slides.py s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/ --whatif
    uv run python scripts/process-slides.py s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/ --local
    uv run python scripts/process-slides.py s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/ --benchmark --whatif
"""

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import boto3
from dotenv import load_dotenv

ENV_PATH = Path(__file__).resolve().parent.parent / "fargate" / ".env"

# Fargate hardware profiles for benchmarking: (vCPU units, memory MB).
# Each entry uses the maximum memory for that vCPU tier.
BENCHMARK_PROFILES = [
    ("2048",   "16384"),    #  2 vCPU,   16 GB
    ("4096",  "30720"),    #  4 vCPU,  30 GB
    ("8192",  "61440"),    #  8 vCPU,  60 GB
    ("16384", "122880"),    # 16 vCPU, 120 GB
]


def env(key: str) -> str:
    """Get a required environment variable."""
    value = os.environ.get(key)
    if not value:
        print(f"ERROR: Missing required env var {key}. Check fargate/.env", file=sys.stderr)
        sys.exit(1)
    return value


@dataclass
class Slide:
    """A single CosMx slide identified by its full S3 path components."""
    bucket: str
    base_path: str  # e.g. CosMx-GBM/.../DecodedFiles/SlideName/ScanId

    @property
    def slide_name(self) -> str:
        parts = self.base_path.rstrip("/").split("/")
        return parts[-2]  # SlideName is parent of ScanId

    @property
    def atomx_run(self) -> str:
        parts = self.base_path.rstrip("/").split("/")
        # .../AtoMxRun/DecodedFiles/SlideName/ScanId
        return parts[-4]

    @property
    def output_prefix(self) -> str:
        """The S3 prefix where process-slide.py uploads results.

        Mirrors SlideContext.output_prefix in process-slide.py — must stay in sync
        so discovery's --skip check looks at the same path the worker uploads to.
        """
        parts = self.base_path.rstrip("/").split("/")
        study = parts[0]
        if len(parts) >= 6:
            experiment = parts[1]
            return f"napari-stitched/{study}/{experiment}/{self.atomx_run}/{self.slide_name}"
        return f"napari-stitched/{study}/{self.atomx_run}/{self.slide_name}"

    def __str__(self) -> str:
        return f"s3://{self.bucket}/{self.base_path}"


def parse_s3_uri(uri: str) -> tuple[str, str]:
    """Parse s3://bucket/prefix into (bucket, prefix)."""
    if not uri.startswith("s3://"):
        raise ValueError(f"Expected s3:// URI, got: {uri}")
    rest = uri[5:]
    bucket, _, prefix = rest.partition("/")
    return bucket, prefix.rstrip("/")


_s3 = None


def _get_s3():
    global _s3
    if _s3 is None:
        _s3 = boto3.client("s3")
    return _s3


def s3_ls(bucket: str, prefix: str) -> list[str]:
    """List immediate children (prefixes) under an S3 path. Returns prefix names only."""
    response = _get_s3().list_objects_v2(
        Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/",
    )
    return [
        p["Prefix"].rstrip("/").rsplit("/", 1)[-1]
        for p in response.get("CommonPrefixes", [])
    ]


def discover_slides(bucket: str, experiment_prefix: str) -> list[Slide]:
    """Walk S3 hierarchy to find all slides under an experiment directory.

    Expected structure:
        experiment_prefix/
            AtoMxRun1/
                DecodedFiles/
                    SlideName1/
                        ScanId1/   <-- this is a slide base path
                    SlideName2/
                        ScanId2/
            AtoMxRun2/
                DecodedFiles/
                    ...
    """
    slides = []
    atomx_runs = s3_ls(bucket, experiment_prefix)

    for run in atomx_runs:
        decoded_path = f"{experiment_prefix}/{run}/DecodedFiles"
        slide_names = s3_ls(bucket, decoded_path)

        for slide_name in slide_names:
            if slide_name == "Logs":
                continue
            slide_path = f"{decoded_path}/{slide_name}"
            scan_ids = [d for d in s3_ls(bucket, slide_path) if d != "Logs"]
            for scan_id in scan_ids:
                base_path = f"{experiment_prefix}/{run}/DecodedFiles/{slide_name}/{scan_id}"
                slides.append(Slide(bucket=bucket, base_path=base_path))

    return slides


def is_already_processed(slide: Slide) -> bool:
    """Check if output already exists for this slide in S3."""
    response = _get_s3().list_objects_v2(
        Bucket=slide.bucket, Prefix=slide.output_prefix + "/", MaxKeys=1,
    )
    return response.get("KeyCount", 0) > 0


_ecs = None


def _get_ecs():
    global _ecs
    if _ecs is None:
        _ecs = boto3.client("ecs", region_name=env("AWS_REGION"))
    return _ecs


def launch_fargate_task(
    slide: Slide,
    whatif: bool,
    cpu: str | None = None,
    memory: str | None = None,
    spot: bool = False,
    segmentation_version: str | None = None,
    channel_names: list[str] | None = None,
) -> str | None:
    """Launch a Fargate task for a single slide. Returns task ID or None for whatif.

    When cpu/memory are provided they override the task definition values,
    allowing the same task definition to be used across different Fargate sizes.
    When spot is True, uses FARGATE_SPOT capacity provider instead of FARGATE.
    """
    command = ["uv", "run", "python", "/app/scripts/process-slide.py"]
    if segmentation_version:
        command += ["--segmentation-version", segmentation_version]
    for pair in channel_names or []:
        command += ["--channel-name", pair]
    command += [slide.bucket, slide.base_path]

    cluster = env("ECS_CLUSTER")
    subnets = env("FARGATE_SUBNETS").split(",")
    security_group = env("FARGATE_SECURITY_GROUP")

    overrides: dict = {
        "containerOverrides": [{
            "name": "process-slide",
            "command": command,
        }],
    }
    if cpu is not None:
        overrides["cpu"] = cpu
    if memory is not None:
        overrides["memory"] = memory

    provider = "FARGATE_SPOT" if spot else "FARGATE"

    if whatif:
        size_info = ""
        if cpu is not None:
            vcpu = int(cpu) // 1024
            mem_gb = int(memory) // 1024
            size_info = f" ({vcpu} vCPU, {mem_gb} GB)"
        print(f"  [whatif] would launch {provider} task{size_info} with command: {' '.join(command)}")
        return None

    kwargs = {
        "cluster": cluster,
        "taskDefinition": "cosmx-process-slide",
        "networkConfiguration": {
            "awsvpcConfiguration": {
                "subnets": subnets,
                "securityGroups": [security_group],
                "assignPublicIp": "DISABLED",
            },
        },
        "overrides": overrides,
    }
    if spot:
        kwargs["capacityProviderStrategy"] = [
            {"capacityProvider": "FARGATE_SPOT", "weight": 1},
        ]
    else:
        kwargs["launchType"] = "FARGATE"

    response = _get_ecs().run_task(**kwargs)
    task_arn = response["tasks"][0]["taskArn"]
    task_id = task_arn.split("/")[-1]
    return task_id


def process_slide_local(
    slide: Slide,
    whatif: bool,
    segmentation_version: str | None = None,
    channel_names: list[str] | None = None,
) -> None:
    """Run process-slide.py locally for a single slide."""
    cmd = ["uv", "run", "python", "scripts/process-slide.py"]
    if whatif:
        cmd.append("--whatif")
    if segmentation_version:
        cmd += ["--segmentation-version", segmentation_version]
    for pair in channel_names or []:
        cmd += ["--channel-name", pair]
    cmd += [slide.bucket, slide.base_path]

    if whatif:
        print(f"  [whatif] would run: {' '.join(cmd)}")
        return

    print(f"  Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover and process all CosMx slides under an S3 experiment directory.",
    )
    parser.add_argument(
        "s3_uri",
        help="S3 URI of the experiment directory (e.g. s3://my-bucket/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/)",
    )
    parser.add_argument(
        "--whatif",
        action="store_true",
        help="Print commands that would be run without executing them.",
    )
    parser.add_argument(
        "--skip",
        action="store_true",
        help="Skip slides that already have output in S3.",
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="Run slides locally and sequentially instead of on Fargate.",
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Launch one Fargate task per hardware profile for the first slide only.",
    )
    parser.add_argument(
        "--spot",
        action="store_true",
        help="Use FARGATE_SPOT capacity provider for cheaper (interruptible) tasks.",
    )
    parser.add_argument(
        "--cpu",
        default=None,
        help="Override Fargate vCPU units (e.g. 16384 for 16 vCPU). "
             "When omitted, uses the task definition default.",
    )
    parser.add_argument(
        "--memory",
        default=None,
        help="Override Fargate memory in MB (e.g. 122880 for 120 GB). "
             "When omitted, uses the task definition default.",
    )
    parser.add_argument(
        "--segmentation-version",
        default=None,
        help="Override segmentation version subdirectory (e.g. Segmentation_uuid_003). "
             "When omitted, auto-detects the highest version from manifest JSONs.",
    )
    parser.add_argument(
        "--channel-name",
        action="append",
        default=[],
        help="Override the kit's channel-to-marker mapping for stitching. "
             "Repeatable. Format: CH=MARKER (CH in B/G/Y/R/U). "
             "Example: --channel-name B=AT8 --channel-name G=6E10. "
             "Use when MorphologyKit metadata in the TIFFs is wrong.",
    )
    args = parser.parse_args()

    if not args.local:
        if not ENV_PATH.exists():
            print(f"ERROR: {ENV_PATH} not found. Copy fargate/.env.example to fargate/.env and fill in your values.", file=sys.stderr)
            sys.exit(1)
        load_dotenv(ENV_PATH)

    bucket, prefix = parse_s3_uri(args.s3_uri)

    print(f"Discovering slides under s3://{bucket}/{prefix}/ ...")
    slides = discover_slides(bucket, prefix)

    if not slides:
        print("No slides found.")
        sys.exit(1)

    print(f"Found {len(slides)} slide(s):\n")
    for slide in slides:
        print(f"  {slide.atomx_run} / {slide.slide_name}")

    if args.skip:
        before = len(slides)
        slides = [s for s in slides if not is_already_processed(s)]
        skipped = before - len(slides)
        if skipped:
            print(f"\nSkipping {skipped} already-processed slide(s).")
        if not slides:
            print("All slides already processed.")
            return

    if args.benchmark:
        slide = slides[0]
        print(f"\nBenchmarking slide: {slide.atomx_run} / {slide.slide_name}")
        print(f"Launching {len(BENCHMARK_PROFILES)} tasks (one per hardware profile):\n")

        task_ids = []
        for cpu, memory in BENCHMARK_PROFILES:
            vcpu = int(cpu) // 1024
            mem_gb = int(memory) // 1024
            label = f"{vcpu} vCPU / {mem_gb} GB"
            print(f"  {label}:")
            task_id = launch_fargate_task(slide, whatif=args.whatif, cpu=cpu, memory=memory, spot=args.spot, segmentation_version=args.segmentation_version, channel_names=args.channel_name)
            if task_id:
                task_ids.append((label, task_id))
                print(f"    Task: {task_id}")
            print()

        if task_ids:
            cluster = env("ECS_CLUSTER")
            region = env("AWS_REGION")
            print("All benchmark tasks launched. Monitor with:")
            print(f"  aws ecs list-tasks --cluster {cluster} --region {region}")
            print(f"\nOr watch a specific task:")
            for label, task_id in task_ids:
                print(f"  {label}: aws ecs describe-tasks --cluster {cluster} --tasks {task_id} --region {region} --query 'tasks[0].lastStatus'")
    else:
        print(f"\n{len(slides)} slide(s) to process:\n")

        if args.local:
            for slide in slides:
                print(f"Processing: {slide.atomx_run} / {slide.slide_name}")
                process_slide_local(slide, whatif=args.whatif, segmentation_version=args.segmentation_version, channel_names=args.channel_name)
                print()
        else:
            task_ids = []
            for slide in slides:
                print(f"Launching: {slide.atomx_run} / {slide.slide_name}")
                task_id = launch_fargate_task(slide, whatif=args.whatif, cpu=args.cpu, memory=args.memory, spot=args.spot, segmentation_version=args.segmentation_version, channel_names=args.channel_name)
                if task_id:
                    task_ids.append((slide, task_id))
                    print(f"  Task: {task_id}")
                print()

            if task_ids:
                cluster = env("ECS_CLUSTER")
                region = env("AWS_REGION")
                print("All tasks launched. Monitor with:")
                print(f"  aws ecs list-tasks --cluster {cluster} --region {region}")
                print(f"\nOr watch a specific task:")
                for slide, task_id in task_ids:
                    print(f"  aws ecs describe-tasks --cluster {cluster} --tasks {task_id} --region {region} --query 'tasks[0].lastStatus'")

    print("\nDone.")


if __name__ == "__main__":
    main()
