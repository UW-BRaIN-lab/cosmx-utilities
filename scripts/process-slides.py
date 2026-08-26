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
import json
import urllib.request
from dotenv import load_dotenv

ENV_PATH = Path(__file__).resolve().parent.parent / "fargate" / ".env"
TASK_DEFINITION_FAMILY = "cosmx-process-slide"
TASK_DEFINITION_TEMPLATE = (
    Path(__file__).resolve().parent.parent
    / "fargate" / "fargate-task-process-slide.json"
)
# Fields pinned by the checked-in template that a run actually depends on.
# cpu/memory are omitted deliberately: --cpu/--memory override them per run.
TASK_DEFINITION_PINNED = ("image", "taskRoleArn", "executionRoleArn",
                          "ephemeralStorage")

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

    The prefix may also point straight at a single AtoMx run (one that contains
    DecodedFiles directly). Two studies of the same slides often need different
    worker flags — a 3D resegmentation takes --input-ndim 3 while the original
    2D run does not — so they have to be launched separately. Slide base paths,
    and therefore output prefixes, are identical either way.
    """
    slides = []
    children = s3_ls(bucket, experiment_prefix)
    single_run = "DecodedFiles" in children
    run_prefixes = [experiment_prefix] if single_run else [
        f"{experiment_prefix}/{run}" for run in children
    ]

    for run_prefix in run_prefixes:
        decoded_path = f"{run_prefix}/DecodedFiles"
        slide_names = s3_ls(bucket, decoded_path)

        for slide_name in slide_names:
            if slide_name == "Logs":
                continue
            slide_path = f"{decoded_path}/{slide_name}"
            scan_ids = [d for d in s3_ls(bucket, slide_path) if d != "Logs"]
            for scan_id in scan_ids:
                base_path = f"{decoded_path}/{slide_name}/{scan_id}"
                slides.append(Slide(bucket=bucket, base_path=base_path))

    return slides


def _registry_digest(image_ref: str) -> str:
    """Resolve a registry image reference to the digest it currently points at.

    The task definition names a moving tag, so the tag alone says nothing about
    which code a task will run. GHCR serves the digest for a public package
    without credentials. Returns "" when it cannot be resolved, which callers
    treat as "cannot tell".
    """
    if "/" not in image_ref:
        return ""
    registry, _, remainder = image_ref.partition("/")
    repository, _, tag = remainder.partition(":")
    tag = tag or "latest"
    if registry != "ghcr.io":
        return ""
    try:
        token_url = (f"https://ghcr.io/token?scope=repository:{repository}:pull"
                     f"&service=ghcr.io")
        with urllib.request.urlopen(token_url, timeout=10) as response:
            token = json.load(response).get("token", "")
        if not token:
            return ""
        request = urllib.request.Request(
            f"https://ghcr.io/v2/{repository}/manifests/{tag}",
            headers={
                "Authorization": f"Bearer {token}",
                "Accept": ", ".join([
                    "application/vnd.oci.image.index.v1+json",
                    "application/vnd.docker.distribution.manifest.list.v2+json",
                    "application/vnd.docker.distribution.manifest.v2+json",
                ]),
            },
            method="HEAD",
        )
        with urllib.request.urlopen(request, timeout=10) as response:
            return response.headers.get("Docker-Content-Digest", "")
    except Exception as e:
        print(f"WARNING: could not resolve digest for {image_ref}: {e}", file=sys.stderr)
        return ""


def expected_code_version() -> str:
    """The code version a task launched now would run, or "" if undeterminable."""
    override = os.environ.get("COSMX_CODE_VERSION")
    if override:
        return override
    try:
        task_def = _get_ecs().describe_task_definition(
            taskDefinition=TASK_DEFINITION_FAMILY)["taskDefinition"]
        image = (task_def.get("containerDefinitions") or [{}])[0].get("image", "")
    except Exception as e:
        print(f"WARNING: could not read task definition: {e}", file=sys.stderr)
        return ""
    return _registry_digest(image)


def read_success_marker(slide: Slide) -> dict | None:
    """The slide's _SUCCESS marker, or None when absent or unreadable."""
    try:
        body = _get_s3().get_object(
            Bucket=slide.bucket, Key=f"{slide.output_prefix}/_SUCCESS",
        )["Body"].read()
    except Exception:
        return None
    if not body.strip():
        return {}          # legacy zero-byte marker: finished, provenance unknown
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return {}


def is_already_processed(slide: Slide, expected_version: str = "") -> bool:
    """Whether this slide's output can be trusted as already done.

    A marker on its own only says some run finished here. That is what let a
    run which completed but wrote corrupt data keep its _SUCCESS, so --skip
    skipped the slide and protected the bad output from being replaced.

    So when the marker records which code produced it, and that differs from
    the code a task would run now, the slide is not treated as done. Where
    either side is unknown -- a legacy marker, or a digest that cannot be
    resolved -- fall back to presence, which is the old behaviour.
    """
    marker = read_success_marker(slide)
    if marker is None:
        return False
    recorded = marker.get("code_version", "")
    if not recorded or not expected_version:
        return True
    if recorded == expected_version:
        return True
    print(f"  {slide.slide_name}: output was produced by different code "
          f"({recorded[:20]}... vs {expected_version[:20]}...), not skipping")
    return False


_ecs = None


def _get_ecs():
    global _ecs
    if _ecs is None:
        _ecs = boto3.client("ecs", region_name=env("AWS_REGION"))
    return _ecs


def worker_flags(args) -> list[str]:
    """Per-slide flags forwarded verbatim to process-slide.py.

    Kept in one place so the Fargate and --local paths cannot drift apart.
    """
    flags: list[str] = []
    if args.segmentation_version:
        flags += ["--segmentation-version", args.segmentation_version]
    for pair in args.channel_name or []:
        flags += ["--channel-name", pair]
    if args.input_ndim is not None:
        flags += ["--input-ndim", str(args.input_ndim)]
    if args.output_ndim is not None:
        flags += ["--output-ndim", str(args.output_ndim)]
    for column in args.column or []:
        flags += ["--column", column]
    if args.annotations_prefix:
        flags += ["--annotations-prefix", args.annotations_prefix]
    if args.fill_from:
        flags += ["--fill-from", args.fill_from]
    return flags


def _template_task_definition() -> dict | None:
    """The checked-in task definition, with ACCOUNT_ID substituted."""
    if not TASK_DEFINITION_TEMPLATE.exists():
        return None
    raw = TASK_DEFINITION_TEMPLATE.read_text().replace(
        "ACCOUNT_ID", env("AWS_ACCOUNT_ID"))
    return json.loads(raw)


def _pinned_fields(task_def: dict) -> dict:
    """Pull the pinned fields out of a task definition, container ones included."""
    container = (task_def.get("containerDefinitions") or [{}])[0]
    values = {}
    for field in TASK_DEFINITION_PINNED:
        if field in container:
            values[field] = container[field]
        elif field in task_def:
            values[field] = task_def[field]
    return values


def task_definition_drift() -> list[tuple[str, object, object]]:
    """Differences between the registered task definition and the repo's template.

    Nothing keeps the two in sync: `register-task-defs.sh` pushes the template
    but nothing re-runs it, so the live definition silently ages. That is not
    hypothetical -- the registered definition kept pointing at a pre-rename
    `ghcr.io/keene-lab/...` image long after the repo moved, and every task
    launched against it died on an image pull it could not authorize. The
    failure surfaces a minute into a run, per task, with an error that says
    nothing about drift.

    Returns [(field, registered, expected)], empty when they agree or when the
    comparison cannot be made.
    """
    template = _template_task_definition()
    if template is None:
        return []
    try:
        registered = _get_ecs().describe_task_definition(
            taskDefinition=TASK_DEFINITION_FAMILY)["taskDefinition"]
    except Exception as e:
        print(f"WARNING: could not read task definition {TASK_DEFINITION_FAMILY}: {e}",
              file=sys.stderr)
        return []

    expected = _pinned_fields(template)
    active = _pinned_fields(registered)
    return [(field, active.get(field), value)
            for field, value in expected.items() if active.get(field) != value]


def launch_fargate_task(
    slide: Slide,
    whatif: bool,
    cpu: str | None = None,
    memory: str | None = None,
    spot: bool = False,
    extra_flags: list[str] | None = None,
) -> str | None:
    """Launch a Fargate task for a single slide. Returns task ID or None for whatif.

    When cpu/memory are provided they override the task definition values,
    allowing the same task definition to be used across different Fargate sizes.
    When spot is True, uses FARGATE_SPOT capacity provider instead of FARGATE.
    """
    command = ["uv", "run", "python", "/app/scripts/process-slide.py"]
    command += extra_flags or []
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
    extra_flags: list[str] | None = None,
) -> None:
    """Run process-slide.py locally for a single slide."""
    cmd = ["uv", "run", "python", "scripts/process-slide.py"]
    if whatif:
        cmd.append("--whatif")
    cmd += extra_flags or []
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
    parser.add_argument(
        "--input-ndim",
        type=int,
        choices=[2, 3],
        default=None,
        help="Dimensionality of the CellLabels TIFFs. Use 3 for a 3D "
             "resegmentation, whose labels are per-z (CellLabels_F###_Z###.tif).",
    )
    parser.add_argument(
        "--output-ndim",
        type=int,
        choices=[2, 3],
        default=None,
        help="Dimensionality of the stitched output. Use 3 to give napari a "
             "z-navigable labels layer; requires --input-ndim 3.",
    )
    parser.add_argument(
        "--column",
        action="append",
        default=[],
        metavar="OUTNAME=SOURCE_HEADER",
        help="Annotation column to carry into _metadata.csv (repeatable). "
             "SOURCE_HEADER may list fallbacks separated by '|' for studies "
             "that renamed a column between AtoMx runs.",
    )
    parser.add_argument(
        "--allow-taskdef-drift",
        action="store_true",
        help="Launch even if the registered task definition differs from "
             "fargate/fargate-task-process-slide.json.",
    )
    parser.add_argument(
        "--annotations-prefix",
        default="",
        metavar="S3URI",
        help="S3 directory of per-FOV annotation sheets named "
             "<slide>_annotations.csv, one per slide. Fills annotation values "
             "AtoMx never captured. Must be readable from the Fargate task.",
    )
    parser.add_argument(
        "--fill-from",
        default="",
        metavar="EXPERIMENT_PREFIX",
        help="Fill annotation values AtoMx left blank from another study of the "
             "same physical slides, joining on FOV.",
    )
    args = parser.parse_args()

    if not args.local:
        if not ENV_PATH.exists():
            print(f"ERROR: {ENV_PATH} not found. Copy fargate/.env.example to fargate/.env and fill in your values.", file=sys.stderr)
            sys.exit(1)
        load_dotenv(ENV_PATH)

    flags = worker_flags(args)

    # Catch a stale registration before spending tasks on it: every task would
    # otherwise fail a minute in, individually, on an error that never mentions
    # the task definition.
    if not args.local and not args.whatif:
        drift = task_definition_drift()
        if drift:
            print(f"ERROR: registered task definition '{TASK_DEFINITION_FAMILY}' "
                  f"does not match fargate/fargate-task-process-slide.json:",
                  file=sys.stderr)
            for field, registered, expected in drift:
                print(f"  {field}:", file=sys.stderr)
                print(f"      registered: {registered}", file=sys.stderr)
                print(f"      expected  : {expected}", file=sys.stderr)
            print("\nRe-register it with:  ./fargate/register-task-defs.sh",
                  file=sys.stderr)
            print("Or pass --allow-taskdef-drift to launch anyway.", file=sys.stderr)
            if not args.allow_taskdef_drift:
                sys.exit(1)
            print("Continuing despite drift (--allow-taskdef-drift).", file=sys.stderr)

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
        expected_version = expected_code_version()
        if not expected_version:
            print("  (could not determine the running code version; skipping on "
                  "marker presence alone)")
        slides = [s for s in slides if not is_already_processed(s, expected_version)]
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
            task_id = launch_fargate_task(slide, whatif=args.whatif, cpu=cpu, memory=memory, spot=args.spot, extra_flags=flags)
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
                process_slide_local(slide, whatif=args.whatif, extra_flags=flags)
                print()
        else:
            task_ids = []
            for slide in slides:
                print(f"Launching: {slide.atomx_run} / {slide.slide_name}")
                task_id = launch_fargate_task(slide, whatif=args.whatif, cpu=args.cpu, memory=args.memory, spot=args.spot, extra_flags=flags)
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
