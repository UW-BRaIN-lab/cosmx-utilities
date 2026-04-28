#!/usr/bin/env python3
"""Launch a CosMx analytics EC2 instance from a custom AMI.

Spins up an r5a.4xlarge (or custom type) with a configurable secondary
EBS data volume that is auto-formatted and mounted at /mnt/cosmx.

Usage:
    # CosMx analytics (memory-optimized, default):
    uv run python ec2/start_ec2.py --name emily-gbm-analysis
    uv run python ec2/start_ec2.py --name emily-gbm-analysis --storage 1024
    uv run python ec2/start_ec2.py --name emily-gbm-analysis --instance-type r5a.8xlarge

    # Napari viewer (GPU + local NVMe for fast slide loading):
    uv run python ec2/start_ec2.py --name emily-napari --napari --raw
    # Auto-selects instance size based on S3 data (g4dn.4xlarge ≤200 GB, g4dn.8xlarge ≤850 GB):
    uv run python ec2/start_ec2.py --name emily-napari --napari --raw --s3 s3://bucket/napari-stitched/study/experiment/

    # Test ami_setup.sh on a raw Ubuntu instance (no AMI needed):
    uv run python ec2/start_ec2.py --name test-ami-setup --raw

    # Test an unmerged branch on the instance (must be pushed to origin):
    uv run python ec2/start_ec2.py --name emily-feature-test --napari --raw --branch my-feature
"""

import argparse
import os
import subprocess
import sys

from botocore.exceptions import ClientError, WaiterError
from dotenv import load_dotenv

from _common import (
    AMI_SETUP_SCRIPT,
    ENV_PATH,
    PROJECT_TAG,
    ROOT_VOLUME_GB,
    boto_session,
    env,
    log,
)

ANALYTICS_INSTANCE_TYPE = "r5a.4xlarge"
DEFAULT_STORAGE_GB = 512

# Napari instance tiers — ordered smallest to largest.
# Each tier maps an instance type to its (total NVMe GB, usable GB) after
# reserving space for OS overhead, filesystem metadata, and temp files.
NAPARI_TIERS = [
    ("g4dn.4xlarge", 225, 200),
    ("g4dn.8xlarge", 900, 850),
]
NAPARI_DEFAULT_INSTANCE_TYPE = NAPARI_TIERS[0][0]

# Bash snippet to find, format, and mount the secondary NVMe data volume.
# On Nitro instances (r5a), /dev/sdf appears as /dev/nvme1n1.
MOUNT_VOLUME_SNIPPET = """\
# ── Mount secondary data volume ──────────────────────────────────────────
DEVICE=""
for d in /dev/nvme*n1; do
    ROOT_DEV=$(findmnt -n -o SOURCE /)
    ROOT_NVME=$(readlink -f "$ROOT_DEV" | sed 's/p[0-9]*$//')
    if [ "$d" != "$ROOT_NVME" ]; then
        DEVICE="$d"
        break
    fi
done

if [ -z "$DEVICE" ]; then
    echo "No secondary NVMe device found" >> /var/log/mount-data.log
else
    if ! blkid "$DEVICE" | grep -q TYPE; then
        mkfs.ext4 "$DEVICE"
    fi
    mkdir -p /mnt/cosmx
    mount "$DEVICE" /mnt/cosmx
    chown ubuntu:ubuntu /mnt/cosmx
    UUID=$(blkid -s UUID -o value "$DEVICE")
    if ! grep -q "$UUID" /etc/fstab; then
        echo "UUID=$UUID /mnt/cosmx ext4 defaults,nofail 0 2" >> /etc/fstab
    fi
    echo "Data volume mounted at /mnt/cosmx" >> /var/log/mount-data.log
fi
"""

# Bash snippet to format and mount the local NVMe instance store on g4dn.
# The instance store is ephemeral — data is lost on stop/terminate.
# ~225 GB fast SSD mounted at /mnt/local for active slide viewing.
# Since Napari mode has no secondary EBS, the only non-root NVMe device
# is the instance store — no need for nvme-cli to distinguish them.
MOUNT_VOLUME_NAPARI_SNIPPET = """\
# ── Mount local NVMe instance store ──────────────────────────────────────
DEVICE=""
for d in /dev/nvme*n1; do
    ROOT_DEV=$(findmnt -n -o SOURCE /)
    ROOT_NVME=$(readlink -f "$ROOT_DEV" | sed 's/p[0-9]*$//')
    if [ "$d" != "$ROOT_NVME" ]; then
        DEVICE="$d"
        break
    fi
done

if [ -z "$DEVICE" ]; then
    echo "No NVMe instance store found" >> /var/log/mount-data.log
else
    mkfs.ext4 -F "$DEVICE"
    mkdir -p /mnt/local
    mount "$DEVICE" /mnt/local
    chown ubuntu:ubuntu /mnt/local
    echo "Instance store $DEVICE mounted at /mnt/local" >> /var/log/mount-data.log
fi
"""

MOUNT_USER_DATA = f"""\
#!/bin/bash
set -euo pipefail
{MOUNT_VOLUME_SNIPPET}"""

MOUNT_USER_DATA_NAPARI = f"""\
#!/bin/bash
set -euo pipefail
{MOUNT_VOLUME_NAPARI_SNIPPET}"""


def s3_total_size_gb(s3_uri: str, region: str) -> float:
    """Return total size in GB of all objects under an S3 prefix."""
    session = boto_session(region)
    s3 = session.client("s3")

    # Parse s3://bucket/prefix
    without_scheme = s3_uri[len("s3://"):]
    bucket, _, prefix = without_scheme.partition("/")

    total_bytes = 0
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            total_bytes += obj["Size"]

    return total_bytes / (1024 ** 3)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Launch a CosMx analytics EC2 instance",
    )
    parser.add_argument(
        "--name",
        required=True,
        help="Instance name tag (e.g. emily-gbm-analysis)",
    )
    parser.add_argument(
        "--storage",
        type=int,
        default=DEFAULT_STORAGE_GB,
        help=f"Secondary EBS data volume size in GB (default: {DEFAULT_STORAGE_GB})",
    )
    parser.add_argument(
        "--instance-type",
        help="EC2 instance type (overrides default for chosen mode)",
    )
    parser.add_argument(
        "--napari",
        action="store_true",
        help="Launch a Napari viewer instance with GPU + local NVMe for "
             "fast slide loading (auto-selects instance size with --s3)",
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Launch a raw Ubuntu instance with ami_setup.sh as user-data "
             "(for testing setup without building an AMI)",
    )
    parser.add_argument(
        "--branch",
        default="main",
        help="Git branch to clone on the instance (default: main). "
             "Must exist on origin — validated before launch.",
    )
    parser.add_argument(
        "--s3",
        metavar="S3_URI",
        help="S3 path to sync to the local NVMe after mount "
             "(e.g. s3://bucket/napari-stitched/study/experiment/). "
             "Requires --raw and --napari.",
    )
    args = parser.parse_args()

    if args.s3 and not (args.raw and args.napari):
        print("ERROR: --s3 requires both --raw and --napari", file=sys.stderr)
        sys.exit(1)

    if args.s3 and not args.s3.startswith("s3://"):
        print("ERROR: --s3 must be an S3 URI (s3://...)", file=sys.stderr)
        sys.exit(1)

    if args.napari and args.storage != DEFAULT_STORAGE_GB:
        print("WARNING: --storage is ignored with --napari (no EBS data volume)", file=sys.stderr)

    load_dotenv(ENV_PATH)
    region = env("AWS_REGION")
    subnet = env("EC2_SUBNET")
    sg = env("EC2_SECURITY_GROUP")
    key_pair = env("EC2_KEY_PAIR")
    profile = env("EC2_INSTANCE_PROFILE")

    dcv_password = os.environ.get("DCV_PASSWORD")

    # For Napari + S3: pick the smallest instance whose NVMe fits the data,
    # unless the user explicitly chose an instance type.
    napari_nvme_total_gb = NAPARI_TIERS[0][1]
    napari_nvme_usable_gb = NAPARI_TIERS[0][2]

    if args.s3:
        log(f"Checking S3 size of {args.s3} ...")
        s3_size = s3_total_size_gb(args.s3, region)

        if s3_size == 0:
            print(f"ERROR: No objects found at {args.s3}", file=sys.stderr)
            sys.exit(1)

        if not args.instance_type:
            # Auto-select the smallest Napari tier that fits
            selected = None
            for instance_type, total_gb, usable_gb in NAPARI_TIERS:
                if s3_size <= usable_gb:
                    selected = (instance_type, total_gb, usable_gb)
                    break
            if selected is None:
                max_tier = NAPARI_TIERS[-1]
                print(
                    f"ERROR: S3 data ({s3_size:.1f} GB) exceeds largest available "
                    f"NVMe capacity (~{max_tier[2]} GB on {max_tier[0]})",
                    file=sys.stderr,
                )
                sys.exit(1)
            args.instance_type, napari_nvme_total_gb, napari_nvme_usable_gb = selected
            log(f"  S3 data size: {s3_size:.1f} GB → auto-selected {args.instance_type} "
                f"(NVMe: ~{napari_nvme_usable_gb} GB usable of {napari_nvme_total_gb} GB)")
        else:
            # User chose a specific instance type — validate against the
            # default tier (we can't know arbitrary instance NVMe sizes).
            log(f"  S3 data size: {s3_size:.1f} GB (user-selected {args.instance_type})")
            log(f"  WARNING: Cannot verify NVMe capacity for custom instance type. "
                "Ensure it has enough local storage.")

    if not args.instance_type:
        args.instance_type = NAPARI_DEFAULT_INSTANCE_TYPE if args.napari else ANALYTICS_INSTANCE_TYPE

    # Validate the branch exists on origin before launch — silent failures
    # inside cloud-init waste 10+ minutes of EC2 time before the user finds out.
    git_branch = args.branch
    try:
        subprocess.check_call(
            ["git", "ls-remote", "--exit-code", "origin", f"refs/heads/{git_branch}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        print("ERROR: git is not installed", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError:
        print(
            f"ERROR: branch '{git_branch}' not found on origin. "
            f"Push the branch first or pass --branch <existing>.",
            file=sys.stderr,
        )
        sys.exit(1)

    mount_snippet = MOUNT_VOLUME_NAPARI_SNIPPET if args.napari else MOUNT_VOLUME_SNIPPET

    if args.raw:
        ami_id = env("UBUNTU_BASE_AMI")
        if not AMI_SETUP_SCRIPT.exists():
            print(f"ERROR: {AMI_SETUP_SCRIPT} not found", file=sys.stderr)
            sys.exit(1)
        setup_script = AMI_SETUP_SCRIPT.read_text()
        user_data = setup_script.rstrip() + "\n\n" + mount_snippet
        if args.s3:
            # rclone with :s3: backend uses the instance's IAM role automatically.
            # --transfers=32 for parallel file downloads, --s3-no-head to skip
            # unnecessary HEAD requests on each object.
            rclone_src = ":s3:" + args.s3[len("s3://"):]
            s3_sync_snippet = f"""
# ── Sync stitched slides from S3 to local NVMe ──────────────────────────
echo "Syncing from {args.s3} to /mnt/local/stitched ..."
sudo -u ubuntu mkdir -p /mnt/local/stitched
sudo -u ubuntu rclone sync '{rclone_src}' /mnt/local/stitched \
    --s3-provider AWS --s3-region {region} --s3-env-auth \
    --transfers 32 --checkers 16 --s3-no-head \
    --progress >> /var/log/ami-setup.log 2>&1
echo "S3 sync complete. Files:" >> /var/log/ami-setup.log
ls -lh /mnt/local/stitched >> /var/log/ami-setup.log 2>&1
"""
            user_data = user_data.rstrip() + "\n" + s3_sync_snippet
    else:
        ami_id = env("EC2_AMI_ID")
        user_data = MOUNT_USER_DATA_NAPARI if args.napari else MOUNT_USER_DATA

    # Inject variables after the set -e line in user-data
    ec2_mode = "napari" if args.napari else "analytics"
    injected = f'export GIT_BRANCH="{git_branch}"\nexport EC2_MODE="{ec2_mode}"\n'
    if dcv_password:
        injected += f'echo "ubuntu:{dcv_password}" | chpasswd\n'
    for marker in ("set -euxo pipefail\n", "set -euo pipefail\n"):
        if marker in user_data:
            user_data = user_data.replace(marker, marker + injected, 1)
            break

    session = boto_session(region)
    ec2 = session.client("ec2")

    mode = "RAW (ami_setup.sh)" if args.raw else "AMI"
    variant = "Napari viewer" if args.napari else "CosMx analytics"
    log(f"Launching {variant} instance [{mode}]")
    log(f"  AMI:            {ami_id}")
    log(f"  Instance type:  {args.instance_type}")
    if args.napari:
        log(f"  Local NVMe:     ~{napari_nvme_total_gb} GB (instance store, mounted at /mnt/local)")
        if args.s3:
            log(f"  S3 sync:        {args.s3} → /mnt/local/stitched")
    else:
        log(f"  Data volume:    {args.storage} GB")
    log(f"  Git branch:     {git_branch}")
    log(f"  Name:           {args.name}")

    resource_tags = [
        {"Key": "Name", "Value": args.name},
        {"Key": "Project", "Value": PROJECT_TAG},
        {"Key": "Purpose", "Value": "napari-viewer" if args.napari else "analytics"},
    ]

    block_devices = [
        {
            "DeviceName": "/dev/sda1",
            "Ebs": {"VolumeSize": ROOT_VOLUME_GB, "VolumeType": "gp3"},
        },
    ]
    if args.napari:
        # Local NVMe instance store only — no EBS data volume needed
        block_devices.append(
            {"DeviceName": "/dev/sdb", "VirtualName": "ephemeral0"}
        )
    else:
        block_devices.append(
            {
                "DeviceName": "/dev/sdf",
                "Ebs": {
                    "VolumeSize": args.storage,
                    "VolumeType": "gp3",
                    "DeleteOnTermination": True,
                },
            },
        )

    try:
        response = ec2.run_instances(
            ImageId=ami_id,
            InstanceType=args.instance_type,
            KeyName=key_pair,
            SubnetId=subnet,
            SecurityGroupIds=[sg],
            IamInstanceProfile={"Name": profile},
            UserData=user_data,
            BlockDeviceMappings=block_devices,
            TagSpecifications=[
                {"ResourceType": "instance", "Tags": resource_tags},
                {"ResourceType": "volume", "Tags": resource_tags},
            ],
            MinCount=1,
            MaxCount=1,
        )

        instance_id = response["Instances"][0]["InstanceId"]
        log(f"Instance launched: {instance_id}")

        log("Waiting for instance to start...")
        ec2.get_waiter("instance_running").wait(InstanceIds=[instance_id])

        desc = ec2.describe_instances(InstanceIds=[instance_id])
        private_ip = desc["Reservations"][0]["Instances"][0].get("PrivateIpAddress", "N/A")

    except (ClientError, WaiterError) as e:
        log(f"AWS error: {e}")
        sys.exit(1)

    log("Instance is running")
    print()
    print(f"  Instance ID:  {instance_id}")
    print(f"  Private IP:   {private_ip}")
    if args.napari:
        print(f"  Local NVMe:   ~{napari_nvme_total_gb} GB (auto-mounted at /mnt/local)")
    else:
        print(f"  Data volume:  {args.storage} GB (auto-mounted at /mnt/cosmx)")
    print()

    if args.raw:
        print("NOTE: --raw mode — ami_setup.sh is running as user-data.")
        print("  Setup takes ~15-30 minutes. Monitor progress:")
        print(f"  aws ssm start-session --target {instance_id} --region {region}")
        print("  Then: tail -f /var/log/ami-setup.log")
        print()

    print("Connect via SSM:")
    print(f"  aws ssm start-session --target {instance_id} --region {region}")
    print()
    print("Connect via DCV (remote desktop on port 8443):")
    if dcv_password:
        print(f"  Open DCV client and connect to {private_ip}")
    else:
        print(f"  Set a password first: aws ssm start-session --target {instance_id}")
        print("  Then run: sudo passwd ubuntu")
        print(f"  Open DCV client and connect to {private_ip}")
    print()

    if args.napari:
        if args.s3:
            print("S3 sync is running automatically as part of instance setup.")
            print("  Monitor progress: tail -f /var/log/ami-setup.log")
            print(f"  Files will appear in /mnt/local/stitched")
            print()
            print("Launch Napari (on the instance via DCV):")
            print("  cd /opt/cosmx-utilities && uv run napari /mnt/local/stitched")
        else:
            print("Copy slide to local NVMe for fast loading:")
            print(f"  # On the instance — sync from S3 to local NVMe:")
            print("  aws s3 sync s3://your-bucket/napari-stitched/study/experiment/ /mnt/local/slide/")
            print()
            print("Launch Napari (on the instance via DCV):")
            print("  cd /opt/cosmx-utilities && uv run napari /mnt/local/slide")
        print()
        print("NOTE: /mnt/local is ephemeral instance storage — data is lost when")
        print("  the instance stops or terminates. Source data lives on S3.")
        print()
    else:
        print("Sync data from S3:")
        print(f"  ./scripts/sync-to-ec2.sh s3://your-bucket/napari-stitched/study/experiment/ {instance_id}")
        print()
        print("Launch Jupyter (on the instance):")
        print("  cd /opt/cosmx-utilities/ec2/analytics && uv run jupyter lab --ip 0.0.0.0 --port 8888")
        print()
        print("Launch Napari (on the instance via DCV):")
        print("  cd /opt/cosmx-utilities && uv run napari /mnt/cosmx")
        print()

    print("When done, stop or terminate the instance:")
    print(f"  aws ec2 stop-instances --instance-ids {instance_id} --region {region}")
    print(f"  aws ec2 terminate-instances --instance-ids {instance_id} --region {region}")
    print()
    print("Find all cosmx-analytics resources by tag:")
    print(f"  aws ec2 describe-instances --region {region} --filters Name=tag:Project,Values={PROJECT_TAG} Name=instance-state-name,Values=running,stopped --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`Name`].Value|[0],State.Name]' --output table")


if __name__ == "__main__":
    main()
