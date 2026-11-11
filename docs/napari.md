# Napari visualization on AWS

Stitch CosMx slides into whole-slide mosaics on AWS Fargate, then view them in Napari on an EC2 GPU instance or your own machine. Back to the [main README](../README.md).

## Install

```bash
# Install CLI tools (no GUI required)
uv sync
uv run stitch-images --help
uv run read-targets --help

# Install with Napari GUI for interactive visualization
uv sync --extra gui
uv run napari /path/to/processed-output
```

## Process slides on Fargate

Discover all slides under an S3 experiment directory and launch one Fargate task per slide:

```bash
# Preview what would be launched
uv run python scripts/process-slides.py s3://bucket/project/study/ --whatif

# Launch Fargate tasks (one per slide, parallel)
uv run python scripts/process-slides.py s3://bucket/project/study/

# Skip already processed slides
uv run python scripts/process-slides.py s3://bucket/project/study/ --skip
```

Behaviour worth knowing:

- **`--skip`** skips a slide only when its `_SUCCESS` marker records the same code version as the image that would run now, so slides re-run automatically after a code change.
- **Uploads replace** the previous output rather than merging into it, so stale zarr chunks never survive. Pass `--no-replace-output` when running locally without delete permission.
- **The S3 URI** may be a whole study or a single AtoMx run. Use one run when two studies cover the same slides with different flags (e.g. a 3D resegmentation).
- **Task-definition drift:** launching refuses if the registered Fargate task definition no longer matches `fargate/fargate-task-process-slide.json`. Re-register with `./fargate/register-task-defs.sh`, or pass `--allow-taskdef-drift`.
- **Task size:** each slide launches at the task definition's default (4 vCPU / 30 GB) unless its transcript file (`flatFiles/<slide>/<slide>_tx_file.csv.gz`) is 2.5 GB or larger, which gets 16 vCPU / 120 GB; at the default, transcript decoding runs out of memory on those slides. `--cpu`/`--memory` always win, `--no-auto-size` turns this off, and an unreadable transcript file keeps the default.

Each Fargate task runs `process-slide.py`, which:
1. Queries segmentation manifests in S3 via DuckDB to find the correct segmentation version
2. Downloads only the needed CellLabels, morphology images, and AnalysisResults
3. Stitches FOV images into whole-slide zarr mosaics (`stitch-images`)
4. Decodes RNA transcript locations (`read-targets`)
5. Generates `_metadata.csv` with metadata annotations such as cell type and assigns deterministic colors
6. Uploads processed results back to S3

### Carry annotations into Napari

`--column OUTNAME=SOURCE_HEADER` (repeatable) selects which flatFiles columns
become colorable annotations in `_metadata.csv`. `SOURCE_HEADER` may list
fallbacks separated by `|`, taking the first header the file actually has — use
this when a study renamed a column between AtoMx runs:

```bash
uv run python scripts/process-slides.py s3://bucket/study/atomx-run/ \
    --column 'Cell Type=RNA_RNA_Cell.Typing.InSituType.2_1_clusters' \
    --column 'Case Specific=Case_specific_SORL1|Case_specific'
```

- **Re-exported typing:** when AtoMx nests a re-run under the original (`<run>/flatFiles/<rerun>/flatFiles/<slide>/`), the export that actually has the requested columns is used.
- **Per-FOV sheets:** when AtoMx never captured the annotations, put one `<slide>_annotations.csv` per slide (one row per FOV; FOV column `FOVs`, `FOV` or `fov`) under a prefix and pass `--annotations-prefix`. A slide column in the sheet is checked against the slide being processed, so the wrong file fails loudly.

  ```bash
  uv run python scripts/process-slides.py s3://bucket/study/atomx-run/ \
      --annotations-prefix s3://bucket/study/annotations \
      --column 'UWA=UWA' --column 'Case Broad=Case_broad'
  ```

- **From another study:** `--fill-from <experiment prefix>` copies FOV-constant columns from another study of the same physical slides, filling blanks only. Per-cell columns (such as cell typing) are refused, because they belong to the segmentation that produced them.

### 3D segmentation

A 3D-resegmented slide has per-z CellLabels (`CellLabels_F#####_Z###.tif`).
Pass `--input-ndim 3 --output-ndim 3` for a z-navigable labels layer in Napari.

3D labels usually sit over a **2D** morphology acquisition (`Morphology2D/`). That is detected automatically and the morphology plane is stored once instead of per z-plane, which keeps a 3D mosaic within a few GiB of the 2D one (copying it to 8 planes would exceed Fargate's 200 GiB storage limit). Existing 2D mosaics need no re-stitching.

> **Use only post-fix 3D exports.** Bruker asked customers to pause AtoMx 3D segmentation until a fixed module was deployed; 3D exports made before the fix should not be used.

### Docker images

The multi-stage Dockerfile produces two image variants:

```bash
# Headless: CLI tools + AWS CLI for Fargate batch processing
docker build --target headless -t cosmx-utilities:headless .

# GUI: Full Napari with Qt for interactive visualization
docker build --target gui -t cosmx-utilities:gui .
```

## View slides on your own machine

If you've received a processed CosMx dataset and want to run Napari locally instead of on the AWS-provisioned EC2 instance, the stitched output is self-contained — there are no AWS dependencies at view time.

**1. Copy the processed data from S3.**

Each slide directory contains `images/` (zarr pyramid), `_metadata.csv`, and `targets.hdf5`. Budget ~50–100 GB per slide; a full study can run several hundred GB. SSD storage is strongly recommended — spinning disks make tile loading painful.

```bash
aws s3 sync s3://<bucket>/napari-stitched/<study>/<experiment>/ /path/to/local/stitched/
# or, if rclone is configured:
rclone sync :s3:<bucket>/napari-stitched/<study>/<experiment>/ /path/to/local/stitched/ \
    --transfers 32 --checkers 16 --progress
```

You'll need AWS credentials with read access to the data bucket — ask the data owner.

**2. Install system dependencies (Linux headless server).**

PyQt6 needs Qt platform libraries that aren't installed by default on bare servers. On Ubuntu/Debian:

```bash
sudo apt install -y libxcb-cursor0 libxkbcommon-x11-0 libegl1 libgl1 libdbus-1-3 \
    libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 \
    libxcb-shape0 libxcb-xinerama0 libxcb-xkb1
```

Napari needs a graphical session. If your server is headless, use VNC, X-forwarding, or NICE DCV. A discrete GPU isn't required — software OpenGL (Mesa) works, just slower on large mosaics.

**3. Clone this repo and install with the GUI extra.**

The `napari-cosmx-fork` is a uv workspace member of this repo, and the top-level wrapper adds dependencies (dask, polars, matplotlib) that the Napari load path uses. Don't install the fork standalone.

```bash
# Install uv if needed: https://docs.astral.sh/uv/getting-started/installation/
git clone https://github.com/UW-BRaIN-lab/cosmx-utilities.git
cd cosmx-utilities
uv sync --extra gui
```

`uv` will provision a Python 3.10 interpreter automatically — the fork pins `>=3.8,<3.11`, so a system `python3.12` won't satisfy it directly.

**4. Launch Napari.**

```bash
uv run napari /path/to/local/stitched
```

**Resource notes.** Loading a stitched slide pulls the zarr pyramid and `targets.hdf5` into memory lazily, but interactive panning at high zoom is RAM-hungry. Plan for at least 64 GB free RAM; 128 GB is comfortable for the typical CosMx slide.

## View slides on an EC2 instance

`ec2/start_ec2.py --napari` starts a GPU instance (`g4dn`) with the stitched slides copied onto its local NVMe drive. Cloning this repo is not enough to use it — you also need AWS access and a local `ec2/.env`.

**One-time setup**

1. **AWS credentials** in the data bucket's account, with the permissions in [`ec2/ec2-operator-policy.json`](../ec2/ec2-operator-policy.json) plus read access to the bucket (an administrator attaches these; commands are in the header of [`ec2/.env.example`](../ec2/.env.example)). The scripts use your AWS CLI credentials (`aws configure sso` / `aws sso login`, or `aws configure`); check with `aws sts get-caller-identity`.
2. **An `ec2/.env` file.** It is gitignored. Copy `ec2/.env.example` to `ec2/.env` and fill in the region, subnet, security group, key pair, instance profile, base Ubuntu AMI, and a `DCV_PASSWORD`. The network and instance profile already exist in the lab's AWS account; ask a lab member for the values. Never commit them.
3. **Tools:** `uv` and the AWS CLI. For the remote desktop you also need the UW VPN and the NICE DCV client, since instances only get private IPs. SSM sessions work without the VPN.

**Launch**

```bash
uv run python ec2/start_ec2.py --name <your-name>-napari --napari --raw \
    --s3 s3://<bucket>/napari-stitched/<study>/<experiment>/
```

The script sums the S3 prefix and picks the smallest instance that fits: `g4dn.4xlarge` for up to 200 GB, `g4dn.8xlarge` for up to 850 GB. Above that it stops with an error, so point `--s3` at a smaller prefix. Pointing `--s3` at a parent prefix keeps the slides nested as `/mnt/local/stitched/<slide>/`; pointing it at one slide puts its files flat in `/mnt/local/stitched/`.

Setup takes about 5 minutes and the copy about 15–25 minutes for several hundred GB. To check progress without logging in, look for the `S3 sync complete` marker in the console output:

```bash
aws ec2 get-console-output --region <region> --instance-id <id> --latest --output text | grep "S3 sync complete"
```

Console output first appears 4–6 minutes after launch. Don't use `ssm describe-instance-information` as a readiness check: many IAM users are denied it, so it returns nothing either way.

**Connect and view**

- DCV: connect to `<private-ip>:8443` on the VPN as user `ubuntu`, with your `DCV_PASSWORD`.
- SSM: `aws ssm start-session --target <id> --region <region>`.

The viewer renders on the CPU (software OpenGL), not the GPU, so open **2–3 slides at a time**, each from its own slide folder. Pointing Napari at the parent `stitched` folder does not work. From the DCV desktop:

```bash
/opt/cosmx-utilities/ec2/open-napari-slides.sh /mnt/local/stitched/<slide>
```

**Cost and cleanup.** Instances bill by the hour until terminated. `/mnt/local` is wiped when an instance stops or terminates, which is safe because the data lives on S3. Terminate when you are done:

```bash
aws ec2 terminate-instances --instance-ids <id> --region <region>
```

## Infrastructure setup

Fargate task definitions, IAM roles, and networking configuration are documented in [`fargate/FARGATE-SETUP.md`](../fargate/FARGATE-SETUP.md). Infrastructure IDs are stored in `fargate/.env` (gitignored) — copy `fargate/.env.example` to get started.
