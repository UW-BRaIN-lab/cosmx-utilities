#!/bin/bash
# Stage 0b: Extract clinical annotations (Case / Block / Region) from each
# Seurat .RDS file into a per-slide sidecar CSV on Kopah. One-time job.
#
# The .RDS files are read directly from AWS S3 source (they're never copied
# to Kopah). Only the resulting clinical CSV lands on Kopah.
#
# Submit (after editing --array=1-N to match manifest length):
#   sbatch pipeline/slurm/0b_extract_clinical.sh
#
# Required env (from pipeline/.env):
#   SOURCE_S3_BUCKET,
#   AWS_SOURCE_ACCESS_KEY_ID, AWS_SOURCE_SECRET_ACCESS_KEY (or AWS_SOURCE_PROFILE),
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_INSITUTYPE (path to R/Seurat container image; see pipeline/containers/)

#SBATCH --job-name=cosmx-extract-clinical
#SBATCH --account=glioblastoma
#SBATCH --partition=compute
#SBATCH --array=1-57
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=pipeline/logs/extract_clinical_%A_%a.out
#SBATCH --error=pipeline/logs/extract_clinical_%A_%a.err

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${MANIFEST:-${PIPELINE_DIR}/manifest.csv}"

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

# Pick this task's slide_id + flat_files_prefix from the manifest
# (header on line 1 → +1 offset; column 2 = slide_id, column 8 = flat_files_prefix).
ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
SLIDE_ID=$(awk -F, -v r="$ROW" 'NR==r {print $2}' "$MANIFEST")
FLAT_PREFIX=$(awk -F, -v r="$ROW" 'NR==r {print $8}' "$MANIFEST")

if [[ -z "${SLIDE_ID:-}" ]]; then
    echo "ERROR: no slide_id at manifest row $ROW" >&2
    exit 1
fi

# The export-batch prefix is everything before "flatFiles/" in the manifest's
# flat_files_prefix column. The Seurat .RDS file lives at that batch root.
BATCH_PREFIX="${FLAT_PREFIX%flatFiles/*}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_clinical_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# Fetch the .RDS straight from AWS S3 (boto3 reads AWS_SOURCE_* from .env).
RDS_FILE="$WORK/${SLIDE_ID}.RDS"
echo "Downloading $SLIDE_ID .RDS from AWS S3..."
uv run python "${PIPELINE_DIR}/python/fetch_seurat_from_s3.py" \
    --slide-id "$SLIDE_ID" \
    --batch-prefix "$BATCH_PREFIX" \
    --output "$RDS_FILE"

OUTPUT_CSV="$WORK/${SLIDE_ID}_clinical.csv"

# Run extraction inside the R/SeuratObject container.
# TODO: build pipeline/containers/insitutype.def and set APPTAINER_INSITUTYPE.
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/r/extract_clinical_annotations.R" \
        --input "$RDS_FILE" \
        --slide-id "$SLIDE_ID" \
        --output "$OUTPUT_CSV"

# Upload only the small CSV sidecar to Kopah; the .RDS is discarded with $WORK.
export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Uploading $OUTPUT_CSV to Kopah..."
s5cmd cp "$OUTPUT_CSV" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/clinical/${SLIDE_ID}_clinical.csv"

echo "Done: $SLIDE_ID"
