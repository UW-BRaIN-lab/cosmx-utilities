#!/bin/bash
# Stage 1: Convert per-slide CosMx flat files into per-slide .h5ad.
# Runs as a Slurm job array (one task per slide).
#
# Submit (after editing --array=1-N to match manifest length):
#   sbatch pipeline/slurm/10_flatfiles_to_anndata.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_SCALESC (path to Python+scanpy container image)

#SBATCH --job-name=cosmx-flatfiles-to-anndata
#SBATCH --account=glioblastoma
#SBATCH --partition=compute
#SBATCH --array=1-57
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/flatfiles_to_anndata_%A_%a.out
#SBATCH --error=pipeline/logs/flatfiles_to_anndata_%A_%a.err

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${MANIFEST:-${PIPELINE_DIR}/manifest.csv}"

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

# Manifest columns (kept in sync with SlideManifestRow in build_manifest.py):
#   1 export_batch
#   2 slide_id
#   3 run_uuid
#   4 run_date
#   5 run_time
#   6 instrument_id
#   7 slot
#   8 flat_files_prefix    (path within SOURCE bucket; we use it for Kopah too)
#   9 decoded_prefix
ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
SLIDE_ID=$(awk -F, -v r="$ROW" 'NR==r {print $2}' "$MANIFEST")
FLAT_PREFIX=$(awk -F, -v r="$ROW" 'NR==r {print $8}' "$MANIFEST")

if [[ -z "${SLIDE_ID:-}" ]]; then
    echo "ERROR: no slide_id at manifest row $ROW" >&2
    exit 1
fi

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_flat_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK/flat"
trap 'rm -rf "$WORK"' EXIT

# s5cmd reads AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / S3_ENDPOINT_URL from env.
export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging flat files for $SLIDE_ID from Kopah..."
S3_FLAT="s3://${KOPAH_BUCKET}/${FLAT_PREFIX}"
for f in "${SLIDE_ID}_exprMat_file.csv.gz" \
         "${SLIDE_ID}_metadata_file.csv.gz" \
         "${SLIDE_ID}_fov_positions_file.csv.gz"; do
    s5cmd cp "${S3_FLAT}${f}" "$WORK/flat/"
done

# Pull the clinical sidecar produced by stage 0b (optional — the Python script
# handles the missing case with a WARN).
s5cmd cp \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/clinical/${SLIDE_ID}_clinical.csv" \
    "$WORK/${SLIDE_ID}_clinical.csv" \
    || echo "WARN: no clinical sidecar yet for $SLIDE_ID"

OUTPUT_H5AD="$WORK/${SLIDE_ID}.h5ad"

# TODO: build pipeline/containers/scalesc.def and set APPTAINER_SCALESC.
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_SCALESC" \
    python "${PIPELINE_DIR}/python/flatfiles_to_anndata.py" \
        --flatfiles-dir "$WORK/flat" \
        --slide-id "$SLIDE_ID" \
        --clinical "$WORK/${SLIDE_ID}_clinical.csv" \
        --manifest "$MANIFEST" \
        --output "$OUTPUT_H5AD"

echo "Uploading $OUTPUT_H5AD to Kopah..."
s5cmd cp "$OUTPUT_H5AD" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/anndata/${SLIDE_ID}.h5ad"

echo "Done: $SLIDE_ID"
