#!/bin/bash
# Two-pass PASS 2 (4c): concatenate the 57 per-slide supervised results and write them
# back into the clustered cohort AnnData -> cosmx_typed.h5ad + UMAP-by-cell_type. CPU.
# (The single-pass equivalent is 90_write_celltypes.sh; this one concatenates a dir.)
#
# Submit (after the 75 array completes):  sbatch pipeline/slurm/76_write_celltypes.sh
#
# Render the marker heatmap BY CELL TYPE afterwards:
#   STAGE3_DIR=stage4 GROUP_KEY=cell_type CLUSTERED_BASENAME=cosmx_typed.h5ad \
#       sbatch pipeline/slurm/50_marker_heatmap.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-write-celltypes
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/write_celltypes_%j.out
#SBATCH --error=pipeline/logs/write_celltypes_%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module load apptainer
export PATH="${HOME}/bin:${PATH}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}/pipeline"
else
    PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
fi

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"
CLUSTERED="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_write_celltypes_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/per_slide_results" "$WORK/qc_plots"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging clustered AnnData + per-slide results from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED}" \
    "$WORK/${CLUSTERED}"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_results/*" \
    "$WORK/per_slide_results/"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/write_celltypes.py" \
        --clustered-h5ad "$WORK/${CLUSTERED}" \
        --insitutype-dir "$WORK/per_slide_results" \
        --output "$WORK/cosmx_typed.h5ad" \
        --qc-plots-dir "$WORK/qc_plots"

echo "Uploading typed AnnData + QC plot to Kopah..."
s5cmd cp "$WORK/cosmx_typed.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/cosmx_typed.h5ad"
s5cmd cp "$WORK/qc_plots/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/qc_plots/" \
    || echo "WARN: no QC plots to upload"

echo "Done: PASS 2 writeback."
