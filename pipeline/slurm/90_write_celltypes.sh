#!/bin/bash
# Stage 4c: write InSituType cell types back into the clustered cohort AnnData.
# CPU only; the rapids-singlecell SIF already carries anndata/scipy/h5py/scanpy.
#
# Submit:
#   sbatch pipeline/slurm/90_write_celltypes.sh
#
# Produces stage4/cosmx_typed.h5ad (cosmx_clustered + obs cell_type/insitutype_prob)
# and a UMAP-by-cell_type QC plot. To render the marker heatmap BY CELL TYPE afterwards:
#   STAGE3_DIR=stage4 GROUP_KEY=cell_type CLUSTERED_BASENAME=cosmx_typed.h5ad \
#       sbatch pipeline/slurm/50_marker_heatmap.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; CPU use, no --nv).

#SBATCH --job-name=cosmx-write-celltypes
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes cosmx_typed.h5ad).
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

# Compute nodes don't expose apptainer/s5cmd on PATH by default (see stage scripts).
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

# Read cosmx_clustered from the stage-3 dir; write the typed AnnData under stage-4 dir.
STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_write_celltypes_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/qc_plots"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging clustered AnnData + InSituType result from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/cosmx_clustered.h5ad" \
    "$WORK/cosmx_clustered.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutype_result.h5" \
    "$WORK/insitutype_result.h5"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/write_celltypes.py" \
        --clustered-h5ad "$WORK/cosmx_clustered.h5ad" \
        --insitutype-h5 "$WORK/insitutype_result.h5" \
        --output "$WORK/cosmx_typed.h5ad" \
        --qc-plots-dir "$WORK/qc_plots"

echo "Uploading typed AnnData + QC plot to Kopah..."
s5cmd cp "$WORK/cosmx_typed.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/cosmx_typed.h5ad"
s5cmd cp "$WORK/qc_plots/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/qc_plots/" \
    || echo "WARN: no QC plots to upload"

echo "Done: stage 4c"
