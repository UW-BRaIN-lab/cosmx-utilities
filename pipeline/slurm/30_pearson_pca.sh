#!/bin/bash
# Stage 3b: batch-corrected quasipoisson Pearson-residual PCA (scPearsonPCA, R).
# Single big-memory CPU job — NOT a GPU job.
#
# Submit:
#   sbatch pipeline/slurm/30_pearson_pca.sh
# Chains into stage 3c:
#   jid=$(sbatch --parsable pipeline/slurm/30_pearson_pca.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/40_cluster.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_SCPEARSON (path to the scpearsonpca.sif R container).

#SBATCH --job-name=cosmx-pearson-pca
#SBATCH --account=glioblastoma
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=350G
#SBATCH --time=08:00:00
#SBATCH --output=pipeline/logs/pearson_pca_%j.out
#SBATCH --error=pipeline/logs/pearson_pca_%j.err

# --mem targets a high-memory node: scPearsonPCA holds the genes x cells counts and
# one working copy as dgCMatrix (~tens of GB at 8M cells). The genes x genes SVD
# itself is small. --cpus-per-task feeds --ncores for the block-wise embedding
# projection. Adjust --mem/--constraint to your allocation's large-memory nodes.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

: "${APPTAINER_SCPEARSON:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_pca_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging the PCA input from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/combined_qc.pca_input.h5" \
    "$WORK/pca_input.h5"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_SCPEARSON" \
    Rscript "${PIPELINE_DIR}/R/pearson_pca.R" \
        --input "$WORK/pca_input.h5" \
        --output "$WORK/embedding.h5" \
        --batch-variable "${BATCH_VARIABLE:-Case}" \
        --npcs "${NPCS:-50}" \
        --scale-max "${SCALE_MAX:-10}" \
        --ncores "${SLURM_CPUS_PER_TASK:-1}"

echo "Uploading embedding to Kopah..."
s5cmd cp "$WORK/embedding.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/embedding.h5"

echo "Done: stage 3b"
