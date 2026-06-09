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
# glioblastoma has no dedicated CPU node → ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes embedding.h5).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/pearson_pca_%j.out
#SBATCH --error=pipeline/logs/pearson_pca_%j.err

# scPearsonPCA's genes x genes SVD is tiny, but the block-wise embedding
# projection runs under parallel::mclapply, which forks one worker per core.
# Fork is copy-on-write, but R's GC touches object headers, so each fork ends up
# duplicating the multi-GB sparse counts — at 2.3M cells x 2000 HVGs, 16 forks
# OOM-killed a 128G node. Keep --cpus-per-task modest (it feeds --ncores) and
# give generous --mem; drop --cpus-per-task to 1 if a larger cohort still OOMs.

set -euo pipefail

# Compute nodes don't expose apptainer/s5cmd on PATH by default. Initialise Lmod
# (a batch shell isn't a login shell, so the module function may be undefined),
# load apptainer, and add ~/bin where the s5cmd binary lives — compute nodes do
# not share the login node's /usr/bin/s5cmd, but home (/mmfs1/home) is shared.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module load apptainer
export PATH="${HOME}/bin:${PATH}"

# Slurm copies the batch script to a spool dir, so "$0" no longer points into the
# repo. Derive the repo from SLURM_SUBMIT_DIR (where sbatch was invoked — submit
# from the repo root); fall back to "$0" for direct, non-Slurm execution.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}/pipeline"
else
    PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
fi

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

# Kopah sub-prefix for Stage 3 outputs. Override (e.g. STAGE3_DIR=stage3_q100)
# to write a parallel run without clobbering the default `stage3/` results.
STAGE3="${STAGE3_DIR:-stage3}"

: "${APPTAINER_SCPEARSON:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_pca_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging the PCA input from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/combined_qc.pca_input.h5" \
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
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/embedding.h5"

echo "Done: stage 3b"
