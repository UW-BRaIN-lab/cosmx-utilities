#!/bin/bash
# Stage 3c: GPU neighbor graph + Leiden clustering + UMAP on the Pearson-PCA
# embedding (rapids-singlecell). Single L40S is sufficient — the embedding is
# small (cells x ~50), so no dask-cuda multi-GPU cluster is needed here.
#
# Submit:
#   sbatch pipeline/slurm/40_cluster.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; run with --nv).

#SBATCH --job-name=cosmx-cluster
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/cluster_%j.out
#SBATCH --error=pipeline/logs/cluster_%j.err

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

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_cluster_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging combined_qc + embedding from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/combined_qc.h5ad" \
    "$WORK/combined_qc.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/embedding.h5" \
    "$WORK/embedding.h5"

apptainer exec --nv \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/cluster_embedding.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --embedding "$WORK/embedding.h5" \
        --output "$WORK/cosmx_clustered.h5ad" \
        --qc-plots-dir "$WORK/qc_plots" \
        --qc-color "${QC_COLOR:-Case,Region,leiden}" \
        --resolution "${RESOLUTION:-1.2}" \
        --n-neighbors "${N_NEIGHBORS:-15}"

echo "Uploading clustered AnnData + QC plots to Kopah..."
s5cmd cp "$WORK/cosmx_clustered.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/cosmx_clustered.h5ad"
# UMAP QC plots (Case / Region / leiden) for reviewing the patient-level
# correction. Non-fatal: the clustered AnnData is already safely uploaded above.
s5cmd cp "$WORK/qc_plots/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/qc_plots/" \
    || echo "WARN: no QC plots to upload"

echo "Done: stage 3c"
