#!/bin/bash
# Stage 3a: concatenate per-slide .h5ad, apply per-cell QC, emit the scPearsonPCA
# input (combined_qc.h5ad + cells.csv). Single big-memory CPU job.
#
# Submit:
#   sbatch pipeline/slurm/20_concat_qc.sh
# Chains into stage 3b:
#   jid=$(sbatch --parsable pipeline/slurm/20_concat_qc.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/30_pearson_pca.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF also carries anndata/scipy for this
#                  CPU-only step).

#SBATCH --job-name=cosmx-concat-qc
# glioblastoma has no dedicated CPU node → ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes its outputs).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/concat_qc_%j.out
#SBATCH --error=pipeline/logs/concat_qc_%j.err

# --mem is generous: concatenating the all-probe cohort holds it in RAM. ckpt
# pools cluster-wide nodes including large-memory ones, so 300G is schedulable
# (it may queue for a big node). Lower it if the cohort is small and you want a
# faster start.

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

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_concat_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/anndata" "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

# s5cmd reads AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / S3_ENDPOINT_URL from env.
export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-slide AnnDatas from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/anndata/*.h5ad" "$WORK/anndata/"

OUT_PREFIX="$WORK/out/combined_qc"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/concat_qc_anndata.py" \
        --anndata-dir "$WORK/anndata" \
        --output "$OUT_PREFIX" \
        --batch-col "${BATCH_COL:-Case}" \
        --cohort "${COHORT:-${PIPELINE_DIR}/cohort_wenyu.csv}" \
        --n-hvg "${N_HVG:-2000}" \
        --min-gene-counts "${MIN_GENE_COUNTS:-50}" \
        --max-area "${MAX_AREA:-30000}" \
        --max-negprobe-prop "${MAX_NEGPROBE_PROP:-0.1}"

echo "Uploading combined_qc artifacts to Kopah..."
# .h5ad: full all-probe record for stage 3c + downstream.
# .pca_input.h5: compact HVG counts + tc + per-batch gene frequency for stage 3b.
s5cmd cp "${OUT_PREFIX}.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/combined_qc.h5ad"
s5cmd cp "${OUT_PREFIX}.pca_input.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/stage3/combined_qc.pca_input.h5"

echo "Done: stage 3a"
