#!/bin/bash
# Stage 4a: emit the compact InSituType input from the stage-3a cohort AnnData.
# CPU only; the rapids-singlecell SIF already carries anndata/scipy/h5py.
#
# Submit:
#   sbatch pipeline/slurm/70_prep_insitutype.sh
# Chains into stage 4b:
#   jid=$(sbatch --parsable pipeline/slurm/70_prep_insitutype.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/80_insitutype.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; CPU use, no --nv).

#SBATCH --job-name=cosmx-prep-insitutype
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes insitutype_input.h5).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=03:00:00
#SBATCH --output=pipeline/logs/prep_insitutype_%j.out
#SBATCH --error=pipeline/logs/prep_insitutype_%j.err

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

# Read combined_qc from the stage-3 dir; write the stage-4 input under stage-4 dir.
STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_prep_insitutype_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging combined_qc from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/combined_qc.h5ad" \
    "$WORK/combined_qc.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/prep_insitutype_inputs.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --output "$WORK/insitutype_input.h5"

echo "Uploading InSituType input to Kopah..."
s5cmd cp "$WORK/insitutype_input.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutype_input.h5"

echo "Done: stage 4a"
