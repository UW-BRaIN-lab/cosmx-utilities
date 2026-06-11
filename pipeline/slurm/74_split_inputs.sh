#!/bin/bash
# Two-pass PASS 2 (prep): split combined_qc.h5ad into per-slide InSituType inputs
# (one <slide_id>.h5 each) for the supervised array. Uses the already-QC'd cohort cells
# so PASS 2 is consistent with PASS 1. CPU; backed read keeps memory modest.
#
# Submit:  sbatch pipeline/slurm/74_split_inputs.sh
# Chain:   jid=$(sbatch --parsable pipeline/slurm/74_split_inputs.sh)
#          sbatch --dependency=afterok:$jid pipeline/slurm/75_insitutype_supervised.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-split-inputs
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/split_inputs_%j.out
#SBATCH --error=pipeline/logs/split_inputs_%j.err

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
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_split_inputs_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/per_slide_inputs"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging combined_qc.h5ad from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/combined_qc.h5ad" \
    "$WORK/combined_qc.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/split_insitutype_inputs.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --output-dir "$WORK/per_slide_inputs" \
        --slide-key "${SLIDE_KEY:-slide_id}"

echo "Uploading per-slide inputs to Kopah..."
s5cmd cp "$WORK/per_slide_inputs/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_inputs/"

n=$(ls "$WORK/per_slide_inputs" | wc -l)
echo "Done: PASS 2 split — ${n} per-slide inputs. Set 75's --array=1-${n}."
