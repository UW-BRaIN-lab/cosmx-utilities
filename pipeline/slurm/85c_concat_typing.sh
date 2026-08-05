#!/bin/bash
# Stage 4 (alt, SHARDED): concatenate the per-slide InSituTree results (85b) into one cohort
# insitutree_result.h5, then hand off to the existing 90 writeback. CPU only.
#
# Submit after the 85b array completes:
#   sbatch --dependency=afterok:<85b_jid> pipeline/slurm/85c_concat_typing.sh
# Then write the cohort labels back onto the clustered AnnData with the UNCHANGED 90:
#   jid=$(sbatch --parsable --dependency=afterok:<85b_jid> pipeline/slurm/85c_concat_typing.sh)
#   STAGE4_DIR=stage4_insitutree RESULT_BASENAME=insitutree_result.h5 \
#       sbatch --dependency=afterok:$jid pipeline/slurm/90_write_celltypes.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_RSC (CPU use; carries h5py/numpy).

#SBATCH --job-name=cosmx-insitutree-concat
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/insitutree_concat_%j.out
#SBATCH --error=pipeline/logs/insitutree_concat_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_insitutree}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutree_concat_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/per_slide_results"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-slide results from Kopah (${STAGE4}/per_slide_results)..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_results/*" \
    "$WORK/per_slide_results/"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/concat_typing_results.py" \
        --results-dir "$WORK/per_slide_results" \
        --output "$WORK/insitutree_result.h5" \
        --summary-csv "$WORK/insitutree_label_counts.csv"

echo "Uploading cohort result + label counts to Kopah (${STAGE4})..."
s5cmd cp "$WORK/insitutree_result.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_result.h5"
s5cmd cp "$WORK/insitutree_label_counts.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_label_counts.csv" \
    || echo "WARN: failed to upload insitutree_label_counts.csv"

echo "Done: cohort InSituTree result assembled."
echo "Next (writeback): STAGE4_DIR=${STAGE4} RESULT_BASENAME=insitutree_result.h5 \\"
echo "    sbatch pipeline/slurm/90_write_celltypes.sh"
