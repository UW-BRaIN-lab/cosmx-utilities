#!/bin/bash
# Stage 5c: aggregate the per-section InSituCNV results and answer the question — do the
# Low_signal cells carry the GBM malignant CNV signature (chr7 gain / chr10 loss)? Pulls
# all persection/*_cnv.h5ad, forms per-group mean CNV profiles, cosine similarity, chr-arm
# summaries, a control-calibrated cnv_score threshold, and the plots. CPU only.
#
# Submit AFTER the stage 5b array has finished:
#   sbatch pipeline/slurm/88_insitucnv_compare.sh
# Outputs (tables + PNGs + SUMMARY.txt) land at stage5_insitucnv/compare/ on Kopah;
# download SUMMARY.txt + the PNGs to read the verdict.
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR   Kopah sub-dir for stage-5 (default stage5_insitucnv)
#   MIN_CELLS    drop groups smaller than this from the comparison (default 200)

#SBATCH --job-name=cosmx-insitucnv-compare
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/insitucnv_compare_%j.out
#SBATCH --error=pipeline/logs/insitucnv_compare_%j.err

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

STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_compare_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/persection" "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-section CNV results from Kopah (${STAGE5}/persection)..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/persection/*" "$WORK/persection/"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/compare_insitucnv_groups.py" \
        --cnv-dir "$WORK/persection" \
        --reference-file "${PIPELINE_DIR}/reference/insitucnv_reference_types.txt" \
        --output-dir "$WORK/out" \
        --min-cells "${MIN_CELLS:-200}"

echo "Uploading comparison tables + plots to Kopah (${STAGE5}/compare)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/compare/"

echo "Done: stage 5c InSituCNV comparison."
