#!/bin/bash
# Stage 5d (diagnostic): edge vs bulk — is the stronger malignant CNV signature in
# infiltrating-edge Low_signal (vs tumor-bulk) biological, or a depth/quality artifact?
# Rebuilds the per-cell signature from the per-section CNV outputs, joins per-cell RNA depth
# from cosmx_typed.h5ad, and does a depth-stratified bulk-vs-edge comparison. CPU only.
#
# Submit AFTER the 5b array + (optionally) 5c compare:
#   sbatch pipeline/slurm/89_insitucnv_edge_vs_bulk.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR   Kopah sub-dir for stage-5 (default stage5_insitucnv)
#   STAGE4_DIR   Kopah sub-dir with the typed AnnData (default stage4_insitutree)
#   TYPED_BASENAME (default cosmx_typed.h5ad)

#SBATCH --job-name=cosmx-insitucnv-edgebulk
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/insitucnv_edgebulk_%j.out
#SBATCH --error=pipeline/logs/insitucnv_edgebulk_%j.err

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
STAGE4="${STAGE4_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_edgebulk_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/persection" "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-section CNV results + typed AnnData from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/persection/*" "$WORK/persection/"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/insitucnv_edge_vs_bulk.py" \
        --cnv-dir "$WORK/persection" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --reference-file "${PIPELINE_DIR}/reference/insitucnv_reference_types.txt" \
        --output-dir "$WORK/out"

echo "Uploading edge-vs-bulk diagnostic to Kopah (${STAGE5}/edge_vs_bulk)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/edge_vs_bulk/"

echo "Done: stage 5d edge-vs-bulk diagnostic."
