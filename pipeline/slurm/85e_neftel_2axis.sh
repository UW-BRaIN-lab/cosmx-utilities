#!/bin/bash
# Step 3 (continuous complement): Neftel 2-axis (butterfly) scoring of the cohort. Scores the four
# Neftel-like modules per cell and places every cell on the malignant-state plane, so gradient /
# between-state tumor cells show as continuous position — the complement to the discrete top-two
# hybrid readout (92b flat posteriors). CPU (scanpy score_genes); big memory for the 7.5M-cell load.
#
# Runs AFTER the 90 writeback (needs cosmx_typed.h5ad). Reads gene_signatures.csv +
# insitutree_hierarchy.json from the repo checkout (bind-mounted).
#
# Submit:
#   STAGE4_DIR=stage4_insitutree sbatch pipeline/slurm/85e_neftel_2axis.sh
# If ckpt preempts it (no mid-run checkpoint), rerun on the dedicated slice:
#   sbatch --account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1 \
#       pipeline/slurm/85e_neftel_2axis.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR      Kopah sub-dir holding cosmx_typed.h5ad (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.

#SBATCH --job-name=cosmx-neftel
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/neftel_%j.out
#SBATCH --error=pipeline/logs/neftel_%j.err

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
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_neftel_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${TYPED_BASENAME} from Kopah (${STAGE4})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" \
    "$WORK/cosmx_typed.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/neftel_2axis.py" \
        --typed-h5ad "$WORK/cosmx_typed.h5ad" \
        --signatures "${PIPELINE_DIR}/reference/gene_signatures.csv" \
        --hierarchy "${PIPELINE_DIR}/reference/insitutree_hierarchy.json" \
        --output-dir "$WORK/out"

echo "Uploading Neftel outputs to Kopah (${STAGE4}/neftel)..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/neftel/"

echo "Done: Neftel 2-axis scoring."
