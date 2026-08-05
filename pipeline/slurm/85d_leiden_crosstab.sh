#!/bin/bash
# Step 3 (Low_signal neighbour diagnosis): cross-tabulate InSituTree cell types against the
# Stage-3c Leiden clustering already carried in cosmx_typed.h5ad. The cheap, already-computed
# route to "which unsupervised neighbourhood does each Low_signal cell sit in, and who are its
# well-typed neighbours?" — no new clustering. CPU only.
#
# Runs AFTER the 90 writeback (needs cosmx_typed.h5ad = clustered + InSituTree cell_type). Number
# 85d groups it with the full-cohort InSituTree family (84/85b/85c) even though it runs post-90.
#
# Submit:
#   STAGE4_DIR=stage4_insitutree sbatch pipeline/slurm/85d_leiden_crosstab.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir holding cosmx_typed.h5ad (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.

#SBATCH --job-name=cosmx-leiden-crosstab
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/leiden_crosstab_%j.out
#SBATCH --error=pipeline/logs/leiden_crosstab_%j.err

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

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_leiden_crosstab_${SLURM_JOB_ID:-local}"
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
    python -u "${PIPELINE_DIR}/python/crosstab_leiden_typing.py" \
        --typed-h5ad "$WORK/cosmx_typed.h5ad" \
        --output-dir "$WORK/out"

echo "Uploading crosstab outputs to Kopah (${STAGE4}/leiden_crosstab)..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/leiden_crosstab/"

echo "Done: Leiden x InSituTree crosstab + Low_signal diagnosis."
