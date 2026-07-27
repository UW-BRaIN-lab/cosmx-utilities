#!/bin/bash
# Stage 5f (Phase 2): per-cell TOP-K posteriors against the fixed InSituTree profile set.
# Flat InSituType::insitutypeML scoring (no clustering, no reference update) to recover each
# cell's top-two profile assignments for the Stage-3 hybrid-vs-continuum readout. Reuses the
# same insitutype_input.h5 Stage 4 consumed; runs in the InSituType SIF on a big-memory CPU node.
#
#   sbatch pipeline/slurm/92_flat_posteriors.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   INPUT_DIR          Kopah sub-dir with insitutype_input.h5 (default stage4, shared).
#   PROFILES_BASENAME  fixed profile matrix in reference/ (default insitutree_profiles.csv).
#   STAGE5_DIR         Kopah sub-dir for stage-5 outputs (default stage5_insitucnv).
#   TOP_K              top assignments per cell (default 3).

#SBATCH --job-name=cosmx-flat-posteriors
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=pipeline/logs/flat_posteriors_%j.out
#SBATCH --error=pipeline/logs/flat_posteriors_%j.err

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

INPUT_DIR="${INPUT_DIR:-stage4}"
PROFILES_BASENAME="${PROFILES_BASENAME:-insitutree_profiles.csv}"
STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
TOP_K="${TOP_K:-3}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_flat_posteriors_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging insitutype_input.h5 (${INPUT_DIR}) + profile matrix from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/insitutype_input.h5" \
    "$WORK/insitutype_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${PROFILES_BASENAME}" \
    "$WORK/profiles.csv"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/flat_posteriors.R" \
        --input "$WORK/insitutype_input.h5" \
        --profiles "$WORK/profiles.csv" \
        --top-k "$TOP_K" \
        --output-csv "$WORK/flat_posteriors.csv"

echo "Uploading flat posteriors to Kopah (${STAGE5}/posteriors)..."
s5cmd cp "$WORK/flat_posteriors.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/posteriors/flat_posteriors.csv"

echo "Done: stage 5f flat posteriors."
