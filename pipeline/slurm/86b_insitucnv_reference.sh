#!/bin/bash
# Stage 5b-ref: build the pooled clean diploid reference for InSituCNV. Runs AFTER 5a prep
# and BEFORE the 5b array. Pools reference-type cells from CONTRALATERAL-uninvolved tissue
# (clean, genuinely diploid), smoothed the same way as the test cells, into one per-gene
# M_log1p reference profile so tumor-bulk sections stop self-referencing a contaminated
# baseline. CPU only.
#
# Submit (after 86_insitucnv_prep.sh):
#   sbatch pipeline/slurm/86b_insitucnv_reference.sh
# Then run the 5b array WITH the reference (87 picks it up automatically if present).
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR    Kopah sub-dir for stage-5 (default stage5_insitucnv)
#   N_NEIGHBORS   spatial smoothing neighbors — MUST match 87 (default 20)

#SBATCH --job-name=cosmx-insitucnv-ref
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/insitucnv_ref_%j.out
#SBATCH --error=pipeline/logs/insitucnv_ref_%j.err

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
# Read the per-section inputs from here (defaults to STAGE5); lets a variant run (e.g. an
# immune-free reference into a separate STAGE5_DIR) reuse the primary run's sections.
SECTIONS="${SECTIONS_DIR:-$STAGE5}"
# Which diploid reference type list to use (default = immune-inclusive).
REF_TYPES="${REF_TYPES_BASENAME:-insitucnv_reference_types.txt}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_ref_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/sections"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "Staging section inputs from Kopah (${SECTIONS}/sections)..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${SECTIONS}/sections/*" "$WORK/sections/"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/build_insitucnv_reference.py" \
        --sections-dir "$WORK/sections" \
        --reference-file "${PIPELINE_DIR}/reference/${REF_TYPES}" \
        --output "$WORK/reference_vector.csv" \
        --n-neighbors "${N_NEIGHBORS:-20}"

echo "Uploading reference vector to Kopah (${STAGE5})..."
s5cmd cp "$WORK/reference_vector.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/reference_vector.csv"

echo "Done: stage 5b-ref (pooled diploid reference)."
