#!/bin/bash
# Stage 5b: InSituCNV copy-number inference, ONE Slurm array task per tissue section.
# Each task stages its section input, runs the spatial-smoothing + infercnvpy recipe, and
# uploads a per-section CNV h5ad. Per-section keeps the spatial neighbor graph inside one
# donor tissue (never bridging the two donors co-mounted on a slide). CPU only.
#
# Submit AFTER stage 5a (needs sections.txt to size the array). Set the array range to the
# section count:
#   s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage5_insitucnv/sections.txt .
#   sbatch --array=0-$(($(wc -l < sections.txt)-1))%32 pipeline/slurm/87_insitucnv.sh
# (The %32 throttle is optional — caps concurrent array tasks so ckpt preemption is gentler.)
# Smoke test on a single section first:  sbatch --array=0 pipeline/slurm/87_insitucnv.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR    Kopah sub-dir for stage-5 (default stage5_insitucnv)
#   N_NEIGHBORS   spatial smoothing neighbors (default 20; MUST match 86b's reference build)
#   WINDOW_SIZE   infercnv window (default 100)   STEP (default 10)
#   DYNAMIC_THRESHOLD (default 1.5)

#SBATCH --job-name=cosmx-insitucnv
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/insitucnv_%A_%a.out
#SBATCH --error=pipeline/logs/insitucnv_%A_%a.err

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
# Read per-section inputs + sections.txt from here (defaults to STAGE5); a variant run
# reuses the primary run's sections while writing results to its own STAGE5_DIR.
SECTIONS="${SECTIONS_DIR:-$STAGE5}"
# Which diploid reference type list to use (default = immune-inclusive).
REF_TYPES="${REF_TYPES_BASENAME:-insitucnv_reference_types.txt}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"
: "${SLURM_ARRAY_TASK_ID:?run this as a Slurm array job (see header)}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Map this array index -> a section id from the prep manifest.
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${SECTIONS}/sections.txt" "$WORK/sections.txt"
N_SECTIONS=$(wc -l < "$WORK/sections.txt")
if (( SLURM_ARRAY_TASK_ID >= N_SECTIONS )); then
    echo "Array index ${SLURM_ARRAY_TASK_ID} >= ${N_SECTIONS} sections; nothing to do."
    exit 0
fi
SECTION=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$WORK/sections.txt")
echo "Array task ${SLURM_ARRAY_TASK_ID}/${N_SECTIONS}: section '${SECTION}'"

s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${SECTIONS}/sections/${SECTION}.h5ad" \
    "$WORK/section.h5ad"

# Use the pooled global diploid reference if 86b built one (avoids per-section reference
# contamination in tumor-bulk sections). Falls back to the per-section reference_cat path
# if no vector is present. Set REFERENCE_VECTOR=none to force the fallback.
REFVEC_ARGS=()
REFVEC_BASENAME="${REFERENCE_VECTOR:-reference_vector.csv}"
if [[ "$REFVEC_BASENAME" != "none" ]] && \
   s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/${REFVEC_BASENAME}" \
       "$WORK/reference_vector.csv" 2>/dev/null; then
    REFVEC_ARGS=(--reference-vector "$WORK/reference_vector.csv")
    echo "Using global diploid reference ${REFVEC_BASENAME}."
else
    echo "No global reference vector; using per-section reference_cat."
fi

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/run_insitucnv.py" \
        --input "$WORK/section.h5ad" \
        --reference-file "${PIPELINE_DIR}/reference/${REF_TYPES}" \
        "${REFVEC_ARGS[@]+"${REFVEC_ARGS[@]}"}" \
        --output "$WORK/${SECTION}_cnv.h5ad" \
        --n-neighbors "${N_NEIGHBORS:-20}" \
        --window-size "${WINDOW_SIZE:-100}" \
        --step "${STEP:-10}" \
        --dynamic-threshold "${DYNAMIC_THRESHOLD:-1.5}" \
        --n-jobs "${SLURM_CPUS_PER_TASK:-8}"

s5cmd cp "$WORK/${SECTION}_cnv.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/persection/${SECTION}_cnv.h5ad"

echo "Done: InSituCNV for section '${SECTION}'."
