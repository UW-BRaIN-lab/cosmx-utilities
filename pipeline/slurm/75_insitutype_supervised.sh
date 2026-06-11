#!/bin/bash
# Two-pass PASS 2 (4b per slide): supervised InSituType assignment of each slide's cells
# to the fixed cohort_profile.csv. Slurm ARRAY — one task per slide. Each task is small
# and fast (supervised ML, no de novo), so the 57 slides type in parallel on ckpt.
#
# Submit (set --array to the slide count printed by 74_split_inputs.sh):
#   sbatch pipeline/slurm/75_insitutype_supervised.sh
# Chain into the writeback after the whole array completes:
#   jid=$(sbatch --parsable pipeline/slurm/75_insitutype_supervised.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/76_write_celltypes.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE.

#SBATCH --job-name=cosmx-insitutype-sup
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --array=1-57
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/insitutype_sup_%A_%a.out
#SBATCH --error=pipeline/logs/insitutype_sup_%A_%a.err

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

STAGE4="${STAGE4_DIR:-stage4}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutype_sup_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

INPUT_PREFIX="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_inputs/"

# This task's slide = the Nth per-slide input (sorted), keyed off the array index.
mapfile -t SLIDES < <(s5cmd ls "$INPUT_PREFIX" | awk '{print $NF}' | grep '\.h5$' | sort)
if (( ${#SLIDES[@]} == 0 )); then
    echo "ERROR: no per-slide inputs at $INPUT_PREFIX" >&2; exit 1
fi
SLIDE_FILE="${SLIDES[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
if [[ -z "$SLIDE_FILE" ]]; then
    echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} > ${#SLIDES[@]} slides; "\
         "set --array=1-${#SLIDES[@]}" >&2
    exit 1
fi
SLIDE="${SLIDE_FILE%.h5}"
echo "Task ${SLURM_ARRAY_TASK_ID}/${#SLIDES[@]}: slide ${SLIDE}"

echo "Staging slide input + cohort profile from Kopah..."
s5cmd cp "${INPUT_PREFIX}${SLIDE_FILE}" "$WORK/${SLIDE_FILE}"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/cohort_profile.csv" \
    "$WORK/cohort_profile.csv"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutype_supervised.R" \
        --input "$WORK/${SLIDE_FILE}" \
        --profile "$WORK/cohort_profile.csv" \
        --output-h5 "$WORK/${SLIDE}.result.h5"

echo "Uploading per-slide result to Kopah..."
s5cmd cp "$WORK/${SLIDE}.result.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_results/${SLIDE}.h5"

echo "Done: PASS 2 slide ${SLIDE}."
