#!/bin/bash
# Stage 5f (SHARDED): per-cell TOP-K posteriors against the fixed InSituTree profile set, one
# task per slide. Full-cohort variant of 92_flat_posteriors.sh. Flat InSituType::insitutypeML
# scoring (NO clustering, NO reference update) recovers each cell's top-two profile assignments
# — the object for the hybrid-vs-continuum readout (adjacent malignant states => continuum;
# malignant + neuronal cross-compartment => candidate hybrid). Reuses the SAME per-slide inputs
# 84 produced, so it shards exactly like the InSituTree array (85b).
#
# Submit (set --array to the count 84 printed):
#   sbatch --array=1-<N> pipeline/slurm/92b_flat_posteriors_sharded.sh
# Chain the concat after the array completes:
#   jid=$(sbatch --parsable --array=1-<N> pipeline/slurm/92b_flat_posteriors_sharded.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/92c_concat_flat_posteriors.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   INPUT_DIR          Kopah sub-dir with per_slide_inputs/ (default stage4, shared with 85b).
#   PROFILES_BASENAME  fixed profile matrix in reference/ (default insitutree_profiles.csv).
#   STAGE5_DIR         Kopah sub-dir for stage-5 outputs (default stage5_insitucnv).
#   TOP_K              top assignments per cell (default 3).

#SBATCH --job-name=cosmx-flat-post-shard
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/flat_post_shard_%A_%a.out
#SBATCH --error=pipeline/logs/flat_post_shard_%A_%a.err

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

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_flat_post_shard_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

INPUT_PREFIX="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/per_slide_inputs/"

mapfile -t SLIDES < <(s5cmd ls "$INPUT_PREFIX" | awk '{print $NF}' | grep '\.h5$' | sort)
if (( ${#SLIDES[@]} == 0 )); then
    echo "ERROR: no per-slide inputs at $INPUT_PREFIX (run 84 first)" >&2; exit 1
fi
SLIDE_FILE="${SLIDES[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
if [[ -z "$SLIDE_FILE" ]]; then
    echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} > ${#SLIDES[@]} slides; "\
         "set --array=1-${#SLIDES[@]}" >&2
    exit 1
fi
SLIDE="${SLIDE_FILE%.h5}"
echo "Task ${SLURM_ARRAY_TASK_ID}/${#SLIDES[@]}: slide ${SLIDE}"

echo "Staging slide input + profile matrix from Kopah..."
s5cmd cp "${INPUT_PREFIX}${SLIDE_FILE}" "$WORK/${SLIDE_FILE}"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${PROFILES_BASENAME}" \
    "$WORK/profiles.csv"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/flat_posteriors.R" \
        --input "$WORK/${SLIDE_FILE}" \
        --profiles "$WORK/profiles.csv" \
        --top-k "$TOP_K" \
        --output-csv "$WORK/${SLIDE}.flat_posteriors.csv"

echo "Uploading per-slide flat posteriors to Kopah (${STAGE5}/posteriors/per_slide)..."
s5cmd cp "$WORK/${SLIDE}.flat_posteriors.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/posteriors/per_slide/${SLIDE}.csv"

echo "Done: sharded flat posteriors slide ${SLIDE}."
