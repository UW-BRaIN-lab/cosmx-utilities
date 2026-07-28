#!/bin/bash
# Stage 5h (Phase 3): summarise the Low_signal rescue — 48 -> X, by compartment, + CNV
# concordance. CPU-only; runs in the InSituCNV SIF.
#
# Submit AFTER 94 (rescue_lowsignal.csv) and 91 (cell_cnv_table.csv.gz) have uploaded:
#   sbatch --dependency=afterok:<94-jobid> pipeline/slurm/95_compare_rescue.sh
# This is a short job, so the ckpt defaults below are fine. If 94 was moved to our dedicated
# allocation, only the dependency changes; to run 95 there too, add:
#   --account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1
# (carrolllab is a collaborator's allocation — do not use it. See memory: hyak-allocations.)
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR          stage-5 sub-dir (default stage5_insitucnv; reads rescue/ + diagnostics/)
#   HIERARCHY_BASENAME  reference/ hierarchy (default insitutree_hierarchy.json)

#SBATCH --job-name=cosmx-compare-rescue
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/compare_rescue_%j.out
#SBATCH --error=pipeline/logs/compare_rescue_%j.err

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
HIERARCHY_BASENAME="${HIERARCHY_BASENAME:-insitutree_hierarchy.json}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_compare_rescue_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging rescue labels + Phase-1 cell table from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/rescue_lowsignal.csv" \
    "$WORK/rescue_lowsignal.csv"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/diagnostics/cell_cnv_table.csv.gz" \
    "$WORK/cell_cnv_table.csv.gz" || echo "WARN: no cell_cnv_table (CNV concordance skipped)"

CELL_TABLE_ARG=()
[[ -f "$WORK/cell_cnv_table.csv.gz" ]] && CELL_TABLE_ARG+=(--cell-table "$WORK/cell_cnv_table.csv.gz")

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/compare_rescue.py" \
        --rescue "$WORK/rescue_lowsignal.csv" \
        --hierarchy "${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}" \
        ${CELL_TABLE_ARG[@]+"${CELL_TABLE_ARG[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading rescue comparison to Kopah (${STAGE5}/rescue)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/"

echo "Done: stage 5h rescue comparison."
