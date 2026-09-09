#!/bin/bash
# Is a de-novo letter's gene programme a cell state, or a slide-level handling artefact?
#
# 75d showed de-novo `t` is separated from every natively-called GBmap class by a heat-shock
# programme (HSPA1A/HSPA1B/HSPB1/DNAJB1/HSP90AA1/HSPH1, mean z +0.3..+1.6 in every t-> column
# against negative in every native). Heat-shock is exactly what warm ischaemia and slow fixation
# produce, so it has to be ruled out before `t` is named a tumour state.
#
# The test is WITHIN a slide: if the programme is a cell state it travels with the letter's cells
# and the letter-vs-rest gap is positive slide after slide; if it is handling, whole sections are
# elevated, the two groups move together (high correlation) and the within-slide gap collapses.
# See python/diagnose_program_by_slide.py for how the summaries read.
#
# Submit (defaults to t + the heat-shock set):
#   sbatch pipeline/slurm/75e_program_slide_confound.sh
# Another letter or another programme:
#   LETTER=b GENES=CES3,A1BG,SNRPA1,GPX1 sbatch pipeline/slurm/75e_program_slide_confound.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/anchor_typing.h5 (default stage4_anchor_pruned).
#   INPUT_DIR   Kopah sub-dir with anchor/anchor_input.h5 + anchor/anchor_cells.csv
#               (default stage4_anchor — 71 wrote the cell/slide sidecar there).
#   LETTER      de-novo letter under test (default t).
#   GENES       comma-separated programme genes (default: the heat-shock set from 75d).

#SBATCH --job-name=cosmx-program-confound
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# Only the programme genes' rows and per-cell column sums are materialised, so this never
# transposes or densifies the counts — far lighter than 75d.
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/program_confound_%j.out
#SBATCH --error=pipeline/logs/program_confound_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_anchor_pruned}"
INPUT="${INPUT_DIR:-stage4_anchor}"
LETTER="${LETTER:-t}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_program_confound_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging anchor counts + labels + the cell/slide sidecar from Kopah..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_cells.csv" "$WORK/anchor_cells.csv"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"

GENES_ARG=()
if [[ -n "${GENES:-}" ]]; then GENES_ARG=(--genes "$GENES"); fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/diagnose_program_by_slide.py" \
        --counts-h5 "$WORK/anchor_input.h5" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --cells-csv "$WORK/anchor_cells.csv" \
        --letter "$LETTER" \
        ${GENES_ARG[@]+"${GENES_ARG[@]}"} \
        --output-csv "$WORK/${LETTER}_program_by_slide.csv"

echo "Uploading the per-slide table to Kopah..."
s5cmd cp "$WORK/${LETTER}_program_by_slide.csv" \
    "${BASE}/${STAGE4}/supervised_gbmap/program_confound/${LETTER}_program_by_slide.csv"

echo "Done. The verdict is in the summary block printed above."
