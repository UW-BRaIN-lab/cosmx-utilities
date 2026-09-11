#!/bin/bash
# PI follow-up: WHY are a de-novo letter's cells a bad fit to the GBmap type they were forced into?
#
# The marker heatmaps (75d) cannot answer this — their matrix is z-scored per gene across the
# displayed groups, which removes amplitude by construction, so a proportionally dimmer copy looks
# identical to a compositionally different population. This job works from raw counts and splits
# the misfit into three measurable causes per letter->destination pair:
#   DEPTH        median counts and genes per cell, letter->D over D [native]
#   FLATNESS     share of a cell's counts in its own top 20 genes (what "no distinctive profile"
#                looks like numerically) — the Low_signal signature
#   COMPOSITION  correlation of depth-corrected mean profiles: same shape, or a different thing
# plus the genes the letter's cells are MISSING relative to their forced type, and what they carry
# instead. See the header of python/denovo_misfit_diagnostics.py.
#
# Submit (defaults to the letters already dissected in 75d):
#   sbatch pipeline/slurm/75g_denovo_misfit.sh
# Other letters, destinations taken from the cross-tab:
#   LETTERS=q,o,v sbatch pipeline/slurm/75g_denovo_misfit.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   INPUT_DIR    Kopah sub-dir holding anchor/anchor_input.h5 (default stage4_anchor).
#   LETTERS      comma-separated letters (default b,t,l,d,e,j).
#   TOP_DEST     destinations per letter, from the cross-tab (default 3).

#SBATCH --job-name=cosmx-denovo-misfit
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# Subsets each group on the CSC before transposing, but b's groups alone are ~800k cells and the
# per-cell top-20 pass walks the CSR rows, so give it room.
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/denovo_misfit_%j.out
#SBATCH --error=pipeline/logs/denovo_misfit_%j.err

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
LETTERS="${LETTERS:-b,t,l,d,e,j}"
TOP_DEST="${TOP_DEST:-3}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

echo "Letters to diagnose: ${LETTERS} (top ${TOP_DEST} destinations each, from the cross-tab)"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_misfit_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging anchor counts + labels from Kopah..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" "$WORK/forced.csv"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_gbmap_crosstab.csv" "$WORK/crosstab.csv"

IFS=',' read -ra LETTER_LIST <<< "$LETTERS"
for letter in "${LETTER_LIST[@]}"; do
    letter="${letter// /}"
    echo
    echo "=== ${letter} ==="
    apptainer exec \
        --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
        --bind "${WORK}:${WORK}" \
        "$APPTAINER_RSC" \
        python "${PIPELINE_DIR}/python/denovo_misfit_diagnostics.py" \
            --counts-h5 "$WORK/anchor_input.h5" \
            --typing-h5 "$WORK/anchor_typing.h5" \
            --forced-csv "$WORK/forced.csv" \
            --letter "$letter" \
            --crosstab "$WORK/crosstab.csv" \
            --top-destinations "$TOP_DEST" \
            --output-dir "$WORK/out/${letter}"

    s5cmd cp "$WORK/out/${letter}/*" \
        "${BASE}/${STAGE4}/supervised_gbmap/misfit/${letter}/"
done

echo
echo "Done. The per-letter verdict tables are printed above; CSVs are on Kopah under"
echo "  ${STAGE4}/supervised_gbmap/misfit/<letter>/"
