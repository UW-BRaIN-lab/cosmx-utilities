#!/bin/bash
# 75l: is the boring explanation the right one -- are the forced cells just shallower?
#
# A cell needs a de-novo bin when no fixed profile fits it. Low RNA depth does that on its own:
# a cell with 80 counts matches everything equally badly, and the fitted de-novo profile, free
# to move, collects it. If the letter's cells were simply thin versions of the natively-called
# ones, every other comparison in the 75 series would be measuring depth, and 75d's attenuated
# markers would have a purely technical cause.
#
# Prints each forced group against ITS OWN native counterpart on four axes -- panel counts,
# genes detected, negprobe fraction, segmented area. (51_qc_by_celltype.sh already does the
# cohort-wide version, by cell type; this is the paired one.)
#
# obs only, backed -- the expression matrix is never read.
#
# Submit:
#   sbatch pipeline/slurm/75l_denovo_origin_qc.sh
#   LETTER=b sbatch pipeline/slurm/75l_denovo_origin_qc.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR    Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   TYPED_DIR     Kopah sub-dir with the fixed-profile typing (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   LETTER        de-novo letter under test (default l).
#   DESTINATIONS  comma-separated GBmap types; default is the letter's own largest.

#SBATCH --job-name=cosmx-origin-qc
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/origin_qc_%j.out
#SBATCH --error=pipeline/logs/origin_qc_%j.err

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
TYPED="${TYPED_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
LETTER="${LETTER:-l}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_origin_qc_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging labels and the fixed-profile typing from Kopah..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" \
    "$WORK/forced_named_posteriors.csv"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_gbmap_crosstab.csv" \
    "$WORK/denovo_vs_gbmap_crosstab.csv"
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

DEST_ARG=(--crosstab "$WORK/denovo_vs_gbmap_crosstab.csv")
if [[ -n "${DESTINATIONS:-}" ]]; then DEST_ARG=(--destinations "$DESTINATIONS"); fi

OUT="$WORK/${LETTER}_qc"
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/denovo_origin_qc.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --forced-csv "$WORK/forced_named_posteriors.csv" \
        --letter "$LETTER" \
        "${DEST_ARG[@]}" \
        --output-dir "$OUT"

echo "Uploading to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/origin_qc"
for f in "$OUT"/*; do
    s5cmd cp "$f" "${DEST}/${LETTER}_$(basename "$f")"
done

echo "Done. The verdict is the paired-medians block printed above."
