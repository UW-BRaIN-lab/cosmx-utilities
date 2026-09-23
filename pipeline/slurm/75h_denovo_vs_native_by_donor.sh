#!/bin/bash
# Do a de-novo letter's cells and the natively-called cells of its destination come from the
# SAME donors, or different ones?
#
# 75d compares those two groups on markers, and its header records why that can only go so far:
# a cell landed in the letter precisely BECAUSE it fit the letter's profile better, so "they
# differ" is definitional. This job asks something the assignment rule does not already answer.
#
# The reason donor is the discriminating axis: the 54 named profiles are FIXED (rescaled
# reference columns that cannot move) while the 27 de-novo profiles are FITTED by EM from this
# cohort's own cells. A fitted profile goes to where the cells are; a fixed one cannot. So a
# de-novo cluster wins wherever this cohort is systematically offset from the reference —
# platform, batch, patient, fixation, segmentation — and the letter would then be absorbing that
# offset rather than marking a cell state. That reading predicts A and B SEPARATE BY DONOR; a
# genuine cell state predicts every donor contributes both, in similar proportion.
#
# The yardstick is the cohort's own: every pair of NAMED types inside the destination's
# compartment. Both members are fixed profiles, so their donor-to-donor ratio carries real
# composition variation but not the fitted-vs-fixed asymmetry.
#
# LABELS ONLY — no counts matrix — so this is seconds of compute, not the 9.6GB read 75d needs.
#
# Submit (l against its own three largest destinations):
#   LETTER=l sbatch pipeline/slurm/75h_denovo_vs_native_by_donor.sh
# Name the destinations by hand:
#   LETTER=l DESTINATIONS=Endo_arterial,Endo_capilar,Pericyte \
#     sbatch pipeline/slurm/75h_denovo_vs_native_by_donor.sh
# The same test as a BATCH test — Export_source is which of the five flat-file exports a slide
# came from:
#   LETTER=l GROUP_BY=export sbatch pipeline/slurm/75h_denovo_vs_native_by_donor.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   LETTER       de-novo letter under test (default l).
#   DESTINATIONS comma-separated GBmap types; default is the letter's own largest.
#   TOP_DEST     how many destinations when DESTINATIONS is unset (default 3).
#   GROUP_BY     case (default), slide, region, block, or export.
#   MIN_CELLS    minimum A+B cells for a group to count (default 25).

#SBATCH --job-name=cosmx-denovo-by-donor
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=pipeline/logs/denovo_by_donor_%j.out
#SBATCH --error=pipeline/logs/denovo_by_donor_%j.err

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
LETTER="${LETTER:-l}"
GROUP="${GROUP_BY:-case}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_by_donor_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging labels from Kopah (no counts matrix needed)..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
# 75c's corrected forced call. NOT 75's re-scored file, which disagrees with the fit on 45% of
# cells — see the 75c self-consistency check.
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" \
    "$WORK/forced_named_posteriors.csv"

DEST_ARG=()
if [[ -n "${DESTINATIONS:-}" ]]; then DEST_ARG=(--destinations "$DESTINATIONS"); fi

OUT="$WORK/${LETTER}_by_${GROUP}"
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/denovo_vs_native_by_group.py" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --forced-csv "$WORK/forced_named_posteriors.csv" \
        --letter "$LETTER" \
        --group-by "$GROUP" \
        --top-dest "${TOP_DEST:-3}" \
        --min-cells "${MIN_CELLS:-25}" \
        ${DEST_ARG[@]+"${DEST_ARG[@]}"} \
        --output-dir "$OUT"

echo "Uploading the two tables to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/denovo_by_group"
for f in "$OUT"/*.csv; do
    s5cmd cp "$f" "${DEST}/${LETTER}_by_${GROUP}_$(basename "$f")"
done

echo "Done. The verdict is in the summary block printed above."
