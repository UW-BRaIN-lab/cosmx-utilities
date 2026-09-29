#!/bin/bash
# The COHORT-WIDE Leiden cluster of every anchor-cohort cell, as a small per-cell CSV.
#
# Not anchor_cells.csv's `slide_cluster`: that is a PER-SLIDE clustering ("<slide>|<n>", 1,075
# strata over 57 slides) used only to stratify the subsample, so the same number means unrelated
# things on different slides. The cohort-wide Stage-3c clustering is in the 7.5M-cell typed
# AnnData's obs; this reads only that obs (backed, no matrix) and keeps the anchor cells.
#
# The result is left in ${OUT_DIR:-$HOME/anchor_napari} so it can be scp'd alongside the other
# Napari inputs, and also uploaded to Kopah for the record.
#
# Submit:
#   sbatch pipeline/slurm/99_anchor_leiden_labels.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   TYPED_DIR       Kopah sub-dir with the typed AnnData (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   LEIDEN_KEY      obs column holding the cohort Leiden (default leiden).
#   OUT_DIR         where the CSV is left on Klone (default $HOME/anchor_napari).

#SBATCH --job-name=cosmx-anchor-leiden
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
# Only obs is read from the 7.5M-cell AnnData (backed mode), plus 2.54M anchor ids.
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=pipeline/logs/anchor_leiden_%j.out
#SBATCH --error=pipeline/logs/anchor_leiden_%j.err

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

TYPED="${TYPED_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
OUT_DIR="${OUT_DIR:-$HOME/anchor_napari}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_leiden_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK" "$OUT_DIR"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging the typed AnnData and the anchor cell list..."
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"
s5cmd cp "${BASE}/stage4_anchor/anchor/anchor_cells.csv" "$WORK/anchor_cells.csv"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/extract_anchor_leiden.py" \
        --h5ad "$WORK/typed.h5ad" \
        --cell-ids-csv "$WORK/anchor_cells.csv" \
        --leiden-key "${LEIDEN_KEY:-leiden}" \
        --output-csv "$WORK/anchor_cohort_leiden.csv"

cp "$WORK/anchor_cohort_leiden.csv" "${OUT_DIR}/anchor_cohort_leiden.csv"
s5cmd cp "$WORK/anchor_cohort_leiden.csv" \
    "${BASE}/stage4_anchor/anchor/anchor_cohort_leiden.csv"

echo "Done. Left at ${OUT_DIR}/anchor_cohort_leiden.csv (and on Kopah)."
echo "The 'Slides spanned per cluster' block above is the check that this is cohort-wide."
