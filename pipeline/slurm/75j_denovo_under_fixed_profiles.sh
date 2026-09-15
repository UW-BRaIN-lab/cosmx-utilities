#!/bin/bash
# If we drop a de-novo letter from the hierarchy, where do its cells actually go?
#
# The triage justified each drop by showing the letter's markers match some GBmap type. That is
# necessary but NOT sufficient: it shows a leaf EXISTS, not that the FIXED profile can capture
# those cells. The sharded InSituTree run is the experiment for exactly that -- typed against
# fixed profiles only, it sent 52.9% of the cohort to Low_signal.
#
# Both runs already exist and both are label tables, so the answer is a join, not a rerun:
#
#   stage4_anchor_pruned/anchor/anchor_typing.h5   semi-supervised, the 27 de-novo letters
#   stage4_insitutree/cosmx_typed.h5ad             the SAME cells, FIXED profiles, Low_signal
#
# A dropped letter whose cells land on its matching leaf was a safe drop. One whose cells go to
# Low_signal well above baseline was not: the leaf existed and still could not hold them, and
# that letter needs its profile kept. This answers the drop question BEFORE the rebuild rather
# than after.
#
# No counts matrix, so this is minutes. The one thing that can sink it is the cell-id join --
# per-cell InSituCNV famously does NOT join to the full cohort because its id format differs --
# so the overlap is checked and reported before anything is computed.
#
# Submit:
#   sbatch pipeline/slurm/75j_denovo_under_fixed_profiles.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR      Kopah sub-dir with anchor/anchor_typing.h5 (default stage4_anchor_pruned).
#   TYPED_DIR       Kopah sub-dir with the fixed-profile typing (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   CELLTYPE_KEY    obs column holding the fixed-profile label (default cell_type).

#SBATCH --job-name=cosmx-denovo-fate
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
# Only obs is read from the 7.5M-cell h5ad (backed mode), plus the 2.54M-row anchor labels.
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=pipeline/logs/denovo_fate_%j.out
#SBATCH --error=pipeline/logs/denovo_fate_%j.err

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
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_fate_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging both label tables from Kopah (no counts matrix)..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/${TYPED_BASENAME}"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/denovo_under_fixed_profiles.py" \
        --anchor-h5 "$WORK/anchor_typing.h5" \
        --typed-h5ad "$WORK/${TYPED_BASENAME}" \
        --celltype-key "${CELLTYPE_KEY:-cell_type}" \
        --output-dir "$WORK/denovo_fate"

echo "Uploading tables to Kopah..."
for f in "$WORK"/denovo_fate/*.csv; do
    s5cmd cp "$f" "${BASE}/${STAGE4}/supervised_gbmap/denovo_fate/$(basename "$f")"
done

echo "Done. The per-letter table printed above is the drop decision."
