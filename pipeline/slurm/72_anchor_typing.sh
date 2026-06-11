#!/bin/bash
# Two-pass PASS 1 (4b-anchor): semi-supervised InSituType on the ANCHOR only — the one
# expensive de-novo EM. Reuses insitutype_typing.R, but on the bounded anchor (not all
# 7.5M cells), with a HIGHER n_clusts range so de novo can resolve both shared and
# patient-private malignant programs. Big-memory CPU job (ckpt).
#
# Submit:  sbatch pipeline/slurm/72_anchor_typing.sh
# Chain:   jid=$(sbatch --parsable pipeline/slurm/72_anchor_typing.sh)
#          sbatch --dependency=afterok:$jid pipeline/slurm/73_build_profile.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE. The panel-restricted
# GBmap reference must be on Kopah at ${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}.

#SBATCH --job-name=cosmx-anchor-typing
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=pipeline/logs/anchor_typing_%j.out
#SBATCH --error=pipeline/logs/anchor_typing_%j.err

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
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_typing_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging anchor input + reference from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_input.h5" \
    "$WORK/anchor_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

# Higher n_clusts than single-pass (default 20:40): the anchor must capture both shared
# and patient-private malignant programs, since PASS 2 only assigns to these types.
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutype_typing.R" \
        --input "$WORK/anchor_input.h5" \
        --reference "$WORK/reference.csv" \
        --output-rds "$WORK/anchor_typing.rds" \
        --output-h5 "$WORK/anchor_typing.h5" \
        --n-clusts "${N_CLUSTS:-20:40}" \
        --update-reference "${UPDATE_REFERENCE:-true}" \
        --rescale "${RESCALE:-true}" \
        --refit "${REFIT:-false}" \
        --refinement "${REFINEMENT:-false}" \
        --n-starts "${N_STARTS:-10}" \
        --max-iters "${MAX_ITERS:-40}"

echo "Uploading anchor typing to Kopah..."
s5cmd cp "$WORK/anchor_typing.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.h5"
s5cmd cp "$WORK/anchor_typing.rds" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.rds" \
    || echo "WARN: failed to upload anchor_typing.rds"

echo "Done: PASS 1 anchor typing."
