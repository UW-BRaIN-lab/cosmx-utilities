#!/bin/bash
# Stage 4b: semi-supervised InSituType cell typing against the GBmap reference (R).
# Single big-memory CPU job — NOT a GPU job.
#
# Submit:
#   sbatch pipeline/slurm/80_insitutype.sh
# Chains into stage 4c:
#   jid=$(sbatch --parsable pipeline/slurm/80_insitutype.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/90_write_celltypes.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_INSITUTYPE (path to the insitutype.sif R container).
# The panel-restricted reference must already be on Kopah at
#   ${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}  (see pipeline/reference/README.md).

#SBATCH --job-name=cosmx-insitutype
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes the result files).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=pipeline/logs/insitutype_%j.out
#SBATCH --error=pipeline/logs/insitutype_%j.err

set -euo pipefail

# Compute nodes don't expose apptainer/s5cmd on PATH by default (see stage scripts).
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
# Where to READ insitutype_input.h5 from. Defaults to the output prefix, but a variant
# run (e.g. refit, a cluster sweep) can read the shared input from the original prefix
# while writing results to its own STAGE4_DIR — avoids duplicating the multi-GB input
# (Kopah server-side copy caps at 5GB anyway). E.g. INPUT_DIR=stage4 STAGE4_DIR=stage4_refit.
INPUT_DIR="${INPUT_DIR:-$STAGE4}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutype_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging InSituType input (from ${INPUT_DIR}) + reference from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/insitutype_input.h5" \
    "$WORK/insitutype_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutype_typing.R" \
        --input "$WORK/insitutype_input.h5" \
        --reference "$WORK/reference.csv" \
        --output-rds "$WORK/insitutype_result.rds" \
        --output-h5 "$WORK/insitutype_result.h5" \
        --n-clusts "${N_CLUSTS:-10:20}" \
        --update-reference "${UPDATE_REFERENCE:-true}" \
        --rescale "${RESCALE:-true}" \
        --refit "${REFIT:-false}" \
        --refinement "${REFINEMENT:-false}" \
        --n-starts "${N_STARTS:-10}" \
        --max-iters "${MAX_ITERS:-40}"

echo "Uploading InSituType result to Kopah..."
s5cmd cp "$WORK/insitutype_result.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutype_result.h5"
# The full .rds (incl. logliks) is the object for a later interactive refineClusters
# pass to condense over-split de novo clusters; non-fatal if upload lags.
s5cmd cp "$WORK/insitutype_result.rds" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutype_result.rds" \
    || echo "WARN: failed to upload insitutype_result.rds"

echo "Done: stage 4b"
