#!/bin/bash
# Stage 4a (sharded): split combined_qc.h5ad into per-slide InSituType/InSituTree inputs
# (one <slide_id>.h5 each) for the full-cohort SHARDED typing arrays. Uses the already-QC'd
# cohort cells so typing is consistent with the rest of the pipeline. CPU; backed read
# keeps memory modest even at ~7.5M cells.
#
# The per-slide inputs land in a NEUTRAL stage-4a location (default stage4/per_slide_inputs)
# and are consumed by BOTH the sharded InSituTree array (85b) and the sharded flat-posteriors
# array (92b) — they are stage-4a artifacts, not InSituTree-specific.
#
# Full-cohort run: isolate outputs with a separate KOPAH_PREFIX (e.g. cosmx-gbm-full) in
# pipeline/.env, exactly as the full-cohort Stage 1-3 run did — do NOT clobber the Wenyu run.
#
# Submit:  sbatch pipeline/slurm/84_split_insitutype_inputs.sh
# Chain (note the array count it prints):
#   jid=$(sbatch --parsable pipeline/slurm/84_split_insitutype_inputs.sh)
#   sbatch --dependency=afterok:$jid --array=1-<N> pipeline/slurm/85b_insitutree_sharded.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_RSC (rapids-singlecell SIF; CPU use).

#SBATCH --job-name=cosmx-split-inputs
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue restarts on
# preemption; the stage is idempotent (re-writes each per-slide <slide_id>.h5).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=pipeline/logs/split_inputs_%j.out
#SBATCH --error=pipeline/logs/split_inputs_%j.err

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

# Read combined_qc from the stage-3 dir; write per-slide inputs under the neutral stage-4a dir.
STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_split_inputs_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/per_slide_inputs"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging combined_qc.h5ad from Kopah (${STAGE3})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/combined_qc.h5ad" \
    "$WORK/combined_qc.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/split_insitutype_inputs.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --output-dir "$WORK/per_slide_inputs" \
        --slide-key "${SLIDE_KEY:-slide_id}"

echo "Uploading per-slide inputs to Kopah (${STAGE4}/per_slide_inputs)..."
s5cmd cp "$WORK/per_slide_inputs/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_inputs/"

n=$(ls "$WORK/per_slide_inputs" | wc -l)
echo "Done: sharded stage-4a split — ${n} per-slide inputs."
echo "Next: sbatch --array=1-${n} pipeline/slurm/85b_insitutree_sharded.sh"
