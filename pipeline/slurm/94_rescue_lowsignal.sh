#!/bin/bash
# Stage 5h (Phase 3): recluster the Low_signal pool + reassign — "is 48% a floor?". One more
# InSituType rescue iteration on the Low_signal cells only (semi-supervised vs the fixed
# InSituTree profiles + de-novo n_clusts range; no rescale/update). InSituType SIF, big-mem CPU.
#
#   sbatch pipeline/slurm/94_rescue_lowsignal.sh
#
# LONG / non-checkpointing run: the ckpt defaults in the #SBATCH block below are PREEMPTIBLE, and
# this ~12 h InSituType job has no mid-run checkpoint — a preemption requeues it from scratch and
# it can thrash. Run it on our dedicated, non-preemptible allocation instead (claims one idle
# L40S so the gpu partition accepts the CPU-only job; the GPU goes unused):
#   sbatch --account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1 \
#       pipeline/slurm/94_rescue_lowsignal.sh
# (carrolllab is a collaborator's allocation — do not use it. See memory: hyak-allocations.)
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   INPUT_DIR          insitutype_input.h5 sub-dir (default stage4, shared).
#   STAGE4_DIR         insitutree_result.h5 sub-dir (labels; default stage4_insitutree).
#   PROFILES_BASENAME  fixed profiles in reference/ (default insitutree_profiles.csv).
#   STAGE5_DIR         stage-5 output sub-dir (default stage5_insitucnv; writes rescue/).
#   N_CLUSTS           de-novo range (default 10:20).

#SBATCH --job-name=cosmx-rescue-lowsignal
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=pipeline/logs/rescue_lowsignal_%j.out
#SBATCH --error=pipeline/logs/rescue_lowsignal_%j.err

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

INPUT_DIR="${INPUT_DIR:-stage4}"
STAGE4="${STAGE4_DIR:-stage4_insitutree}"
PROFILES_BASENAME="${PROFILES_BASENAME:-insitutree_profiles.csv}"
STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
N_CLUSTS="${N_CLUSTS:-10:20}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_rescue_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging counts + InSituTree labels + profiles from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/insitutype_input.h5" \
    "$WORK/insitutype_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_result.h5" \
    "$WORK/insitutree_result.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${PROFILES_BASENAME}" \
    "$WORK/profiles.csv"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/rescue_lowsignal.R" \
        --input "$WORK/insitutype_input.h5" \
        --labels-h5 "$WORK/insitutree_result.h5" \
        --profiles "$WORK/profiles.csv" \
        --n-clusts "$N_CLUSTS" \
        --output-csv "$WORK/rescue_lowsignal.csv" \
        --output-profiles "$WORK/rescued_profiles.csv" \
        --output-rds "$WORK/rescue_lowsignal.rds"

echo "Uploading rescue results to Kopah (${STAGE5}/rescue)..."
s5cmd cp "$WORK/rescue_lowsignal.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/rescue_lowsignal.csv"
s5cmd cp "$WORK/rescued_profiles.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/rescued_profiles.csv"
s5cmd cp "$WORK/rescue_lowsignal.rds" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/rescue_lowsignal.rds" \
    || echo "WARN: failed to upload rescue_lowsignal.rds"

echo "Done: stage 5h Low_signal rescue."
