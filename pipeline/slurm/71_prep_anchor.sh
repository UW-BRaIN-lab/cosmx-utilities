#!/bin/bash
# Two-pass PASS 1 (4a-anchor): build the stratified anchor input from the clustered
# cohort AnnData. Needs a GPU: the anchor is stratified on PER-SLIDE Leiden clusters
# computed on the scPearsonPCA embedding via rapids-singlecell (same idiom as stage 3c
# 40_cluster.sh). Single L40S is sufficient — each slide's embedding is small.
#
# Submit:  sbatch pipeline/slurm/71_prep_anchor.sh
# Chain:   jid=$(sbatch --parsable pipeline/slurm/71_prep_anchor.sh)
#          sbatch --dependency=afterok:$jid pipeline/slurm/72_anchor_typing.sh
#
# Tune the anchor size with CAP_PER_STRATUM (cells per slide x cluster stratum); tune
# the per-slide clustering with RESOLUTION / N_NEIGHBORS.
# Required env (pipeline/.env): KOPAH_*, APPTAINER_RSC (run with --nv).

#SBATCH --job-name=cosmx-prep-anchor
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/prep_anchor_%j.out
#SBATCH --error=pipeline/logs/prep_anchor_%j.err

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

STAGE3="${STAGE3_DIR:-stage3}"
# Isolate the reference-rebuild outputs under their own stage dir (anchor input + fit).
STAGE4="${STAGE4_DIR:-stage4_anchor}"
CLUSTERED="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_prep_anchor_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${CLUSTERED} from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED}" \
    "$WORK/${CLUSTERED}"

apptainer exec --nv \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/prep_insitutype_anchor.py" \
        --clustered-h5ad "$WORK/${CLUSTERED}" \
        --output "$WORK/anchor_input.h5" \
        --cluster-mode "${CLUSTER_MODE:-per-slide}" \
        --resolution "${RESOLUTION:-1.2}" \
        --n-neighbors "${N_NEIGHBORS:-15}" \
        --cluster-key "${CLUSTER_KEY:-leiden}" \
        --slide-key "${SLIDE_KEY:-slide_id}" \
        --cap-per-stratum "${CAP_PER_STRATUM:-2000}" \
        --seed "${SEED:-0}"

echo "Uploading anchor input to Kopah..."
s5cmd cp "$WORK/anchor_input.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_input.h5"
s5cmd cp "$WORK/anchor_cells.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_cells.csv"

echo "Done: PASS 1 anchor input."
