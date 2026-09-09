#!/bin/bash
# Diagnostic for the 79.3% named control in 75/75b: decompose it from STORED data.
#
# 75b reported that only 79.3% of already-named anchor cells keep their name under the forced
# supervised run. Per the InSituType source that should be 100%: insitutype()'s final phase is
# itself an insitutypeML() call on all cells with the profiles it returns, same
# estimateBackground() path, one "all" cohort, reference_sds NULL for RNA — and dropping 27
# LOSING columns cannot change an argmax. So a premise is wrong, and this finds out which.
#
# insitutype() returns $logliks (cells x 81) and 72 saved the full result, so no re-scoring is
# needed — the answer is already in anchor_typing.rds. See the header of
# R/diagnose_forced_selfconsistency.R for the A/B/C decomposition this prints.
#
# Submit (defaults target the pruned k=27 run, after 75 has produced its posteriors):
#   sbatch pipeline/slurm/75c_diagnose_selfconsistency.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).

#SBATCH --job-name=cosmx-forced-selfconsist
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# The stored loglik matrix is ~2.5M x 81 doubles (~1.6GB) and the script holds a couple of
# working copies for the max.col passes, on top of whatever else the rds carries.
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/forced_selfconsist_%j.out
#SBATCH --error=pipeline/logs/forced_selfconsist_%j.err

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
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_forced_selfconsist_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}"

echo "Staging the FULL anchor typing result (.rds, has \$logliks) + 75's posteriors..."
s5cmd cp "${BASE}/anchor/anchor_typing.rds" "$WORK/anchor_typing.rds"
s5cmd cp "${BASE}/supervised_gbmap/supervised_gbmap_posteriors.csv" "$WORK/posteriors.csv"
ls -la "$WORK"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/diagnose_forced_selfconsistency.R" \
        --typing-rds "$WORK/anchor_typing.rds" \
        --posteriors "$WORK/posteriors.csv" \
        --output-csv "$WORK/forced_selfconsistency_by_label.csv" \
        --output-margins-csv "$WORK/forced_selfconsistency_margins.csv"

echo "Uploading diagnostic tables to Kopah (${STAGE4}/supervised_gbmap)..."
s5cmd cp "$WORK/forced_selfconsistency_by_label.csv" \
    "${BASE}/supervised_gbmap/forced_selfconsistency_by_label.csv"
s5cmd cp "$WORK/forced_selfconsistency_margins.csv" \
    "${BASE}/supervised_gbmap/forced_selfconsistency_margins.csv"

echo "Done. The A/B/C decomposition is printed above in this log."
