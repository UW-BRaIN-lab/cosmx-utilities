#!/bin/bash
# Compare two Stage 3 runs (e.g. a 50- vs 100-count QC pass). CPU only.
# Stages both runs' cosmx_clustered.h5ad from Kopah and diffs them into small
# comparison CSVs (see pipeline/python/compare_runs.py).
#
# Submit (defaults compare stage3 vs stage3_q100):
#   sbatch pipeline/slurm/60_compare_runs.sh
#   STAGE3_DIR_A=stage3 STAGE3_DIR_B=stage3_q100 sbatch pipeline/slurm/60_compare_runs.sh
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-compare-runs
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/compare_runs_%j.out
#SBATCH --error=pipeline/logs/compare_runs_%j.err

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

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

A="${STAGE3_DIR_A:-stage3}"
B="${STAGE3_DIR_B:-stage3_q100}"
# Which AnnData basename + obs key to compare. Defaults diff the Stage-3c Leiden runs;
# for two Stage-4 typed runs set CLUSTERED_BASENAME=cosmx_typed.h5ad GROUP_KEY=cell_type
# (e.g. STAGE3_DIR_A=stage4 STAGE3_DIR_B=stage4_refit).
CLUSTERED_BASENAME="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"
GROUP_KEY="${GROUP_KEY:-leiden}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_compare_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${CLUSTERED_BASENAME} for '$A' and '$B' from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${A}/${CLUSTERED_BASENAME}" "$WORK/a.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${B}/${CLUSTERED_BASENAME}" "$WORK/b.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/compare_runs.py" \
        --run-a "$WORK/a.h5ad" --label-a "$A" \
        --run-b "$WORK/b.h5ad" --label-b "$B" \
        --group-key "$GROUP_KEY" \
        --output-dir "$WORK/out"

echo "Uploading comparison CSVs to Kopah..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/comparisons/${A}_vs_${B}/"

echo "Done: comparison ${A} vs ${B}."
