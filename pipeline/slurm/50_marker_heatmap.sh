#!/bin/bash
# Marker-heatmap compute: per-Leiden-cluster markers + pseudobulk z-score matrix
# from cosmx_clustered.h5ad. CPU only (no GPU); the rapids-singlecell SIF already
# carries anndata/scipy/pandas. Render the PDF afterwards with
# pipeline/R/marker_heatmap.R (locally, or in an R container).
#
# Submit:
#   sbatch pipeline/slurm/50_marker_heatmap.sh
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-marker-heatmap
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes its CSVs).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/marker_heatmap_%j.out
#SBATCH --error=pipeline/logs/marker_heatmap_%j.err

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

# Kopah sub-prefix for Stage 3 outputs. Override (e.g. STAGE3_DIR=stage3_q100)
# to write a parallel run without clobbering the default `stage3/` results.
STAGE3="${STAGE3_DIR:-stage3}"
# Which clustered AnnData to read. Defaults to the stage-3c output; for the cell-type
# heatmap point it at the stage-4c typed file (with GROUP_KEY=cell_type STAGE3_DIR=stage4).
CLUSTERED_BASENAME="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_marker_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${CLUSTERED_BASENAME} from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED_BASENAME}" \
    "$WORK/clustered.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/marker_pseudobulk.py" \
        --clustered-h5ad "$WORK/clustered.h5ad" \
        --output-dir "$WORK/out" \
        --group-key "${GROUP_KEY:-leiden}" \
        --top-n "${TOP_N:-5}" \
        --min-group-n "${MIN_GROUP_N:-10}"

echo "Uploading marker-heatmap CSVs to Kopah..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/marker_heatmap/"

echo "Done: marker-heatmap compute. Render with pipeline/R/marker_heatmap.R."
