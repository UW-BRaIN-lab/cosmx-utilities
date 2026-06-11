#!/bin/bash
# Stage 3 QC: per-cluster batch-mixing report for a clustered AnnData. Quantifies
# whether each Leiden cluster draws broadly across donors (batch correction held)
# or is dominated by one/two patients (residual patient structure / patient-private
# biology). CPU only; the rapids-singlecell SIF already carries anndata/pandas.
#
# Submit:
#   sbatch pipeline/slurm/45_cluster_batch_entropy.sh
# Override the batch/cluster columns or the clustered object via env:
#   BATCH_COL=Case CLUSTER_KEY=leiden CLUSTERED_BASENAME=cosmx_clustered.h5ad
#
# Runs as a Slurm job (NOT on the login node): reading a multi-million-row obs in
# backed mode still materialises the full index + columns, which the login-node
# memory arbiter kills. ckpt big-mem is plenty.
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-cluster-entropy
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes its CSV).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=pipeline/logs/cluster_entropy_%j.out
#SBATCH --error=pipeline/logs/cluster_entropy_%j.err

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

# Read the clustered object from the stage-3 dir; write the report alongside it.
STAGE3="${STAGE3_DIR:-stage3}"
CLUSTERED="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_entropy_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${CLUSTERED} from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED}" \
    "$WORK/${CLUSTERED}"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/cluster_batch_entropy.py" \
        --h5ad "$WORK/${CLUSTERED}" \
        --cluster-key "${CLUSTER_KEY:-leiden}" \
        --batch-key "${BATCH_COL:-Case}" \
        --csv "$WORK/cluster_batch_entropy.csv"

echo "Uploading entropy report to Kopah..."
s5cmd cp "$WORK/cluster_batch_entropy.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/qc/cluster_batch_entropy.csv"

echo "Done: per-cluster batch-mixing report."
