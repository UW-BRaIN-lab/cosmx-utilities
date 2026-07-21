#!/bin/bash
# Per-cell-type QC distributions from a typed AnnData (qc_by_celltype.py): is a given
# cell type (default the InSituTree Low_signal sink) low-RNA-quality or real-but-flat?
# CPU only; the rapids SIF carries anndata/pandas/matplotlib.
#
# Submit (defaults = InSituTree typed run, highlight Low_signal):
#   sbatch pipeline/slurm/51_qc_by_celltype.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir holding the typed AnnData (default stage4_insitutree)
#   TYPED_BASENAME (default cosmx_typed.h5ad)
#   GROUP_KEY    (default cell_type)   HIGHLIGHT (default Low_signal)
#   OUT_SUBDIR   Kopah sub-dir under STAGE4_DIR for outputs (default qc_by_celltype)

#SBATCH --job-name=cosmx-qc-celltype
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/qc_by_celltype_%j.out
#SBATCH --error=pipeline/logs/qc_by_celltype_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
OUT_SUBDIR="${OUT_SUBDIR:-qc_by_celltype}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_qc_celltype_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${TYPED_BASENAME} from Kopah (${STAGE4})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/qc_by_celltype.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --group-key "${GROUP_KEY:-cell_type}" \
        --highlight "${HIGHLIGHT:-Low_signal}" \
        --output-dir "$WORK/out"

echo "Uploading QC-by-celltype outputs to Kopah..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${OUT_SUBDIR}/"

echo "Done: QC by cell type."
