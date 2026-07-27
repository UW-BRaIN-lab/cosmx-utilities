#!/bin/bash
# Stage 5g (Phase 2): hybrid state vs continuum among the CNV-high Low_signal cells. Scores the
# malignant + Neuronal_NLGN3 gene modules, joins the Phase-1 per-cell table + the flat
# posteriors, and produces the top-two structure / co-expression / spatial-null readouts.
# CPU-only; runs in the InSituCNV SIF (has scanpy + sklearn).
#
# Submit AFTER 91 (cell_cnv_table.csv.gz) and 92 (flat_posteriors.csv) have uploaded:
#   sbatch pipeline/slurm/93_insitucnv_hybrid_continuum.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR         stage-5 sub-dir (default stage5_insitucnv; reads diagnostics/ + posteriors/)
#   STAGE4_DIR         typed-AnnData sub-dir (default stage4_insitutree)
#   TYPED_BASENAME     (default cosmx_typed.h5ad)
#   SIGNATURES_BASENAME  reference/ gene modules (default gene_signatures.csv)
#   HIERARCHY_BASENAME   reference/ hierarchy (default insitutree_hierarchy.json)

#SBATCH --job-name=cosmx-hybrid-continuum
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/hybrid_continuum_%j.out
#SBATCH --error=pipeline/logs/hybrid_continuum_%j.err

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

STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
STAGE4="${STAGE4_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
SIGNATURES_BASENAME="${SIGNATURES_BASENAME:-gene_signatures.csv}"
HIERARCHY_BASENAME="${HIERARCHY_BASENAME:-insitutree_hierarchy.json}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_hybrid_continuum_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging typed AnnData + Phase-1 cell table + flat posteriors from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/diagnostics/cell_cnv_table.csv.gz" \
    "$WORK/cell_cnv_table.csv.gz"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/posteriors/flat_posteriors.csv" \
    "$WORK/flat_posteriors.csv"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/insitucnv_hybrid_continuum.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --cell-table "$WORK/cell_cnv_table.csv.gz" \
        --posteriors "$WORK/flat_posteriors.csv" \
        --signatures "${PIPELINE_DIR}/reference/${SIGNATURES_BASENAME}" \
        --hierarchy "${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}" \
        --output-dir "$WORK/out"

echo "Uploading hybrid-vs-continuum results to Kopah (${STAGE5}/hybrid)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/hybrid/"

echo "Done: stage 5g hybrid-vs-continuum."
