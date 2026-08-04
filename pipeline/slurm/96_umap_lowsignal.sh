#!/bin/bash
# Stage 5i: UMAP of the Low_signal pool coloured by the de-novo (unsupervised) clusters, region,
# donor, and CNV call. GPU via rapids-singlecell (rapids-singlecell SIF, --nv), on our dedicated
# non-preemptible gpu-l40s allocation.
#
#   sbatch pipeline/slurm/96_umap_lowsignal.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   typed-AnnData sub-dir (default stage4_insitutree)
#   STAGE5_DIR   stage-5 sub-dir for rescue/ + diagnostics/ (default stage5_insitucnv; set to
#                stage5_insitucnv_n20w100 to colour by the corrected CNV call)
#   TYPED_BASENAME (default cosmx_typed.h5ad)

#SBATCH --job-name=cosmx-umap-lowsignal
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/umap_lowsignal_%j.out
#SBATCH --error=pipeline/logs/umap_lowsignal_%j.err

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
STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_umap_lowsignal_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging typed AnnData + rescue labels (+ CNV table) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/rescue/rescue_lowsignal.csv" \
    "$WORK/rescue_lowsignal.csv"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/diagnostics/cell_cnv_table.csv.gz" \
    "$WORK/cell_cnv_table.csv.gz" || echo "WARN: no cell_cnv_table (CNV colouring skipped)"

OPT_ARGS=()
[[ -f "$WORK/cell_cnv_table.csv.gz" ]] && OPT_ARGS+=(--cell-table "$WORK/cell_cnv_table.csv.gz")

apptainer exec --nv --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/umap_lowsignal.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --rescue "$WORK/rescue_lowsignal.csv" \
        ${OPT_ARGS[@]+"${OPT_ARGS[@]}"} \
        --output-h5ad "$WORK/out/lowsignal_umap.h5ad" \
        --output-dir "$WORK/out"

echo "Uploading UMAP to Kopah (${STAGE5}/umap)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/umap/"

echo "Done: stage 5i Low_signal UMAP."
