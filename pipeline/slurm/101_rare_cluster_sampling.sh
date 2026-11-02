#!/bin/bash
# Figure: how one small cohort-wide Leiden cluster reaches the anchor-cohort subsample.
# Panel A cohort UMAP with the cluster highlighted, B one slide in tissue (sampled cells filled,
# unsampled hollow), C the share of every cohort Leiden cluster that is in the anchor cohort. CPU only.
# See pipeline/python/plot_rare_cluster_sampling.py.
#
# Submit:
#   CLUSTER=21 sbatch pipeline/slurm/101_rare_cluster_sampling.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   CLUSTER         cohort-wide Leiden cluster to highlight (required).
#   SLIDE           slide for panel B (default: the slide with most cells of the cluster).
#   TYPED_DIR       Kopah sub-dir with the typed AnnData (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   STAGE4_DIR      Kopah sub-dir holding anchor/anchor_cells.csv (default stage4_anchor).
#   OUT_DIR         where the PNG and CSV are left on Klone (default $HOME/rare_cluster_sampling).

#SBATCH --job-name=cosmx-rare-cluster
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# obs of the 7.5M-cell AnnData (backed, no matrix) plus the 2.5M anchor ids.
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/rare_cluster_%j.out
#SBATCH --error=pipeline/logs/rare_cluster_%j.err

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

: "${CLUSTER:?set CLUSTER to the cohort-wide Leiden cluster, e.g. CLUSTER=21 sbatch ...}"
TYPED="${TYPED_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
STAGE4="${STAGE4_DIR:-stage4_anchor}"
OUT_DIR="${OUT_DIR:-$HOME/rare_cluster_sampling}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_rare_cluster_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK" "$OUT_DIR"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging the typed AnnData and the anchor cohort's cell list..."
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_cells.csv" "$WORK/anchor_cells.csv"

SLIDE_ARGS=()
if [[ -n "${SLIDE:-}" ]]; then
    SLIDE_ARGS=(--slide "$SLIDE")
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/plot_rare_cluster_sampling.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --anchor-cells-csv "$WORK/anchor_cells.csv" \
        --cluster "$CLUSTER" \
        "${SLIDE_ARGS[@]}" \
        --output-stem "$WORK/rare_cluster_${CLUSTER}"

cp "$WORK"/rare_cluster_"${CLUSTER}"*.png "$WORK"/rare_cluster_"${CLUSTER}"_sampling_by_cluster.csv "$OUT_DIR/"
s5cmd cp "$WORK/rare_cluster_${CLUSTER}*" "${BASE}/${STAGE4}/anchor/"

echo "Done. Left in ${OUT_DIR} (and on Kopah under ${STAGE4}/anchor/)."
