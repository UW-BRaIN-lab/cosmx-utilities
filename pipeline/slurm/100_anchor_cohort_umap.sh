#!/bin/bash
# Figure: cohort-wide UMAP with the anchor cohort coloured by its InSituType (k=27) call and
# every cell outside the anchor cohort in grey. CPU only.
# Positions come from the typed AnnData's cohort-wide Stage-3c UMAP; labels come from the
# anchor cohort's own fit (anchor_typing.h5). See pipeline/python/plot_anchor_cohort_umap.py.
#
# Submit:
#   sbatch pipeline/slurm/100_anchor_cohort_umap.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   TYPED_DIR       Kopah sub-dir with the typed AnnData (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   STAGE4_DIR      Kopah sub-dir holding anchor/anchor_typing.h5 (default stage4_anchor_pruned, the
#                   k=27 fit on the pruned gene panel; stage4_anchor is the earlier full-panel fit).
#   OUT_DIR         where the PNGs are left on Klone (default $HOME/anchor_umap).

#SBATCH --job-name=cosmx-anchor-umap
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# ~7.5M points drawn twice (all labels, de novo letters) plus the AnnData obs index.
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=pipeline/logs/anchor_umap_%j.out
#SBATCH --error=pipeline/logs/anchor_umap_%j.err

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

TYPED="${TYPED_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
STAGE4="${STAGE4_DIR:-stage4_anchor_pruned}"
OUT_DIR="${OUT_DIR:-$HOME/anchor_umap}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_umap_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK" "$OUT_DIR"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging the typed AnnData and the anchor cohort's typing..."
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/plot_anchor_cohort_umap.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --anchor-typing-h5 "$WORK/anchor_typing.h5" \
        --label-colors-csv "${PIPELINE_DIR}/reference/anchor_label_colors.csv" \
        --output-stem "$WORK/umap_anchor_cohort"

cp "$WORK"/umap_anchor_cohort_*.png "$OUT_DIR/"
s5cmd cp "$WORK/umap_anchor_cohort_*.png" "${BASE}/${STAGE4}/anchor/"

echo "Done. Left in ${OUT_DIR} (and on Kopah under ${STAGE4}/anchor/)."
