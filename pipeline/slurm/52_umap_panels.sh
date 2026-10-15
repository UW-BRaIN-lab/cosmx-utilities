#!/bin/bash
# Per-group UMAP panels (plot_umap_panels.py): one panel per slide / donor against a grey
# all-cells backdrop, plus the group-by-cluster share table. The all-in-one UMAP coloured by
# batch is overplotted and cannot show mixing; this can. CPU only, reads obs + obsm only.
#
# Submit (defaults: panels by slide_id and by Region):
#   sbatch pipeline/slurm/52_umap_panels.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE3_DIR     Kopah sub-dir holding the clustered AnnData (default stage3)
#   PANEL_KEYS     space-separated obs columns, one figure each (default "slide_id Region")
#   PANEL_NCOLS    panels per row (default 5; 8 reads better for 40 donors)
#   OUT_SUBDIR     Kopah sub-dir under STAGE3_DIR for the outputs (default umap_panels)

#SBATCH --job-name=cosmx-umap-panels
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/umap_panels_%j.out
#SBATCH --error=pipeline/logs/umap_panels_%j.err

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
OUT_SUBDIR="${OUT_SUBDIR:-umap_panels}"
PANEL_KEYS="${PANEL_KEYS:-slide_id Region}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_umap_panels_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging cosmx_clustered.h5ad from Kopah (${STAGE3})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/cosmx_clustered.h5ad" "$WORK/clustered.h5ad"

for key in $PANEL_KEYS; do
    apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
        python -u "${PIPELINE_DIR}/python/plot_umap_panels.py" \
            --h5ad "$WORK/clustered.h5ad" \
            --out-dir "$WORK/out" \
            --key "$key" \
            --ncols "${PANEL_NCOLS:-5}"
done

echo "Uploading panel figures to Kopah..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${OUT_SUBDIR}/"

echo "Done: UMAP panels."
