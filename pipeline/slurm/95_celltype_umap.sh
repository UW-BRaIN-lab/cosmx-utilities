#!/bin/bash
# Stage 4 figure: annotated cell-type UMAP for one InSituType run. CPU only.
# Annotates the de-novo letters on the fly from the run's denovo_annotations table and
# draws a readable, count-sorted legend (pipeline/python/plot_celltype_umap.py).
#
# Submit one per run (RUN = the run's Kopah sub-prefix; a matching annotation table must
# exist at pipeline/reference/denovo_annotations/<RUN>.csv):
#   RUN=stage4_extl3_rescale sbatch pipeline/slurm/95_celltype_umap.sh
#   for r in stage4 stage4_refit stage4_ext_l3 stage4_extl3_rescale; do
#     RUN=$r sbatch pipeline/slurm/95_celltype_umap.sh; done
#
# ANNOTATION_TABLE overrides the RUN-derived table name. Needed when two studies share a
# sub-prefix: the retina cohort's typing also lives under stage4/, so RUN=stage4 alone would
# label its clusters with the GBM annotations. The retina run is:
#   RUN=stage4 ANNOTATION_TABLE=retina_stage4 sbatch pipeline/slurm/95_celltype_umap.sh
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-celltype-umap
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/celltype_umap_%j.out
#SBATCH --error=pipeline/logs/celltype_umap_%j.err

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

RUN="${RUN:?set RUN to the Kopah sub-prefix for the run, e.g. stage4_extl3_rescale}"
MAPPING="${PIPELINE_DIR}/reference/denovo_annotations/${ANNOTATION_TABLE:-$RUN}.csv"
[[ -f "$MAPPING" ]] || { echo "ERROR: no annotation table $MAPPING" >&2; exit 1; }

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_umap_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${RUN}/cosmx_typed.h5ad from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/cosmx_typed.h5ad" "$WORK/typed.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/plot_celltype_umap.py" \
        --h5ad "$WORK/typed.h5ad" \
        --mapping "$MAPPING" \
        --output "$WORK/umap_cell_type_annotated.png" \
        --title "${RUN}"

echo "Uploading annotated UMAP to Kopah..."
s5cmd cp "$WORK/umap_cell_type_annotated.png" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/qc_plots/umap_cell_type_annotated.png"

echo "Done: annotated cell-type UMAP for ${RUN}"
