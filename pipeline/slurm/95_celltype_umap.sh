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
# For the InSituTree typing (no de-novo letters, so no annotation table):
#   RUN=stage4 TYPED_BASENAME=cosmx_typed_tree.h5ad \
#       OUTPUT_BASENAME=umap_insitutree_75.png TITLE="InSituTree, 75 leaves" \
#       sbatch pipeline/slurm/95_celltype_umap.sh
#
# Knobs: TYPED_BASENAME, COLOR_KEY (obs column, default cell_type), OUTPUT_BASENAME, TITLE,
# NAPARI_COLORS (default true -- the same hash colours as the Napari columns),
# PANELS_OUTPUT (a second PNG, one small panel per label).
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

# Loads .env WITHOUT letting the ambient shell override the study's identity; see
# slurm/_load_env.sh.
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/slurm/_load_env.sh"

RUN="${RUN:?set RUN to the Kopah sub-prefix for the run, e.g. stage4_extl3_rescale}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
COLOR_KEY="${COLOR_KEY:-cell_type}"
OUTPUT_BASENAME="${OUTPUT_BASENAME:-umap_cell_type_annotated.png}"
TITLE="${TITLE:-$RUN}"

# An InSituTree typing has no de-novo letters to rename, so the annotation table is
# OPTIONAL. It is still required when one is named, so a typo cannot silently draw an
# unannotated plot.
PLOT_ARGS=()
if [[ -n "${ANNOTATION_TABLE:-}" || -f "${PIPELINE_DIR}/reference/denovo_annotations/${RUN}.csv" ]]; then
    MAPPING="${PIPELINE_DIR}/reference/denovo_annotations/${ANNOTATION_TABLE:-$RUN}.csv"
    [[ -f "$MAPPING" ]] || { echo "ERROR: no annotation table $MAPPING" >&2; exit 1; }
    PLOT_ARGS+=(--mapping "$MAPPING")
    echo "Renaming de-novo letters from ${MAPPING}"
else
    echo "No annotation table for ${RUN}; colouring by ${COLOR_KEY} as stored"
fi
# Same hash colours the Napari columns use, so a label is the same colour in both.
# NB `[[ ... ]] && cmd` as a bare statement returns 1 when the test fails, and under
# `set -e` that kills the job. Always use an if-block here.
if [[ "${NAPARI_COLORS:-true}" == "true" ]]; then
    PLOT_ARGS+=(--napari-colors)
fi

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_umap_${SLURM_JOB_ID:-local}"
if [[ -n "${PANELS_OUTPUT:-}" ]]; then
    PLOT_ARGS+=(--panels-output "$WORK/${PANELS_OUTPUT}")
fi
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${RUN}/${TYPED_BASENAME} from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/plot_celltype_umap.py" \
        --h5ad "$WORK/typed.h5ad" \
        --color-key "$COLOR_KEY" \
        --output "$WORK/${OUTPUT_BASENAME}" \
        --title "${TITLE}" \
        ${PLOT_ARGS[@]+"${PLOT_ARGS[@]}"}

echo "Uploading annotated UMAP to Kopah..."
s5cmd cp "$WORK/${OUTPUT_BASENAME}" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/qc_plots/${OUTPUT_BASENAME}"
if [[ -n "${PANELS_OUTPUT:-}" ]]; then
    s5cmd cp "$WORK/${PANELS_OUTPUT}" \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/qc_plots/${PANELS_OUTPUT}"
fi

echo "Done: annotated cell-type UMAP for ${RUN}"
