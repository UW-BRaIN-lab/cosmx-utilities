#!/bin/bash
# Were this cohort's slides all segmented the same way?
#
# AtoMx changed its compartment output mid-cohort: slides segmented 2026-04-16
# carry a `cytoplasm` compartment, slides from 2026-05-07 on do not, from an
# IDENTICAL config. If that update also moved cell BOUNDARIES then cell area
# shifts with the segmentation date, and every cross-slide comparison in the
# cohort carries a batch effect keyed to when segmentation ran, not to biology.
#
# Cell geometry is already in the flat-file metadata, so this needs no images.
#
# Submit from the repo root:
#   MANIFEST=$HOME/cosmx-utilities-fullcohort/pipeline/manifest.csv \
#       sbatch pipeline/slurm/13_segmentation_homogeneity.sh
#
# Optional env:
#   SLIDE_GROUPS  CSV of slide_id plus grouping columns -- use it to supply the
#                 segmentation dates read off the AtoMx UI, which are not in the
#                 flat files. Pair with GROUP_BY.
#   GROUP_BY      column to compare across (default: whichever of version /
#                 cellSegmentationSetId actually separates the slides)
#   OUT_DIR       default <submit dir>/segmentation_homogeneity

#SBATCH --job-name=cosmx-segmentation-homogeneity
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/segmentation_homogeneity_%j.out
#SBATCH --error=pipeline/logs/segmentation_homogeneity_%j.err

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
MANIFEST="${MANIFEST:-${PIPELINE_DIR}/manifest.csv}"

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_morph_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/flat"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

# Manifest columns: 2 = slide_id, 8 = flat_files_prefix (see 10_flatfiles_to_anndata.sh).
echo "Staging metadata files listed in $MANIFEST ..."
STAGED=0
while IFS=, read -r _ SLIDE_ID _ _ _ _ _ FLAT_PREFIX _; do
    [[ -z "${SLIDE_ID:-}" || "$SLIDE_ID" == "slide_id" ]] && continue
    SRC="s3://${KOPAH_BUCKET}/${FLAT_PREFIX}${SLIDE_ID}_metadata_file.csv.gz"
    # </dev/null: s5cmd would otherwise read the manifest off this loop's stdin
    # and silently eat the remaining slides.
    if s5cmd cp "$SRC" "$WORK/flat/" </dev/null; then
        STAGED=$((STAGED + 1))
    else
        echo "WARN: could not stage $SRC" >&2
    fi
done < "$MANIFEST"

if [[ "$STAGED" -eq 0 ]]; then
    echo "ERROR: staged no metadata files from $MANIFEST" >&2
    exit 1
fi
echo "Staged $STAGED metadata file(s)"

OUT_DIR="${OUT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}/morphology_contrast}"
mkdir -p "$OUT_DIR"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    --bind "${OUT_DIR}:${OUT_DIR}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/segmentation_homogeneity.py" \
        --flatfiles-dir "$WORK/flat" \
        --manifest "$MANIFEST" \
        ${SLIDE_GROUPS:+--slide-groups "$SLIDE_GROUPS"} \
        ${GROUP_BY:+--group-by "$GROUP_BY"} \
        --output "${OUT_DIR}/segmentation_by_slide.csv" \
        --figure "${OUT_DIR}/segmentation_area.png"

echo "Done. Results in $OUT_DIR"
