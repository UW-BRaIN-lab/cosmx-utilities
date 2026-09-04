#!/bin/bash
# Diagnostic: why does DAPI read hazy where the Histone marker reads crisp?
#
# Reads only the per-cell morphology intensities already in the flat-file
# metadata, so it needs no exprMat and no new acquisition. One job streams the
# whole cohort (metadata files are small next to exprMat), rather than a job array.
#
# Submit from the repo root:
#   sbatch pipeline/slurm/11_morphology_contrast.sh
#
# Against another cohort:
#   MANIFEST=pipeline/manifest_retina.csv MAX_AREA=50000 \
#       sbatch pipeline/slurm/11_morphology_contrast.sh
#
# Optional env:
#   MAX_AREA     drop cells above this area; merged-cell blobs otherwise inflate the
#                intensity-vs-area correlation. Retina uses 50000, GBM QC uses 30000.
#   OUT_DIR      where results land (default <submit dir>/morphology_contrast)
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (path to the Python container image)

#SBATCH --job-name=cosmx-morphology-contrast
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/morphology_contrast_%j.out
#SBATCH --error=pipeline/logs/morphology_contrast_%j.err

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
    python "${PIPELINE_DIR}/python/morphology_channel_contrast.py" \
        --flatfiles-dir "$WORK/flat" \
        --manifest "$MANIFEST" \
        ${MAX_AREA:+--max-area "$MAX_AREA"} \
        --output "${OUT_DIR}/morphology_contrast_by_slide.csv" \
        --fov-output "${OUT_DIR}/morphology_contrast_by_fov.csv" \
        --figure "${OUT_DIR}/morphology_contrast.png"

echo "Done. Results in $OUT_DIR"
