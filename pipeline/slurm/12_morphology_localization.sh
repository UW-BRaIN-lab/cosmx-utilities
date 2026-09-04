#!/bin/bash
# Does a morphology channel mark nuclei, or wash across the section?
#
# Stage 11 answered what the per-cell metadata can answer. It cannot see an
# extracellular wash: Mean.<channel> is measured only inside cell masks. This
# stage scores the morphology TIFFs against CompartmentLabels instead, which does.
#
# Morphology TIFFs are ~150MB per FOV, so this samples FOVs rather than taking a
# whole slide. Four to six per slide is ample for a background measurement.
#
# Submit from the repo root:
#   CELLSTATS_URI=s3://brainlabkg/CosMx-GBM/<run>/<slide>/<scan>/CellStatsDir \
#   SLIDE_ID=7583G27583G7 FOVS=1,50,100,150 \
#       sbatch pipeline/slurm/12_morphology_localization.sh
#
# Required env:
#   CELLSTATS_URI  s3:// prefix of the slide's CellStatsDir (no trailing slash)
#   SLIDE_ID       slide identifier, used to name the output
#   FOVS           comma-separated FOV numbers to sample
#
# Optional env:
#   SEGMENTATION_DIR  a Segmentation_<uuid>_<nnn> subdir to take CompartmentLabels
#                     from. Set this whenever the slide was resegmented, or the
#                     compartments will come from the ORIGINAL segmentation and
#                     will not match the flat files.
#   USE_SOURCE=1      read from the AWS source bucket instead of Kopah
#   CHANNEL_NAMES     e.g. "B=AT8 G=6E10", passed through as --channel-name
#   OUT_DIR           default <submit dir>/morphology_localization

#SBATCH --job-name=cosmx-morphology-localization
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/morphology_localization_%j.out
#SBATCH --error=pipeline/logs/morphology_localization_%j.err

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

: "${CELLSTATS_URI:?set CELLSTATS_URI to the CellStatsDir s3:// prefix}"
: "${SLIDE_ID:?set SLIDE_ID}"
: "${FOVS:?set FOVS, e.g. FOVS=1,50,100,150}"

if [[ "${USE_SOURCE:-0}" == "1" ]]; then
    export AWS_ACCESS_KEY_ID="${AWS_SOURCE_ACCESS_KEY_ID:?}"
    export AWS_SECRET_ACCESS_KEY="${AWS_SOURCE_SECRET_ACCESS_KEY:?}"
    unset S3_ENDPOINT_URL
else
    export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
    export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
    export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"
fi

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_localization_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/fovs"
trap 'rm -rf "$WORK"' EXIT

# CompartmentLabels live under CellStatsDir/FOV<nnnnn>/ for the original
# segmentation, and under CellStatsDir/<Segmentation_dir>/FOV<nnnnn>/ for a
# resegmentation. Morphology2D is never duplicated per segmentation.
LABEL_BASE="$CELLSTATS_URI"
if [[ -n "${SEGMENTATION_DIR:-}" ]]; then
    LABEL_BASE="${CELLSTATS_URI}/${SEGMENTATION_DIR}"
    echo "Taking CompartmentLabels from $SEGMENTATION_DIR"
else
    echo "WARN: SEGMENTATION_DIR unset - using the ORIGINAL segmentation's" >&2
    echo "      compartments, which will not match a resegmented slide" >&2
fi

STAGED=0
IFS=',' read -ra FOV_LIST <<< "$FOVS"
for FOV in "${FOV_LIST[@]}"; do
    PADDED=$(printf "%05d" "$FOV")
    echo "Staging FOV $PADDED ..."
    # The morphology filename carries a run-specific timestamp prefix, so glob it.
    if ! s5cmd cp "${CELLSTATS_URI}/Morphology2D/*_F${PADDED}.TIF" "$WORK/fovs/" </dev/null; then
        echo "WARN: no morphology TIFF for FOV $PADDED" >&2
        continue
    fi
    if ! s5cmd cp "${LABEL_BASE}/FOV${PADDED}/CompartmentLabels_F${PADDED}.tif" \
            "$WORK/fovs/" </dev/null; then
        echo "WARN: no CompartmentLabels for FOV $PADDED" >&2
        continue
    fi
    STAGED=$((STAGED + 1))
done

if [[ "$STAGED" -eq 0 ]]; then
    echo "ERROR: staged no FOVs from $CELLSTATS_URI" >&2
    exit 1
fi
echo "Staged $STAGED FOV(s)"

OUT_DIR="${OUT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}/morphology_localization}"
mkdir -p "$OUT_DIR"

CHANNEL_ARGS=()
for pair in ${CHANNEL_NAMES:-}; do
    CHANNEL_ARGS+=(--channel-name "$pair")
done

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    --bind "${OUT_DIR}:${OUT_DIR}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/morphology_localization.py" \
        --fov-dir "$WORK/fovs" \
        --slide-id "$SLIDE_ID" \
        ${CHANNEL_ARGS[@]+"${CHANNEL_ARGS[@]}"} \
        --output "${OUT_DIR}/${SLIDE_ID}_localization.csv"

echo "Done. Results in $OUT_DIR"
