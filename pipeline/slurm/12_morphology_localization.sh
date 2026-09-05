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
# Submit from the repo root. The manifest already records where each slide's raw
# data lives, so normally you only name the slide:
#   MANIFEST=$HOME/cosmx-utilities-fullcohort/pipeline/manifest.csv \
#   SLIDE_ID=7329A67822A7 FOVS=1,50,100,150 \
#       sbatch pipeline/slurm/12_morphology_localization.sh
#
# Required env:
#   SLIDE_ID       slide identifier, matching manifest column 2
#   FOVS           comma-separated FOV numbers to sample
#   MANIFEST       pipeline manifest; CELLSTATS_URI is derived from its
#                  decoded_prefix column (9), which points at
#                  <export>/DecodedFiles/<slide>/<scan>/ -- CellStatsDir sits
#                  inside that, NOT beside DecodedFiles
#   CELLSTATS_URI  set this instead to point somewhere the manifest does not cover
#
# Optional env:
#   SEGMENTATION_DIR  a Segmentation_<uuid>_<nnn> subdir to take CompartmentLabels
#                     from. Set this whenever the slide was resegmented, or the
#                     compartments will come from the ORIGINAL segmentation and
#                     will not match the flat files.
#   USE_SOURCE=0      read from Kopah instead of the AWS source bucket. Defaults
#                     to source: only flatFiles were migrated to Kopah, so the
#                     morphology TIFFs exist ONLY in source S3.
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

# lmod's init references unset vars (LD_LIBRARY_PATH), which set -u turns into
# noise in the error log that looks like a real failure. Relax u just for it.
if ! command -v module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
    set -u
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

: "${SLIDE_ID:?set SLIDE_ID}"
: "${FOVS:?set FOVS, e.g. FOVS=1,50,100,150}"

# Derive the CellStatsDir from the manifest unless it was named explicitly. The
# decoded_prefix column is the slide base path process-slides.py discovers, and
# CellStatsDir is its child.
if [[ -z "${CELLSTATS_URI:-}" ]]; then
    : "${MANIFEST:?set MANIFEST or CELLSTATS_URI}"
    : "${SOURCE_S3_BUCKET:?SOURCE_S3_BUCKET missing from .env}"
    # build_manifest.py writes with csv.DictWriter, whose default line terminator
    # is CRLF, so the LAST column carries a trailing \r. Unstripped it lands
    # mid-URI and s5cmd silently matches nothing. Stage 11 reads column 8 and
    # never saw this; column 9 is the last one.
    DECODED=$(awk -F, -v s="$SLIDE_ID" '$2==s {sub(/\r$/, "", $9); print $9; exit}' "$MANIFEST")
    if [[ -z "$DECODED" ]]; then
        echo "ERROR: $SLIDE_ID has no decoded_prefix in $MANIFEST" >&2
        exit 1
    fi
    CELLSTATS_URI="s3://${SOURCE_S3_BUCKET}/${DECODED%/}/CellStatsDir"
    echo "Resolved CellStatsDir from manifest:"
    echo "  $CELLSTATS_URI"
fi

if [[ "${USE_SOURCE:-1}" == "1" ]]; then
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
    echo "WARN: SEGMENTATION_DIR unset - using the ORIGINAL segmentation," >&2
    echo "      whose compartments will not match a resegmented slide." >&2
    echo "      Segmentation directories present on this slide:" >&2
    s5cmd ls "${CELLSTATS_URI}/" </dev/null 2>/dev/null \
        | grep -o 'Segmentation_[^/]*' | sort -u | sed 's/^/        /' >&2 || true
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
