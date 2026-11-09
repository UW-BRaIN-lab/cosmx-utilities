#!/bin/bash
# Reference triage 3 of 3: do the same cell type's profiles from different atlases carry
# REGION, or atlas signature?
#
# The combined reference holds three astrocytes, three microglia, two endothelia, two
# pericytes, two smooth muscles, two fibroblasts and two pigmented epithelia. If each copy
# wins in its own tissue they are biology and all stay; if they compete for the same tissue
# they are three noisy copies of one profile and only one belongs. See
# atlas_duplicates_by_region.py and reference/atlas_duplicate_classes.csv.
#
# Region comes from the committed FOV sheet joined on (slide_id, fov), NOT from obs['Region']
# -- a typed AnnData built before the latest annotation round carries the stale export value.
# Reads obs only (backed), so it is minutes at cohort scale. CPU only.
#
# Submit:
#   sbatch pipeline/slurm/79_atlas_duplicates.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR        Kopah sub-dir holding the typed AnnData (default stage4)
#   TYPED_BASENAME    the typed AnnData (default cosmx_typed.h5ad). Must be an INSITUTYPE
#                     typing -- the tree's leaves are a different vocabulary.
#   LABEL_KEY         obs column holding the call (default cell_type)
#   FOV_ANNOTATIONS_BASENAME  under pipeline/reference/ (default fov_annotations.csv)
#   CLASSES_BASENAME  under pipeline/reference/ (default atlas_duplicate_classes.csv)
#   INCLUDE_EXCLUDED  true keeps FOVs the SME ruled out (default false, as the refit will)
#   OUT_SUBDIR        Kopah sub-dir under STAGE4_DIR (default reference_triage)

#SBATCH --job-name=cosmx-atlas-duplicates
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/atlas_duplicates_%j.out
#SBATCH --error=pipeline/logs/atlas_duplicates_%j.err

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

STAGE4="${STAGE4_DIR:-stage4}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
OUT_SUBDIR="${OUT_SUBDIR:-reference_triage}"

# Both tables live in the repo, so the result reflects the committed rulings and a reader can
# tell which version produced it from the branch alone.
FOV_CSV="${PIPELINE_DIR}/reference/${FOV_ANNOTATIONS_BASENAME:-fov_annotations.csv}"
CLASSES_CSV="${PIPELINE_DIR}/reference/${CLASSES_BASENAME:-atlas_duplicate_classes.csv}"
[[ -f "$FOV_CSV" ]] || { echo "ERROR: no such FOV annotations: $FOV_CSV" >&2; exit 1; }
[[ -f "$CLASSES_CSV" ]] || { echo "ERROR: no such class table: $CLASSES_CSV" >&2; exit 1; }

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_atlas_duplicates_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${TYPED_BASENAME} from Kopah (${STAGE4})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" \
    "$WORK/typed.h5ad"

# NB `[[ ... ]] && cmd` as a bare statement returns 1 when the test fails, and under `set -e`
# that kills the job. Always use an if-block here.
EXCLUDED_ARG=()
if [[ "${INCLUDE_EXCLUDED:-false}" == "true" ]]; then
    EXCLUDED_ARG=(--include-excluded)
    echo "Keeping the FOVs the SME excluded (INCLUDE_EXCLUDED=true)"
fi

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/atlas_duplicates_by_region.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --fov-annotations "$FOV_CSV" \
        --classes "$CLASSES_CSV" \
        --label-key "${LABEL_KEY:-cell_type}" \
        ${EXCLUDED_ARG[@]+"${EXCLUDED_ARG[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading the atlas-duplicate tables to Kopah..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${OUT_SUBDIR}/atlas_duplicates/"

echo "Done: atlas duplicates by region. The tables are printed above."
