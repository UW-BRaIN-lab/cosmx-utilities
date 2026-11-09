#!/bin/bash
# Per-(donor, region) QC from the stage-3 cohort (qc_by_region.py): is ONE region of ONE
# donor usable, when the rest of that donor is not? Reads obs only (backed), so it is
# minutes even at cohort scale. CPU only; the rapids SIF carries anndata/pandas.
#
# The slide-level QC table cannot answer this -- a slide carrying both eye and brain
# averages them. Region comes from the FOV annotation sheet, joined on (slide_id, fov),
# so it reflects the SME's rulings rather than whatever the upstream export said.
#
# Submit (whole cohort, no focus donor):
#   sbatch pipeline/slurm/21_qc_by_region.sh
#
# One donor against the other donors' SAME region, which is the comparison that means
# something (optic nerve is sparser than retina in every donor):
#   FOCUS_SLIDE=uwa7761eyes sbatch pipeline/slurm/21_qc_by_region.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE3_DIR    Kopah sub-dir holding the cohort (default stage3)
#   COMBINED_BASENAME  (default combined_qc.h5ad)
#   FOV_ANNOTATIONS_BASENAME  under pipeline/reference/ (default fov_annotations.csv)
#   FOCUS_SLIDE   a slide_id to compare against the rest of the cohort (default none)
#   FLOORS        comma-separated panel-count floors (default 50,100,150,200)
#   OUT_SUBDIR    Kopah sub-dir under STAGE3_DIR for outputs (default qc_by_region)

#SBATCH --job-name=cosmx-qc-region
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/qc_by_region_%j.out
#SBATCH --error=pipeline/logs/qc_by_region_%j.err

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

STAGE3="${STAGE3_DIR:-stage3}"
COMBINED_BASENAME="${COMBINED_BASENAME:-combined_qc.h5ad}"
OUT_SUBDIR="${OUT_SUBDIR:-qc_by_region}"

# The annotations live in the repo, not on Kopah, so the table reflects the committed
# rulings and a reader can tell which version produced it from the branch alone.
FOV_CSV="${PIPELINE_DIR}/reference/${FOV_ANNOTATIONS_BASENAME:-fov_annotations.csv}"
[[ -f "$FOV_CSV" ]] || { echo "ERROR: no such FOV annotations: $FOV_CSV" >&2; exit 1; }

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_qc_region_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${COMBINED_BASENAME} from Kopah (${STAGE3})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${COMBINED_BASENAME}" \
    "$WORK/combined_qc.h5ad"

FOCUS_ARG=()
if [[ -n "${FOCUS_SLIDE:-}" ]]; then
    FOCUS_ARG=(--focus-slide "$FOCUS_SLIDE")
    echo "Focus donor: ${FOCUS_SLIDE}"
fi

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/qc_by_region.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --fov-annotations "$FOV_CSV" \
        --floors "${FLOORS:-50,100,150,200}" \
        ${FOCUS_ARG[@]+"${FOCUS_ARG[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading QC-by-region tables to Kopah..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${OUT_SUBDIR}/"

echo "Done: QC by region. The tables are also printed above, so the log alone answers it."
