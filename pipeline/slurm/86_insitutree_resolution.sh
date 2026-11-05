#!/bin/bash
# Stage 4d′: read InSituTree's posteriors to find where down the tree it stopped resolving.
#
# Reads the summaryAnnotation CSV from 85 and the hierarchy it ran against. Cheap: one CSV,
# a few minutes. Run it after every InSituTree fit, because a leaf holding an implausibly
# large population and a leaf holding a real one look identical in the count table that 85
# prints -- only the posterior read against chance for its contest tells them apart.
#
# Submit:
#   HIERARCHY_BASENAME=retina_hierarchy.json sbatch pipeline/slurm/86_insitutree_resolution.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY, APPTAINER_RSC.
# Required args:
#   HIERARCHY_BASENAME           the hierarchy JSON the fit used
# Optional:
#   TREE_STAGE_DIR=stage4_insitutree      where 85 put its outputs
#   SUMMARY_BASENAME=insitutree_summary.csv
#   FLAT_RESULT_BASENAME=insitutype_result.h5   the FLAT InSituType run, in STAGE4_DIR.
#                                Set it: the posterior saturates at ~1.000 and says nothing,
#                                so how a leaf fared in OPEN competition is the real test.
#   DENOVO_ANNOTATIONS=          the keep/new_name CSV 84 used, so kept de novo leaves are
#                                matched to the letter they came from rather than looking
#                                like profiles that never won anything
#   TOP=30                                labels to print (all are written to CSV)

#SBATCH --job-name=cosmx-insitutree-resolution
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=pipeline/logs/insitutree_resolution_%j.out
#SBATCH --error=pipeline/logs/insitutree_resolution_%j.err

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

# Loads .env WITHOUT overwriting what the caller passed; see slurm/_load_env.sh.
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/slurm/_load_env.sh"

TREE_STAGE_DIR="${TREE_STAGE_DIR:-stage4_insitutree}"
SUMMARY_BASENAME="${SUMMARY_BASENAME:-insitutree_summary.csv}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"
: "${HIERARCHY_BASENAME:?set it to the hierarchy JSON the fit used}"

HIERARCHY_JSON="${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}"
[[ -f "$HIERARCHY_JSON" ]] || { echo "ERROR: no such hierarchy: $HIERARCHY_JSON" >&2; exit 1; }

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutree_resolution_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${SUMMARY_BASENAME} (${TREE_STAGE_DIR}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${TREE_STAGE_DIR}/${SUMMARY_BASENAME}" \
    "$WORK/summary.csv"

FLAT_ARGS=()
if [[ -n "${FLAT_RESULT_BASENAME:-}" ]]; then
    STAGE4="${STAGE4_DIR:-stage4}"
    echo "Staging the flat run ${FLAT_RESULT_BASENAME} (${STAGE4}) for comparison..."
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${FLAT_RESULT_BASENAME}" \
        "$WORK/flat_result.h5"
    FLAT_ARGS=(--flat-typing-h5 "$WORK/flat_result.h5")
    if [[ -n "${DENOVO_ANNOTATIONS:-}" ]]; then
        ANNOTATIONS_CSV="${PIPELINE_DIR}/reference/denovo_annotations/${DENOVO_ANNOTATIONS}"
        [[ -f "$ANNOTATIONS_CSV" ]] \
            || { echo "ERROR: no such annotations: $ANNOTATIONS_CSV" >&2; exit 1; }
        FLAT_ARGS+=(--denovo-annotations "$ANNOTATIONS_CSV")
    fi
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/diagnose_insitutree_resolution.py" \
        --summary-csv "$WORK/summary.csv" \
        --hierarchy "$HIERARCHY_JSON" \
        --top "${TOP:-30}" \
        ${FLAT_ARGS[@]+"${FLAT_ARGS[@]}"} \
        --out-dir "$WORK/out"

echo "Uploading resolution tables to Kopah (${TREE_STAGE_DIR}/resolution)..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${TREE_STAGE_DIR}/resolution/"

echo "Done: InSituTree resolution diagnostic."
