#!/bin/bash
# Stage 4b'': decompose de novo InSituType clusters against the Stage-3c Leiden partition,
# with the per-component QC depth that separates "population the reference lacks" from
# "pool of low-signal cells". Reads stage-3c + stage-4 directly, so it runs as soon as
# 80_insitutype.sh has finished; it needs neither the 4c writeback nor InSituTree.
#
# Submit:
#   sbatch pipeline/slurm/82_denovo_components.sh
# Decompose something other than the default letters:
#   DENOVO_CLUSTERS=o,a,m sbatch pipeline/slurm/82_denovo_components.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; CPU use, no --nv).
# Optional:
#   DENOVO_CLUSTERS=o            comma-separated de novo labels to decompose
#   PARTITION_COL=leiden         obs column holding the unsupervised partition
#   DEPTH_COLS=total_counts,qc_genes_detected
#   MIN_COMPONENT_PCT=1.0        roll components below this % of the cluster into __other__
#   REFERENCE_BASENAME=          reference profile CSV in pipeline/reference/; adds each
#                                component's nearest reference types (the column that says
#                                which hierarchy branch it belongs under). Unset = omit.
#   TOP_REFERENCE=2              nearest reference types to report per component
#   FOV_ANNOTATIONS_BASENAME=    FOV annotation CSV in pipeline/reference/; drops cells in
#                                FOVs the SME marked exclude=1 BEFORE anything is measured.
#                                Set it whenever the result feeds a reference profile.
#   CLUSTERED_BASENAME=cosmx_clustered.h5ad
#   RESULT_BASENAME=insitutype_result.h5

#SBATCH --job-name=cosmx-denovo-components
# Short and idempotent (re-writes its CSVs), so ckpt's preemption costs nothing here --
# unlike the InSituType fits, which must go to a dedicated partition.
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=00:45:00
# Without REFERENCE_BASENAME the AnnData is read BACKED and only obs is materialized, where
# 64G would do. With it, the cells of the clusters being decomposed are pulled into memory
# for the pseudobulk -- still a subset, never the whole matrix, but 128G leaves room for a
# letter the size of `o` or `a` (~100k cells x 6k probes).
#SBATCH --output=pipeline/logs/denovo_components_%j.out
#SBATCH --error=pipeline/logs/denovo_components_%j.err

set -euo pipefail

# Compute nodes don't expose apptainer/s5cmd on PATH by default (see stage scripts).
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

STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"
CLUSTERED="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"
RESULT_BASENAME="${RESULT_BASENAME:-insitutype_result.h5}"
DENOVO_CLUSTERS="${DENOVO_CLUSTERS:-o}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_components_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${CLUSTERED} (${STAGE3}) + ${RESULT_BASENAME} (${STAGE4}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED}" \
    "$WORK/clustered.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${RESULT_BASENAME}" \
    "$WORK/typing_result.h5"

# One --cluster per label, so the script can emit a table per letter in a single pass.
CLUSTER_ARGS=()
IFS=',' read -ra LABELS <<< "$DENOVO_CLUSTERS"
for label in "${LABELS[@]}"; do
    # Trim the ends only. `tr -d` would delete INTERNAL spaces too, turning a named type
    # like "Mona_T cell" into "Mona_Tcell", which matches nothing -- and the script would
    # then warn about a label that was never wrong.
    label="$(echo "$label" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    [[ -n "$label" ]] && CLUSTER_ARGS+=(--cluster "$label")
done
echo "Decomposing: ${LABELS[*]}"

# The reference lives in the repo (as 85_insitutree.sh reads the hierarchy), not on Kopah.
REFERENCE_ARG=()
if [[ -n "${REFERENCE_BASENAME:-}" ]]; then
    REFERENCE_CSV="${PIPELINE_DIR}/reference/${REFERENCE_BASENAME}"
    [[ -f "$REFERENCE_CSV" ]] || { echo "ERROR: no such reference: $REFERENCE_CSV" >&2; exit 1; }
    REFERENCE_ARG=(--reference-csv "$REFERENCE_CSV")
    echo "Scoring components against ${REFERENCE_BASENAME}"
fi

FOV_ARG=()
if [[ -n "${FOV_ANNOTATIONS_BASENAME:-}" ]]; then
    FOV_CSV="${PIPELINE_DIR}/reference/${FOV_ANNOTATIONS_BASENAME}"
    [[ -f "$FOV_CSV" ]] || { echo "ERROR: no such FOV annotations: $FOV_CSV" >&2; exit 1; }
    FOV_ARG=(--fov-annotations "$FOV_CSV")
    echo "Dropping SME-excluded FOVs per ${FOV_ANNOTATIONS_BASENAME}"
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/diagnose_denovo_components.py" \
        --clustered-h5ad "$WORK/clustered.h5ad" \
        --typing-h5 "$WORK/typing_result.h5" \
        --partition-col "${PARTITION_COL:-leiden}" \
        --depth-cols "${DEPTH_COLS:-total_counts,qc_genes_detected}" \
        --min-component-pct "${MIN_COMPONENT_PCT:-1.0}" \
        --top-reference "${TOP_REFERENCE:-2}" \
        --out-dir "$WORK/out" \
        ${REFERENCE_ARG[@]+"${REFERENCE_ARG[@]}"} \
        ${FOV_ARG[@]+"${FOV_ARG[@]}"} \
        "${CLUSTER_ARGS[@]}"

echo "Uploading tables to Kopah (${STAGE4}/denovo_components)..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/denovo_components/"

echo "Done: de novo component tables. Tables above; CSVs in ${STAGE4}/denovo_components/."
