#!/bin/bash
# Stage 4b′: post-typing diagnostics — spatial concordance + sibling adjudication.
# Reads the stage-4b result and the stage-3a cohort AnnData directly, so it can run as
# soon as 80_insitutype.sh finishes; it does NOT need stage 3b/3c or the 4c writeback.
#
# Submit:
#   sbatch pipeline/slurm/81_diagnose_typing.sh
# Chain it straight off the typing job:
#   jid=$(sbatch --parsable pipeline/slurm/80_insitutype.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/81_diagnose_typing.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; CPU use, no --nv).
# Optional:
#   REGION_COL=Region            obs column holding the anatomical label
#   HIERARCHY_BASENAME=          hierarchy JSON in pipeline/reference/; the terminal
#                                competition nodes become the sibling groups to
#                                adjudicate. Unset = concordance only.
#   DIAG_MARKERS=NEAT1           genes to report per sibling (assay-offset probes)
#   RESULT_BASENAME=insitutype_result.h5   set insitutree_result.h5 to diagnose 85's output

#SBATCH --job-name=cosmx-diagnose-typing
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue restarts on
# preemption; the stage is idempotent (re-writes its CSVs).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
# 64G is generous: the AnnData is read BACKED, so only obs plus a couple of gene columns
# are materialized. The floor is set by the obs table, not the count matrix.
#SBATCH --output=pipeline/logs/diagnose_typing_%j.out
#SBATCH --error=pipeline/logs/diagnose_typing_%j.err

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

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

STAGE3="${STAGE3_DIR:-stage3}"
STAGE4="${STAGE4_DIR:-stage4}"
RESULT_BASENAME="${RESULT_BASENAME:-insitutype_result.h5}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_diagnose_typing_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging combined_qc (${STAGE3}) + ${RESULT_BASENAME} (${STAGE4}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/combined_qc.h5ad" \
    "$WORK/combined_qc.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${RESULT_BASENAME}" \
    "$WORK/typing_result.h5"

# Hierarchy is read from the repo (as 85_insitutree.sh does), not from Kopah.
HIER_ARG=()
if [[ -n "${HIERARCHY_BASENAME:-}" ]]; then
    HIER_ARG=(--hierarchy "${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}")
    echo "Adjudicating sibling groups from ${HIERARCHY_BASENAME}"
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/diagnose_typing_concordance.py" \
        --combined-h5ad "$WORK/combined_qc.h5ad" \
        --typing-h5 "$WORK/typing_result.h5" \
        --region-col "${REGION_COL:-Region}" \
        --markers "${DIAG_MARKERS-NEAT1}" \
        --out-dir "$WORK/out" \
        ${HIER_ARG[@]+"${HIER_ARG[@]}"}

echo "Uploading diagnostics to Kopah (${STAGE4}/typing_diagnostics)..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/typing_diagnostics/"

echo "Done: typing diagnostics. Tables above; CSVs in ${STAGE4}/typing_diagnostics/."
