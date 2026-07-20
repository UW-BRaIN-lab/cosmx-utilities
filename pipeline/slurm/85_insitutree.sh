#!/bin/bash
# Stage 4 (alt): SUPERVISED hierarchical cell typing with InSituTree (R).
# Single big-memory CPU job — NOT a GPU job. Reuses the SAME insitutype_input.h5 that
# InSituType uses (no new stage-4a prep); reads a Core-L4-GBmap + de-novo-tumor profile
# matrix and a committed cell-type hierarchy JSON.
#
# Submit:
#   sbatch pipeline/slurm/85_insitutree.sh
# Chains into the writeback (reuses 90 with the InSituTree result + its own STAGE4_DIR):
#   jid=$(sbatch --parsable pipeline/slurm/85_insitutree.sh)
#   STAGE4_DIR=stage4_insitutree RESULT_BASENAME=insitutree_result.h5 \
#       sbatch --dependency=afterok:$jid pipeline/slurm/90_write_celltypes.sh
#
# Prerequisite: build the profile matrix once and put it on Kopah reference/:
#   s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage4/insitutype_result.h5 .   # Core-L4 rescale run
#   uv run python pipeline/python/prep_insitutree_profiles.py \
#       --profiles-h5 insitutype_result.h5 \
#       --annotations pipeline/reference/denovo_annotations/stage4.csv \
#       --output pipeline/reference/insitutree_profiles.csv
#   s5cmd cp pipeline/reference/insitutree_profiles.csv \
#       s3://$KOPAH_BUCKET/$KOPAH_PREFIX/reference/insitutree_profiles.csv
# The hierarchy JSON is read straight from the repo checkout (bind-mounted), not Kopah.
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_INSITUTYPE (path to the insitutype.sif R container; carries InSituTree).

#SBATCH --job-name=cosmx-insitutree
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue
# restarts on preemption; the stage is idempotent (re-writes the result files).
# InSituTree is supervised (no de-novo phase) so this is faster than the InSituType run;
# if ckpt preemption thrashes it, rerun on the non-preemptible slice with --gpus=0:
#   sbatch --account=glioblastoma --partition=gpu-l40s --qos=glioblastoma-gpu-l40s \
#       --gpus=0 --cpus-per-task=16 pipeline/slurm/85_insitutree.sh
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=pipeline/logs/insitutree_%j.out
#SBATCH --error=pipeline/logs/insitutree_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_insitutree}"
# Where to READ insitutype_input.h5 from — defaults to the shared `stage4` prefix so the
# multi-GB input is not duplicated (same pattern as 80_insitutype.sh's INPUT_DIR).
INPUT_DIR="${INPUT_DIR:-stage4}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-insitutree_profiles.csv}"
HIERARCHY_BASENAME="${HIERARCHY_BASENAME:-insitutree_hierarchy.json}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutree_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging InSituTree input (from ${INPUT_DIR}) + profile matrix from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/insitutype_input.h5" \
    "$WORK/insitutype_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutree_typing.R" \
        --input "$WORK/insitutype_input.h5" \
        --reference "$WORK/reference.csv" \
        --hierarchy "${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}" \
        --output-rds "$WORK/insitutree_result.rds" \
        --output-h5 "$WORK/insitutree_result.h5" \
        --output-csv "$WORK/insitutree_summary.csv" \
        --quantile-absolute "${QUANTILE_ABS:-0.5}" \
        --quantile-percent "${QUANTILE_PCT:-0.5}" \
        --excluded-genes "${EXCLUDED_GENES:-}"

echo "Uploading InSituTree result to Kopah..."
s5cmd cp "$WORK/insitutree_result.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_result.h5"
s5cmd cp "$WORK/insitutree_summary.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_summary.csv" \
    || echo "WARN: failed to upload insitutree_summary.csv"
# The full .rds (nested per-node result) is the object for re-summarizing at a coarser
# level without re-running; non-fatal if the upload lags.
s5cmd cp "$WORK/insitutree_result.rds" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutree_result.rds" \
    || echo "WARN: failed to upload insitutree_result.rds"

echo "Done: stage 4 InSituTree"
