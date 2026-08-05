#!/bin/bash
# Stage 4 (alt, SHARDED): SUPERVISED hierarchical cell typing with InSituTree, one task per
# slide. Full-cohort (~7.5M cell) variant of 85_insitutree.sh — instead of one big-memory
# monolith, a Slurm ARRAY types each slide's cells against the SAME fixed profile matrix +
# hierarchy, then 85c concatenates the per-slide results.
#
# Why sharding is valid here: InSituTree is supervised (no cohort-level de-novo/EM phase) and
# its per-node gene selection is driven by the FIXED profiles, not the data — so a cell's node
# assignments depend only on its own counts + the fixed reference. Cells are independent; typing
# one slide at a time is semantically equivalent to the monolith EXCEPT that InSituType's
# per-cell background (estimateBackground) is estimated within each slide rather than cohort-wide
# (the slide is the natural batch — a defensible, arguably better difference). Sharding keeps each
# task small (fits ckpt, requeue-safe) and removes the ~7.5M-cell / 377G-node memory risk.
#
# Reference: reuse the validated pilot profiles (Core-L4 GBmap + 6 de-novo tumor states) —
# insitutree_profiles.csv + insitutree_hierarchy.json, UNCHANGED. InSituTree's branch-local
# competition is exactly what makes the 54 fine Core-L4 types safe at scale (no attractor blowup).
#
# Full-cohort run: isolate outputs with a separate KOPAH_PREFIX (e.g. cosmx-gbm-full) in
# pipeline/.env (as Stage 1-3 did). Per-slide inputs come from 84's neutral stage-4a location.
#
# Submit (set --array to the count 84 printed):
#   sbatch --array=1-<N> pipeline/slurm/85b_insitutree_sharded.sh
# Chain into the concat + writeback after the whole array completes:
#   jid=$(sbatch --parsable --array=1-<N> pipeline/slurm/85b_insitutree_sharded.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/85c_concat_typing.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE (R SIF; carries InSituTree).

#SBATCH --job-name=cosmx-insitutree-shard
# Supervised + per-slide => small and fast, so the array types slides in parallel on ckpt.
# --requeue restarts on preemption; each task is idempotent (re-writes its per-slide result).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
# 128G gives margin on the largest 2-tissue slides (~0.5M cells); the 2.33M monolith used
# 256G, so ~0.5M/slide fits comfortably. Bump with `--mem=` if a dense slide ever OOMs.
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/insitutree_shard_%A_%a.out
#SBATCH --error=pipeline/logs/insitutree_shard_%A_%a.err

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

# Read per-slide inputs from the neutral stage-4a location; write per-slide results under the
# InSituTree stage dir (mirrors 85_insitutree.sh's INPUT_DIR / STAGE4_DIR split).
INPUT_DIR="${INPUT_DIR:-stage4}"
STAGE4="${STAGE4_DIR:-stage4_insitutree}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-insitutree_profiles.csv}"
HIERARCHY_BASENAME="${HIERARCHY_BASENAME:-insitutree_hierarchy.json}"
UPLOAD_RDS="${UPLOAD_RDS:-0}"   # 57 nested per-node .rds are heavy; off by default.

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitutree_shard_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

INPUT_PREFIX="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT_DIR}/per_slide_inputs/"

# This task's slide = the Nth per-slide input (sorted), keyed off the array index.
mapfile -t SLIDES < <(s5cmd ls "$INPUT_PREFIX" | awk '{print $NF}' | grep '\.h5$' | sort)
if (( ${#SLIDES[@]} == 0 )); then
    echo "ERROR: no per-slide inputs at $INPUT_PREFIX (run 84 first)" >&2; exit 1
fi
SLIDE_FILE="${SLIDES[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
if [[ -z "$SLIDE_FILE" ]]; then
    echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} > ${#SLIDES[@]} slides; "\
         "set --array=1-${#SLIDES[@]}" >&2
    exit 1
fi
SLIDE="${SLIDE_FILE%.h5}"
echo "Task ${SLURM_ARRAY_TASK_ID}/${#SLIDES[@]}: slide ${SLIDE}"

echo "Staging slide input + profile matrix from Kopah..."
s5cmd cp "${INPUT_PREFIX}${SLIDE_FILE}" "$WORK/${SLIDE_FILE}"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutree_typing.R" \
        --input "$WORK/${SLIDE_FILE}" \
        --reference "$WORK/reference.csv" \
        --hierarchy "${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}" \
        --output-rds "$WORK/${SLIDE}.result.rds" \
        --output-h5 "$WORK/${SLIDE}.result.h5" \
        --output-csv "$WORK/${SLIDE}.summary.csv" \
        --quantile-absolute "${QUANTILE_ABS:-0.5}" \
        --quantile-percent "${QUANTILE_PCT:-0}" \
        --excluded-genes "${EXCLUDED_GENES:-}"

echo "Uploading per-slide result to Kopah (${STAGE4})..."
s5cmd cp "$WORK/${SLIDE}.result.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_results/${SLIDE}.h5"
# Per-slide multi-level summaryAnnotation (annotLevel_k / probs_annotLevel_k) for optional
# coarser re-summary without re-running; non-fatal if the upload lags.
s5cmd cp "$WORK/${SLIDE}.summary.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_summaries/${SLIDE}.csv" \
    || echo "WARN: failed to upload ${SLIDE}.summary.csv"
if [[ "$UPLOAD_RDS" == "1" ]]; then
    s5cmd cp "$WORK/${SLIDE}.result.rds" \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/per_slide_rds/${SLIDE}.rds" \
        || echo "WARN: failed to upload ${SLIDE}.result.rds"
fi

echo "Done: sharded InSituTree slide ${SLIDE}."
