#!/bin/bash
# InSituTree REFERENCE REBUILD, fit step: semi-supervised InSituType (RESCALE=true) on the
# stratified ANCHOR (71) — the one expensive de-novo EM. Its post-rescale $profiles (54 Core-L4
# named types re-scaled to THIS 38-donor cohort + fresh de-novo clusters capturing patient-private
# programs) become the new InSituTree reference via prep_insitutree_profiles.py. Reuses
# insitutype_typing.R on the bounded anchor (~2.3M, not all 7.5M — the full cohort's EM is
# intractable, but the Wenyu ~2.33M rescale completed, so anchor-scale is proven-feasible).
#
# Submit (after 71):  sbatch pipeline/slurm/72_anchor_typing.sh
# Then (Phase B, interactive): marker-inspect the anchor de-novo clusters, annotate keepers, and
#   run prep_insitutree_profiles.py on anchor_typing.h5 -> new insitutree_profiles.csv.
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE. The panel-restricted Core-L4
# GBmap reference must be on Kopah at ${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}.

#SBATCH --job-name=cosmx-anchor-typing
# The ~2.3M-cell rescale EM is LONG and does NOT checkpoint — a ckpt requeue restarts it
# from scratch and thrashes (cf. the Stage-5 rescue). Run on the dedicated, non-preemptible
# gpu-l40s slice. --gres=gpu:1 claims a GPU the CPU-only InSituType fit won't use, purely to
# schedule on the dedicated partition (Emily-approved). Adjust --qos if your allocation differs.
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
# 256G matches the pilot's 2.33M InSituType rescale; keeps the request under the
# glioblastoma/gpu-l40s association memory cap (320G + a co-running job hit AssocGrpMemLimit).
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=pipeline/logs/anchor_typing_%j.out
#SBATCH --error=pipeline/logs/anchor_typing_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_anchor}"
# Core GBmap level-4 (54 fine types) — the reference the InSituTree pilot validated. NOT the
# Extended-L3 set the two-pass used (InSituTree's branch-local competition makes the fine types
# safe). Swap this for a Neutrophil-augmented panel only if the PI wants neutrophils resolvable.
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_typing_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging anchor input + reference from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_input.h5" \
    "$WORK/anchor_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

# Optional FAQ pruned panel (73_select_genes.sh output). KEEP_GENES is a Kopah key under
# ${KOPAH_PREFIX}/ to a kept_genes.txt; unset = full panel (pilot behavior).
KEEP_ARG=()
if [[ -n "${KEEP_GENES:-}" ]]; then
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${KEEP_GENES}" "$WORK/kept_genes.txt"
    KEEP_ARG=(--keep-genes "$WORK/kept_genes.txt")
    echo "Restricting anchor fit to pruned panel from ${KEEP_GENES}"
fi

# n_clusts default 10:20 (pilot); raise via N_CLUSTS to the de-novo K the 74 sweep found. The
# de-novo EM discovers the tumor programs — shared and 38-donor patient-private — that become
# the Malignant leaves of the rebuilt reference.
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutype_typing.R" \
        --input "$WORK/anchor_input.h5" \
        --reference "$WORK/reference.csv" \
        --output-rds "$WORK/anchor_typing.rds" \
        --output-h5 "$WORK/anchor_typing.h5" \
        --n-clusts "${N_CLUSTS:-10:20}" \
        --update-reference "${UPDATE_REFERENCE:-true}" \
        --rescale "${RESCALE:-true}" \
        --refit "${REFIT:-false}" \
        --refinement "${REFINEMENT:-false}" \
        --n-starts "${N_STARTS:-10}" \
        --max-iters "${MAX_ITERS:-40}" \
        ${KEEP_ARG[@]+"${KEEP_ARG[@]}"}

echo "Uploading anchor typing to Kopah..."
s5cmd cp "$WORK/anchor_typing.h5" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.h5"
s5cmd cp "$WORK/anchor_typing.rds" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.rds" \
    || echo "WARN: failed to upload anchor_typing.rds"

echo "Done: PASS 1 anchor typing."
