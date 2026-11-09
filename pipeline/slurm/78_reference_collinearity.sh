#!/bin/bash
# Reference triage 2 of 3: how far apart are the reference cell types, on the genes the fit
# actually sees? Two profiles that are nearly identical split a cell's evidence between them
# and BOTH lose, which is the mechanism behind an unassigned fraction that no K can fix.
#
# Needs no fit -- reference plus kept-gene list only, so it can be run BEFORE a refit to
# choose the type list, and after 73 to check what the new gene selection did to it.
#
# Submit (full panel):
#   sbatch pipeline/slurm/78_reference_collinearity.sh
# On the gene set a fit used, which is the number that decides a merge:
#   KEEP_GENES=stage4_anchor/gene_selection/kept_genes.txt \
#       sbatch pipeline/slurm/78_reference_collinearity.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR          Kopah sub-dir for outputs (default stage4)
#   REFERENCE_BASENAME  the reference under ${KOPAH_PREFIX}/reference/ -- the retina
#                       checkout's .env sets retina_combined_panel.csv
#   KEEP_GENES          Kopah key of a kept_genes.txt (73_select_genes.sh); unset = full panel
#   THRESHOLD           merge-candidate cosine bar (default 0.99)
#   REPORT_FLOOR        list pairs at or above this (default 0.95)
#   OUT_SUBDIR          Kopah sub-dir under STAGE4_DIR (default reference_triage)

#SBATCH --job-name=cosmx-ref-collinearity
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=pipeline/logs/reference_collinearity_%j.out
#SBATCH --error=pipeline/logs/reference_collinearity_%j.err

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
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"
OUT_SUBDIR="${OUT_SUBDIR:-reference_triage}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_ref_collinearity_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging reference ${REFERENCE_BASENAME}..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

# The pruned panel is what the fit sees, so a merge decision taken on the full panel can be
# wrong in either direction. Unset = full panel, and the report says which it used.
KEEP_ARG=()
if [[ -n "${KEEP_GENES:-}" ]]; then
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${KEEP_GENES}" "$WORK/kept_genes.txt"
    KEEP_ARG=(--keep-genes "$WORK/kept_genes.txt")
    echo "Scoring on the pruned panel from ${KEEP_GENES}"
else
    echo "No KEEP_GENES: scoring on the FULL panel (see the stage header)"
fi

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/reference_collinearity.py" \
        --reference "$WORK/reference.csv" \
        --threshold "${THRESHOLD:-0.99}" \
        --report-floor "${REPORT_FLOOR:-0.95}" \
        ${KEEP_ARG[@]+"${KEEP_ARG[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading the collinearity tables to Kopah..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${OUT_SUBDIR}/collinearity/"

echo "Done: reference collinearity."
