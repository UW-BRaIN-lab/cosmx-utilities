#!/bin/bash
# Stage 4b': anchor report — how many reference cell types FAIL TO ANCHOR, and why.
#
# insitutype() with update_reference_profiles=TRUE rescales each reference profile from its
# own ANCHOR CELLS, and silently discards any type left with <= insufficient_anchors_thresh
# (default 20) anchors — that type's profile is never adapted to the platform, so it rarely
# wins a cell and its cells land in the de novo sink instead. Nothing in the typing output
# records which types this happened to. This stage reproduces the anchor selection and
# reports it per type.
#
# CHEAP and DIAGNOSTIC-ONLY: it runs InSituType's anchor functions, not the EM, so it needs
# no typing run and can be pointed at a stage-4a input to explain a run that already
# happened. Start with SUBSAMPLE to get the picture in minutes.
#
# Submit:
#   SUBSAMPLE=200000 sbatch pipeline/slurm/82_anchor_report.sh      # fast preview
#   sbatch pipeline/slurm/82_anchor_report.sh                       # full cohort
# Diagnose a pruned-panel run (pass the SAME kept_genes.txt that run used):
#   KEEP_GENES=stage4_anchor/gene_selection/kept_genes.txt \
#       sbatch pipeline/slurm/82_anchor_report.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_INSITUTYPE (the insitutype.sif R container; carries InSituType).
# Optional:
#   STAGE4_DIR=stage4              where insitutype_input.h5 lives, and where the report
#                                  is written (under <STAGE4_DIR>/anchor_report/)
#   REFERENCE_BASENAME=gbmap_level4_panel.csv   reference under ${KOPAH_PREFIX}/reference/
#   KEEP_GENES=                    Kopah key of a kept_genes.txt (73_select_genes.sh)
#   SUBSAMPLE=0                    cells to sample (0 = all); a preview, unbiased because
#                                  anchor selection is per-type and threshold-based
#   REPORT_TAG=                    suffix for the output dir, to keep variants side by side
#                                  (e.g. REPORT_TAG=pruned -> anchor_report_pruned/)
#   The four anchor knobs, defaulting to insitutype()'s own values:
#   N_ANCHOR_CELLS=2000  MIN_ANCHOR_COSINE=0.3  MIN_ANCHOR_LLR=0.03
#   INSUFFICIENT_ANCHORS_THRESH=20

#SBATCH --job-name=cosmx-anchor-report
# glioblastoma has no dedicated CPU node -> ckpt (free, preemptible). --requeue restarts on
# preemption; the stage is idempotent (re-writes its CSV + summary).
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=04:00:00
# get_anchor_stats materializes TWO dense cells x types matrices (cos and llr). At 2.5M
# cells x ~75 types that is ~1.5GB each, on top of the cells x genes counts — hence 256G,
# the same envelope as the typing stage. SUBSAMPLE shrinks it proportionally.
#SBATCH --output=pipeline/logs/anchor_report_%j.out
#SBATCH --error=pipeline/logs/anchor_report_%j.err

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

STAGE4="${STAGE4_DIR:-stage4}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"
OUT_DIR_NAME="anchor_report${REPORT_TAG:+_${REPORT_TAG}}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_report_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging InSituType input (${STAGE4}) + reference (${REFERENCE_BASENAME}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/insitutype_input.h5" \
    "$WORK/insitutype_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

# Same optional pruned panel as 80_insitutype.sh / 72_anchor_typing.sh. Pass the SAME file
# the run being diagnosed used: pruning changes which genes carry the cosine, so it changes
# which types anchor.
KEEP_ARG=()
if [[ -n "${KEEP_GENES:-}" ]]; then
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${KEEP_GENES}" "$WORK/kept_genes.txt"
    KEEP_ARG=(--keep-genes "$WORK/kept_genes.txt")
    echo "Restricting the anchor report to the pruned panel from ${KEEP_GENES}"
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/anchor_report.R" \
        --input "$WORK/insitutype_input.h5" \
        --reference "$WORK/reference.csv" \
        --output-csv "$WORK/out/anchor_report.csv" \
        --output-summary "$WORK/out/anchor_report_summary.txt" \
        --n-anchor-cells "${N_ANCHOR_CELLS:-2000}" \
        --min-anchor-cosine "${MIN_ANCHOR_COSINE:-0.3}" \
        --min-anchor-llr "${MIN_ANCHOR_LLR:-0.03}" \
        --insufficient-anchors-thresh "${INSUFFICIENT_ANCHORS_THRESH:-20}" \
        --subsample "${SUBSAMPLE:-0}" \
        --seed "${SEED:-0}" \
        ${KEEP_ARG[@]+"${KEEP_ARG[@]}"}

echo "Uploading anchor report to Kopah (${STAGE4}/${OUT_DIR_NAME})..."
s5cmd cp "$WORK/out/*" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${OUT_DIR_NAME}/"

echo "Done: anchor report. Summary above; CSV in ${STAGE4}/${OUT_DIR_NAME}/."
