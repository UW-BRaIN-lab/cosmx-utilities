#!/bin/bash
# Why did a cell choose a de-novo letter over a named GBmap type, and is the boundary a real gap?
#
# The end of the chain that started with the PI's question. 75d compared the two groups on
# markers (definitional, per its own header); 75h ruled out donor and batch and found the
# question inverts -- ~91% of vascular cells go to `l` in EVERY one of the 38 donors, so the
# native cells are the exception, not the rule. This job goes at the decision itself.
#
# PART 1 (always runs, model-free): the assignment is argmax over the 81 stored log-likelihoods,
# so margin = loglik(letter) - loglik(destination) is literally why each cell went where it did.
# A and B sit either side of zero BY CONSTRUCTION, so the question is not whether they separate
# but where zero falls in the DENSITY -- through the peak of one population (arbitrary cut) or in
# a trough (two populations).
#
# PART 2 (gated): the likelihood is a sum over genes, so the margin decomposes exactly, per gene.
# That means recomputing InSituType's likelihood, which is the one step that can be wrong, so it
# is validated against the stored margin and REFUSES to emit a decomposition that does not
# reproduce it. A first-run validation failure is a normal outcome, not a bug: the log prints the
# lldist source so the model can be corrected. Part 1's outputs are written either way.
#
# Submit:
#   LETTER=l DESTINATION=Endo_capilar sbatch pipeline/slurm/75i_denovo_margin_decomposition.sh
# Part 1 only (much lighter -- no counts read):
#   LETTER=l DESTINATION=Endo_capilar SKIP_GENES=1 \
#     sbatch pipeline/slurm/75i_denovo_margin_decomposition.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir with anchor/ (default stage4_anchor_pruned).
#   INPUT_DIR    Kopah sub-dir with anchor/anchor_input.h5 (default stage4_anchor).
#   LETTER       de-novo letter (default l).
#   DESTINATION  named GBmap type to compare against (default Endo_capilar).
#   SUBSAMPLE    cells per group for Part 2 (default 20000).
#   NB_SIZE      negative-binomial size for the recomputed likelihood (default 10).
#   SKIP_GENES   set to 1 to run Part 1 only.

#SBATCH --job-name=cosmx-margin-decomp
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# The stored loglik matrix is ~2.5M x 81 doubles (~1.6GB) and Part 2 additionally reads the
# 9.6GB counts and densifies a subsample-by-gene block three times over.
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/margin_decomp_%j.out
#SBATCH --error=pipeline/logs/margin_decomp_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_anchor_pruned}"
INPUT="${INPUT_DIR:-stage4_anchor}"
LETTER="${LETTER:-l}"
DESTINATION="${DESTINATION:-Endo_capilar}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_margin_decomp_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

# The .rds, not the h5 — the compact h5 drops $logliks, and the whole job is about them.
echo "Staging the full insitutype result from Kopah..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.rds" "$WORK/anchor_typing.rds"

COUNTS_ARG=()
if [[ "${SKIP_GENES:-0}" != "1" ]]; then
    echo "Staging the anchor counts (Part 2)..."
    s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
    COUNTS_ARG=(--counts-h5 "$WORK/anchor_input.h5")
fi

OUT="$WORK/${LETTER}_vs_${DESTINATION}"
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/denovo_margin_decomposition.R" \
        --typing-rds "$WORK/anchor_typing.rds" \
        --letter "$LETTER" \
        --destination "$DESTINATION" \
        --subsample "${SUBSAMPLE:-20000}" \
        --nb-size "${NB_SIZE:-10}" \
        ${COUNTS_ARG[@]+"${COUNTS_ARG[@]}"} \
        --output-dir "$OUT"

echo "Uploading tables to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/margin_decomposition"
for f in "$OUT"/*.csv; do
    s5cmd cp "$f" "${DEST}/${LETTER}_vs_${DESTINATION}_$(basename "$f")"
done

echo "Done. Part 1's verdict is the density_at_zero block printed above."
