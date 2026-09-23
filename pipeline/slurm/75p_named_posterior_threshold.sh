#!/bin/bash
# Would the FAQ's 80% posterior threshold catch b's cells if the rebuild had no flat leaf?
#
# InSituTree has no Low_signal outcome of its own, so dropping b's Low_signal_denovo leaf means
# every b cell is forced onto a named leaf, and "Low_signal" becomes whatever we flag afterwards.
# The InSituType FAQ calls cells below a posterior of 0.8 "unclassified", but b's cells are
# near-ties (median margin 4.1 vs 25.1), which can still read as posteriors near 0.98. This job
# measures it from the fit's own stored logliks: drop the de-novo columns, take the posterior
# over the named ones, and report per label what share each posterior cut and margin cut flags.
#
# The posterior formula (plain softmax, or with a cluster-frequency prior) is checked against
# the fit's stored `prob` first; see R/named_posterior_threshold.R. Reads only the .rds -- no
# counts -- so it is light.
#
# Submit:
#   sbatch pipeline/slurm/75p_named_posterior_threshold.sh
# Another focus letter or cut-offs:
#   FOCUS=t THRESHOLDS=0.8,0.9 MARGINS=2,5,10 sbatch pipeline/slurm/75p_named_posterior_threshold.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/anchor_typing.rds (default stage4_anchor_pruned).
#   FOCUS       letter to cross-tabulate posterior against margin (default b).
#   THRESHOLDS  posterior cut-offs (default 0.8,0.95,0.99).
#   MARGINS     top1-vs-top2 log-likelihood cut-offs (default 2,5,10,25).

#SBATCH --job-name=cosmx-named-posterior
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# The stored loglik matrix is ~2.5M x 81 doubles (~1.6GB); blocks of 250k rows keep the
# softmax working copies small, but readRDS briefly holds the whole result.
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/named_posterior_%j.out
#SBATCH --error=pipeline/logs/named_posterior_%j.err

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
FOCUS="${FOCUS:-b}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_named_posterior_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

# The .rds, not the h5 -- the compact h5 drops $logliks.
echo "Staging the full insitutype result from Kopah..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.rds" "$WORK/anchor_typing.rds"

OUT="$WORK/named_posterior"
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/named_posterior_threshold.R" \
        --typing-rds "$WORK/anchor_typing.rds" \
        --focus "$FOCUS" \
        --thresholds "${THRESHOLDS:-0.8,0.95,0.99}" \
        --margins "${MARGINS:-2,5,10,25}" \
        --output-dir "$OUT"

echo "Uploading tables to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/named_posterior"
for f in "$OUT"/*.csv "$OUT"/*.csv.gz; do
    [[ -e "$f" ]] || continue
    s5cmd cp "$f" "${DEST}/$(basename "$f")"
done

echo "Done. Read posterior_formula_check first: if formula_validated is FALSE, the posteriors"
echo "are indicative only. Then the ${FOCUS} grid and its READ line in the log above."
