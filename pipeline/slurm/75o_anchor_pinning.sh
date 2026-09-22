#!/bin/bash
# Are 75c's label-vs-argmax mismatches just InSituType's auto-selected ANCHOR cells?
#
# insitutype() is called here with anchors=NULL, but update_reference_profiles=TRUE makes
# updateReferenceProfiles() DERIVE them, and the last thing insitutype() does is
#     out$clust[!is.na(anchors)] <- anchors[...]
# a hard overwrite that ignores the likelihoods. That predicts every feature of the mismatch
# set: named-labelled only, zero de-novo, arbitrarily large deficits, 99.5% losing to a letter.
#
# The fit saved $anchors, so this is a join — no counts, minutes. It decides whether the
# "native groups were not assigned by the likelihoods" reading stands or inverts.
#
# Submit:
#   sbatch pipeline/slurm/75o_anchor_pinning.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/anchor_typing.rds (default stage4_anchor_pruned).

#SBATCH --job-name=cosmx-anchor-pinning
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/anchor_pinning_%j.out
#SBATCH --error=pipeline/logs/anchor_pinning_%j.err

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
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_anchor_pinning_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}"

# The .rds, not the h5 — the compact h5 drops $logliks and $anchors.
echo "Staging the full insitutype result from Kopah..."
s5cmd cp "${BASE}/anchor/anchor_typing.rds" "$WORK/anchor_typing.rds"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/check_anchor_pinning.R" \
        "$WORK/anchor_typing.rds" "$WORK/anchor_pinning_by_label.csv"

echo "Uploading..."
s5cmd cp "$WORK/anchor_pinning_by_label.csv" \
    "${BASE}/supervised_gbmap/anchor_pinning_by_label.csv"

echo "Done. The verdict is the anchor x mismatch table printed above."
