#!/bin/bash
# FULLY SUPERVISED GBmap re-score of the anchor cells (companion to the semi-supervised 72).
#
# The k=27 pruned anchor fit (72) is SEMI-supervised: each cell picks either a named GBmap
# Core-L4 type or one of 27 de-novo letters. The PI's question is what those letters would be
# CALLED if the de-novo option were removed. This job answers it: score the SAME anchor cells
# against the NAMED GBmap profiles only, with InSituType::insitutypeML — supervised ML
# assignment, no de-novo clustering, no reference update. Every cell is forced onto a GBmap
# type. 75b then cross-tabs letters vs forced calls and draws the Sankey.
#
# Two steps, two containers:
#  1. (RSC) prep_supervised_profiles.py drops the 27 de-novo columns from the anchor run's
#     converged $profiles, keeping the named ones. Sourcing the profiles from the SAME
#     anchor_typing.h5 is what makes the comparison fair: those columns are already
#     CosMx-rescaled (72 ran rescale=TRUE) and already on the pruned gene panel, whereas the
#     raw scRNA-seq gbmap_level4_panel.csv is not — and insitutypeML does NOT rescale, so
#     feeding it the raw panel would fold a platform-scale correction into the answer.
#  2. (INSITUTYPE) flat_posteriors.R scores anchor_input.h5 against those profiles and emits
#     per-cell top-K posteriors. No --keep-genes is needed: flat_posteriors.R intersects the
#     panel with the profile genes, so the pruned gene space comes along with the profiles.
#
# Submit (defaults target the pruned k=27 run):
#   sbatch pipeline/slurm/75_supervised_gbmap.sh
# Then:
#   sbatch --dependency=afterok:<jobid> pipeline/slurm/75b_denovo_vs_supervised.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir holding anchor/anchor_typing.h5 (default stage4_anchor_pruned)
#                and receiving supervised_gbmap/ outputs.
#   INPUT_DIR    Kopah sub-dir holding anchor/anchor_input.h5 (default stage4_anchor — 72 wrote
#                the pruned typing to a fresh prefix but reuses the shared 9.6GB input).
#   TOP_K        top assignments per cell (default 3; top2 shows the runner-up GBmap type).
#   PROFILES_KEY Kopah key (under KOPAH_PREFIX/) of a ready-made named-profile CSV to score
#                against, INSTEAD of extracting from the anchor fit. Unset = extract (the
#                default, and the fair comparison). Set it to reference/gbmap_level4_panel.csv
#                for the raw, un-rescaled scRNA-seq GBmap variant — see the scale caveat above
#                before reading anything into that run.

#SBATCH --job-name=cosmx-supervised-gbmap
# One long insitutypeML pass over ~2.5M cells that does NOT checkpoint — a ckpt requeue would
# restart it from scratch (and re-stage the 9.6GB input). Run on the dedicated, non-preemptible
# gpu-l40s slice, as 72 does. --gres=gpu:1 claims a GPU this CPU-only fit won't use, purely to
# schedule on the dedicated partition. carrolllab/cpu-g2-mem2x is a collaborator's allocation —
# do not use it.
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
# Well under 72's 256G: supervised ML is one scoring pass (counts + a cells x 53 loglik matrix),
# not the multi-start de-novo EM. Leaves headroom under the association memory cap.
#SBATCH --mem=192G
#SBATCH --time=08:00:00
#SBATCH --output=pipeline/logs/supervised_gbmap_%j.out
#SBATCH --error=pipeline/logs/supervised_gbmap_%j.err

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
TOP_K="${TOP_K:-3}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_supervised_gbmap_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

OUT_PREFIX="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/supervised_gbmap"

echo "Staging anchor typing (${STAGE4}) + anchor input (${INPUT}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.h5" \
    "$WORK/anchor_typing.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${INPUT}/anchor/anchor_input.h5" \
    "$WORK/anchor_input.h5"

if [[ -n "${PROFILES_KEY:-}" ]]; then
    echo "Step 1/2: staging ready-made profiles from ${PROFILES_KEY} (SKIPPING extraction)..."
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${PROFILES_KEY}" \
        "$WORK/gbmap_named_profiles.csv"
else
    echo "Step 1/2: extracting the named-only GBmap profiles from the anchor fit..."
    apptainer exec \
        --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
        --bind "${WORK}:${WORK}" \
        "$APPTAINER_RSC" \
        python "${PIPELINE_DIR}/python/prep_supervised_profiles.py" \
            --profiles-h5 "$WORK/anchor_typing.h5" \
            --output "$WORK/gbmap_named_profiles.csv"
fi

echo "Step 2/2: supervised insitutypeML scoring against the named profiles..."
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/flat_posteriors.R" \
        --input "$WORK/anchor_input.h5" \
        --profiles "$WORK/gbmap_named_profiles.csv" \
        --top-k "$TOP_K" \
        --output-csv "$WORK/supervised_gbmap_posteriors.csv"

echo "Uploading supervised GBmap results to Kopah (${STAGE4}/supervised_gbmap)..."
s5cmd cp "$WORK/gbmap_named_profiles.csv" "${OUT_PREFIX}/gbmap_named_profiles.csv"
s5cmd cp "$WORK/supervised_gbmap_posteriors.csv" \
    "${OUT_PREFIX}/supervised_gbmap_posteriors.csv"

echo "Done: fully supervised GBmap re-score. Next: 75b_denovo_vs_supervised.sh"
