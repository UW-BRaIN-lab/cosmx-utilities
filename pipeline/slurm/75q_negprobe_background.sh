#!/bin/bash
# Is letter `b`'s "carries instead" signature a cell state, or probe background?
#
# 75g RETRACTED the "Low_signal = expression spread thin" framing: at matched depth `b`'s cells
# are 1.44-1.72x MORE concentrated in their top genes than the types they are forced into. So
# they have a profile. The open question is what that profile IS, because the genes carrying it
# (SNRPA1, A1BG, SNRPD1, VPS29, NDUFA3, PTEN, GPX1, F2RL3, B9D2) recur in `l`'s and `e`'s lists
# wherever real signal is weak, while `b` is MISSING the panel's most abundant transcripts.
#
# `/neg` in anchor_input.h5 is the per-cell mean negprobe count, so a gene whose signal is probe
# background should track it across cells once depth is held fixed. See the header of
# python/negprobe_background_test.py for the two hypotheses and what separates them.
#
# Submit (defaults target `b` on the pruned k=27 run):
#   sbatch pipeline/slurm/75q_negprobe_background.sh
# Another letter, or another signature:
#   LETTER=l sbatch pipeline/slurm/75q_negprobe_background.sh
#   LETTER=e SIGNATURE=SNRPA1,A1BG,GPX1 sbatch pipeline/slurm/75q_negprobe_background.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/anchor_typing.h5 (default stage4_anchor_pruned).
#   INPUT_DIR   Kopah sub-dir with anchor/anchor_input.h5 (default stage4_anchor).
#   LETTER      cell type the per-gene affinity is computed within (default b).
#   SIGNATURE   comma-separated set under test (default: b's recurring "carries instead" genes).
#   CONTRAST    comma-separated ambient-RNA check (default: the abundant transcripts b misses).
#   STRATA      depth quantile bins (default 10). Depth drives expression AND negprobe counts, so
#               an unstratified correlation calls the whole panel background — see the test
#               test_depth_coupling_alone_is_NOT_called_background.
#   NULL_SETS   matched-null sets (default 200). Affinity is computed once for every gene, so the
#               null is nearly free here.

#SBATCH --job-name=cosmx-negprobe-bg
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# Holds the full CSC plus one depth stratum's nonzeros at a time; `b` alone is ~778k cells.
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=pipeline/logs/negprobe_bg_%j.out
#SBATCH --error=pipeline/logs/negprobe_bg_%j.err

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
LETTER="${LETTER:-b}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

echo "Letter under test: ${LETTER}   strata: ${STRATA:-10}   null sets: ${NULL_SETS:-200}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_negprobe_bg_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging anchor counts (/neg lives here) + labels from Kopah..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"

SIG_ARG=()
if [[ -n "${SIGNATURE:-}" ]]; then SIG_ARG=(--signature-genes "$SIGNATURE"); fi
CONTRAST_ARG=()
if [[ -n "${CONTRAST:-}" ]]; then CONTRAST_ARG=(--contrast-genes "$CONTRAST"); fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/negprobe_background_test.py" \
        --counts-h5 "$WORK/anchor_input.h5" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --letter "$LETTER" \
        --depth-strata "${STRATA:-10}" \
        --null-sets "${NULL_SETS:-200}" \
        ${SIG_ARG[@]+"${SIG_ARG[@]}"} \
        ${CONTRAST_ARG[@]+"${CONTRAST_ARG[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading the tables to Kopah (${STAGE4}/supervised_gbmap/negprobe_bg/${LETTER})..."
s5cmd cp "$WORK/out/*" \
    "${BASE}/${STAGE4}/supervised_gbmap/negprobe_bg/${LETTER}/"

echo "Done. The verdict is in the three blocks printed above."
