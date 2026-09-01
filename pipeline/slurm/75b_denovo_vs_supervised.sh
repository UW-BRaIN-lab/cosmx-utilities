#!/bin/bash
# De-novo letters vs their forced GBmap names: cross-tab + Sankey (consumes 75).
#
# Joins the semi-supervised anchor calls (72's anchor_typing.h5: 27 de-novo letters + named
# GBmap types) against the fully supervised re-score of the same cells (75's
# supervised_gbmap_posteriors.csv), then draws the Sankey the PI asked for. The named rows of
# the all-labels cross-tab are the control: cells the semi-supervised fit already named should
# mostly keep that name under the forced run.
#
# Draws TWO Sankeys. The compartment one is the figure to read: GBmap's fine leaves are
# mutually collinear on this panel, so the named control re-calls only ~79% of cells onto the
# same leaf but ~97% onto the same compartment. The leaf Sankey is kept for the detail behind
# it. COLLAPSE_MAP defines the roll-up and is a judgement call — edit that CSV and re-run this
# job (the 52-minute re-score in 75 does NOT need to repeat) to redraw under a different one.
#
# Cheap and fast — re-run it freely to re-cut the figures (e.g. MIN_PROB for a
# confidence-gated view). Everything runs in the RSC container (pandas + matplotlib).
#
# Submit (after 75):
#   sbatch --dependency=afterok:<75 jobid> pipeline/slurm/75b_denovo_vs_supervised.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR    Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   ANNOTATIONS   repo path to the de-novo annotation CSV, used for readable Sankey labels
#                 (default reference/denovo_annotations/fullcohort_pruned_k27.csv).
#   MIN_PROB      drop forced calls below this top1 posterior before cross-tabbing (default 0
#                 = keep every cell; insitutypeML always assigns one).
#   COLLAPSE_MAP  repo path to the GBmap-type -> compartment map (default
#                 reference/gbmap_compartments.csv). Set to the empty string to skip the
#                 compartment roll-up and draw only the leaf-level Sankey.
#   MIN_FRAC      Sankey ribbon draw threshold as a fraction of all cells (default 0.0015).

#SBATCH --job-name=cosmx-denovo-vs-gbmap
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/denovo_vs_gbmap_%j.out
#SBATCH --error=pipeline/logs/denovo_vs_gbmap_%j.err

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
ANNOTATIONS="${ANNOTATIONS:-reference/denovo_annotations/fullcohort_pruned_k27.csv}"
MIN_PROB="${MIN_PROB:-0}"
MIN_FRAC="${MIN_FRAC:-0.0015}"
COLLAPSE_MAP="${COLLAPSE_MAP-reference/gbmap_compartments.csv}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_vs_gbmap_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}"

echo "Staging anchor typing + supervised posteriors from Kopah..."
s5cmd cp "${BASE}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/supervised_gbmap/supervised_gbmap_posteriors.csv" "$WORK/posteriors.csv"

COLLAPSE_ARG=()
if [[ -n "$COLLAPSE_MAP" ]]; then
    COLLAPSE_ARG=(--collapse-map "${PIPELINE_DIR}/${COLLAPSE_MAP}")
fi

echo "Cross-tabbing de-novo letters against their forced GBmap calls..."
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/crosstab_denovo_vs_supervised.py" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --posteriors "$WORK/posteriors.csv" \
        --annotations "${PIPELINE_DIR}/${ANNOTATIONS}" \
        --min-prob "$MIN_PROB" \
        ${COLLAPSE_ARG[@]+"${COLLAPSE_ARG[@]}"} \
        --output-dir "$WORK/out"

draw_sankey() {
    apptainer exec \
        --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
        --bind "${WORK}:${WORK}" \
        "$APPTAINER_RSC" \
        python "${PIPELINE_DIR}/python/plot_crosstab_sankey.py" \
            --crosstab "$WORK/out/$1" \
            --label-left "de novo (k=27, semi-supervised)" \
            --label-right "$2" \
            --min-frac "$MIN_FRAC" \
            --output "$WORK/out/$3"
}

if [[ -n "$COLLAPSE_MAP" ]]; then
    echo "Drawing the compartment Sankey (the one to read)..."
    draw_sankey denovo_vs_compartment_crosstab.csv \
        "forced GBmap compartment" denovo_vs_compartment_sankey.png
fi

echo "Drawing the leaf-level Sankey (detail)..."
draw_sankey denovo_vs_gbmap_crosstab.csv \
    "forced GBmap call (supervised)" denovo_vs_gbmap_sankey.png

echo "Uploading cross-tabs + Sankey to Kopah (${STAGE4}/supervised_gbmap)..."
s5cmd cp "$WORK/out/*" "${BASE}/supervised_gbmap/"

echo "Done. Fetch the small outputs to the Mac with:"
echo "  s5cmd cp '${BASE}/supervised_gbmap/*.csv' ."
echo "  s5cmd cp '${BASE}/supervised_gbmap/*.png' ."
