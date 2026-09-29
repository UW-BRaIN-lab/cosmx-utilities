#!/bin/bash
# PI request: heatmaps comparing a de-novo letter's cells — split by the GBmap class the forced
# supervised run puts them in — against the cells that got that same GBmap class NATIVELY.
#
# "For the cells that were different enough to be put into a de-novo cell type when they were
# allowed to be, how do they compare to the cells that received a GBmap call without needing a
# de-novo bin?" Two comparisons, both in one figure per letter:
#   down each pair    t->OPC-like  vs  OPC-like [native]     (WHICH genes differ)
#   across the row    t->OPC-like  vs  t->AC-like  vs ...    (is the forced split real?)
# See the header of python/denovo_vs_native_pseudobulk.py for why the second axis is usually the
# more informative one, and why "they differ" is definitional rather than a finding.
#
# Computes only; RENDER ON THE MAC — the InSituType/RSC containers have no ComplexHeatmap
# (see the heatmap render gotcha in the reference rebuild notes). This job emits the three small
# CSVs marker_heatmap.R reads; scp them over and run it there.
#
# Submit (defaults to the two comparisons the PI asked for):
#   sbatch pipeline/slurm/75d_denovo_vs_native_markers.sh
# Sweep several letters, taking each one's own largest destinations automatically — the right
# mode for the letters that scatter across many named types (b 11 leaves, u 9, d 8, e 7, t 6):
#   LETTERS=b,d,e,j,z sbatch pipeline/slurm/75d_denovo_vs_native_markers.sh
# One ad-hoc comparison with destinations named by hand instead:
#   LETTER=z DESTINATIONS=MES-like_hypoxia_MHC,Astrocyte,AC-like \
#     sbatch pipeline/slurm/75d_denovo_vs_native_markers.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR   Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   INPUT_DIR    Kopah sub-dir holding anchor/anchor_input.h5 (default stage4_anchor).
#   LETTERS      comma-separated letters to sweep; destinations come from the cross-tab.
#   TOP_DEST     with LETTERS, destinations per letter (default 3).
#   TOP_N        markers selected per group (default 8).
#   MIN_GROUP_N  drop groups smaller than this (default 50).

#SBATCH --job-name=cosmx-denovo-vs-native
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# Reads the 9.6GB anchor counts, but subsets to the compared cells (~200k) on the CSC before
# transposing, so the working set is far smaller than the full 2.5M-cell matrix.
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=pipeline/logs/denovo_vs_native_%j.out
#SBATCH --error=pipeline/logs/denovo_vs_native_%j.err

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
TOP_N="${TOP_N:-8}"
MIN_GROUP_N="${MIN_GROUP_N:-50}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

# Default comparisons, from the corrected Sankey. `t` also includes OPC because the PI named it,
# though the corrected run sends only 0.7% of t there (the old, buggy Sankey showed 9.4%).
COMPARISONS=(
    "t:OPC-like,AC-like,MES-like_hypoxia_MHC,OPC"
    "l:Endo_capilar,Endo_arterial,Pericyte"
)
if [[ -n "${LETTER:-}" ]]; then
    : "${DESTINATIONS:?set DESTINATIONS alongside LETTER}"
    COMPARISONS=("${LETTER}:${DESTINATIONS}")
elif [[ -n "${LETTERS:-}" ]]; then
    # Empty destination list => the script derives them from the cross-tab.
    COMPARISONS=()
    IFS=',' read -ra _letters <<< "$LETTERS"
    for _l in "${_letters[@]}"; do COMPARISONS+=("${_l// /}:"); done
fi
TOP_DEST="${TOP_DEST:-3}"

# State the resolved plan up front. An older checkout silently ignores LETTERS and runs the
# defaults instead, which is invisible until you go looking for outputs that were never made.
echo "Comparisons to run (${#COMPARISONS[@]}): ${COMPARISONS[*]}"
if [[ -n "${LETTERS:-}" ]]; then
    echo "  mode: LETTERS sweep, destinations from the cross-tab (TOP_DEST=${TOP_DEST})"
elif [[ -n "${LETTER:-}" ]]; then
    echo "  mode: single ad-hoc comparison"
else
    echo "  mode: built-in defaults (set LETTERS=... to sweep other letters)"
fi

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_denovo_vs_native_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging anchor counts + labels from Kopah..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" "$WORK/forced.csv"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_gbmap_crosstab.csv" "$WORK/crosstab.csv"

for spec in "${COMPARISONS[@]}"; do
    letter="${spec%%:*}"
    dests="${spec#*:}"
    outdir="$WORK/out/${letter}_vs_native"
    echo
    if [[ -n "$dests" ]]; then
        DEST_ARG=(--destinations "$dests")
        echo "=== ${letter} vs native [${dests}] ==="
    else
        DEST_ARG=(--crosstab "$WORK/crosstab.csv" --top-destinations "$TOP_DEST")
        echo "=== ${letter} vs native [top ${TOP_DEST} destinations from the cross-tab] ==="
    fi
    apptainer exec \
        --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
        --bind "${WORK}:${WORK}" \
        "$APPTAINER_RSC" \
        python "${PIPELINE_DIR}/python/denovo_vs_native_pseudobulk.py" \
            --counts-h5 "$WORK/anchor_input.h5" \
            --typing-h5 "$WORK/anchor_typing.h5" \
            --forced-csv "$WORK/forced.csv" \
            --letter "$letter" \
            "${DEST_ARG[@]}" \
            --top-n "$TOP_N" \
            --min-group-n "$MIN_GROUP_N" \
            --output-dir "$outdir"

    echo "Uploading ${letter} marker inputs to Kopah..."
    s5cmd cp "$outdir/*" \
        "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_native/${letter}/"
done

echo
echo "Done. Fetch to the Mac and render there (no ComplexHeatmap in the container)."
echo "NOTE: every letter writes the SAME three filenames, so they must stay in per-letter"
echo "subdirectories — do NOT pass --flatten here or one letter overwrites the other."
echo "  s5cmd cp '${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_native/*' ~/denovo_vs_native/"
echo "  # then, from a Mac terminal:"
echo "  scp -r emilyek@klone.hyak.uw.edu:denovo_vs_native ~/keene-lab/cosmx-utilities/"
echo "  Rscript pipeline/R/marker_heatmap.R denovo_vs_native/t t_heatmap \"\" \"t forced vs native\""
