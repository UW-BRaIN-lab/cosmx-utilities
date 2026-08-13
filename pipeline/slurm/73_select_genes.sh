#!/bin/bash
# InSituType 6k-panel gene selection (FAQ union) — the CHEAP tuning step before the anchor
# re-fit. Computes per-gene means over the stage-4a anchor + the reference, retains genes
# that are above background OR moderate-to-high in >=1 profile, and writes kept_genes.txt.
# Sweep BG_MULT / REF_QUANTILE until the union lands in the 3000-5000 band (see summary).
#
# Submit (reads the same anchor_input.h5 as 72):  sbatch pipeline/slurm/73_select_genes.sh
# Then feed kept_genes.txt to 74_choose_k.sh (K-sweep) and the anchor re-fit.
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE.

#SBATCH --job-name=cosmx-gene-select
# Pure column-means over the anchor — CPU-only, small and fast. ckpt is fine (short job).
#SBATCH --account=glioblastoma
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/gene_select_%j.out
#SBATCH --error=pipeline/logs/gene_select_%j.err

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
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_gene_select_${SLURM_JOB_ID:-local}"
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

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/select_informative_genes.R" \
        --input "$WORK/anchor_input.h5" \
        --reference "$WORK/reference.csv" \
        --out-dir "$WORK/gene_selection" \
        --bg-mult "${BG_MULT:-3}" \
        --ref-quantile "${REF_QUANTILE:-0.5}" \
        --target-min "${TARGET_MIN:-3000}" \
        --target-max "${TARGET_MAX:-5000}" \
        --max-cells "${GENE_SELECT_MAX_CELLS:-500000}"

echo "Uploading gene-selection outputs to Kopah (${STAGE4}/gene_selection)..."
for f in kept_genes.txt excluded_genes.csv gene_selection.csv gene_selection_summary.txt; do
    s5cmd cp "$WORK/gene_selection/$f" \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/gene_selection/$f"
done

echo "=== gene selection summary ==="
cat "$WORK/gene_selection/gene_selection_summary.txt"
echo "Done: gene selection."
