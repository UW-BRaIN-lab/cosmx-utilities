#!/bin/bash
# 75m: co-embed a de-novo letter's cells with the natively-called cells of its destinations.
#
# Do the forced cells intermingle with the native cluster they were assigned to, or sit as a
# satellite island beside it? A co-embedding over the compared cells alone is the right
# instrument: a cohort-wide UMAP has no resolution inside a compartment that is 4% of the cells.
#
# The figure answers it by eye; the number beside it is read off the kNN graph the embedding
# already built (PCA space, not the distorted 2-D projection) as a group x group
# observed/expected matrix. Read a forced row against ITS OWN native column, and calibrate with
# the native-vs-native entries -- those are two genuinely distinct populations measured the same
# way, so they set the scale for what "different" looks like here. Separation alone is not
# evidence of shoehorning: the assignment split these cells BECAUSE they differ.
#
# This one reads the expression matrix, so it is the only job in the set that needs a GPU -- and
# it falls back to scanpy on CPU if rapids is unavailable.
#
# Submit:
#   sbatch pipeline/slurm/75m_denovo_origin_umap.sh
#   LETTER=c DESTINATIONS=Pericyte,SMC,SMC_COL sbatch pipeline/slurm/75m_denovo_origin_umap.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR    Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   TYPED_DIR     Kopah sub-dir with the fixed-profile typing (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   LETTER        de-novo letter under test (default l).
#   DESTINATIONS  comma-separated GBmap types; default is the letter's own largest.
#   N_HVG (2000)  N_PCS (50)  N_NEIGHBORS (15)  SEED (0)

#SBATCH --job-name=cosmx-origin-umap
#SBATCH --account=glioblastoma
#SBATCH --partition=gpu-l40s
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=pipeline/logs/origin_umap_%j.out
#SBATCH --error=pipeline/logs/origin_umap_%j.err

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
TYPED="${TYPED_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
LETTER="${LETTER:-l}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_origin_umap_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging labels and the fixed-profile typing from Kopah..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" \
    "$WORK/forced_named_posteriors.csv"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_gbmap_crosstab.csv" \
    "$WORK/denovo_vs_gbmap_crosstab.csv"
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

DEST_ARG=(--crosstab "$WORK/denovo_vs_gbmap_crosstab.csv")
if [[ -n "${DESTINATIONS:-}" ]]; then DEST_ARG=(--destinations "$DESTINATIONS"); fi

OUT="$WORK/${LETTER}_umap"
apptainer exec --nv \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/denovo_origin_umap.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --forced-csv "$WORK/forced_named_posteriors.csv" \
        --letter "$LETTER" \
        --n-hvg "${N_HVG:-2000}" \
        --n-pcs "${N_PCS:-50}" \
        --n-neighbors "${N_NEIGHBORS:-15}" \
        --seed "${SEED:-0}" \
        "${DEST_ARG[@]}" \
        --output-dir "$OUT" \
        --output-h5ad "$OUT/origin_umap.h5ad"

# Every name is prefixed with the letter, as 75k and 75l do. Without it a second letter
# OVERWRITES the first's figures and enrichment matrix in place -- the same collision that
# destroyed 75e's per-FOV CSVs. Only the .h5ad escaped it, by carrying the letter in its name.
echo "Uploading to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/origin_umap"
for f in "$OUT"/*; do
    s5cmd cp "$f" "${DEST}/${LETTER}_$(basename "$f")"
done

echo "Done. The figures are umap_group.png / umap_origin.png; the number is above."
