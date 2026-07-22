#!/bin/bash
# Stage 5a: build per-tissue-section InSituCNV inputs from the typed cohort AnnData.
# Reads stage4_insitutree/cosmx_typed.h5ad, adds obsm['spatial'] from the CosMx centroids,
# annotates gene positions (autosomes) from the BioMart table baked into insitucnv.sif,
# and splits into one gene-only raw-count h5ad per (slide, donor) tissue section. CPU only.
#
# Submit:
#   sbatch pipeline/slurm/86_insitucnv_prep.sh
# Then read the section count it prints (sections.txt) to size the array in stage 5b:
#   s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage5_insitucnv/sections.txt .
#   sbatch --array=0-$(($(wc -l < sections.txt)-1)) pipeline/slurm/87_insitucnv.sh
#
# Prerequisite: build insitucnv.sif once (see pipeline/containers/insitucnv.def), point
# APPTAINER_INSITUCNV in pipeline/.env at it.
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE4_DIR      Kopah sub-dir with the typed AnnData (default stage4_insitutree)
#   TYPED_BASENAME  (default cosmx_typed.h5ad)
#   STAGE5_DIR      Kopah sub-dir for stage-5 outputs (default stage5_insitucnv)
#   MIN_CELLS       skip sections smaller than this (default 500)
#   MIN_REF_CELLS   skip sections with fewer diploid-reference cells (default 50)

#SBATCH --job-name=cosmx-insitucnv-prep
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=pipeline/logs/insitucnv_prep_%j.out
#SBATCH --error=pipeline/logs/insitucnv_prep_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
STAGE5="${STAGE5_DIR:-stage5_insitucnv}"

: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_prep_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${TYPED_BASENAME} from Kopah (${STAGE4})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/prep_insitucnv_input.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --reference-file "${PIPELINE_DIR}/reference/insitucnv_reference_types.txt" \
        --output-dir "$WORK/out" \
        --min-cells "${MIN_CELLS:-500}" \
        --min-reference-cells "${MIN_REF_CELLS:-50}"

echo "Uploading section inputs + manifest to Kopah (${STAGE5})..."
s5cmd cp "$WORK/out/sections/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/sections/"
s5cmd cp "$WORK/out/sections_manifest.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/sections_manifest.csv"
s5cmd cp "$WORK/out/sections.txt" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/sections.txt"

echo "Done: stage 5a prep. Section count = $(wc -l < "$WORK/out/sections.txt")."
