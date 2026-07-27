#!/bin/bash
# Stage 5e (diagnostic, Phase 1): characterise the Low_signal population beyond the CNV
# malignant/normal split — reference/negative-control integrity (Stage 0), edge dilution
# (Stage 1), the area x count doublet screen (Stage 2), and the spatial margin map. Rebuilds
# the per-cell malignant-signature from the per-section CNV outputs and joins per-cell Area /
# depth / typing prob from cosmx_typed.h5ad. Writes a reusable per-cell master table
# (cell_cnv_table.csv.gz) for later phases. CPU only; NO new CNV/typing run required.
#
# Submit AFTER the 5b array (persection/ present on Kopah):
#   sbatch pipeline/slurm/91_insitucnv_lowsignal_diagnostics.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUCNV from pipeline/.env):
#   STAGE5_DIR           Kopah sub-dir for stage-5 (default stage5_insitucnv)
#   STAGE4_DIR           Kopah sub-dir with the typed AnnData (default stage4_insitutree)
#   TYPED_BASENAME       (default cosmx_typed.h5ad)
#   SENSITIVE_THRESHOLD  2nd, more sensitive malignant cutoff to report (default 0.45; the
#                        validated bimodal trough); set to "none" to skip.
#   SIG_THRESHOLD        override the strict cutoff (default: 95th pct of neg controls).
#   MAP_SECTIONS         sections to render in the spatial map (default 6).

#SBATCH --job-name=cosmx-insitucnv-diag
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=pipeline/logs/insitucnv_diag_%j.out
#SBATCH --error=pipeline/logs/insitucnv_diag_%j.err

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

STAGE5="${STAGE5_DIR:-stage5_insitucnv}"
STAGE4="${STAGE4_DIR:-stage4_insitutree}"
TYPED_BASENAME="${TYPED_BASENAME:-cosmx_typed.h5ad}"
SENSITIVE_THRESHOLD="${SENSITIVE_THRESHOLD:-0.45}"
MAP_SECTIONS="${MAP_SECTIONS:-6}"
: "${APPTAINER_INSITUCNV:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_insitucnv_diag_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/persection" "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-section CNV results + typed AnnData from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/persection/*" "$WORK/persection/"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

# optional args (kept out of the array so an unset override doesn't pass an empty flag)
EXTRA_ARGS=()
[[ -n "${SIG_THRESHOLD:-}" ]] && EXTRA_ARGS+=(--sig-threshold "$SIG_THRESHOLD")
[[ "${SENSITIVE_THRESHOLD}" != "none" ]] && EXTRA_ARGS+=(--sensitive-threshold "$SENSITIVE_THRESHOLD")

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUCNV" \
    python -u "${PIPELINE_DIR}/python/insitucnv_lowsignal_diagnostics.py" \
        --cnv-dir "$WORK/persection" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --reference-file "${PIPELINE_DIR}/reference/insitucnv_reference_types.txt" \
        --map-sections "$MAP_SECTIONS" \
        ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
        --output-dir "$WORK/out"

echo "Uploading Phase-1 diagnostics to Kopah (${STAGE5}/diagnostics)..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/diagnostics/"

echo "Done: stage 5e Low_signal Phase-1 diagnostics."
