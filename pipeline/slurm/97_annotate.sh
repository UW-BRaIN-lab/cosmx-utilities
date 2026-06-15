#!/bin/bash
# Stage 4: write a run's de-novo annotations into its typed AnnData. CPU only.
# Relabels obs.cell_type so de-novo letters read "a - MES/AC-like tumor" (named GBmap
# types unchanged) via annotate_denovo.py + the run's denovo_annotations table.
#
# Submit (RUN = the run's Kopah sub-prefix; a matching table must exist at
# pipeline/reference/denovo_annotations/<RUN>.csv):
#   RUN=stage4_extl3_rescale sbatch pipeline/slurm/97_annotate.sh
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_RSC.

#SBATCH --job-name=cosmx-annotate
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/annotate_%j.out
#SBATCH --error=pipeline/logs/annotate_%j.err

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

RUN="${RUN:?set RUN to the Kopah sub-prefix for the run, e.g. stage4_extl3_rescale}"
MAPPING="${PIPELINE_DIR}/reference/denovo_annotations/${RUN}.csv"
[[ -f "$MAPPING" ]] || { echo "ERROR: no annotation table $MAPPING" >&2; exit 1; }

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_annotate_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${RUN}/cosmx_typed.h5ad from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/cosmx_typed.h5ad" "$WORK/typed.h5ad"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/annotate_denovo.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --mapping "$MAPPING" \
        --output "$WORK/cosmx_typed_annotated.h5ad"

echo "Uploading annotated AnnData to Kopah..."
s5cmd cp "$WORK/cosmx_typed_annotated.h5ad" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/cosmx_typed_annotated.h5ad"

echo "Done: annotated AnnData for ${RUN} -> ${RUN}/cosmx_typed_annotated.h5ad"
