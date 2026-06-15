#!/bin/bash
# Stage 4 figure: InSituType flightpath plot for one run (R). CPU only.
# Loads the run's full insitutype_result.rds, subsamples cells, and renders the
# flightpath (pipeline/R/plot_flightpath.R).
#
# Submit one per run (RUN = the run's Kopah sub-prefix):
#   RUN=stage4_extl3_rescale sbatch pipeline/slurm/96_flightpath.sh
#   for r in stage4 stage4_refit stage4_ext_l3 stage4_extl3_rescale; do
#     RUN=$r sbatch pipeline/slurm/96_flightpath.sh; done
#
# Required env (from pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE.

#SBATCH --job-name=cosmx-flightpath
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/flightpath_%j.out
#SBATCH --error=pipeline/logs/flightpath_%j.err

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
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_flightpath_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${RUN}/insitutype_result.rds from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/insitutype_result.rds" "$WORK/result.rds"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/plot_flightpath.R" \
        --input "$WORK/result.rds" \
        --output "$WORK/flightpath.png"

echo "Uploading flightpath to Kopah..."
s5cmd cp "$WORK/flightpath.png" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RUN}/qc_plots/flightpath.png"

echo "Done: flightpath for ${RUN}"
