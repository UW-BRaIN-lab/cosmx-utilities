#!/bin/bash
# Two-pass PASS 1 (profile): build the fixed cohort profile (genes x cell types, incl.
# de novo) from the anchor's raw counts + labels via InSituType::getRNAprofiles. CPU,
# quick. Output cohort_profile.csv is what every slide is typed against in PASS 2.
#
# Submit:  sbatch pipeline/slurm/73_build_profile.sh
# Chain:   jid=$(sbatch --parsable pipeline/slurm/73_build_profile.sh)
#          sbatch --dependency=afterok:$jid pipeline/slurm/75_insitutype_supervised.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE.

#SBATCH --job-name=cosmx-build-profile
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/build_profile_%j.out
#SBATCH --error=pipeline/logs/build_profile_%j.err

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

STAGE4="${STAGE4_DIR:-stage4}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_build_profile_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging anchor input + anchor typing from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_input.h5" \
    "$WORK/anchor_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_typing.h5" \
    "$WORK/anchor_typing.h5"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/insitutype_profile.R" \
        --anchor-input "$WORK/anchor_input.h5" \
        --anchor-typing "$WORK/anchor_typing.h5" \
        --output-csv "$WORK/cohort_profile.csv"

echo "Uploading cohort profile to Kopah..."
s5cmd cp "$WORK/cohort_profile.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/cohort_profile.csv"

echo "Done: PASS 1 cohort profile."
