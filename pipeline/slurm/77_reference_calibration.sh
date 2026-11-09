#!/bin/bash
# Reference triage 1 of 3: which reference cell types did InSituType actually CALIBRATE?
#
# With update_reference_profiles=TRUE, updateReferenceProfiles() pins the cells that best
# match each reference profile and rescales that profile toward them. A type that attracts no
# cell past those bars is never calibrated; its fitted profile can collapse to a few expected
# counts, after which it cannot win a cell and so appears nowhere in the output. This is the
# only place such a type is visible. See reference_calibration.R.
#
# Reads the FULL .rds that 80 saves (clust, profiles, anchors) -- a join, no counts, minutes.
#
# Submit (after a fit exists; run 1 lives in stage4/):
#   sbatch pipeline/slurm/77_reference_calibration.sh
#
# Env knobs (KOPAH_*, APPTAINER_INSITUTYPE from pipeline/.env):
#   STAGE4_DIR          Kopah sub-dir holding the fit (default stage4)
#   RESULT_RDS_BASENAME the fit to read (default insitutype_result.rds)
#   REFERENCE_BASENAME  the reference the fit used, under ${KOPAH_PREFIX}/reference/
#                       -- the retina checkout's .env sets retina_combined_panel.csv
#   OUT_SUBDIR          Kopah sub-dir under STAGE4_DIR (default reference_triage)

#SBATCH --job-name=cosmx-ref-calibration
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/reference_calibration_%j.out
#SBATCH --error=pipeline/logs/reference_calibration_%j.err

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

# Loads .env WITHOUT letting the ambient shell override the study's identity; see
# slurm/_load_env.sh.
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/slurm/_load_env.sh"

STAGE4="${STAGE4_DIR:-stage4}"
RESULT_RDS_BASENAME="${RESULT_RDS_BASENAME:-insitutype_result.rds}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"
OUT_SUBDIR="${OUT_SUBDIR:-reference_triage}"

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_ref_calibration_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${RESULT_RDS_BASENAME} (${STAGE4}) and reference ${REFERENCE_BASENAME}..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${RESULT_RDS_BASENAME}" \
    "$WORK/fit.rds"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/reference_calibration.R" \
        --result-rds "$WORK/fit.rds" \
        --reference "$WORK/reference.csv" \
        --output-csv "$WORK/reference_calibration.csv"

echo "Uploading the calibration table to Kopah..."
s5cmd cp "$WORK/reference_calibration.csv" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${OUT_SUBDIR}/reference_calibration.csv"

echo "Done: reference calibration. The table is printed above, so the log alone answers it."
