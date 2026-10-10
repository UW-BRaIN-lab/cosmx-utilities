#!/bin/bash
# Per-donor RNA depth table (donor_depth_table.py): depth metrics and the share of each
# donor's cells in the shallow clusters, next to slide, run and diagnosis. Tests whether a
# group comparison (e.g. SORL1 carriers) is confounded with sequencing depth. CPU only,
# reads obs only.
#
# Submit:
#   sbatch pipeline/slurm/53_donor_depth.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC, MANIFEST from pipeline/.env):
#   STAGE3_DIR           Kopah sub-dir holding the clustered AnnData (default stage3)
#   ANNOTATIONS          repo path of the FOV table (default pipeline/reference/fov_annotations_sorl1.csv)
#   SHALLOW_MAX_MEDIAN   cluster median total_counts below which a cluster is "shallow" (default 560)
#   OUT_SUBDIR           Kopah sub-dir under STAGE3_DIR for the outputs (default depth)

#SBATCH --job-name=cosmx-donor-depth
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=pipeline/logs/donor_depth_%j.out
#SBATCH --error=pipeline/logs/donor_depth_%j.err

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

STAGE3="${STAGE3_DIR:-stage3}"
OUT_SUBDIR="${OUT_SUBDIR:-depth}"
ANNOTATIONS="${ANNOTATIONS:-pipeline/reference/fov_annotations_sorl1.csv}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"
: "${MANIFEST:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_donor_depth_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging cosmx_clustered.h5ad from Kopah (${STAGE3})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/cosmx_clustered.h5ad" "$WORK/clustered.h5ad"

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" \
    --bind "$(dirname "$MANIFEST"):$(dirname "$MANIFEST")" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/donor_depth_table.py" \
        --h5ad "$WORK/clustered.h5ad" \
        --annotations "${SLURM_SUBMIT_DIR:-$PWD}/${ANNOTATIONS}" \
        --manifest "$MANIFEST" \
        --shallow-max-median "${SHALLOW_MAX_MEDIAN:-560}" \
        --output "$WORK/out/donor_depth.csv"

echo "Uploading donor depth table to Kopah..."
s5cmd cp "$WORK/out/donor_depth.csv" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${OUT_SUBDIR}/donor_depth.csv"

echo "Done: donor depth table."
