#!/bin/bash
# Per-cell table with the cell type from every InSituType run (export_cell_labels.py): slide,
# donor, diagnosis, Leiden, UMAP, depth, and celltype_<run> + prob_<run> per run. Built so the
# typing can be used offline without the multi-GB AnnData. CPU only, reads obs + obsm only.
#
# Submit (RUNS = space-separated name:kopah-subdir[:result-file] entries; the file defaults to
# insitutype_result.h5, so an InSituTree run is  tree:stage4_insitutree_k21:insitutree_result.h5):
#   RUNS="k10:stage4_k10 k21:stage4_k21 k30:stage4_k30" sbatch pipeline/slurm/54_export_cell_labels.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC, from pipeline/.env):
#   RUNS          required; each subdir must hold insitutype_result.h5
#   STAGE3_DIR    Kopah sub-dir holding cosmx_clustered.h5ad (default stage3)
#   ANNOTATIONS   repo path of the FOV table (default pipeline/reference/fov_annotations_sorl1.csv)
#   OUT_SUBDIR    Kopah sub-dir under STAGE3_DIR for the output (default labels)

#SBATCH --job-name=cosmx-export-labels
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=pipeline/logs/export_labels_%j.out
#SBATCH --error=pipeline/logs/export_labels_%j.err

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

: "${RUNS:?set RUNS='name:subdir name:subdir ...'}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"
STAGE3="${STAGE3_DIR:-stage3}"
OUT_SUBDIR="${OUT_SUBDIR:-labels}"
ANNOTATIONS="${ANNOTATIONS:-pipeline/reference/fov_annotations_sorl1.csv}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_export_labels_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging cosmx_clustered.h5ad from Kopah (${STAGE3})..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/cosmx_clustered.h5ad" "$WORK/clustered.h5ad"

RUN_ARGS=()
for spec in $RUNS; do
    IFS=: read -r name subdir resultfile <<< "$spec"
    resultfile="${resultfile:-insitutype_result.h5}"
    echo "Staging ${subdir}/${resultfile} as '${name}'..."
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${subdir}/${resultfile}" "$WORK/${name}.h5"
    RUN_ARGS+=(--run "${name}=$WORK/${name}.h5")
done

apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/export_cell_labels.py" \
        --clustered-h5ad "$WORK/clustered.h5ad" \
        --annotations "${SLURM_SUBMIT_DIR:-$PWD}/${ANNOTATIONS}" \
        --output "$WORK/out/cell_labels.csv.gz" \
        "${RUN_ARGS[@]}"

echo "Uploading cell_labels.csv.gz to Kopah..."
s5cmd cp "$WORK/out/cell_labels.csv.gz" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${OUT_SUBDIR}/cell_labels.csv.gz"
echo "Done: cell labels."
