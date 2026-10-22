#!/bin/bash
# SEA-AD MTG "Whole Taxonomy" (1,378,211 cells, ~51 GB H5AD) -> InSituType reference profiles.
# Downloads the H5AD from CELLxGENE (resumable), then averages raw counts per cell type.
# CELLxGENE collection 1ca90a2d-2943-483d-b678-b809bf464c30, dataset c2876b1b-06d8-4d96-a56b-5304f815b99a.
#
# Submit:  sbatch pipeline/slurm/76_seaad_mtg_reference.sh
# Override (env at submit time):
#   SEAAD_DIR     work dir (default /gscratch/scrubbed/$USER/seaad_mtg; scrubbed = transient, 21-day purge)
#   LABEL_COL     obs column with the cell-type labels (default Supertype; Subclass = coarser)
#   OBS_FILTER    e.g. "disease=normal" to build from healthy donors only (default: all cells)
# Output: ${SEAAD_DIR}/SEAAD_MTG_${LABEL_COL}_pseudobulk_profiles.csv  (genes x types, mean raw counts)
# Then: pipeline/reference_builder/examples/seaad_mtg.R rescales it into the final InSituType profile.

#SBATCH --job-name=cosmx-seaad-ref
# One sequential pass over 51 GB; I/O-bound, memory-light (sparse 50k-cell chunks). ckpt is fine
# because the download resumes and the pass restarts cleanly.
#SBATCH --account=glioblastoma
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=pipeline/logs/seaad_ref_%j.out
#SBATCH --error=pipeline/logs/seaad_ref_%j.err

set -euo pipefail

SEAAD_H5AD_URL="https://datasets.cellxgene.cziscience.com/732b4aa7-e441-4196-b05e-127433425613.h5ad"
SEAAD_H5AD_BYTES=50775398656
SEAAD_DIR="${SEAAD_DIR:-/gscratch/scrubbed/${USER}/seaad_mtg}"
LABEL_COL="${LABEL_COL:-Supertype}"
OBS_FILTER="${OBS_FILTER:-}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_DIR="${SLURM_SUBMIT_DIR}"
else
    REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
fi

mkdir -p "${SEAAD_DIR}" "${REPO_DIR}/pipeline/logs"
H5AD="${SEAAD_DIR}/SEAAD_MTG_whole_taxonomy.h5ad"
OUTPUT="${SEAAD_DIR}/SEAAD_MTG_${LABEL_COL}_pseudobulk_profiles.csv"

# -C - resumes a partial download after a preemption.
if [[ "$(stat -c %s "${H5AD}" 2>/dev/null || echo 0)" -ne "${SEAAD_H5AD_BYTES}" ]]; then
    curl --fail --location -C - --retry 5 --output "${H5AD}" "${SEAAD_H5AD_URL}"
fi
DOWNLOADED_BYTES="$(stat -c %s "${H5AD}")"
if [[ "${DOWNLOADED_BYTES}" -ne "${SEAAD_H5AD_BYTES}" ]]; then
    echo "ERROR: ${H5AD} is ${DOWNLOADED_BYTES} bytes, expected ${SEAAD_H5AD_BYTES}" >&2
    exit 1
fi

FILTER_ARGS=()
if [[ -n "${OBS_FILTER}" ]]; then
    FILTER_ARGS=(--obs-filter "${OBS_FILTER}")
fi

cd "${REPO_DIR}"
uv run --no-project --with anndata --with scipy --with pandas --with numpy \
    python -I pipeline/reference_builder/preprocess_h5ad_profiles.py \
    --h5ad "${H5AD}" --cell-type-col "${LABEL_COL}" \
    ${FILTER_ARGS[@]+"${FILTER_ARGS[@]}"} \
    --prefix SEAAD_ --output "${OUTPUT}"
echo "Wrote ${OUTPUT}"
