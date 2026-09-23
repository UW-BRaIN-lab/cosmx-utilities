#!/bin/bash
# 75k: do a de-novo letter's cells sit where their assigned cell type lives?
#
# The one diagnostic in this set that reads nothing the typing produced. For each cell it
# measures the fraction of its k nearest SPATIAL neighbours that are mural (pericyte / SMC /
# perivascular fibroblast), because vessels are multicellular tubes: real endothelium is wrapped
# by mural cells, and a tumour cell shoehorned into an endothelial profile is not. Mural rather
# than endothelial is deliberate -- asking whether cells called endothelial neighbour cells
# called endothelial would just read the typing back out.
#
# Neighbourhoods come from the FIXED-PROFILE run's 7.5M cells, not the 2.54M anchor: the anchor
# is a Leiden-stratified subsample with a cap per stratum, so its local density is a function of
# how abundant each cluster was. Each group is then compared with size-matched random cells from
# the SAME FOVs, which holds tissue density and FOV composition fixed.
#
# Labels + obs only -- the expression matrix is never read -- so this is minutes.
#
# Submit (l against its own three largest destinations):
#   sbatch pipeline/slurm/75k_denovo_origin_spatial.sh
# Another letter, destinations named by hand:
#   LETTER=c DESTINATIONS=Pericyte,SMC,SMC_COL \
#     sbatch pipeline/slurm/75k_denovo_origin_spatial.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR    Kopah sub-dir with anchor/ + supervised_gbmap/ (default stage4_anchor_pruned).
#   TYPED_DIR     Kopah sub-dir with the fixed-profile typing (default stage4_insitutree).
#   TYPED_BASENAME  default cosmx_typed.h5ad.
#   LETTER        de-novo letter under test (default l).
#   DESTINATIONS  comma-separated GBmap types; default is the letter's own largest.
#   MURAL_TYPES   types defining the vessel wall (default pericyte/SMC/perivascular
#                 fibroblast). NOT InSituType's `anchors`, and not our anchor cohort --
#                 see the terminology note in CLAUDE.md.
#   K             spatial neighbours per cell (default 15).
#   N_PERM        size-matched draws for the null (default 200).
#   N_FOVS        example FOV panels to draw (default 6). They are chosen for LEGIBILITY --
#                 a FOV must hold MIN_GROUP_CELLS cells of BOTH groups (default 10) and its mural
#                 share must fall inside MURAL_QUANTILES (default 0.25,0.75) of the eligible
#                 FOVs. Ranking on forced-cell count instead, as this used to, picks the most
#                 vessel-dense FOVs in the cohort -- mural wall-to-wall, no visible
#                 architecture, and often only a handful of native cells to compare against.
#   MIN_GROUP_CELLS / MURAL_QUANTILES  the two knobs above.
#   EXCLUDE_COMPARED  set to 1 to hold the compared cells out of the neighbourhood reference.
#                 REQUIRED for a letter whose own cells carry an anchor type under the
#                 fixed-profile run (c is ~72% Pericyte), or each group partly supplies its own
#                 evidence. The job reports that overlap per group either way.

#SBATCH --job-name=cosmx-origin-spatial
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/origin_spatial_%j.out
#SBATCH --error=pipeline/logs/origin_spatial_%j.err

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

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_origin_spatial_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging labels and the fixed-profile typing from Kopah..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
# 75c's corrected forced call. NOT 75's re-scored file, which disagrees with the fit on 45%
# of cells -- see the 75c self-consistency check.
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" \
    "$WORK/forced_named_posteriors.csv"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/denovo_vs_gbmap_crosstab.csv" \
    "$WORK/denovo_vs_gbmap_crosstab.csv"
s5cmd cp "${BASE}/${TYPED}/${TYPED_BASENAME}" "$WORK/typed.h5ad"

DEST_ARG=(--crosstab "$WORK/denovo_vs_gbmap_crosstab.csv")
if [[ -n "${DESTINATIONS:-}" ]]; then DEST_ARG=(--destinations "$DESTINATIONS"); fi

OUT="$WORK/${LETTER}_spatial"
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/denovo_origin_spatial.py" \
        --typed-h5ad "$WORK/typed.h5ad" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --forced-csv "$WORK/forced_named_posteriors.csv" \
        --letter "$LETTER" \
        --mural-types "${MURAL_TYPES:-Pericyte,SMC,SMC_COL,Perivascular_fibroblast,Scavenging_pericyte}" \
        --k "${K:-15}" \
        --n-permutations "${N_PERM:-200}" \
        --n-example-fovs "${N_FOVS:-6}" \
        --min-group-cells "${MIN_GROUP_CELLS:-10}" \
        --mural-quantiles "${MURAL_QUANTILES:-0.25,0.75}" \
        ${EXCLUDE_COMPARED:+--exclude-compared} \
        "${DEST_ARG[@]}" \
        --output-dir "$OUT"

echo "Uploading to Kopah..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/origin_spatial"
for f in "$OUT"/*; do
    s5cmd cp "$f" "${DEST}/${LETTER}_$(basename "$f")"
done

echo "Done. The verdict is the observed-vs-null table printed above."
