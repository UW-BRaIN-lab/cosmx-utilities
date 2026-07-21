#!/bin/bash
# Compare two TYPED runs cell-by-cell: adjusted Rand index, a per-type agreement table,
# a row-normalized cross-tab heatmap (compare_external_typing.py), and Sankey diagrams
# in both directions (plot_crosstab_sankey.py). CPU only; the rapids SIF carries
# anndata/pandas/scikit-learn/matplotlib.
#
# Both runs must share cell IDs (they do when both derive from the same
# cosmx_clustered.h5ad — same Stage 1/3 cells), so we align on the shared
# <slide>_F<fov>_C<cell> id with the same regex on both sides.
#
# Defaults compare the Extended-L3 rescale keeper (LEFT / Sankey source) against the
# InSituTree run (RIGHT / Sankey destination). Override the env knobs to diff any two
# typed runs.
#
# Submit (defaults = keeper vs InSituTree):
#   sbatch pipeline/slurm/98_compare_typing.sh
# Or point it at other runs / the keeper's annotated file:
#   LEFT_BASENAME=cosmx_typed_annotated.h5ad sbatch pipeline/slurm/98_compare_typing.sh
#
# Env knobs (from pipeline/.env for KOPAH_*, APPTAINER_RSC):
#   LEFT_DIR/LEFT_BASENAME/LEFT_KEY/LEFT_LABEL     source run (Sankey left)
#   RIGHT_DIR/RIGHT_BASENAME/RIGHT_KEY/RIGHT_LABEL  dest run   (Sankey right)
#   OUT_SUBDIR   Kopah sub-dir under RIGHT_DIR for the outputs.

#SBATCH --job-name=cosmx-compare-typing
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/compare_typing_%j.out
#SBATCH --error=pipeline/logs/compare_typing_%j.err

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

LEFT_DIR="${LEFT_DIR:-stage4_extl3_rescale}"
LEFT_BASENAME="${LEFT_BASENAME:-cosmx_typed.h5ad}"
LEFT_KEY="${LEFT_KEY:-cell_type}"
LEFT_LABEL="${LEFT_LABEL:-keeper (Ext-L3 rescale)}"

RIGHT_DIR="${RIGHT_DIR:-stage4_insitutree}"
RIGHT_BASENAME="${RIGHT_BASENAME:-cosmx_typed.h5ad}"
RIGHT_KEY="${RIGHT_KEY:-cell_type}"
RIGHT_LABEL="${RIGHT_LABEL:-InSituTree}"

OUT_SUBDIR="${OUT_SUBDIR:-comparison_vs_keeper}"
ID_RE='(.+)_F(\d+)_C(\d+)$'   # our cell-id format, shared by both runs

: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_compare_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging typed AnnDatas from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${LEFT_DIR}/${LEFT_BASENAME}"  "$WORK/left.h5ad"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RIGHT_DIR}/${RIGHT_BASENAME}" "$WORK/right.h5ad"

# compare_external_typing.py takes the "ours" side as a light per-cell CSV; dump the
# RIGHT run's labels (obs read backed — no full load) so both sides align on cell id.
apptainer exec --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -c "import anndata; a=anndata.read_h5ad('$WORK/right.h5ad', backed='r'); \
        a.obs[['$RIGHT_KEY']].rename(columns={'$RIGHT_KEY':'cell_type'}).to_csv('$WORK/right.csv')"

# LEFT = external (Sankey source / rows); RIGHT = ours (Sankey dest / cols).
apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/compare_external_typing.py" \
        --external-h5ad "$WORK/left.h5ad" --external-key "$LEFT_KEY" --label-external "$LEFT_LABEL" \
        --ours-csv "$WORK/right.csv" --label-ours "$RIGHT_LABEL" \
        --our-id-regex "$ID_RE" --external-id-regex "$ID_RE" \
        --output-dir "$WORK/out"

# Sankey source->dest (keeper -> InSituTree) straight from the cross-tab...
apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/plot_crosstab_sankey.py" \
        --crosstab "$WORK/out/external_vs_ours_crosstab.csv" \
        --label-left "$LEFT_LABEL" --label-right "$RIGHT_LABEL" \
        --output "$WORK/out/sankey_left_to_right.png"

# ...and the reverse direction (InSituTree -> keeper) from the transposed cross-tab.
apptainer exec --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -c "import pandas as pd; \
        pd.read_csv('$WORK/out/external_vs_ours_crosstab.csv', index_col=0).T \
          .to_csv('$WORK/out/crosstab_transposed.csv')"
apptainer exec --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" --bind "${WORK}:${WORK}" "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/plot_crosstab_sankey.py" \
        --crosstab "$WORK/out/crosstab_transposed.csv" \
        --label-left "$RIGHT_LABEL" --label-right "$LEFT_LABEL" \
        --output "$WORK/out/sankey_right_to_left.png"

echo "Uploading comparison outputs to Kopah..."
s5cmd cp "$WORK/out/*" "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${RIGHT_DIR}/${OUT_SUBDIR}/"

echo "Done: typing comparison (${LEFT_LABEL} vs ${RIGHT_LABEL})."
