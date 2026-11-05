#!/bin/bash
# Stage 5: add typing columns to the stitched Napari `_metadata.csv` files, one pass per
# column, each merged into the output of the last so ONE new file per slide carries them all.
#
# Why a job rather than running the script by hand: the inputs now live on Kopah -- the
# clustered AnnData (for obs['leiden'], which the component splits key on) and a different
# typing result per column. Staging those by hand for each pass is where mistakes happen.
#
# The passes come from a committed spec CSV rather than flags, so what a reviewer was shown
# is recoverable afterwards:
#   column,stage_dir,result_basename,split_spec,note
# `split_spec` is optional and names a CSV in pipeline/reference/ (see --split-spec).
#
# Submit (the existing merged files must already be on Klone -- scp them up first):
#   NAPARI_PASSES=napari_passes_retina.csv \
#     MERGE_INTO_DIR=$HOME/retina_napari_metadata_merged \
#     sbatch pipeline/slurm/87_napari_metadata.sh
#
# Then scp $HOME/napari_metadata_out back to the Mac and upload to S3 from there -- Klone has
# Kopah credentials, not AWS ones.
#
# Runs under `uv run`, NOT apptainer, unlike every other stage. The colour function is
# DuckDB's hash(), shared by six scripts so a label keeps one colour whichever wrote it, and
# duckdb is a project dependency that the analysis containers do not carry. This script needs
# no RAPIDS -- it reads obs backed and writes CSVs.
#
# It uses --no-project with an explicit dependency list rather than the project environment:
# `uv sync` would also build napari-cosmx-fork, a GUI package irrelevant here, and a job
# should not depend on a project env being in any particular state. Warm the cache once on a
# login node (same command, it is cheap) so the compute node does not need to resolve.
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY.
# Set UV_CACHE_DIR and UV_PROJECT_ENVIRONMENT to /gscratch paths in .env: the default puts
# both under $HOME, whose 10 GB quota the uv cache fills quickly.
# Required args:
#   NAPARI_PASSES                spec CSV in pipeline/reference/
#   MERGE_INTO_DIR               directory of existing <slide>_metadata.csv to add to
# Optional:
#   STAGE3_DIR=stage3            where cosmx_clustered.h5ad lives
#   CLUSTERED_BASENAME=cosmx_clustered.h5ad
#   NAPARI_OUT_DIR=$HOME/napari_metadata_out    where the finished files are left

#SBATCH --job-name=cosmx-napari-metadata
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=00:40:00
# The AnnData is read BACKED (obs only) once per pass; the floor is the obs table plus the
# metadata CSVs being rewritten, not the count matrix.
#SBATCH --output=pipeline/logs/napari_metadata_%j.out
#SBATCH --error=pipeline/logs/napari_metadata_%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
# s5cmd lives in ~/bin; no apptainer module needed, this stage runs under uv.
export PATH="${HOME}/bin:${PATH}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_DIR="${SLURM_SUBMIT_DIR}"
else
    REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
fi
PIPELINE_DIR="${REPO_DIR}/pipeline"

# Loads .env WITHOUT overwriting what the caller passed; see slurm/_load_env.sh.
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/slurm/_load_env.sh"

STAGE3="${STAGE3_DIR:-stage3}"
# TWO AnnData files, because the pieces live apart: the Napari join key
# (obs['flatfile_cell_id']) is in the stage-3a combined file, and obs['leiden'], which the
# component splits key on, is in the stage-3c clustered one -- 3c does not carry the key
# forward. Only passes with a split_spec need the second.
COMBINED="${COMBINED_BASENAME:-combined_qc.h5ad}"
CLUSTERED="${CLUSTERED_BASENAME:-cosmx_clustered.h5ad}"
NAPARI_OUT_DIR="${NAPARI_OUT_DIR:-${HOME}/napari_metadata_out}"

# What celltypes-to-napari-metadata.py imports. Deliberately unquoted below so each --with
# is a separate word. duckdb is the one that is NOT in the analysis containers.
UV_DEPS="--with duckdb --with pandas --with numpy --with h5py --with anndata"

: "${NAPARI_PASSES:?set it to a spec CSV in pipeline/reference/}"
: "${MERGE_INTO_DIR:?set it to the directory of existing <slide>_metadata.csv files}"

PASSES_CSV="${PIPELINE_DIR}/reference/${NAPARI_PASSES}"
[[ -f "$PASSES_CSV" ]] || { echo "ERROR: no such passes spec: $PASSES_CSV" >&2; exit 1; }
[[ -d "$MERGE_INTO_DIR" ]] || { echo "ERROR: no such directory: $MERGE_INTO_DIR" >&2; exit 1; }

existing_count=$(find "$MERGE_INTO_DIR" -name '*_metadata.csv' | wc -l)
(( existing_count > 0 )) || { echo "ERROR: no *_metadata.csv in $MERGE_INTO_DIR" >&2; exit 1; }
echo "Merging into ${existing_count} existing metadata file(s) from ${MERGE_INTO_DIR}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_napari_metadata_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${COMBINED} (${STAGE3}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${COMBINED}" "$WORK/combined.h5ad"

# Only fetched if some pass actually needs a split column.
if tail -n +2 "$PASSES_CSV" | awk -F, '{print $4}' | grep -q '[^[:space:]]'; then
    echo "Staging ${CLUSTERED} (${STAGE3}) for the split column(s)..."
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE3}/${CLUSTERED}" "$WORK/clustered.h5ad"
fi

# Each pass merges into the previous pass's output, so the last one holds every column.
cp -r "$MERGE_INTO_DIR" "$WORK/metadata"

pass_number=0
# Skip the header; ignore blank lines and comments.
while IFS=, read -r column stage_dir result_basename split_spec _note; do
    [[ -z "${column// }" || "${column:0:1}" == "#" ]] && continue
    pass_number=$((pass_number + 1))
    echo ""
    echo "=== pass ${pass_number}: column '${column}' from ${stage_dir}/${result_basename} ==="

    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${stage_dir}/${result_basename}" \
        "$WORK/typing_${pass_number}.h5"

    SPLIT_ARGS=()
    if [[ -n "${split_spec// }" ]]; then
        SPLIT_CSV="${PIPELINE_DIR}/reference/${split_spec}"
        [[ -f "$SPLIT_CSV" ]] || { echo "ERROR: no such split spec: $SPLIT_CSV" >&2; exit 1; }
        SPLIT_ARGS=(--split-spec "$SPLIT_CSV" --split-obs-h5ad "$WORK/clustered.h5ad")
        echo "    applying splits from ${split_spec}"
    fi

    uv run --no-project ${UV_DEPS} \
        python -u "${REPO_DIR}/scripts/celltypes-to-napari-metadata.py" \
            --typing-h5 "$WORK/typing_${pass_number}.h5" \
            --combined-h5ad "$WORK/combined.h5ad" \
            --column "$column" \
            --merge-into "$WORK/metadata" \
            --out-dir "$WORK/metadata" \
            ${SPLIT_ARGS[@]+"${SPLIT_ARGS[@]}"}
done < <(tail -n +2 "$PASSES_CSV")

(( pass_number > 0 )) || { echo "ERROR: $PASSES_CSV defined no passes" >&2; exit 1; }

mkdir -p "$NAPARI_OUT_DIR"
cp "$WORK"/metadata/*_metadata.csv "$NAPARI_OUT_DIR/"
echo ""
echo "Done: ${pass_number} pass(es) over ${existing_count} slide(s)."
echo "Files are in ${NAPARI_OUT_DIR} -- scp them to the Mac and upload to S3 from there."
echo "Columns now present:"
head -1 "$(find "$NAPARI_OUT_DIR" -name '*_metadata.csv' | head -1)" | tr ',' '\n' | grep -v '_color$'
