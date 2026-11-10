#!/bin/bash
# Stage 4c′: build the InSituTree profile matrix and put it on Kopah for 85.
#
# This was a hand-typed prerequisite in 85's header. It is a stage because getting it wrong
# is silent: InSituTree does no platform rescaling of its own, so every column of this matrix
# must come from ONE InSituType run's post-rescale $profiles, and every hierarchy leaf must
# be a column of it and vice versa. A mismatch surfaces as cells vanishing into branches
# rather than as an error, so the lockstep check runs here, before the fit.
#
# Submit:
#   DENOVO_ANNOTATIONS=retina_stage4.csv HIERARCHY_BASENAME=retina_hierarchy.json \
#     TREE_PROFILES_BASENAME=retina_insitutree_profiles.csv \
#     sbatch pipeline/slurm/83_prep_insitutree_profiles.sh
# Then chain the fit off it:
#   jid=$(sbatch --parsable pipeline/slurm/83_prep_insitutree_profiles.sh)
#   sbatch --dependency=afterok:$jid pipeline/slurm/85_insitutree.sh
#
# Required env (from pipeline/.env):
#   KOPAH_ENDPOINT_URL, KOPAH_BUCKET, KOPAH_PREFIX,
#   KOPAH_ACCESS_KEY_ID, KOPAH_SECRET_ACCESS_KEY,
#   APPTAINER_RSC (the rapids-singlecell SIF; CPU use, no --nv).
# Required args:
#   DENOVO_ANNOTATIONS           CSV in pipeline/reference/denovo_annotations/ with
#                                keep + new_name columns (the data-driven schema; a CSV
#                                without them falls back to the GBM pilot rename map)
#   HIERARCHY_BASENAME           hierarchy JSON in pipeline/reference/
#   TREE_PROFILES_BASENAME       output name, uploaded to Kopah reference/
# Optional:
#   SOURCE_RESULT_BASENAME=insitutype_result.h5   the InSituType run to source profiles from
#   DROP_NAMED=                  comma-separated named types to drop from the matrix

#SBATCH --job-name=cosmx-prep-insitutree-profiles
# Minutes of pandas on one h5; ckpt preemption costs nothing and it is idempotent.
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=pipeline/logs/prep_insitutree_profiles_%j.out
#SBATCH --error=pipeline/logs/prep_insitutree_profiles_%j.err

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
SOURCE_RESULT_BASENAME="${SOURCE_RESULT_BASENAME:-insitutype_result.h5}"

: "${APPTAINER_RSC:?must be set in pipeline/.env}"
: "${DENOVO_ANNOTATIONS:?set it to a CSV in pipeline/reference/denovo_annotations/}"
: "${HIERARCHY_BASENAME:?set it to a hierarchy JSON in pipeline/reference/}"
: "${TREE_PROFILES_BASENAME:?set it to the output profile matrix name}"

ANNOTATIONS_CSV="${PIPELINE_DIR}/reference/denovo_annotations/${DENOVO_ANNOTATIONS}"
HIERARCHY_JSON="${PIPELINE_DIR}/reference/${HIERARCHY_BASENAME}"
for required_file in "$ANNOTATIONS_CSV" "$HIERARCHY_JSON"; do
    [[ -f "$required_file" ]] || { echo "ERROR: no such file: $required_file" >&2; exit 1; }
done

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_prep_insitutree_profiles_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging ${SOURCE_RESULT_BASENAME} (${STAGE4}) from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${SOURCE_RESULT_BASENAME}" \
    "$WORK/source_result.h5"

DROP_ARGS=()
if [[ -n "${DROP_NAMED:-}" ]]; then
    IFS=',' read -ra DROP_TYPES <<< "$DROP_NAMED"
    for drop_type in "${DROP_TYPES[@]}"; do
        drop_type="$(echo "$drop_type" | sed 's/^ *//; s/ *$//')"
        [[ -n "$drop_type" ]] && DROP_ARGS+=(--drop-named "$drop_type")
    done
    echo "Dropping named types: ${DROP_TYPES[*]}"
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u "${PIPELINE_DIR}/python/prep_insitutree_profiles.py" \
        --profiles-h5 "$WORK/source_result.h5" \
        --annotations "$ANNOTATIONS_CSV" \
        --output "$WORK/${TREE_PROFILES_BASENAME}" \
        ${DROP_ARGS[@]+"${DROP_ARGS[@]}"}

# Lockstep: every leaf a column and every column a leaf. InSituTree will not complain --
# a leaf with no profile simply never wins, and a column with no leaf is never offered, so
# cells quietly go somewhere else. Fail here instead.
echo "Checking hierarchy leaves against profile columns..."
apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python -u -c "
import json, sys
import pandas as pd

columns = set(pd.read_csv('$WORK/${TREE_PROFILES_BASENAME}', index_col=0, nrows=1).columns)
leaves = []
def walk(node):
    if isinstance(node, dict):
        for child in node.values(): walk(child)
    elif isinstance(node, list): leaves.extend(node)
walk(json.load(open('$HIERARCHY_JSON')))

duplicated = sorted({leaf for leaf in leaves if leaves.count(leaf) > 1})
no_leaf = sorted(columns - set(leaves))
no_column = sorted(set(leaves) - columns)
print(f'{len(leaves)} leaves ({len(set(leaves))} unique), {len(columns)} profile columns')
if duplicated: print('LEAVES APPEARING MORE THAN ONCE:', duplicated)
if no_leaf:    print('PROFILE COLUMNS WITH NO LEAF:', no_leaf)
if no_column:  print('LEAVES WITH NO PROFILE COLUMN:', no_column)
if duplicated or no_leaf or no_column:
    sys.exit('ERROR: hierarchy and profile matrix are not in lockstep')
print('LOCKSTEP OK')
"

echo "Uploading ${TREE_PROFILES_BASENAME} to Kopah reference/..."
s5cmd cp "$WORK/${TREE_PROFILES_BASENAME}" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${TREE_PROFILES_BASENAME}"

# NOT REFERENCE_BASENAME: .env defines that for the atlas reference and would overwrite it.
echo "Done. Run 85 with:"
echo "  TREE_PROFILES_BASENAME=${TREE_PROFILES_BASENAME} \\"
echo "  HIERARCHY_BASENAME=${HIERARCHY_BASENAME}"
