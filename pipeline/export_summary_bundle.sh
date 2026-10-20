#!/bin/bash
# Pull every small output needed to write up a study offline into ONE tarball on Klone scratch:
# figures, CSV/TXT tables, compact InSituType results (not the .rds or the typed AnnData), and
# the run logs. Run on a login node (not sbatch). Kopah keeps everything; this is only a
# convenience copy for when Hyak is unreachable.
#
# Usage:  bash pipeline/export_summary_bundle.sh [bundle-name]
# Env (pipeline/.env): KOPAH_*. Extra Kopah sub-dirs under ${KOPAH_PREFIX}/ via EXTRA_DIRS="a b".

set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

NAME="${1:-sorl1_summary_bundle}"
DEST="/gscratch/scrubbed/${USER}/${NAME}"
mkdir -p "$DEST"

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"
BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

# Whole sub-dirs of small files (figures, tables, gene selection, K sweep).
SMALL_DIRS="stage3/qc_plots stage3/umap_panels stage3/depth stage3/qc stage3/labels stage3/marker_heatmap stage4/gene_selection stage4/k_sweep_pruned reference ${EXTRA_DIRS:-}"
for d in $SMALL_DIRS; do
    mkdir -p "$DEST/$d"
    s5cmd cp "${BASE}/${d}/*" "$DEST/$d/" 2>/dev/null || echo "  (nothing under ${d})"
done

# Typing runs: everything EXCEPT the full .rds and the multi-GB typed AnnData.
for d in $(s5cmd ls "${BASE}/" | awk '$NF ~ /^stage4_/ {print $NF}'); do
    d="${d%/}"
    mkdir -p "$DEST/$d"
    s5cmd cp --exclude "*.rds" --exclude "cosmx_typed.h5ad" "${BASE}/${d}/*" "$DEST/$d/"
done

mkdir -p "$DEST/logs"
cp "${PIPELINE_DIR}"/logs/*.out "${PIPELINE_DIR}"/logs/*.err "$DEST/logs/" 2>/dev/null || true
cp "${PIPELINE_DIR}/.env.example" "$DEST/"

tar -C "$(dirname "$DEST")" -czf "${DEST}.tar.gz" "$NAME"
ls -lh "${DEST}.tar.gz"
du -sh "$DEST"
echo "Done: ${DEST}.tar.gz"
