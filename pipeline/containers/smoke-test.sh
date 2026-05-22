#!/bin/bash
# Smoke test for the scalesc.sif Apptainer container.
#
# Runs three checks against $APPTAINER_SCALESC:
#   1. Imports — scanpy / anndata / rapids_singlecell / scalesc all load.
#   2. GPU    — cupy sees a device + can allocate. Auto-skipped when --nv
#               is not available (e.g. on a login node or compute partition).
#   3. Stage 1 end-to-end — pulls the migrated slide from Kopah and runs
#               flatfiles_to_anndata.py inside the container, then asserts
#               the post-merge invariants (obs.index unique + slide-prefixed,
#               FOV column present). Auto-skipped when KOPAH creds are unset.
#
# Source pipeline/.env first so APPTAINER_SCALESC and KOPAH_* are populated:
#     set -a; source pipeline/.env; set +a
#     module load apptainer
#     bash pipeline/containers/smoke-test.sh

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

TEST_SLIDE_ID="${TEST_SLIDE_ID:-7134A77439A6}"
TEST_FLAT_PREFIX="${TEST_FLAT_PREFIX:-CosMx-GBM/Wenyu-cohort/Wenyupart142926_04_05_2026_11_38_05_67/flatFiles/${TEST_SLIDE_ID}}"

: "${APPTAINER_SCALESC:?must be set (source pipeline/.env first)}"

if ! command -v apptainer >/dev/null 2>&1; then
    echo "ERROR: apptainer not on PATH. Run 'module load apptainer' first." >&2
    exit 1
fi
if [[ ! -f "$APPTAINER_SCALESC" ]]; then
    echo "ERROR: APPTAINER_SCALESC=$APPTAINER_SCALESC does not exist." >&2
    exit 1
fi

echo "==> Tier 1: Python imports inside $APPTAINER_SCALESC"
apptainer exec "$APPTAINER_SCALESC" python -c "
import anndata, pandas, scanpy, rapids_singlecell, scalesc
print('  anndata          ', anndata.__version__)
print('  pandas           ', pandas.__version__)
print('  scanpy           ', scanpy.__version__)
print('  rapids_singlecell', rapids_singlecell.__version__)
print('  scalesc          ', getattr(scalesc, '__version__', 'unknown'))
"

if apptainer exec --nv "$APPTAINER_SCALESC" nvidia-smi -L >/dev/null 2>&1; then
    echo "==> Tier 2: GPU visible — running cupy + rapids smoke"
    apptainer exec --nv "$APPTAINER_SCALESC" python -c "
import cupy as cp
n = cp.cuda.runtime.getDeviceCount()
assert n > 0, 'no CUDA devices visible inside container'
print('  device count:', n)
x = cp.arange(10**6, dtype='float32')
print('  cupy sum:', float(x.sum()), 'on device', x.device)
import rapids_singlecell as rsc
print('  rapids_singlecell loaded:', rsc.__version__)
"
else
    echo "==> Tier 2: SKIPPED — no GPU visible. For GPU coverage, run inside"
    echo "    'salloc --account=glioblastoma --partition=gpu-l40s --gpus=1 ...'"
fi

if [[ -z "${KOPAH_ACCESS_KEY_ID:-}" ]] || [[ -z "${KOPAH_SECRET_ACCESS_KEY:-}" ]]; then
    echo "==> Tier 3: SKIPPED — KOPAH creds not set (source pipeline/.env)."
    exit 0
fi
if ! command -v s5cmd >/dev/null 2>&1; then
    echo "==> Tier 3: SKIPPED — s5cmd not on PATH." >&2
    exit 0
fi
if [[ ! -f "$PIPELINE_DIR/manifest.csv" ]]; then
    echo "==> Tier 3: SKIPPED — $PIPELINE_DIR/manifest.csv missing; run build_manifest.py first."
    exit 0
fi

WORK="${TMPDIR:-/tmp}/scalesc_smoke_$$"
mkdir -p "$WORK/flat"
trap 'rm -rf "$WORK"' EXIT

echo "==> Tier 3: staging $TEST_SLIDE_ID from Kopah to $WORK/flat"
AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID" \
AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY" \
S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL" \
    s5cmd cp "s3://${KOPAH_BUCKET}/${TEST_FLAT_PREFIX}/*" "$WORK/flat/"

echo "==> Tier 3: running flatfiles_to_anndata.py inside container"
apptainer exec \
    --bind "$WORK:$WORK" \
    --bind "$PIPELINE_DIR:$PIPELINE_DIR" \
    "$APPTAINER_SCALESC" \
    python "$PIPELINE_DIR/python/flatfiles_to_anndata.py" \
        --flatfiles-dir "$WORK/flat" \
        --slide-id "$TEST_SLIDE_ID" \
        --manifest "$PIPELINE_DIR/manifest.csv" \
        --output "$WORK/${TEST_SLIDE_ID}.h5ad"

echo "==> Tier 3: asserting invariants on the output AnnData"
apptainer exec --bind "$WORK:$WORK" "$APPTAINER_SCALESC" python - <<PYEOF
import anndata as ad
a = ad.read_h5ad("$WORK/${TEST_SLIDE_ID}.h5ad")
print("  shape:", a.shape)
print("  obs columns:", list(a.obs.columns))
assert a.obs.index.is_unique, "obs.index not unique"
assert "FOV" in a.obs.columns, "FOV column missing"
first = a.obs.index[0]
assert first.startswith("${TEST_SLIDE_ID}_F"), f"unexpected obs index: {first!r}"
fov = a.obs["FOV"].iloc[0]
assert fov.startswith("${TEST_SLIDE_ID}_F"), f"unexpected FOV value: {fov!r}"
print("  PASS — obs unique, FOV present, slide-prefixed IDs")
PYEOF

echo
echo "All tiers PASS."
