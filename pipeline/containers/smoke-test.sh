#!/bin/bash
# Smoke test for the rapids-singlecell.sif Apptainer container.
#
# Runs four checks against $APPTAINER_RSC:
#   1. Imports     — scanpy / anndata / rapids_singlecell / dask_cuda load.
#   2. GPU         — cupy sees a device + can allocate. Auto-skipped when
#                    --nv is not available (e.g. on a login node).
#   3. Multi-GPU   — dask-cuda LocalCUDACluster spins up one worker per GPU.
#                    Auto-skipped when fewer than 2 GPUs are visible.
#   4. Stage 1     — pulls the migrated slide from Kopah and runs
#                    flatfiles_to_anndata.py inside the container, then
#                    asserts the post-merge invariants (obs.index unique +
#                    slide-prefixed, FOV column present). Auto-skipped
#                    when KOPAH creds are unset.
#
# Source pipeline/.env first so APPTAINER_RSC and KOPAH_* are populated:
#     set -a; source pipeline/.env; set +a
#     module load apptainer
#     bash pipeline/containers/smoke-test.sh

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

TEST_SLIDE_ID="${TEST_SLIDE_ID:-7134A77439A6}"
TEST_FLAT_PREFIX="${TEST_FLAT_PREFIX:-CosMx-GBM/Wenyu-cohort/Wenyupart142926_04_05_2026_11_38_05_67/flatFiles/${TEST_SLIDE_ID}}"

: "${APPTAINER_RSC:?must be set (source pipeline/.env first)}"

if ! command -v apptainer >/dev/null 2>&1; then
    echo "ERROR: apptainer not on PATH. Run 'module load apptainer' first." >&2
    exit 1
fi
if [[ ! -f "$APPTAINER_RSC" ]]; then
    echo "ERROR: APPTAINER_RSC=$APPTAINER_RSC does not exist." >&2
    exit 1
fi

echo "==> Tier 1: Python imports inside $APPTAINER_RSC"
apptainer exec "$APPTAINER_RSC" python -c "
import anndata, pandas, scanpy, rapids_singlecell, dask_cuda
print('  anndata          ', anndata.__version__)
print('  pandas           ', pandas.__version__)
print('  scanpy           ', scanpy.__version__)
print('  rapids_singlecell', rapids_singlecell.__version__)
print('  dask_cuda        ', dask_cuda.__version__)
"

GPU_COUNT=0
if apptainer exec --nv "$APPTAINER_RSC" nvidia-smi -L >/dev/null 2>&1; then
    GPU_COUNT=$(apptainer exec --nv "$APPTAINER_RSC" nvidia-smi -L | wc -l)
fi

if (( GPU_COUNT > 0 )); then
    echo "==> Tier 2: $GPU_COUNT GPU(s) visible — running cupy + rapids smoke"
    apptainer exec --nv "$APPTAINER_RSC" python -c "
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
    echo "    'salloc --account=glioblastoma --partition=gpu-l40s --gpus=2 ...'"
fi

if (( GPU_COUNT >= 2 )); then
    echo "==> Tier 3: $GPU_COUNT GPUs visible — verifying dask-cuda multi-GPU cluster"
    apptainer exec --nv "$APPTAINER_RSC" python - <<'PYEOF'
from dask_cuda import LocalCUDACluster
from dask.distributed import Client

with LocalCUDACluster() as cluster, Client(cluster) as client:
    workers = client.scheduler_info()["workers"]
    print(f"  workers up: {len(workers)}")
    assert len(workers) >= 2, f"expected >=2 dask-cuda workers, got {len(workers)}"
    devices = sorted(w["name"] for w in workers.values())
    print(f"  device names: {devices}")
    import cupy as cp
    def device_id():
        return cp.cuda.runtime.getDevice()
    seen = set(client.run(device_id).values())
    print(f"  distinct CUDA device IDs across workers: {sorted(seen)}")
    assert len(seen) >= 2, f"workers landed on the same GPU: {seen}"
print("  PASS — dask-cuda spans all visible GPUs")
PYEOF
else
    if (( GPU_COUNT == 1 )); then
        echo "==> Tier 3: SKIPPED — only 1 GPU visible. Multi-GPU path needs --gpus=2."
    else
        echo "==> Tier 3: SKIPPED — no GPUs visible."
    fi
fi

if [[ -z "${KOPAH_ACCESS_KEY_ID:-}" ]] || [[ -z "${KOPAH_SECRET_ACCESS_KEY:-}" ]]; then
    echo "==> Tier 4: SKIPPED — KOPAH creds not set (source pipeline/.env)."
    exit 0
fi
if ! command -v s5cmd >/dev/null 2>&1; then
    echo "==> Tier 4: SKIPPED — s5cmd not on PATH." >&2
    exit 0
fi
if [[ ! -f "$PIPELINE_DIR/manifest.csv" ]]; then
    echo "==> Tier 4: SKIPPED — $PIPELINE_DIR/manifest.csv missing; run build_manifest.py first."
    exit 0
fi

WORK="${TMPDIR:-/tmp}/rsc_smoke_$$"
mkdir -p "$WORK/flat"
trap 'rm -rf "$WORK"' EXIT

echo "==> Tier 4: staging $TEST_SLIDE_ID from Kopah to $WORK/flat"
AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID" \
AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY" \
S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL" \
    s5cmd cp "s3://${KOPAH_BUCKET}/${TEST_FLAT_PREFIX}/*" "$WORK/flat/"

echo "==> Tier 4: running flatfiles_to_anndata.py inside container"
apptainer exec \
    --bind "$WORK:$WORK" \
    --bind "$PIPELINE_DIR:$PIPELINE_DIR" \
    "$APPTAINER_RSC" \
    python "$PIPELINE_DIR/python/flatfiles_to_anndata.py" \
        --flatfiles-dir "$WORK/flat" \
        --slide-id "$TEST_SLIDE_ID" \
        --manifest "$PIPELINE_DIR/manifest.csv" \
        --output "$WORK/${TEST_SLIDE_ID}.h5ad"

echo "==> Tier 4: asserting invariants on the output AnnData"
apptainer exec --bind "$WORK:$WORK" "$APPTAINER_RSC" python - <<PYEOF
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
