# pipeline/containers

Apptainer recipes for the CosMx GBM pipeline on Hyak. Each `.def` file
produces a `.sif` that the Slurm jobs in `pipeline/slurm/` invoke with
`apptainer exec`.

## rapids-singlecell.sif

Runtime for stages 1–6 of the pipeline:

- **Stage 1** (`compute` partition): scanpy / anndata / pandas / boto3 →
  per-slide flat files → `.h5ad`. GPU not used but the same SIF works.
- **Stage 3** (`gpu-l40s` partition): rapids-singlecell + Harmony for
  GPU-accelerated QC / batch correction on the subsampled cohort
  AnnData, distributed across 2x L40S via dask-cuda.
- **Stage 6** (`gpu-l40s`): GPU embedding on the full cohort.

Base: `nvidia/cuda:12.4.1-cudnn-devel-ubuntu22.04` (devel, not runtime —
CuPy's NVRTC JIT needs the CUDA headers). Python 3.10 venv at `/opt/venv`
(on `PATH` via `%environment`), built with `uv` from a pinned lock file.

### Why rapids-singlecell, not ScaleSC

Earlier iterations used ScaleSC on top of rapids-singlecell. ScaleSC's
chunked PCA/Harmony was designed around a single A100 80GB; on the
gpu-l40s partition (2x L40S, ~46GB VRAM each, no NVLink) it ran into
GPU memory allocation failures because its chunk sizes assume more
headroom per device. dask-cuda's `LocalCUDACluster` spreads cell-axis
chunks across both L40S cards, giving ~92GB aggregate VRAM and earning
both GPUs' keep. rapids-singlecell drives the underlying GPU primitives
either way — ScaleSC was only adding scale-out logic we no longer need.

### Build

On any Linux host with Apptainer + fakeroot. Klone login nodes work.
**Build from the repo root** — the recipe's `%files` section copies
`pipeline/containers/rapids-singlecell.lock` using a path relative to the
build's working directory, so the CWD must be `~/cosmx-utilities`. Keep
the cache + tmpdir on `/gscratch` (home is only 10 GB and Apptainer's
layer cache will overflow it); the output SIF can go there too.

```bash
ssh klone.hyak.uw.edu
module load apptainer

# Build + cache space (community scrubbed has no quota, 21-day purge).
mkdir -p /gscratch/scrubbed/$USER/{apptainer-cache,apptainer-tmp,containers}
export APPTAINER_CACHEDIR=/gscratch/scrubbed/$USER/apptainer-cache
export APPTAINER_TMPDIR=/gscratch/scrubbed/$USER/apptainer-tmp

cd ~/cosmx-utilities                 # build context = repo root (for %files)
apptainer build --fakeroot \
    /gscratch/scrubbed/$USER/containers/rapids-singlecell.sif \
    pipeline/containers/rapids-singlecell.def
```

Do **not** build from `/gscratch/scrubbed/$USER` with a path to the `.def`
elsewhere — `%files` would look for the lock under the wrong directory and
the build fails.

Expect a SIF in the ~11–15 GB range, dominated by the CUDA + cuDNN
**devel** base (needed for CuPy's NVRTC headers) and the RAPIDS wheels.
The build installs with `uv pip sync` against `rapids-singlecell.lock`, so
it reproduces the exact validated version set rather than re-resolving.

If `--fakeroot` is not available, build on any Linux box with Apptainer
installed (e.g. a local VM) and `scp` the SIF to klone.

### Dependencies and the lock file

`rapids-singlecell.lock` is a full `pip freeze` of a validated build,
including transitive pins (cupy, cudf, cuml, rmm). The `.def` installs it
with `uv pip sync ... --index-strategy unsafe-best-match` (RAPIDS wheels
are split across PyPI and `pypi.nvidia.com`). To intentionally bump
versions: edit/regenerate the lock, rebuild, and re-run the smoke test —
then commit the new lock alongside. Pinning the whole tree is deliberate:
an unpinned cupy silently jumped 13.x → 14.1.0 between builds and broke
GPU kernel compilation.

### Install

The SIF lives where you built it; just point `pipeline/.env` at it:

```
APPTAINER_RSC=/gscratch/scrubbed/<username>/containers/rapids-singlecell.sif
```

Scrubbed is purged after 21 days of no access. The slurm jobs `apptainer
exec` it on every run, so as long as the pipeline runs at least every
three weeks, the SIF stays put. If it goes idle, rebuild — the recipe
is the source of truth.

### Smoke test

`smoke-test.sh` bundles four checks: imports, single-GPU sanity, a
dask-cuda multi-GPU cluster check, and an end-to-end stage 1 run
against the migrated test slide. Auto-skips sections that can't run in
the current environment (no `--nv`, only 1 GPU, no Kopah creds).

```bash
set -a; source pipeline/.env; set +a
module load apptainer

# Imports only — runs anywhere, including the login node:
bash pipeline/containers/smoke-test.sh

# Imports + GPU + multi-GPU + end-to-end — submit inside a 2-GPU allocation:
salloc --account=glioblastoma --partition=gpu-l40s \
    --time=30:00 --cpus-per-task=8 --mem=32G --gpus=2
set -a; source pipeline/.env; set +a
bash pipeline/containers/smoke-test.sh
```

Override the test slide with `TEST_SLIDE_ID=...` if you've migrated a
different one. `--nv` is what exposes the host's NVIDIA driver to the
container; the script adds it automatically when a GPU is available.

### Multi-GPU memory tips

The L40S cards have no NVLink, so cross-GPU traffic goes over PCIe
Gen4. Practical implications when wiring up Stage 3:

- Initialize RMM with managed memory + a pool before any CuPy/cuDF
  imports, so a chunk that briefly exceeds VRAM spills to host RAM
  instead of OOM'ing:

  ```python
  import rmm
  rmm.reinitialize(managed_memory=True, pool_allocator=True)
  ```

- Spin up the cluster with `LocalCUDACluster()` (no args picks up
  `CUDA_VISIBLE_DEVICES` from Slurm).

- Build the kNN graph once on a single GPU after PCA/HVG — the
  post-reduction data is small and avoids PCIe-bound communication for
  Leiden / UMAP.
