# pipeline/containers

Apptainer recipes for the CosMx GBM pipeline on Hyak. Each `.def` file
produces a `.sif` that the Slurm jobs in `pipeline/slurm/` invoke with
`apptainer exec`.

## rapids-singlecell.sif

Runtime for the Python stages of the pipeline:

- **Stage 1** (`compute` partition): scanpy / anndata / pandas / boto3 →
  per-slide flat files → `.h5ad`. GPU not used but the same SIF works.
- **Stage 3a** (`compute` partition): concatenate per-slide `.h5ad`,
  per-cell QC, restrict to the study cohort (`pipeline/cohort_wenyu.csv`,
  matching `(Case, Block)` to drop the extraneous tissue co-mounted on each
  slide), emit the scPearsonPCA input. CPU-only (anndata/scipy).
- **Stage 3c** (`gpu-l40s` partition): neighbor graph + Leiden + UMAP on
  the Pearson-PCA embedding, plus UMAP QC plots colored by `Case` / `Region` /
  `leiden` for reviewing the patient-level batch correction. Runs on a
  **single** L40S — it operates on the small cells × ~50 embedding, so no
  dask-cuda cluster is needed here.

The actual Pearson-residual PCA (stage 3b) does **not** run in this
container — see `scpearsonpca.sif` below for why.

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

# Imports + GPU + multi-GPU + end-to-end — submit inside a 2-GPU allocation.
# Use --nodes=1 --gpus-per-node=2 (NOT plain --gpus=2): the latter can place
# 1 GPU each on two nodes, and dask-cuda's LocalCUDACluster is single-node, so
# Tier 3 would only see one GPU.
salloc --account=glioblastoma --partition=gpu-l40s \
    --nodes=1 --gpus-per-node=2 --cpus-per-task=8 --mem=32G --time=45:00
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

## scpearsonpca.sif

R runtime for **stage 3b**: the batch-corrected quasipoisson Pearson-residual
PCA (`scPearsonPCA::sparse_quasipoisson_pca_seurat_batch`). Batch correction is
at the **patient level** (the `Case` column), per Patrick Danaher's guidance that
per-slide effects are typically minor — matching the patient-batch example in the
[scPearsonPCA post](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/pearsonpca/).
A CPU/R stage that runs on the `compute` partition, **not** gpu-l40s.

Following that example, stage 3a hands this stage an **HVG subset** (~2000 genes
for the 6k panel) plus **precomputed all-gene** `totalcounts` and per-batch
`gene_frequency`; the R stage never sees the full panel.

Base: `rocker/r-ver:4.4.2`. Packages: `data.table`, `Matrix`, `RSpectra`,
`hdf5r` from a dated Posit Package Manager snapshot, plus `scPearsonPCA` pinned
to a commit SHA. No Seurat — the driver calls with `return_seurat_reduction =
FALSE` and gets back plain `cell.embeddings` + `feature.loadings`.

### Why a separate R stage, not rapids-singlecell's Pearson residuals

The vignette's dimension reduction is quasipoisson Pearson-residual PCA. The
obvious GPU route, `rsc.pp.normalize_pearson_residuals`, **densifies** an
`n_cells × n_genes` residual matrix on a single GPU (it has no dask path in
0.13.1). At this cohort's ~8M cells that is ~64 GB for 2000 genes — it does not
fit on a 46 GB L40S, and even both cards combined is too tight once PCA scratch
is added.

`scPearsonPCA` is purpose-built for exactly this: it computes the PCA via SVD of
the **genes × genes** cross-product (≈6000²) and projects cell embeddings in
blocks, so the dense `cells × genes` residual matrix is never formed. Peak
memory is the sparse counts (`genes × cells`) plus one working copy — tens of GB
on a big-memory CPU node, independent of how the residuals are defined. It is
the faithful *and* the only memory-feasible option here, not a compromise.

> **dgCMatrix nonzero ceiling.** R's `dgCMatrix` index slots are 32-bit, so the
> counts matrix must stay under 2³¹ (~2.1B) nonzeros. This is exactly why stage 3a
> hands R only the **HVG subset** rather than the full panel — the all-gene
> `totalcounts` and per-batch `gene_frequency` are computed in Python (scipy uses
> 64-bit indices, so the full cohort is fine) and passed in precomputed, so the
> dense residual matrix is never formed and R's matrix stays small. Stage 3a still
> errors clearly if even the HVG-subset nonzeros would cross the ceiling (reduce
> `--n-hvg` or split the cohort).

### Build

Same flow as `rapids-singlecell.sif`, but this recipe has no `%files`, so the
build CWD doesn't matter (the repo-root convention still works):

```bash
ssh klone.hyak.uw.edu
module load apptainer
export APPTAINER_CACHEDIR=/gscratch/scrubbed/$USER/apptainer-cache
export APPTAINER_TMPDIR=/gscratch/scrubbed/$USER/apptainer-tmp

cd ~/cosmx-utilities
apptainer build --fakeroot \
    /gscratch/scrubbed/$USER/containers/scpearsonpca.sif \
    pipeline/containers/scpearsonpca.def
```

Point `pipeline/.env` at it:

```
APPTAINER_SCPEARSON=/gscratch/scrubbed/<username>/containers/scpearsonpca.sif
```

### Reproducibility

R is pinned by the `rocker/r-ver:4.4.2` base; CRAN packages by the dated PPM
snapshot URL in the recipe; `scPearsonPCA` by commit SHA. To move forward, bump
the snapshot date and the SHA together and rebuild.
