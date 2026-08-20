# cosmx-utilities

A scalable, end-to-end pipeline for processing [CosMx Spatial Molecular Imager](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/single-cell-imaging-overview/) data from the [BRaIN Lab](https://dlmp.uw.edu/research-labs/brainlab) at the University of Washington. 

This pipeline automates the workflow for converting raw CosMx slide exports to interactive spatial visualization in [Napari](https://napari.org), using AWS cloud infrastructure for on-demand, cost-effective compute and a fork of [napari-cosmx](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/napari-cosmx-intro/) for headless image processing and data exploration.

## Pipeline overview

After exporting a study from the AtoMx Spatial Informatics Platform, raw slide data is staged on AWS S3. Each slide is processed in parallel on AWS Fargate: FOV images are stitched into whole-slide mosaics, RNA transcript locations are decoded, and cell metadata is generated, all formatted for viewing in Napari. Processed data is uploaded back to S3, where it can be loaded into Napari for visualization across slides.

```mermaid
flowchart TD
    classDef action fill:#4b2e83,stroke:none,color:#ffffff
    classDef data fill:#b7a57a,stroke:none,color:#32006e
    classDef viewer fill:#ffc700,stroke:none,color:#32006e
    classDef futureAction fill:#4b2e83,stroke:#b7a57a,stroke-width:2px,stroke-dasharray: 5 5,color:#ffffff
    classDef done fill:#4b2e83,stroke:none,color:#ffffff
    classDef running fill:#ffc700,stroke:#4b2e83,stroke-width:3px,color:#32006e
    classDef todo fill:#4b2e83,stroke:#b7a57a,stroke-width:2px,stroke-dasharray: 5 5,color:#ffffff

    style fargate fill:#c5b4e3,stroke:none
    style ec2 fill:#c5b4e3,stroke:none
    style hyak fill:#e8e1cc,stroke:#85754d

    atomx["AtoMx study export"]:::data
    sftp["Download slides from BSB SFTP endpoint"]:::action
    s3raw["Raw slide data<br>(AWS S3)"]:::data
    s3out["Processed slide data<br>(AWS S3)"]:::data
    napari["Explore in Napari with napari-cosmx-fork"]:::viewer

    dockerfile["Multi-stage Dockerfile"]:::data
    headless["Headless image<br>(napari-cosmx-fork CLI + DuckDB + AWS CLI)"]:::action

    dockerfile --> headless
    headless --> stitch
    headless --> targets
    headless --> meta

    subgraph fargate["AWS Fargate"]
        direction TB
        stitch["Stitch FOV images into whole-slide mosaics with napari-cosmx-fork"]:::action
        targets["Decode RNA transcript locations with napari-cosmx-fork"]:::action
        meta["Generate cell metadata and labels using DuckDB SQL"]:::action
        stitch --> targets --> meta
    end

    subgraph ec2["AWS EC2 (auto-provisioned)"]
        direction TB
        ami["Build custom AMI<br>(create_ami.py)"]:::action
        launch["Launch GPU instance<br>(start_ec2.py --napari)"]:::action
        sync["Sync processed data from S3 to NVMe"]:::action
        dcv["Connect via DCV remote desktop"]:::action
        ami --> launch --> sync --> dcv
    end

    migrate["Migrate flat files S3 → Kopah<br>migrate_s3_to_kopah.py"]:::done
    kopahflat["Kopah flat files<br>(brainlabkg/CosMx-GBM)"]:::data

    subgraph hyak["UW Hyak (Klone) — cell-typing pipeline · color = progress"]
        direction TB
        stage1["Stage 1 · flat files → per-slide .h5ad<br>10_flatfiles_to_anndata.sh · ckpt array · rapids-singlecell.sif"]:::done
        anndata["…/anndata/ · 28 per-slide .h5ad"]:::data
        stage3a["Stage 3a · concat + QC + cohort filter (12 donors)<br>2000 HVGs + per-Case gene frequency<br>20_concat_qc.sh · ckpt · rapids-singlecell.sif"]:::done
        combined["…/stage3/ · combined_qc.h5ad + pca_input.h5<br>2.33M cohort cells × 6519 probes"]:::data
        stage3b["Stage 3b · quasipoisson Pearson-residual PCA<br>batch-corrected by patient (Case) · scPearsonPCA<br>30_pearson_pca.sh · ckpt · scpearsonpca.sif"]:::done
        embedding["…/stage3/embedding.h5 · 2.33M × 50 PCs"]:::data
        stage3c["Stage 3c · neighbors → Leiden (1.2) → UMAP<br>+ Case / Region QC plots<br>40_cluster.sh · gpu-l40s · rapids-singlecell.sif"]:::done
        clustered["…/stage3/ · cosmx_clustered.h5ad + qc_plots/"]:::data
        stage4a["Stage 4a · gene counts + per-cell negprobe mean<br>70_prep_insitutype.sh · ckpt · rapids-singlecell.sif"]:::done
        itinput["…/stage4/ · insitutype_input.h5"]:::data
        stage4b["Stage 4b · semi-supervised InSituType vs GBmap (level 4)<br>rescale, n_clusts 10:20 · 80_insitutype.sh · ckpt · insitutype.sif"]:::done
        itresult["…/stage4/ · insitutype_result.h5 + .rds"]:::data
        stage4c["Stage 4c · write cell_type back to obs + UMAP<br>90_write_celltypes.sh · ckpt · rapids-singlecell.sif"]:::done
        typed["…/stage4/ · cosmx_typed.h5ad + qc_plots/"]:::data
        stage1 --> anndata --> stage3a --> combined --> stage3b --> embedding --> stage3c --> clustered --> stage4a --> itinput --> stage4b --> itresult --> stage4c --> typed
    end

    atomx --> sftp
    sftp --> s3raw
    s3raw -- per slide (parallel) --> stitch
    meta --> s3out
    s3raw -- flat files --> migrate
    migrate --> kopahflat
    kopahflat -- per slide (Slurm array) --> stage1
    clustered -.->|cell-type labels| napari
    s3out --> launch
    dcv --> napari
```

**Progress legend (Hyak cell-typing pipeline):** solid purple = done · yellow = running now · dashed border = not yet run / planned · gold = data artifact on Kopah. Batch correction is at the patient (`Case`) level; tissue regions (tumor bulk / infiltrating edge / contralateral) are preserved as biological signal.

## Key design decisions

**Ephemeral Fargate compute** — Using Fargate's ephemeral storage rather than manually provisioning EC2 instances and EBS volumes for each study avoids the accumulation of long-lived AWS resources that are common with ad hoc processing workflows and costly to maintain.

**DuckDB for metadata streaming** — Cell metadata is streamed from S3 with DuckDB SQL queries and transformed into Napari metadata layers, enabling interactive overlay of cell annotations without storing large metadata files locally.

## Repository structure

```
cosmx-utilities/
├── napari-cosmx-fork/       # Fork of napari-cosmx with optional GUI dependencies
├── scripts/
│   ├── process-slide.py     # Single slide: detect segmentation version, download, stitch, upload
│   ├── process-slides.py    # Discover all slides in S3 and launch Fargate tasks
│   ├── generate-slide-metadata.py  # Generate _metadata.csv with cell metadata and deterministic colors
│   ├── compare-cell-typing.py  # Compare two cell-type columns per run: abundance bars + Sankey
│   └── diagnose-cell-typing.py  # Characterize de-novo/low-signal clusters + gene-pruning + marker heatmap
├── fargate/                 # Fargate task definitions and IAM configuration
├── ec2/                     # EC2 auto-provisioning for analytics and Napari viewer instances
│   ├── start_ec2.py         # Launch analytics (r5a) or Napari (g4dn GPU) instances
│   ├── create_ami.py        # Build custom AMI with all dependencies pre-installed
│   └── ami_setup.sh         # User-data setup script (deps, volumes, DCV)
├── Dockerfile               # Multi-stage: headless (Fargate) and GUI (desktop) targets
└── pyproject.toml           # uv workspace with optional [gui] extra
```

## Quick start

### Local development

```bash
# Install CLI tools (no GUI required)
uv sync
uv run stitch-images --help
uv run read-targets --help

# Install with Napari GUI for interactive visualization
uv sync --extra gui
uv run napari /path/to/processed-output
```

### Processing slides on Fargate

Discover all slides under an S3 experiment directory and launch one Fargate task per slide:

```bash
# Preview what would be launched
uv run python scripts/process-slides.py s3://bucket/project/study/ --whatif

# Launch Fargate tasks (one per slide, parallel)
uv run python scripts/process-slides.py s3://bucket/project/study/

# Skip already processed slides
uv run python scripts/process-slides.py s3://bucket/project/study/ --skip
```

The S3 URI may point at a whole study or at a single AtoMx run. Point it at one
run when two studies cover the same slides but need different flags — a 3D
resegmentation takes `--input-ndim 3` while the original 2D run does not.

Each Fargate task runs `process-slide.py`, which:
1. Queries segmentation manifests in S3 via DuckDB to find the correct segmentation version
2. Downloads only the needed CellLabels, morphology images, and AnalysisResults
3. Stitches FOV images into whole-slide zarr mosaics (`stitch-images`)
4. Decodes RNA transcript locations (`read-targets`)
5. Generates `_metadata.csv` with metadata annotations such as cell type and assigns deterministic colors
6. Uploads processed results back to S3

#### Carrying annotations into Napari

`--column OUTNAME=SOURCE_HEADER` (repeatable) selects which flatFiles columns
become colorable annotations in `_metadata.csv`. `SOURCE_HEADER` may list
fallbacks separated by `|`, taking the first header the file actually has — use
this when a study renamed a column between AtoMx runs:

```bash
uv run python scripts/process-slides.py s3://bucket/study/atomx-run/ \
    --column 'Cell Type=RNA_RNA_Cell.Typing.InSituType.2_1_clusters' \
    --column 'Case Specific=Case_specific_SORL1|Case_specific'
```

When cell typing is re-run, AtoMx nests the new export under the original run
(`<run>/flatFiles/<rerun>/flatFiles/<slide>/`). Both are discovered, and the one
that actually has the requested columns is used — the segmentation ID cannot
tell them apart, since a re-export carries the same segmentation.

When AtoMx never captured the annotations at all, supply them as per-FOV
sheets — one row per FOV, one file per slide named `<slide>_annotations.csv` —
and point `--annotations-prefix` at the directory holding them:

```bash
uv run python scripts/process-slides.py s3://bucket/study/atomx-run/ \
    --annotations-prefix s3://bucket/study/annotations \
    --column 'UWA=UWA' --column 'Case Broad=Case_broad'
```

The sheet's FOV column may be named `FOVs`, `FOV`, or `fov`. If it also names
the slide (`Flow Cells`, `Run_Tissue_name`, ...), that is checked against the
slide being processed, so pointing at the wrong file fails loudly instead of
labelling cells with another case.

If instead the annotations exist on another study covering the same physical
slides, `--fill-from <experiment prefix>` transfers them by joining on FOV. Only blanks are filled, and only for columns constant within an
FOV; per-cell columns such as cell typing belong to the segmentation that
produced them and are refused, with a line saying which were skipped.

#### 3D segmentation

A 3D-resegmented slide has per-z CellLabels (`CellLabels_F#####_Z###.tif`).
Pass `--input-ndim 3 --output-ndim 3` for a z-navigable labels layer in Napari.

Such a slide is usually 3D labels over a **2D** morphology acquisition
(`Morphology2D/`, no `_Z###` in the filenames). That is detected automatically,
and the single morphology plane is stored once rather than copied to every z
plane. This matters for more than tidiness: morphology dominates mosaic size, so
duplicating it across 8 planes would take a ~30 GiB slide to ~230 GiB, past
Fargate's 200 GiB ephemeral storage ceiling (which is also the Fargate maximum,
so it cannot simply be raised). Stored once, a 3D mosaic costs only a few GiB
more than the 2D one.

Existing 2D mosaics are unaffected and need no re-stitching — they carry no
`ndim` attribute, and the reader defaults to 2D.

### Docker images

The multi-stage Dockerfile produces two image variants:

```bash
# Headless: CLI tools + AWS CLI for Fargate batch processing
docker build --target headless -t cosmx-utilities:headless .

# GUI: Full Napari with Qt for interactive visualization
docker build --target gui -t cosmx-utilities:gui .
```

## Viewing processed slides on your own server

If you've received a processed CosMx dataset and want to run Napari locally instead of on the AWS-provisioned EC2 instance, the stitched output is self-contained — there are no AWS dependencies at view time.

**1. Copy the processed data from S3.**

Each slide directory contains `images/` (zarr pyramid), `_metadata.csv`, and `targets.hdf5`. Budget ~50–100 GB per slide; a full study can run several hundred GB. SSD storage is strongly recommended — spinning disks make tile loading painful.

```bash
aws s3 sync s3://<bucket>/napari-stitched/<study>/<experiment>/ /path/to/local/stitched/
# or, if rclone is configured:
rclone sync :s3:<bucket>/napari-stitched/<study>/<experiment>/ /path/to/local/stitched/ \
    --transfers 32 --checkers 16 --progress
```

You'll need AWS credentials with read access to the data bucket — ask the data owner.

**2. Install system dependencies (Linux headless server).**

PyQt6 needs Qt platform libraries that aren't installed by default on bare servers. On Ubuntu/Debian:

```bash
sudo apt install -y libxcb-cursor0 libxkbcommon-x11-0 libegl1 libgl1 libdbus-1-3 \
    libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 \
    libxcb-shape0 libxcb-xinerama0 libxcb-xkb1
```

Napari needs a graphical session. If your server is headless, use VNC, X-forwarding, or NICE DCV. A discrete GPU isn't required — software OpenGL (Mesa) works, just slower on large mosaics.

**3. Clone this repo and install with the GUI extra.**

The `napari-cosmx-fork` is a uv workspace member of this repo, and the top-level wrapper adds dependencies (dask, polars, matplotlib) that the Napari load path uses. Don't install the fork standalone.

```bash
# Install uv if needed: https://docs.astral.sh/uv/getting-started/installation/
git clone https://github.com/UW-BRaIN-lab/cosmx-utilities.git
cd cosmx-utilities
uv sync --extra gui
```

`uv` will provision a Python 3.10 interpreter automatically — the fork pins `>=3.8,<3.11`, so a system `python3.12` won't satisfy it directly.

**4. Launch Napari.**

```bash
uv run napari /path/to/local/stitched
```

**Resource notes.** Loading a stitched slide pulls the zarr pyramid and `targets.hdf5` into memory lazily, but interactive panning at high zoom is RAM-hungry. Plan for at least 64 GB free RAM; 128 GB is comfortable for the typical CosMx slide.

## Infrastructure setup

Fargate task definitions, IAM roles, and networking configuration are documented in [`fargate/FARGATE-SETUP.md`](./fargate/FARGATE-SETUP.md). Infrastructure IDs are stored in `fargate/.env` (gitignored) — copy `fargate/.env.example` to get started.

## Hyak cell-typing pipeline (in progress)

We are extending the pipeline onto the University of Washington's Hyak HPC cluster ([Klone](https://hyak.uw.edu/docs/)) for GPU-accelerated cell typing and batch correction, scheduled with [Slurm](https://slurm.schedmd.com/overview.html) and run from [Apptainer](https://apptainer.org/docs/user/main/) images (converted from the Docker build), with working storage on UW Kopah. The stages live in [`pipeline/`](./pipeline) and are shown in the diagram above:

- **Stage 1** converts each slide's CosMx flat files into an AnnData `.h5ad`.
- **Stage 3a** concatenates the cohort, applies per-cell QC, restricts to the study donors, and selects highly variable genes.
- **Stage 3b** computes a patient-level, batch-corrected quasi-Poisson Pearson-residual PCA with [scPearsonPCA](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/pearsonpca/), following the Bruker CosMx analysis vignette (the SVD is taken over the genes × genes cross-product, so the dense residual matrix is never formed).
- **Stage 3c** builds the neighbor graph, Leiden clustering, and UMAP on the GPU with [rapids-singlecell](https://rapids-singlecell.readthedocs.io/) (`gpu-l40s` partition), and emits UMAP QC plots for reviewing the batch correction.
- **Stage 4** types cells with [InSituType](https://github.com/Nanostring-Biostats/InSituType) (R) in semi-supervised mode against the CZI [GBmap](https://www.gbmap.org/) reference (level-4 annotation, restricted to the panel): cells are matched to named GBmap types while de novo clusters (letters a, b, c…) are discovered for off-reference tumour populations. Stage 4a emits the gene counts + per-cell negprobe mean, 4b runs the typing (gentle `rescale`-only reference update, `n_clusts` 10:20), and 4c writes `cell_type` back into the clustered AnnData — so the marker-heatmap/QC tooling re-renders by cell type with a one-flag flip.
- **Stage 5** infers copy-number variation with [InSituCNV](https://github.com/Moldia/InSituCNV-manuscript) (Moldia lab; [infercnvpy](https://infercnvpy.readthedocs.io/) + spatial neighbor smoothing, Python) to test whether the large transcriptionally-flat `Low_signal` fraction is actually tumour, via the GBM copy-number signature (chr7 gain / chr10 loss). CNV is orthogonal to expression typing. The spatial neighbor graph is built **per tissue block** (`slide × donor × block`) — each CosMx slide co-mounts two tissue blocks that are different regions (e.g. tumor bulk + contralateral uninvolved), often from the same donor, so a cohort-, slide-, or even donor-level graph would smooth counts (and subtract the diploid reference) across physically distinct tissues (and, for two-donor slides, across patients). Stage 5a (`86_insitucnv_prep.sh`) builds per-section inputs (raw gene counts + `obsm['spatial']` + autosomal gene positions), 5b (`87_insitucnv.sh`, a Slurm array) runs inferCNV per section against a conservative diploid reference (oligodendrocyte/neuron/immune/vascular), and 5c (`88_insitucnv_compare.sh`) aggregates per-group mean CNV profiles + cosine similarity + a control-calibrated `cnv_score`, with the known-malignant states as a chr7+/chr10− positive control and contralateral `Low_signal` as a negative control. Two CPU-only diagnostics reuse the per-section CNV outputs (no rerun): 5d (`89_insitucnv_edge_vs_bulk.sh`) checks whether the stronger edge-vs-bulk `Low_signal` signature is biological or a depth artifact, and 5e (`91_insitucnv_lowsignal_diagnostics.sh`) characterises `Low_signal` beyond the malignant/normal split — reference/negative-control integrity, edge-dilution at the margin, an area×count doublet screen, and a spatial margin map — and writes a reusable per-cell master table (`cell_cnv_table.csv.gz`). Two further stages resolve the malignant `Low_signal` cells: 5f/5g (`92_flat_posteriors.sh` → `93_insitucnv_hybrid_continuum.sh`) recover each cell's top-two profile posteriors (flat `insitutypeML` re-score) and score the Neftel-like malignant modules + a neuronal/NLGN3 program (`reference/gene_signatures.csv`) to separate a genuine cross-compartment **hybrid** state from a malignant **continuum**; and 5h (`94_rescue_lowsignal.sh` → `95_compare_rescue.sh`) reclusters the `Low_signal` pool and reassigns it (one more CosMx-native rescue iteration) to test whether the untyped fraction is still reference-limited (collapses) or an irreducible core.

> **Run note — long jobs and `ckpt`.** The Slurm scripts default to the preemptible `glioblastoma-ckpt`/`ckpt` allocation, which is fine for the short stages. The one exception is the `94` rescue (~12 h InSituType with no mid-run checkpoint): on `ckpt` it gets preempted and requeued from scratch. Run it (and, if you like, `95`) on the lab's dedicated, non-preemptible allocation instead — claiming one otherwise-idle L40S so the GPU partition accepts the CPU-only job:
> ```bash
> sbatch --account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1 \
>     pipeline/slurm/94_rescue_lowsignal.sh
> ```

Still planned: [scvi-tools](https://scvi-tools.org) integration and interactive Napari sessions via [Open OnDemand](https://www.openondemand.org).

Pipeline tools, Fargate infrastructure templates, and napari-cosmx-fork are publicly available in this repository and on [GHCR](https://github.com/UW-BRaIN-lab/cosmx-utilities/pkgs/container/cosmx-utilities).
