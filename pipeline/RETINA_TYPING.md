# Retina / brain / optic-nerve cell typing — run book

Semi-supervised InSituType typing of the retina study on Hyak Klone, using the same
stage-1 → stage-4 machinery as the GBM cohort. This file is the run sequence; the
reference itself (and an important caveat about it) is documented in
`pipeline/reference/README.md`.

## Study design

12 slides, **one donor per slide**. Every slide carries a cross-section of **retina +
optic nerve + brain**, so the reference is a *combined* ocular/brain atlas
(`retina_combined_panel.csv`, 5,993 genes × 68 types) rather than a single-tissue one.

The flat-file metadata has a per-cell **`Region`** column — `Retina`, `Optic nerve`,
`Gray matter`, `White matter`, `Adjacent soft tissue` — which stage 1 carries into `obs`
for free. Two consequences worth holding onto:

- **`Region` is biology, not batch.** The PCA/QC batch is the donor, and since one donor
  == one slide, `BATCH_COL=slide_id`.
- **`Region` is the typing validation axis.** Retinal types must land in `Retina`, Allen
  cortical types in `Gray`/`White matter`, and Mona optic-nerve types in `Optic nerve`.
  A retinal rod call in cortex is a typing failure you can see without any ground truth.

Per-region sequencing depth differs a lot, and the optic nerve is the shallowest —
median counts / genes per cell on `eyes7517`:

| Region | cells | median counts | median genes |
|---|---|---|---|
| White matter | 29,834 | 420 | 319 |
| Optic nerve | 21,002 | **318** | 248 |
| Retina | 8,428 | 630 | 460 |
| Gray matter | 6,040 | 693 | 499 |
| Adjacent soft tissue | 1,032 | 153 | 122 |

Scale: ~66k cells/slide → **~900k cells**, ~840k after QC. That is well under the 2.33M
Wenyu cohort that `80_insitutype.sh` typed in one job, so this study needs **no anchor
stage and no sharding** — a single flat fit is fine.

## Setup (once)

Run from a **separate Klone checkout with its own `pipeline/.env`** — the slurm scripts
`source pipeline/.env` with `set -a`, so anything defined there beats a submit-time
`VAR=... sbatch ...`. This is the same isolation the full-cohort run used.

```bash
cd ~
git clone git@github.com:UW-BRaIN-lab/cosmx-utilities.git cosmx-utilities-retina
cd cosmx-utilities-retina
cp pipeline/.env.example pipeline/.env
```

Then set in `pipeline/.env` (on top of the Kopah keys and `APPTAINER_*` paths):

```ini
SOURCE_S3_PREFIX=CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26
KOPAH_PREFIX=cosmx-retina
MANIFEST=/mmfs1/home/emilyek/cosmx-utilities-retina/pipeline/manifest_retina.csv
BATCH_COL=slide_id
COHORT=
MIN_GENE_COUNTS=50
MAX_AREA=50000
MAX_NEGPROBE_PROP=0.1
REFERENCE_BASENAME=retina_combined_panel.csv
QC_COLOR=slide_id,Region,leiden
```

`MAX_AREA=50000` is the one QC threshold that differs from GBM: retinal/optic-nerve glia
are genuinely large (median 10,086 px), and GBM's `30000` discards 3.6% of cells against
`50000`'s 0.9%, while still catching merged-segment doublets.

Kopah layout, following the existing convention (raw flat files sit outside
`KOPAH_PREFIX` because `migrate_s3_to_kopah.py` mirrors the source key):

| what | Kopah location |
|---|---|
| raw flat files | `s3://brainlabkg/CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26/…` |
| pipeline outputs | `s3://brainlabkg/cosmx-retina/{anndata,stage3,stage4}/` |
| reference | `s3://brainlabkg/cosmx-retina/reference/retina_combined_panel.csv` |

## Stage 0 — migrate + stage the reference

Build the manifest (12 rows). Manifests are generated, not committed — note the
`--export-batch` filter, without which the prefix yields each slide **twice** (it holds an
initial segmentation and a resegmentation side by side, plus a norefit export):

```bash
uv run python pipeline/python/build_manifest.py \
    --source-bucket keene-cosmx-data \
    --source-prefix CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26 \
    --export-batch Resegmentationcosmxretinabrain22626_01_04_2026_15_10_18_504 \
    --output pipeline/manifest_retina.csv
```

Mirror the flat files S3 → Kopah (skips the ~640 MB/slide `_tx_file` by default):

```bash
uv run python pipeline/python/migrate_s3_to_kopah.py --dry-run
```

```bash
uv run python pipeline/python/migrate_s3_to_kopah.py
```

Stage the reference:

```bash
set -a; source pipeline/.env; set +a
AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID" \
AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY" \
S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL" \
    s5cmd cp pipeline/reference/retina_combined_panel.csv \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/retina_combined_panel.csv"
```

## Stage 1 — flat files → per-slide AnnData

12 slides, so override the array (the `#SBATCH --array=1-57` in the file is for GBM):

```bash
sbatch --array=1-12 pipeline/slurm/10_flatfiles_to_anndata.sh
```

## Stage 3a — concatenate + QC

```bash
sbatch pipeline/slurm/20_concat_qc.sh
```

Check in the log: `Cohort filter DISABLED`, ~840k of ~900k cells kept, and **12
`slide_id` batches**.

Stages 3b/3c (`30_pearson_pca.sh`, `40_cluster.sh`) are **not** needed for typing. Run
them if you want the Leiden clusters for the post-typing crosstab and the `Region` QC
UMAPs — that is also the cheapest early read on whether the three tissues separate.

## Stage 4a — InSituType input

```bash
sbatch pipeline/slurm/70_prep_insitutype.sh
```

Produces `cosmx-retina/stage4/insitutype_input.h5` (gene counts CSC + per-cell negprobe
mean). Everything below reads that one file.

## Stage 4a′ — gene selection (FAQ pruning)

The 6k FAQ recommends typing on a well-chosen 3,000–5,000 gene subset. `INPUT_KEY` is
what points 73/74 at a plain stage-4a input instead of the full-cohort anchor:

```bash
INPUT_KEY=stage4/insitutype_input.h5 sbatch pipeline/slurm/73_select_genes.sh
```

Read `gene_selection_summary.txt` in the job output. If the kept union is outside
3,000–5,000, sweep and resubmit — e.g. to keep fewer genes:

```bash
INPUT_KEY=stage4/insitutype_input.h5 BG_MULT=5 REF_QUANTILE=0.6 \
    sbatch pipeline/slurm/73_select_genes.sh
```

Expect a *different* number than the ~4,375 from the earlier flat-file diagnostic: that
used max-enrichment across data clusters, this uses the FAQ union (above-background OR
informative in ≥1 reference profile).

## Stage 4a″ — de novo K sweep

`insitutype()` picks K by minimum AIC on only a **2,000-cell** subset, which is how the
GBM pilot's `10:20` returned 20 — the ceiling, not an optimum. Sweep it explicitly, and
run it **twice** to see whether pruning moves K (it did for GBM: 22 → 27, because the
per-cluster AIC penalty is `2 × n_genes`):

```bash
INPUT_KEY=stage4/insitutype_input.h5 sbatch pipeline/slurm/74_choose_k.sh
```

```bash
INPUT_KEY=stage4/insitutype_input.h5 \
    KEEP_GENES=stage4/gene_selection/kept_genes.txt \
    K_SWEEP_OUT=k_sweep_pruned \
    sbatch pipeline/slurm/74_choose_k.sh
```

Confirm in `k_sweep_summary.txt` that the minimum is **interior** (not `CENSORED` at the
range edge); widen `K_SWEEP_RANGE` if it isn't. Defaults are the validated fast ones
(`K_SWEEP_SUBSET=10000`, `K_SWEEP_REPS=2`) — the 20000/3 combination takes ~11 h and
dies at walltime without writing anything.

## Stage 4b — semi-supervised InSituType

Pass the chosen K as a **single value, not a range**: a range makes `insitutype()` re-run
its own 2,000-cell K selection, which is the flaw the sweep exists to avoid.

```bash
N_CLUSTS=<K from the pruned sweep> \
    KEEP_GENES=stage4/gene_selection/kept_genes.txt \
    sbatch pipeline/slurm/80_insitutype.sh
```

`RESCALE=true` / `REFIT=false` are the defaults and are what we want: rescale gives the
gentle platform correction (and the `$profiles` that InSituTree will later need), refit
is the aggressive variant we avoid.

## After the fit — triage before anything else

1. **Region concordance** — cross-tabulate `cell_type` × `Region`. This is the strongest
   available check, and it needs no ground truth.
2. **De novo characterisation** — `pipeline/python/inspect_denovo_profiles.py` on the
   result's `/profiles`: per de novo cluster, top specific markers, nearest named type by
   cosine, and size. Expect at least one flat low-signal sink; the question is how much
   smaller it is than the 66–74% seen in the AtoMx Refit=FALSE calls.
3. **Marker heatmap** — `profile_marker_inputs.py` (compute on Klone) then
   `marker_heatmap.R` (render on the Mac — the InSituType SIF has no ComplexHeatmap).

Then, and only then, design the InSituTree hierarchy. It is the real endpoint for this
study: the combined reference has the same cell type once per source atlas at cosine
0.96–0.99, plus fine types that are unresolvable at 6k depth, and branch-local
competition is what makes those tractable. See the collinearity section of
`pipeline/reference/README.md`.

## Gotchas

- **`--export-batch` is mandatory** when building this manifest, or you get 24 rows.
- **`obs` is wide.** This study's flat-file metadata carries ~190 columns, including
  ~100 AtoMx neighborhood-composition floats and AtoMx's own
  `RNA_RNA_Cell.Typing.InSituType.1_1_clusters` call. All of it lands in `obs`. That is a
  feature — the AtoMx call is what lets `scripts/compare-cell-typing.py` diff our R typing
  against theirs — but it does inflate the per-slide `.h5ad` and stage-3a memory.
- **`.env` beats the submit environment** for any variable it defines. `INPUT_KEY`,
  `KEEP_GENES`, `N_CLUSTS` and `K_SWEEP_OUT` are deliberately left unset in `.env` so
  they can be set per submission.
- **`MANIFEST` must be an absolute path** — slurm copies the batch script out of the
  repo, so a relative path resolves against the spool dir.
- Long non-checkpointing fits belong on the dedicated `gpu-l40s` slice
  (`--account=glioblastoma --partition=gpu-l40s --qos=normal`), not `ckpt`, or preemption
  requeues them from scratch. At ~840k cells the `80` job should be short enough for
  `ckpt`; escalate if it thrashes.
