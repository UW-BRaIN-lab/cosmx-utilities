# Retina / brain / optic-nerve cell typing — run book

Semi-supervised InSituType typing of the retina study on Hyak Klone, using the same
stage-1 → stage-4 machinery as the GBM cohort. This file is the run sequence; the
reference itself (and an important caveat about it) is documented in
`pipeline/reference/README.md`.

## Study design

12 slides, **one donor per slide**. Every slide carries a cross-section of **retina +
optic nerve + brain**, so the reference is a *combined* ocular/brain atlas
(`retina_combined_panel.csv`, 5,993 genes × 64 types) rather than a single-tissue one.
Duplicate source types were reviewed by a count-discrimination test: oligodendrocyte and
OPC were merged (only 10 and 7 of 400 expected counts distinguish the copies), the rest
compete, and the Allen `Splatter`/`Miscellaneous` QC categories were dropped as
generic-neuron attractors. See `pipeline/reference/README.md`.

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
BATCH_VARIABLE=slide_id
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

Mirror the flat files S3 → Kopah (skips the ~640 MB/slide `_tx_file` by default). No
`--manifest` flag needed — `build_manifest.py --output` and `migrate_s3_to_kopah.py
--manifest` both default to `$MANIFEST`, which `.env` sets:

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

The hierarchy JSON is read from the repo, not Kopah, so it needs no upload.

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

## Stage 3b/3c — batch-corrected PCA + Leiden

**Run these.** InSituType does not read the PCA, so it is tempting to skip them, but the
cost is larger than it looks:

- `40_cluster.sh` (3c) stages 3b's `embedding.h5`, so skipping 3b skips 3c as well — no
  Leiden, no UMAP, no `cosmx_clustered.h5ad`.
- `90_write_celltypes.sh` stages `cosmx_clustered.h5ad`, so without 3b/3c there is **no
  `cosmx_typed.h5ad`** at the end — nothing with `cell_type` in `obs` for Napari or any
  spatial work. (`81_diagnose_typing.sh` deliberately reads `combined_qc.h5ad` plus the
  result h5 instead, so the diagnostics do not depend on this — but that routes around the
  gap rather than closing it.)
- What they provide is the **reference-free** evidence this study most needs: a Leiden ×
  `cell_type` crosstab (`85d`) checks that the typing respects unsupervised structure,
  which matters precisely because a 64-type collinear reference leaves the EM unconfident.
  A UMAP coloured by `Region` is also the most direct read on whether retina / optic nerve
  / brain separate at all, and `QC_COLOR=slide_id` shows whether the 12 donors mix.

`BATCH_VARIABLE` (3b, R) must match `BATCH_COL` (3a, Python) — both are set in `.env`
above so they cannot drift:

```bash
sbatch pipeline/slurm/30_pearson_pca.sh
```

```bash
sbatch pipeline/slurm/40_cluster.sh
```

⚠️ 3b needs `APPTAINER_SCPEARSON`, and `scpearsonpca.sif` was **purged** from
`/gscratch/glioblastoma/$USER/containers/` in the August 2026 cleanup. Check first, and
rebuild if absent — `scpearsonpca.def` is pinned (dated Posit snapshot + commit SHA), so
the rebuild is reproducible; the build command is in its header and takes ~20–40 min on a
login node.

Aside on method: scPearsonPCA exists rather than rapids-singlecell's
`normalize_pearson_residuals` because at ~8M cells the rapids path densifies an
n_cells × n_genes residual matrix that will not fit a 46GB L40S. At retina's ~840k × 2000
HVGs (~6.7GB dense) that constraint does not bind — but swapping methods would need new
code and would make retina-vs-GBM comparison murkier, so rebuild the pinned container.

## Stage 4a — InSituType input

```bash
sbatch pipeline/slurm/70_prep_insitutype.sh
```

Produces `cosmx-retina/stage4/insitutype_input.h5` (gene counts CSC + per-cell negprobe
mean). Everything below reads that one file.

## Stage 4a′ — gene selection (FAQ pruning)

The 6k FAQ recommends typing on a well-chosen 3,000–5,000 gene subset. `INPUT_KEY` is
what points 73/74 at a plain stage-4a input instead of the full-cohort anchor.

⚠️ **`73`/`74` default `STAGE4_DIR` to `stage4_anchor`**, not `stage4` like `70`/`80` — they
were written for the full-cohort anchor workflow. `INPUT_KEY` redirects only the *input*, so
their **outputs** still land under `stage4_anchor/`. Either set `STAGE4_DIR=stage4` in `.env`
so all four stages agree, or take the `stage4_anchor/` path when passing `KEEP_GENES` on to
`74` and `80`. Check with `s5cmd ls "s3://$KOPAH_BUCKET/$KOPAH_PREFIX/"` rather than
assuming.

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

**1–2 are one job.** `81_diagnose_typing.sh` reads the stage-4b result and the stage-3a
AnnData directly, so it runs the moment `80` finishes — no 3b/3c, no writeback needed:

```bash
HIERARCHY_BASENAME=retina_hierarchy.json sbatch pipeline/slurm/81_diagnose_typing.sh
```

1. **Spatial concordance** — modal `Region` and its share per type, plus the same
   aggregated by source atlas (HRCA types should concentrate in `Retina`, Allen in
   `Gray`/`White matter`, Mona in `Optic nerve`). It flags **spatially diffuse** types —
   modal share under 40% with ≥1000 cells — which is what a low-signal catch-all looks
   like, and the quantitative form of the neuropathologist's "cells inappropriately
   everywhere".
2. **Sibling adjudication** — for every terminal competition node in the hierarchy, the
   normalized mutual information of the sibling choice against `Region` versus against
   `qc_area` / `total_counts` / `qc_genes_detected` (quantile-binned to the same
   cardinality so the numbers are comparable). `Region` higher → real regional biology,
   keep the split. Technical higher → the assay offset sorting cells by segmentation
   quality, collapse that node to its parent. It also reports NEAT1 fold-change across
   siblings and which sibling carries it, since NEAT1 is the direct assay-offset probe.

   **Read both numbers, not just the region one.** A node whose siblings sort by cell area
   still yields a confident, spatially-structured-looking call, because segmentation
   quality itself varies by tissue. That is the failure mode that looks like success.
3. **De novo characterisation** — `pipeline/python/inspect_denovo_profiles.py` on the
   result's `/profiles`: per de novo cluster, top specific markers, nearest named type by
   cosine, and size. Expect at least one flat low-signal sink; the question is how much
   smaller it is than the 66–74% seen in the AtoMx Refit=FALSE calls.
4. **Marker heatmap** — `profile_marker_inputs.py` (compute on Klone) then
   `marker_heatmap.R` (render on the Mac — the InSituType SIF has no ComplexHeatmap).

## Stage 4c — write cell types back

```bash
sbatch pipeline/slurm/90_write_celltypes.sh
```

Produces `cosmx-retina/stage4/cosmx_typed.h5ad` — `combined_qc` obs plus `cell_type` and
`insitutype_prob` — which is the artifact downstream spatial work and Napari re-colouring
read. Requires 3b/3c (it stages `cosmx_clustered.h5ad`).

## Stage 4d — InSituTree

`pipeline/reference/retina_hierarchy.json` is already built: 11 top-level branches over
the 64 named types, depth 4. Two things must happen before `85_insitutree.sh` can run:

1. **CosMx-scale profiles.** InSituTree does no platform rescaling, so its `--reference`
   must come from the `RESCALE=true` InSituType run above, via
   `python/prep_insitutree_profiles.py` — not from `retina_combined_panel.csv`.
2. **De novo leaves added in lockstep.** Every hierarchy leaf must be a column of that
   profile matrix and vice versa, so the kept de novo clusters from triage need leaves
   (a `Low_signal` branch for the sink, and whichever real programs survive).

Then:

```bash
TREE_PROFILES_BASENAME=retina_insitutree_profiles.csv     HIERARCHY_BASENAME=retina_hierarchy.json     sbatch --account=glioblastoma --partition=gpu-l40s --qos=gpu-l40s     pipeline/slurm/85_insitutree.sh
```

What the hierarchy is *for* is worth being precise about: it does not decide that the
three astrocyte copies are one thing. It makes the broad call first (robust at 400 counts),
then defers the hard call to its own node where you can read the posterior — and a
near-uniform posterior there is itself the answer. See
`pipeline/reference/README.md`.

Note the **regional astrocyte question is not answered by the typing at all**. Call
astrocytes as one population, then test differential expression across `Region` within
them: that uses our own tissue instead of three labs' protocols, has real spatial ground
truth, and can discover regional programs rather than being limited to the atlases.

## ⚠ KNOWN CONTAMINATION — must be revisited before anything is published

**`denovo_fibroblast`'s profile carries ~9% ciliary body and cornea.** Accepted deliberately
on 2026-11-04 to get an exploratory InSituTree run out; it is not a finished result.

The SME excluded 57 FOVs on `uwa7761eyes` (2–15, 24–33, 36–68) as ciliary body / cornea,
tissue no type in the HRCA + Monavarfeshani + Allen reference covers. Stage 1 predates that
exclusion, so those cells went through QC and typing carrying the stale AtoMx `Region`
(which labels them `Retina`). Measured against the stitched Napari metadata:

| | |
|---|---|
| cells in those 57 FOVs | 15,384 (14,433 typed) |
| of them called `a` | 64.5% |
| `a` cells inside those FOVs | 9,304 — **9.1% of `a` cohort-wide**, 72% of `a` on that slide |

`denovo_fibroblast` is derived from `a`, so roughly a tenth of its profile is that tissue.
Because the profile is frozen into the InSituTree reference, the contamination then competes
at the TOP node against `/Fibroblast`, `/Pigmented` and `/Vascular/Mural` for the whole
cohort — not only on that donor.

**Why no flag fixes it.** `prep_insitutree_profiles.py` never reads cells; it takes the
already-aggregated `/profiles` matrix out of the InSituType result, where the mean over `a`
was computed when InSituType ran. `--fov-annotations` on stage 70 keeps those cells out of
the InSituTree *run*, which is worth doing on its own — InSituTree has no de novo escape
hatch and would otherwise force all 14,433 into some leaf, most likely the Fibroblast and
Pigmented leaves whose counts are the agreed test of whether `denovo_fibroblast` is right.
But the *profile* is upstream of both.

**To clear it, one of:**

1. Re-run stages 1 → 4 with `--fov-annotations` + `--exclude-regions EXCLUDE_not_in_atlas`.
   Clean, and necessary anyway if `uwa7761eyes` is dropped entirely. Note that de novo
   letters are arbitrary per run — every letter-keyed decision and table has to be re-mapped
   by profile cosine afterwards.
2. Recompute the `denovo_fibroblast` column from cells, filtered, instead of from
   `/profiles`. Gate it on the scale check this file's sibling script warns about:
   recompute `a` UNFILTERED first and confirm it reproduces `/profiles['a']`, since every
   column of the InSituTree matrix must share one scaling or the per-branch likelihoods are
   incomparable.

### What the 2026-11-04 re-run with the exclusion showed

Stage 70 dropped **14,433 cells (1.7% of 849,510), leaving 835,077**, and `82` was re-run
against the filtered cohort. Two independent confirmations that the exclusion hit the right
cells: `a` fell 102,468 → 93,164, a drop of **9,304**, matching the count measured separately
from the stitched Napari metadata; and the cohort's dominant region flipped from Retina 23.3%
to Gray matter 23.5%, which is what must happen when 14,433 cells wrongly labelled Retina go.

**The ciliary cells were spread across `a`, not concentrated.** Leiden component 11 (3,960
cells, 581 median counts — the deepest) vanished entirely; component 5 lost 3,224; every other
component lost a few hundred. Component 22 lost only 48 of 6,338 and is **not** ciliary.

**`a` is a low-signal population almost throughout.** Against a cohort median of 369 counts,
after the exclusion exactly ONE of its nine components sits above cohort depth:

| component | % of `a` | median counts | vs cohort | nearest (cosine) | runner-up gap |
|---|---|---|---|---|---|
| 5 | 51.2 | 190 | 51% | Mona_fibroblast (0.670) | 0.012 |
| 4 | 9.9 | 218 | 59% | Mona_fibroblast (0.771) | 0.004 |
| 0 | 9.7 | 101 | 27% | Mona_fibroblast (0.715) | 0.006 |
| 3 | 7.1 | 236 | 64% | Mona_pericyte (0.666) | 0.007 |
| **22** | 6.8 | **447** | **121%** | Mona_fibroblast (0.635) | 0.004 |
| 6 | 6.4 | 234 | 63% | HRCA_microglial (0.580) | 0.000 |
| 2 | 3.1 | 209 | 57% | HRCA_microglial (0.688) | 0.009 |
| 1 | 2.9 | 256 | 69% | Mona_endothelium (0.689) | 0.006 |
| 7 | 1.0 | 124 | 34% | Mona_fibroblast (0.708) | 0.006 |

**Do not read the `nearest` column as a branch assignment for `a`.** Every runner-up gap above
is ≤0.012, and in each case the two types sit under DIFFERENT top-level branches (fibroblast vs
pericyte vs microglia). That is the top of a flat distribution, not a discrimination — which is
what cells at 27–69% of cohort depth should produce. A tie WITHIN a branch is fine and expected,
since the reference is collinear there by design (component 3's pericyte-vs-smooth-muscle,
component 12 of `o`'s H1-vs-H2 horizontal).

**Consequence: `denovo_fibroblast` was NOT built this round.** A leaf built from `a`'s
fibroblast-nearest components would be dominated by component 5 — half of `a`, at half cohort
depth — on a 0.012 margin over pericyte. `a` gets no leaf until this is resolved; see the
component 22 question below.

### Component 22 — a candidate the decomposition surfaced

6,290 cells, 6.8% of `a`, and the only component at or above cohort depth (447 counts, 331
genes vs 369/288). Region-coherent at 84.8% retina, and 91.8% of its Leiden cluster is `a`, so
it is a tight population rather than a slice of a larger one. Its best reference match is only
0.635, with pericyte 0.004 behind.

**Colouring it in Napari.** `pipeline/reference/napari_splits_retina.csv` breaks it out as its
own value so the neuropathologist can see it rather than read about it:

```
label,column,values,new_name
a,leiden,22,a_deep_retinal
```

Pass it to the metadata writer with `--split-spec`, and point `--combined-h5ad` at
**`cosmx_clustered.h5ad`**, not `combined_qc.h5ad` — the split needs `obs['leiden']`, which
only the clustered file carries:

```
uv run python scripts/celltypes-to-napari-metadata.py \
    --typing-h5 insitutree_result.h5 \
    --combined-h5ad cosmx_clustered.h5ad \
    --split-spec pipeline/reference/napari_splits_retina.csv \
    --column insitutree --out-dir retina_napari_insitutree \
    --merge-into retina_napari_metadata_merged
```

This rides along with the single deferred `_metadata.csv` regeneration that will also carry
the InSituTree calls and the corrected `uwa7761eyes` Region — one new file, not three.

A split that matches no cell is a FATAL error, not a warning. The whole point is to make a
population visible; one that silently moves nothing leaves it indistinguishable from the
cluster around it, and the output looks entirely normal to whoever opens it.

Deep AND matching nothing is the combination that separates a population the atlases lack from
a flat low-signal pool: a faint cluster matches nothing BECAUSE it is faint, and this one is not
faint. It is the strongest candidate in `a` for the sclera / choroidal stroma that §02 of the
review page originally asked the neuropathologist about. Not named or given a leaf without her
— flagged for the next review.

**Read any InSituTree result from this run with the 57 FOVs filtered out at analysis time**
— one slide, a known FOV list — or the Fibroblast and Pigmented leaf counts will look
inflated for a reason that has nothing to do with the hierarchy working.

Related and also open: correcting `Region` in the pipeline at all requires re-running stage 1,
and `--exclude-regions` drops those cells *before* QC stats, HVG and PCA, which changes the
clustering and the fit. So it has to happen before InSituTree, not after.

## Proposed: de novo leaves in the hierarchy (PREPARED, awaiting neuropathology sign-off)

The 2026-11-04 sign-off dissolved ten de novo letters on the premise that the hierarchy would
absorb them into the named types they resembled. The first InSituTree run shows it does not,
for a reason that was measurable beforehand: **those branches' named profiles never fit this
data.**

| branch | its named leaves won | the de novo letter won | ratio |
|---|---|---|---|
| Retinal_neuron/Ganglion/Other_ganglion | 132 | `n` 15,039 | **113.9x** |
| Cortical_neuron/Excitatory | 244 | `c,i,k` 23,352 | **95.7x** |
| Cortical_neuron/Inhibitory | 222 | `r,q` 8,483 | **38.2x** |
| Vascular/Endothelial | 2,224 | `j` 44,513 | **20.0x** |
| Myeloid | 2,821 | `f` 48,058 | **17.0x** |
| Oligodendrocyte_lineage | 15,723 | `h,d` 169,620 | **10.8x** |
| *Astrocyte — already has de novo leaves* | *6,052* | *`m,s` 147,183* | *24.3x* |

Counts are from the flat InSituType run, where all 83 types competed freely. Forcing cells to
choose among profiles that lost by these margins does not produce a cautious answer — it
produces a confident wrong one: `Allen_Deep-layer near-projecting`, a rare L5/L6 subset, took
**99.5%** of its branch at a posterior of 1.000.

**Do not use the posterior to catch this.** Every contest in that run returned a median
posterior of 1.000 at the theoretical maximum against chance. A multinomial likelihood over a
few hundred counts and thousands of genes saturates, so the posterior is 1.0 whether the call
is right or wrong — including through any `--min-prob` filter. The flat-run count is the test
that works; `86_insitutree_resolution.sh` with `FLAT_RESULT_BASENAME` computes it.

**The control is in the same run.** `Astrocyte` is this proposal already implemented, and it is
the healthiest node in the tree: 55.5 / 35.5 between its two de novo leaves, both populations
surviving. At 24.3x it sits in the middle of the table, so there is no line that includes
Excitatory and excludes Oligodendrocyte_lineage.

**Prepared, not applied.** The current `retina_stage4.csv` and `retina_hierarchy.json` are left
alone so the 66-leaf run stays reproducible. The proposal is a parallel pair:

```
DENOVO_ANNOTATIONS=retina_stage4_denovo_leaves.csv \
  HIERARCHY_BASENAME=retina_hierarchy_denovo_leaves.json \
  TREE_PROFILES_BASENAME=retina_insitutree_profiles_denovo_leaves.csv \
  sbatch pipeline/slurm/84_prep_insitutree_profiles.sh
```

then 85 with the same `HIERARCHY_BASENAME` and that `TREE_PROFILES_BASENAME`. **76 leaves,
76 columns, lockstep verified locally.** 12 letters kept, 7 dropped — `a` (low signal), `o`
(components span several branches), `g`/`b`/`l` (stress state), `e` (incoherent), `p` (9 cells).

**Naming:** annotation first, number last, matching `denovo_astrocyte_1`; numbered only where a
branch takes more than one, so a unique leaf is not `_1` with no sibling. The finer identity
(`r` is SST+, `d` is newly-formed oligodendrocyte) stays in the annotations CSV rather than in
a type name that every downstream table would inherit as a claim.

**This reverses part of her decision, so it is her call.** The justification is new information
— the tree demonstrably mis-assigned those cells — which is what warrants reopening it. One
caveat to put to her honestly: a de novo profile wins back the cells it was fit to, so this
re-derives the InSituType partition rather than testing it. The evidence that these are real
populations remains the markers and regional distribution, not their size in the tree.

## Gotchas

**`pipeline/.env` used to OVERWRITE variables you passed on the command line — fixed
2026-11-04, but know the shape of it.** Every stage did `set -a; source .env; set +a`, so
`FOO=bar sbatch stage.sh` lost to a `FOO=` line in `.env`. It cost a run: `.env` sets
`REFERENCE_BASENAME=retina_combined_panel.csv` for stage 4b, stage 85 read the same name for
its InSituTree profile matrix, and the fit started against the **atlas** reference. It did not
fail on a missing file — it ran, and died later with `hierarchy leaves absent from reference
profile columns`, an error pointing at the hierarchy, which was correct. The tell was in the
log: `reference: 5993 genes x 64 cell types` where 84 had written `3,349 genes x 66`.

All 60 stages now load `.env` through `slurm/_load_env.sh`, which keeps the caller's value and
**echoes every override it kept**, so a log records the values a run actually used rather than
the ones requested. 85's profile matrix is also `TREE_PROFILES_BASENAME` now, not
`REFERENCE_BASENAME`, because one name standing for both the atlas file and the InSituTree
profile matrix is what made the collision possible.

**The habit that would have caught it earlier:** `pipeline/.env.example` is committed and lists
every name `.env` defines. Before relying on an override, check whether that name is already
spoken for there — and read back the dimensions or counts a stage reports instead of assuming
it took your value.



- **R progress goes to `.err`, Python progress to `.out`.** The R stages (`30`, `73`, `74`,
  `80`, `85`) log with `message()`, which writes to stderr — `pearson_pca.R` for instance has
  four `message()` calls and no `print`/`cat` at all, so its `.out` holds only the wrapper's
  bash `echo`s. Grep `pipeline/logs/<stage>_<jobid>.err` for anything R said. (`73`/`74` also
  `cat` their summary files from bash at the end, so those land in `.out`.) The Python stages
  (`10`, `20`, `70`, `81`) use `print()` and behave as expected.
- **A silent R log is usually not a hang.** These scripts log at phase boundaries only, so
  `30` prints nothing between `counts: ...` and `writing embedding` while it computes
  residuals and the SVD. Compare elapsed time against the job's `--time` before worrying.
- **`--export-batch` is mandatory** when building this manifest, or you get 24 rows.
- **Raw flat files land OUTSIDE `KOPAH_PREFIX`.** `migrate_s3_to_kopah.py` sets
  `dst_key = src_key`, so they mirror the source path at
  `s3://brainlabkg/CosMx-retina/…` while only pipeline outputs go under
  `cosmx-retina/`. This looks wrong in the transfer log but is deliberate — it is
  what makes the flat files shareable between runs.
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
