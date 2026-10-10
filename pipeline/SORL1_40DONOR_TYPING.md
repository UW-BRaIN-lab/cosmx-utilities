# Maddie SORL1 3D cohort (40 donors) — run book

Semi-supervised InSituType on the Hyak pipeline against the CosMx-scale human brain
reference (`Brain_6k`, Bruker/Nanostring `CosMx-Cell-Profiles`, healthy frontal cortex),
then InSituTree only if the diagnostics call for it. Same machinery and `.env` pattern as
`RETINA_TYPING.md`; this file records what differs.

## Study

- Source: `s3://keene-cosmx-data/CosMx-Maddie/20261002_SORL1_3D_40donors_05_10_2026_17_07_58_724/`
  — **10 slides × 4 donors = 40 donors**, 3D segmentation, **862,476 cells** (one flat-file
  row per cell). No anchor stage and no sharding needed (retina fit ~900k in one job).
- 4 slides re-use the pilot's names (`20260708_UWA_{599_657_710_741, 6766_6966_7068_7137,
  7163_7181_7665_7796, 787_795_6589_6745}`, 16 donors, identical FOV grids). 6 slides are new
  (`6562A7…`, `6760A7…`, `699A7…`, `728A7…`, `7678A7…`, `7999A7…`, 24 donors).
- Panel: the 6,175-probe 6k Discovery panel, identical to `cosmx_6k_panel_genes.txt`.
- **Annotations are blank in the flat files for 9 of 10 slides** (only `7163_…_7796` is
  annotated; the exported column is spelled `Case_braod`). Donor therefore comes from the
  per-FOV table `reference/fov_annotations_sorl1.csv` via `FOV_ANNOTATIONS`, which sets
  `Region` = donor id. It is built from Maddie's sheet
  (`reference/sorl1_annotations_combined.csv`) and covers all 1,868 FOVs / 10 slides / 40
  donors; the builder fails if any slide or FOV is missing:
  `build_fov_annotations_sorl1.py --annotations-csv reference/sorl1_annotations_combined.csv
  --manifest <manifest> --fov-counts reference/sorl1_fov_counts.csv --output reference/fov_annotations_sorl1.csv`.
- Donors 741 (now AD+LATE) and 7068 FOV 100 (SORL1 R953C) are correct in her sheet.
- **Confound to remember downstream:** every SORL1 carrier (599, 6589, 7068, 7137, 7181,
  7796) is on one of the four pilot-named slides; the six new slides have none. SORL1 status
  is therefore partly confounded with slide / run.

## Reference

`reference/brain6k_panel.csv` = `Brain_6k.profiles.csv` restricted to the panel
(6,077 of 6,175 genes, 98.4%; 13 types). Regenerate with
`prep_insitutype_reference.py --reference-csv Brain_6k.profiles.csv
--panel-genes pipeline/reference/cosmx_6k_panel_genes.txt --output pipeline/reference/brain6k_panel.csv`.

Brain_6k has **no** fibroblast / pericyte / T-cell / reactive-astrocyte / disease-associated
microglia types, so AD pathology states must be recovered as de novo clusters. That is why K
is chosen from the data (stage 4a′) rather than fixed.

## `.env` (separate Klone checkout, like retina)

```ini
SOURCE_S3_PREFIX=CosMx-Maddie
KOPAH_PREFIX=cosmx-maddie-sorl1
MANIFEST=<checkout>/pipeline/manifest_sorl1_40donors.csv
BATCH_COL=Region            # = donor id (FOV table); no `Case` column in this study
BATCH_VARIABLE=Region       # keep in sync with BATCH_COL (3b, R)
COHORT=
MAX_AREA=50000              # see QC note
REFERENCE_BASENAME=brain6k_panel.csv
FOV_ANNOTATIONS=reference/fov_annotations_sorl1.csv
QC_COLOR=slide_id,Region,leiden
REGION_KEY=none             # marker heatmap (50): Region is the donor id, so no Region split
STAGE4_DIR=stage4
INPUT_KEY=stage4/insitutype_input.h5
```

QC note (measured on the raw metadata): `Area > 30000` drops 3.9% of cells (3D cells are
large; 99th percentile 39.9k px), as GBM's cut did on retina. `MAX_AREA=50000` is the same
choice made for retina — confirm the 3a log before accepting. Median counts per cell range
**207 – 1,372 across slides** (6.6×), so check slide 7 (`…787_795_6589_6745`, 207) after QC.

## Stage 0 — manifest + migration

```bash
uv run python pipeline/python/build_manifest.py \
    --source-bucket keene-cosmx-data --source-prefix CosMx-Maddie \
    --export-batch 20261002_SORL1_3D_40donors_05_10_2026_17_07_58_724 \
    --output pipeline/manifest_sorl1_40donors.csv
uv run python pipeline/python/migrate_s3_to_kopah.py --dry-run
uv run python pipeline/python/migrate_s3_to_kopah.py
```

Needs `SOURCE_S3_PREFIX` + the Kopah keys in `.env`, so run on Klone. Copies exprMat,
metadata, fov_positions and the 2D polygons per slide (~0.3 GB/slide); the `_Z###`
polygons and the ~2 GB `_tx_file` are skipped (stage 1 reads neither). Then stage the
reference to `s3://$KOPAH_BUCKET/$KOPAH_PREFIX/reference/` (command in `RETINA_TYPING.md`).

## Stages

```bash
sbatch --array=1-10 pipeline/slurm/10_flatfiles_to_anndata.sh
sbatch pipeline/slurm/20_concat_qc.sh      # expect: cohort filter DISABLED, 40 `Region` batches
sbatch pipeline/slurm/30_pearson_pca.sh
sbatch pipeline/slurm/40_cluster.sh
sbatch pipeline/slurm/70_prep_insitutype.sh
sbatch pipeline/slurm/73_select_genes.sh   # FAQ 3-5k gene band
```

## Batch-mixing figures (52)

`sbatch pipeline/slurm/52_umap_panels.sh` makes one UMAP panel per slide and per donor
(highlighted over all other cells in grey) plus `<key>_by_leiden_share.csv`. Use
`PANEL_NCOLS=8` for the 40 donor panels. The all-in-one UMAPs are overplotted and cannot show
mixing; these can.

## Depth confound table (53)

`sbatch pipeline/slurm/53_donor_depth.sh` writes `stage3/depth/donor_depth.csv`: per donor,
median counts / genes, the share of cells in the shallow clusters (median `total_counts` below
560, i.e. clusters 1, 3, 4, 5, 8, 13), and slide, run date, diagnosis and carrier status. The
SORL1 carriers all sit on the earlier-run pilot slides and carry the largest shallow share
(80% vs 47% for new-slide non-carriers), so a carrier-vs-non-carrier comparison needs depth
accounted for. The table is post-QC (50-count floor), so it understates how shallow the
lowest donors are.

## Choosing the number of de novo clusters (74)

`74_choose_k.sh` runs InSituType's AIC sweep; the pilot at 10:20 censored at the ceiling,
so sweep wide and check the curve **turns over inside the range**:

```bash
KEEP_GENES=stage4/gene_selection/kept_genes.txt K_SWEEP_RANGE=5:40 K_SWEEP_OUT=k_sweep_pruned \
    sbatch pipeline/slurm/74_choose_k.sh
```

The job does not checkpoint, so if ckpt preempts and requeues it, rerun on a non-preempted
partition (see the header of `74_choose_k.sh`; do not use the collaborator allocation).
AIC alone over-splits, so treat its optimum as an upper bound and confirm with a second,
reference-free criterion: run 80 at the AIC K and at a smaller K, then compare the Leiden ×
`cell_type` crosstab (85d) — keep the smallest K at which no further de novo cluster is a
distinct Leiden population. Note from retina: K was *not* the lever there (the anchoring
was), so do not over-invest in the sweep.

## Stage 4 + InSituTree decision

**`UPDATE_REFERENCE=false` for this study.** `Brain_6k` is already CosMx-derived, so there is no
scRNA-seq-to-CosMx platform shift to correct. InSituType's own documentation says
`update_reference_profiles` is TRUE by default, "however, if the reference profiles are from the
same platform as the study, then FALSE could be better". TRUE turns on anchor selection, and
`rescale`/`refit` only apply under it; the package then overwrites the anchor cells' labels with
the anchor assignments, ignoring the likelihoods (the "pinned cells" in CLAUDE.md). So run:
`UPDATE_REFERENCE=false RESCALE=false REFIT=false`, and check the log line
`insitutype: ... update_reference=FALSE rescale=FALSE refit=FALSE` before trusting a run.
The first fits (K=21, K=10, `stage4_k21`/`stage4_k10`, `rescale=TRUE`) and the second pair
(`*_norescale`, `rescale=FALSE` but still `update_reference=TRUE`) are superseded by
`stage4_k21_asis` / `stage4_k10_asis`. The first two were run with the retina defaults without
checking where the reference came from.

## Resolved items

- Slide names `6562 A7 576 …` and `…6802 A6` carry a typo made when the slides were
  entered; the real donors are **526** and **6803**, as in the annotation sheet. The slide
  names are export directory names and stay as they are, so never infer donors from them.
