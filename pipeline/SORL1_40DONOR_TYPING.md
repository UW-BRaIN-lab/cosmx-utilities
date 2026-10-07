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
- **Annotations are blank in the flat files for 9 of 10 slides** (only `7163_…_7796` carries
  `UWA`/`Case_braod`/`Case_specific`/`SORL1_mutation`; note the exported column is spelled
  `Case_braod`). Donor therefore comes from the per-FOV table
  `reference/fov_annotations_sorl1.csv` via `FOV_ANNOTATIONS`, which sets `Region` = donor id.
  It currently covers the 4 pilot-named slides only; **the 6 new slides need Maddie's
  FOV → donor → diagnosis table before stage 1** (see Open items).
- Donor 741 was mislabelled as a SORL1 carrier; corrected to AD+LATE (38 FOVs, slide
  `…599_657_710_741` FOV 1–38) in `fov_annotations_sorl1.csv`.

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

```bash
sbatch pipeline/slurm/80_insitutype.sh     # N_CLUSTS=<chosen K>, RESCALE=true, REFIT=false
sbatch pipeline/slurm/81_diagnose_typing.sh
sbatch pipeline/slurm/90_write_celltypes.sh
```

Go to InSituTree (85) only if the diagnostics show **(a)** a large flat/low-confidence sink,
or **(b)** unresolved sibling confusion (Astrocyte A/B, Microglia A/B, Inhibitory A/B/C,
L2/3-L4-L6) — its value is deferring those calls with a readable posterior. It needs a
hierarchy JSON; Brain_6k ships `Brain_6k.celltypeslist.R` as the starting tree.

## Open items

1. Maddie: FOV → donor → diagnosis (case_broad, case_group, SORL1 variant) for the 6 new slides.
2. Maddie: confirm her source annotation CSVs for the pilot slides carry the 741 fix (the
   committed table is corrected, but regenerating from an old CSV would undo it).
3. `7068` FOV 100 has `case_group=AD+LATE` but `sorl1_mutation=AD+LATE SORL1 R953C`
   (rest of the donor says SORL1) — one-FOV inconsistency in the pilot table; confirm.
