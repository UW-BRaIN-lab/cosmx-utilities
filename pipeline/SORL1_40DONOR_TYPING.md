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

## Offline copies (54, export_summary_bundle.sh)

`RUNS="k10:stage4_k10 k21:stage4_k21 k30:stage4_k30" sbatch pipeline/slurm/54_export_cell_labels.sh`
writes `stage3/labels/cell_labels.csv.gz`: one row per cell (slide, donor, diagnosis, Leiden,
UMAP, depth) with `celltype_<run>` and `prob_<run>` for every run. `bash pipeline/export_summary_bundle.sh`
packs every small output (figures, tables, compact results, logs; not the `.rds` or typed AnnData)
into one tarball on Klone scratch.

## Marker lists (top_markers_from_profiles.py)

`python/top_markers_from_profiles.py --profiles-h5 insitutype_result.h5 --names names.csv` writes the top 20
genes for every type (named and de novo) from the fitted profiles, ranked by log2 fold-change against the other
types, with an expression floor (`--min-profile 0.2`). Fold-enrichment with no floor lets near-zero genes top the
list (the inspector's "markers" for the flat clusters were of that kind). A cluster with no gene above 2-fold is
either flat (every fold-change negative, K=21 cluster `i`) or shares its genes with sibling clusters (`a`, `q`);
the most-expressed genes tell them apart.

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

```bash
UPDATE_REFERENCE=false RESCALE=false REFIT=false sbatch pipeline/slurm/80_insitutype.sh   # + N_CLUSTS=<K>, STAGE4_DIR, INPUT_DIR, KEEP_GENES
sbatch pipeline/slurm/90_write_celltypes.sh      # then 85d (Leiden x type) and 51 (QC by type)
```

`81_diagnose_typing.sh` tests the call against `Region`, which is the donor id here, so its
spatial-concordance test is uninformative for this study; use the 85d cross-tab instead.

Go to InSituTree (85) only if the diagnostics show **(a)** a large flat/low-confidence sink,
or **(b)** unresolved sibling confusion (Astrocyte A/B, Microglia A/B, Inhibitory A/B/C,
L2/3-L4-L6) — its value is deferring those calls with a readable posterior. It needs a
hierarchy JSON; Brain_6k ships `Brain_6k.celltypeslist.R` as the starting tree.

## InSituTree: profiles and hierarchy (K=21)

InSituTree is fully supervised, so the new populations found by InSituType have to be turned into a reference
first: the supported route is semi-supervised InSituType, then a reference built from its fitted profiles, then
InSituTree with those profiles and a hierarchy. The profiles must come from ONE `rescale=TRUE` fit; the K=21
calibrated run (`stage4_k21`) is that.

- `reference/denovo_annotations/stage4_k21.csv`: `keep`/`new_name` per cluster. 8 kept: **6 new cell types**
  (`c`, `g`, `s`, `m`, `p` excitatory; `l` VIP interneuron) plus the two flat clusters (`i`, `u`). 13 not kept.
- **Rule for keeping a de novo leaf:** it must differ from the leaves it competes with in its branch, at least 10
  genes at 2-fold over its siblings. Marker lists against all 34 types are misleading here, because a cluster's
  siblings share its genes and a no-floor fold-change lets near-zero genes top the list. Failed the test and
  dropped: `q`, `h`, `r`, `d` (0, 0, 1, 4 genes versus the other excitatory leaves); `n` (0 versus reference OPC);
  `o` (5 genes, all neuronal: ambient RNA); `t` (15 genes, stress and blood genes: HSPA1A/B, NEAT1, HBA1/2).
  `l` has only 4 (CALB2, VIP, CNR1, RGS12) but they are the established VIP markers; kept by judgment.
  Dropped clusters' cells fall to the reference type or nearest leaf. `a`, `j`, `f`, `b`, `k`, `e` were excluded
  earlier (hybrid profile, donor-skewed, or donor-specific stress states).
- `reference/sorl1_hierarchy.json`: Brain_6k's tree with `Glial` moved to the top level (the source file nests it inside
  `Neuronal`; the K=21 profiles separate glia from neurons). The excitatory leaves are one flat branch (cosine
  does not support sub-groups: within 0.941 vs between 0.924). `Low_signal_neuronal_denovo` (cosine 0.93 to
  an excitatory leaf) sits in the Neuron branch so it competes with the neurons it resembles; the generic flat
  cluster is the top-level `Low_signal` sink. Checked with InSituTree's own `plotProfileSimilarity()` (plots in
  `sorl1_summary_for_maddie/insitutree/`).
- `reference/sorl1_insitutree_profiles.csv`: the matching 3,323-gene x 21-type matrix.
- `python/tests/test_sorl1_hierarchy.py` checks that leaves and matrix columns match both ways and that only the
  distinct de novo leaves are kept.

**How InSituTree picks genes (read from `supervisedSubcluster.R`).** At each node it takes the children's profiles
(an internal child is the mean of its leaves), normalises each to sum 100, and keeps genes whose range across
the children exceeds the 0.5 quantile on both absolute and percentage difference. That is a loose filter
(about 1,000 to 1,350 of the 3,323 genes pass at each node here), not a per-leaf "maximum contrast". The marker genes
of every kept leaf pass at their node (PCP4, RORB, BCL11B, RIT2, SEMA3E, HS3ST4, SYNPR, NR4A2; VIP, CALB2, CNR1, RGS12).
Near-identical siblings are the case to avoid: at an [Oligodendrocyte, Oligodendrocyte_denovo] node the
selection includes the real oligodendrocyte genes and the neuronal ambient genes, so the split between the two
would reflect expression level and contamination, not a cell type.

```bash
# on Klone: rebuild the matrix from the fit (only if the K=21 result changes), upload it
s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage4_k21/insitutype_result.h5 .
uv run --no-project --python 3.11 --with h5py --with numpy --with pandas python pipeline/python/prep_insitutree_profiles.py \
    --profiles-h5 insitutype_result.h5 --annotations pipeline/reference/denovo_annotations/stage4_k21.csv \
    --output pipeline/reference/sorl1_insitutree_profiles.csv
s5cmd cp pipeline/reference/sorl1_insitutree_profiles.csv s3://$KOPAH_BUCKET/$KOPAH_PREFIX/reference/sorl1_insitutree_profiles.csv
sed -i '/^REFERENCE_BASENAME=/d' pipeline/.env   # .env beats submit-time values; it points at the raw 13-type reference
INPUT_DIR=stage4 STAGE4_DIR=stage4_insitutree_k21 REFERENCE_BASENAME=sorl1_insitutree_profiles.csv \
    HIERARCHY_BASENAME=sorl1_hierarchy.json sbatch pipeline/slurm/85_insitutree.sh
```

The R job prints the tree (`printTree`) and stops if a leaf is missing from the matrix. Check its log for
`hierarchy: 21 leaves over the reference's 21 types`. Calls are read from the finest resolved level; low
`prob` marks the cells InSituTree could not place within a branch. The `Low_signal` sink should be compared
across carriers and non-carriers before interpreting it (the flat clusters are carrier-enriched).

## Resolved items

- Slide names `6562 A7 576 …` and `…6802 A6` carry a typo made when the slides were
  entered; the real donors are **526** and **6803**, as in the annotation sheet. The slide
  names are export directory names and stay as they are, so never infer donors from them.
