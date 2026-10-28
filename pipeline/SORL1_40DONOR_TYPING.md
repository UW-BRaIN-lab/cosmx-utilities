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

**`UPDATE_REFERENCE=true RESCALE=true REFIT=false` for this study** (the "calibrated" setting;
the K=21 run is `stage4_k21`). `Brain_6k` is CosMx-derived, so we first tried using its profiles
as they are (`UPDATE_REFERENCE=false RESCALE=false REFIT=false`, `stage4_k21_asis` /
`stage4_k10_asis`). That is what InSituType's documentation suggests when the reference is from
the same platform as the study, but it failed here: about 1% of cells got a Brain_6k name at
K=21 and nearly everything went to de novo clusters. Do not use `_asis` as the baseline.
`UPDATE_REFERENCE=true` turns on anchor selection, so some cells are pinned to a reference type
and their labels ignore the likelihoods (the "pinned cells" in CLAUDE.md). `rescale` and
`refit` only apply under it, and `update=TRUE` with both `rescale` and `refit` FALSE is invalid
(the package stops), which is why the `*_norescale` pair is not a usable run either. Check the
log line `insitutype: ... update_reference=TRUE rescale=TRUE refit=FALSE` before trusting a run.

```bash
UPDATE_REFERENCE=true RESCALE=true REFIT=false sbatch pipeline/slurm/80_insitutype.sh   # + N_CLUSTS=<K>, STAGE4_DIR, INPUT_DIR, KEEP_GENES
sbatch pipeline/slurm/90_write_celltypes.sh      # then 85d (Leiden x type) and 51 (QC by type)
```

InSituType fits are slow and do not checkpoint: submit `80` on `gpu-l40s`, not `ckpt`
(`--account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1`), or preemption
restarts it from zero.

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

- `reference/denovo_annotations/stage4_k21.csv`: `keep`/`new_name` per cluster. 9 kept: **7 new cell types**
  (`c`, `g`, `s`, `m`, `p` excitatory; `l` VIP interneuron; `o` oligodendrocyte) plus the two flat
  clusters (`i`, `u`). 12 not kept.
- **v1 -> v2 -> v3, and the final choice (v2).** v1 (`stage4_insitutree_k21`) had neither `o` nor `n`: the OPC leaf held
  62,887 cells, only 22% real OPC (Leiden 0), 24% oligodendrocytes (Leiden 2), 38% very shallow Leiden-8 cells; 57% of
  Leiden 2 reached Oligodendrocyte. **v2** (`stage4_insitutree_k21_v2`) restored `o` and nested the glia so the first decision
  is OPC versus the oligodendrocyte lineage: oligodendrocytes (Leiden 2+22) called oligodendrocyte lineage 66% -> 76%, OPC
  18% -> 8%. The OPC leaf stayed 63,186 cells, 23% real OPC. **v3** (`stage4_insitutree_k21_v3`) added `n` as
  `OPC_low_count_denovo` beside reference OPC. **It did not help and is not used:** the reference OPC leaf fell to 9% real
  OPC while the "low count" leaf took the real OPC cells (35% Leiden 0, higher counts), total OPC calls (62.4k vs 63.2k) and
  real-OPC recall (76% vs 75%) were unchanged, and both leaves had 28-30% of cells below posterior 0.9 (posterior <0.99 rose
  21.2% -> 23.2%). `n` has 0 genes at 2-fold over reference OPC, so InSituTree cannot separate them by anything but noise.
  **The repo's `sorl1_hierarchy.json` and matrix are v2, byte-identical to what the v2 run used.**
- **OPC calls are unreliable for shallow cells, and the posterior does not warn.** In v2's OPC leaf the share that is real OPC
  (Leiden 0) rises with depth: <150 counts 6%, 150-299 9%, 300-599 16%, 600-1,199 38%, 1,200+ 78%, while the median
  posterior stays 0.99-1.00. Use `opc_confident_<run>` (called OPC AND (Leiden 0 OR total_counts >= 600)):
  22,289 cells (2.8% of all), 65% in Leiden 0. Add it with
  `OPC_FLAG="tree_v2=OPC" RUNS="... tree_v2:stage4_insitutree_k21_v2:insitutree_result.h5" sbatch pipeline/slurm/54_export_cell_labels.sh`.
- **Rule for keeping a de novo leaf (and the exception, `o`):** it must differ from the leaves it competes with in its branch, at least 10
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
- `reference/sorl1_insitutree_profiles.csv`: the matching 3,323-gene x 22-type matrix.
- `python/tests/test_sorl1_hierarchy.py` checks that leaves and matrix columns match both ways and that only the
  distinct de novo leaves are kept.

**How InSituTree picks genes (read from `supervisedSubcluster.R`).** At each node it takes the children's profiles
(an internal child is the mean of its leaves), normalises each to sum 100, and keeps genes whose range across
the children exceeds the 0.5 quantile on both absolute and percentage difference. That is a loose filter,
not a per-leaf "maximum contrast". Note `85_insitutree.sh` defaults `QUANTILE_PCT=0` (package default 0.5), so
the job selects on absolute difference alone: expect exactly half the genes (1,661 of 3,323) at every node in the
log (`These are the selected N genes`); with the package default it would be about 1,000 to 1,350. The marker genes
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
INPUT_DIR=stage4 STAGE4_DIR=stage4_insitutree_k21_v2 REFERENCE_BASENAME=sorl1_insitutree_profiles_v2.csv \
    HIERARCHY_BASENAME=sorl1_hierarchy.json sbatch pipeline/slurm/85_insitutree.sh
```

The R job prints the tree (`printTree`) and stops if a leaf is missing from the matrix. Check its log for
`hierarchy: 22 leaves over the reference's 22 types`. Calls are read from the finest resolved level; low
`prob` marks the cells InSituTree could not place within a branch. The `Low_signal` sink should be compared
across carriers and non-carriers before interpreting it (the flat clusters are carrier-enriched).

## InSituTree v4 (after Maddie's spatial review, 2026-10) — run on Klone (job 42090220); superseded by v5

Maddie read K=21 and tree v2 in Napari (sheet `SORL1_v3_cell_categories.xlsx`) and found `L6 neuron` cells in L2/3.
The labels table explains it: the v2 `L6 neuron` leaf (14,461 cells) is 32% K=21 L6 cells but **31% `a` and 12% `o`**
(glia, median 574 counts vs 2,006), so it is a shallow-glia attractor; the v2 `L2/3` leaf is only 10.7% K=21 L2/3
(h 17.7, d 15.4, b 12, j 10.1, f 8.8, a 8, q 7.9%). The Brain_6k L6 profile is 0.986 cosine to L4 (0.54 to Oligodendrocyte).
(Her sheet calls this "v3"; that name is taken by the rejected OPC_low_count run, so this is **v4**.)

What v4 changes (files: `reference/sorl1_hierarchy_v4.json`, `denovo_annotations/stage4_k21_v4.csv`, `sorl1_insitutree_profiles_v4.csv`,
`sorl1_hierarchy_v4_collapse.csv`; 29 leaves = 12 named + 17 de novo):
- **Drop the reference `L6 neuron`** (`prep_insitutree_profiles.py --drop-named "L6 neuron"`). `p` and `g` were already leaves.
- **Add `h`, `r`** to a new `Excitatory_deep` branch with c, m, p, g; `Excitatory_upper` = L2/3; `Excitatory_L4` = L4 + s.
- **Add `t`** as `Astrocyte_NEAT1_HSPA1A_denovo` (neutral name: it differs from Astrocyte A/B by stress/blood genes).
- **Sink leaves `d, j, q, f, b`** beside `i` in `Neuron_low_signal`, so homeless cells stop landing in L2/3/L4. They are not cell types; after
  the run merge them with `sorl1_hierarchy_v4_collapse.csv`. A branch with ONE leaf labels cells by the branch name, so the L2/3 leaf
  is reported as `Excitatory_upper`, and all six sink leaves keep their own names until collapsed.
- `n` stays out (tested as v3); `a` stays out (Maddie: 0/20 enriched genes); `k`, `e` are flagged by their K=21 label, not by the tree.
- v2 emitted no internal-node labels, so a deep cell does NOT stop at `Excitatory_deep`; do not expect it to.

```bash
# on Klone, after git pull (check .git/index.lock first if the last pull was interrupted)
s5cmd cp pipeline/reference/sorl1_insitutree_profiles_v4.csv s3://$KOPAH_BUCKET/$KOPAH_PREFIX/reference/sorl1_insitutree_profiles_v4.csv
INPUT_DIR=stage4 STAGE4_DIR=stage4_insitutree_k21_v4 REFERENCE_BASENAME=sorl1_insitutree_profiles_v4.csv \
    HIERARCHY_BASENAME=sorl1_hierarchy_v4.json sbatch pipeline/slurm/85_insitutree.sh
```

Acceptance (all from the labels table; export with `54`, adding `tree_v4:stage4_insitutree_k21_v4:insitutree_result.h5` to RUNS):
1. Cross-tab each excitatory leaf by K=21 source: the glia-derived share (a, o, n, t) of the deep and L2/3 leaves should be near zero.
2. Where the K=21 L6 cluster, h, r and the sink clusters went (expect h, r to the deep branch).
3. **Sink share by group** (Maddie): carriers are shallower and may land in the sink more often, which would shift every neuron comparison.
   `celltype_by_group.py --label celltype_tree_v4 --collapse-csv pipeline/reference/sorl1_hierarchy_v4_collapse.csv`, then read the
   `Neuron_low_signal` row; also report it within slide (July slides only) and against median counts per donor.
4. Re-test t and h trends on Control / AD Pure / AD+LATE donors only (`--exclude-slide-prefix 20260708`).
5. Cortical depth per leaf (pia to white matter) needs cell coordinates plus a boundary annotation; not available locally.

## InSituTree v5 (after the SEA-AD marker review and Maddie's spatial check)

The SEA-AD MTG atlas did not work as the reference (2-6% of cells named vs 28% for Brain_6k; its profiles sit at cosine
0.25-0.42 from the CosMx clusters). Its marker genes did: `marker_overlap_vs_reference.py` gives 9 K=21 clusters a strong
SEA-AD identity (`denovo_annotations/stage4_k21_seaad_markers.csv`), and Maddie confirmed the layer position of the deep
excitatory ones in Napari on 2026-10-27. v5 changes only what that review settled; the tree is otherwise v4.

| change | why |
|---|---|
| `p` -> `Excitatory_L6b_like` | SEA-AD L6b (L6 CT nearly tied); Maddie: marks L6 well |
| `g` -> `Excitatory_L6_IT_Car3_like` | SEA-AD L6 IT Car3; marks L6 well |
| `c` -> `Excitatory_L5_6_NP_like` | SEA-AD L5/6 NP; marks L5/6 well |
| `m` -> `Excitatory_L5_ET_like` | SEA-AD L5 ET; marks L5 well, sparse |
| add `a` as leaf `Astrocyte_low_signal` | with no leaf its 31,238 cells scatter in v4: 21% Astrocyte B, 8% A, 11% OPC, **13% `Excitatory_upper`**, 7% Low_signal |

Files: `denovo_annotations/stage4_k21_v5.csv` (v4 plus the five changes, each with a `v5_note`), `sorl1_hierarchy_v5.json`
(30 leaves = 12 named + 18 de novo), `sorl1_insitutree_profiles_v5.csv` (built by `prep_insitutree_profiles.py` from the
K=21 fit; the 29 shared columns are identical to v4 apart from the renames), and two report-level maps for
`celltype_by_group.py --collapse-csv`:
- `sorl1_hierarchy_v5_collapse.csv` (analysis level): the neuron sink leaves -> `Neuron_low_signal` (as v4), `L4 neuron` ->
  `Excitatory_L4_5_IT_like`, `r` and `s` -> `Excitatory_L5_IT_like`. The L4 neuron's markers match L5 IT best (L4 IT second)
  and it sits above L5/6 NP, hence "L4/5". `r` and `s` overlap spatially and share an SEA-AD call; `s` is absent from about
  half the donors (19 of 40 have >= 20 cells), so they are reported together.
- `sorl1_hierarchy_v5_parent.csv`: each leaf's parent branch (`Astrocyte`, `Microglia`, `Excitatory_deep`, ...), so a split
  such as Astrocyte A vs B is reported at both levels. Both maps only relabel; they never move a cell.

**`Astrocyte_low_signal` is a sink, not a sub-state.** `a` had 0/20 enriched genes against Astrocyte A/B, and in v4's own
notes its profile correlates more with neurons (0.87) than with astrocytes (0.78). It is a leaf so that its cells have
somewhere to go; do not read it as a cell type. It only catches the cells that the tree routes to `Glia` first (see the
measured result below).

What v5 does NOT fix: `Excitatory_upper` still receives cells from the unresolved clusters (in v4, 38% of `h`, 24% of
`j`, 23% of `f`, 22% of `d`, 12% of `q`). They look L2/3-like to the profiles (cosine ~0.96) and the sink leaves do not
hold them. Report that leaf as a mixed excitatory group, not L2/3; for L2/3 use the flat K=21 `L2/3 neuron` call (3,094
cells), which Maddie found spatially accurate. No cluster is a confirmed L4 population yet (`j` leans L4 but is diffuse).

```bash
# on Klone, from the Maddie checkout
git fetch origin && git checkout sorl1-insitutree-v5 && git pull
# .env beats submit-time values: an earlier SEA-AD test left REFERENCE_BASENAME=seaad_mtg_panel.csv in it, which would
# silently replace the v5 profiles with the SEA-AD file
sed -i '/^REFERENCE_BASENAME=/d' pipeline/.env
set -a && source pipeline/.env && set +a
s5cmd cp pipeline/reference/sorl1_insitutree_profiles_v5.csv s3://$KOPAH_BUCKET/$KOPAH_PREFIX/reference/sorl1_insitutree_profiles_v5.csv
INPUT_DIR=stage4 STAGE4_DIR=stage4_insitutree_k21_v5 REFERENCE_BASENAME=sorl1_insitutree_profiles_v5.csv \
    HIERARCHY_BASENAME=sorl1_hierarchy_v5.json \
    sbatch --account=glioblastoma --partition=gpu-l40s --qos=normal --gres=gpu:1 pipeline/slurm/85_insitutree.sh
```

Acceptance (export the labels table with `54`, adding `tree_v4:stage4_insitutree_k21_v4:insitutree_result.h5` and
`tree_v5:stage4_insitutree_k21_v5:insitutree_result.h5` to RUNS, with `OUT_SUBDIR=labels_v5` so `labels/` is not overwritten):
1. **Composition of `Astrocyte_low_signal` by K=21 source:** mostly `a` plus the astrocyte clusters. If neuron clusters contribute
   a meaningful share, `a`'s neuron-like profile is pulling neurons and the leaf should be reconsidered.
2. Where `a` goes now: the share landing in `Excitatory_upper` (13% in v4).
3. The renamed leaves (`p`, `g`, `c`, `m`) should be near their v4 sizes; only small shifts from the new leaf are expected.
4. `Excitatory_upper`: the glia-derived share should fall; the leaks from `h`, `d`, `j`, `f` will remain (above).
5. Re-run the disease comparisons at both levels with the two collapse maps, and within one run before calling anything biology
   (see the depth and run confounds above).

**Measured (run 42538943, 793,380 cells):**
- **Check 1 passes, but the leaf is mixed.** `Astrocyte_low_signal` is 25,144 cells; no named neuron cluster feeds it. 58% is `a`
  (14,521); the rest is `i` (2,916, the low-signal neuronal sink), `t` (1,967, reactive astrocyte), `n` (1,903, OPC, which has no leaf),
  `o` (922), `l` (544) and `k` (281). About 40% is other clusters landing there. Median posterior 1.000, 10.7% below 0.8.
- **Check 2 failed: the leak did not fall.** Of `a`'s 31,238 cells, 46.5% reach `Astrocyte_low_signal` and 11.7% still land in
  `Excitatory_upper` (13% in v4); the rest scatter (Low_signal 5.0%, OPC 3.8%, Microglia A 3.3%, Inhibitory C 3.3%, ...).
- **Why:** the tree decides top-down. At the first level 37.4% of `a` cells go to `Neuron` (55.5% `Glia`, 5.0% `Low_signal`, 2.1%
  `Vascular`), and a cell routed to `Neuron` cannot reach a leaf under `Glia`. So `a`'s hybrid, neuron-like profile (0.87 with
  neurons vs 0.78 with astrocytes) is not fixed by giving it a leaf. `a` is 3.9% of cells; v5 keeps the leaf and reports it as a sink,
  and any further change (for example moving `a` under `Low_signal`) is a new run with its own checks.
- **Check 3 passes.** The renamed leaves kept 96-98% of their v4 cells: `p` 4,362 -> 4,215, `g` 4,688 -> 4,590, `c` 2,268 -> 2,190,
  `m` 2,401 -> 2,362. The renames are label-only; the small losses are cells that went to `Astrocyte_low_signal`.
- **Check 4: a small fall only.** The glia-derived share (`a`/`o`/`n`/`t`) of `Excitatory_upper` went 15.3% -> 13.8% (27,719 -> 26,976
  cells), all of it from `a` (4,112 -> 3,664). The flat K=21 `L2/3 neuron` is only 2,867 cells (10.6%) of the leaf; the rest is unresolved
  clusters (`h` 4,357, `i` 4,082, `a` 3,664, `d` 2,609, `q` 2,558, `j` 1,987, `f` 1,286). Report `Excitatory_upper` as a mixed excitatory
  group, never as L2/3.

## Resolved items

- Slide names `6562 A7 576 …` and `…6802 A6` carry a typo made when the slides were
  entered; the real donors are **526** and **6803**, as in the annotation sheet. The slide
  names are export directory names and stay as they are, so never infer donors from them.
