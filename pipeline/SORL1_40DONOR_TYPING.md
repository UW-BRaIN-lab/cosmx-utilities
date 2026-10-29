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
   (see the depth and run confounds above). Done: see "Disease comparisons on v5" below.

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

**v4 -> v5: where the Vascular and Low_signal cells went.** Vascular fell 111,411 -> 108,734 (3,257 left, 580 joined) and Low_signal
138,518 -> 133,416 (6,008 left, 906 joined). The new leaf is not the main destination: 13% of the Vascular leavers and 11% of the
Low_signal leavers went to `Astrocyte_low_signal`; the rest went to OPC (Vascular 1,030, Low_signal 2,228) and microglia (1,531, 1,254).
Only 316 of the 3,257 Vascular leavers (9.7%) are K=21 Endothelial; the sources are `i` (1,004), `n` (524), Microglia A (348), `o` (284)
and `a` (277). No donor dominates (largest 10.3%). The leavers were borderline in v4: median v4 posterior 0.635 (Vascular) and 0.751
(Low_signal) against 1.000 for the cells that stayed (lower quartile 0.51-0.59 vs 1.0). So the drops are reassignments of
low-confidence cells at the edges. Why adding `a` moved them is a hypothesis (an internal node's profile is the mean of its leaves, so
the Astrocyte and Glia profiles changed); it has not been tested. The Vascular disease result did not change (q 0.40).

### Disease comparisons on v5 (check 5)

`celltype_by_group.py --label celltype_tree_v5`, per-donor % of cells averaged by group, Kruskal-Wallis across groups with BH q within
each run, exploratory. Three designs because every SORL1 donor is on the July-run slides (`20260708`): all donors (12 Control / 11 AD Pure
/ 11 AD+LATE / 6 SORL1), without the July slides (10 / 7 / 7 / 0), and July slides only (2 / 4 / 4 / 6). The July-only run keeps the
run constant but has 2 Control donors.

| Type | Design | Control | AD Pure | AD+LATE | SORL1 | q | trend q |
|---|---|---|---|---|---|---|---|
| Microglia B | all donors | 1.86 | 3.17 | 3.73 | 6.69 | 0.0145 | 0.199 |
| | no July | 1.93 | 2.77 | 4.29 | n/a | 0.81 | 0.087 |
| | July only | 1.52 | 3.88 | 2.75 | 6.69 | 0.086 | 0.97 |
| OPC | all donors | 5.23 | 6.81 | 7.11 | 14.8 | 0.030 | 0.426 |
| | no July | 4.46 | 5.36 | 5.26 | n/a | 0.81 | 0.65 |
| | July only | 9.11 | 9.34 | 10.3 | 14.8 | 0.27 | 0.97 |
| Neuron_low_signal (sink) | all donors | 9.6 | 10.4 | 10.6 | 9.4 | 0.987 | 0.859 |
| | no July | 7.23 | 6.55 | 5.38 | n/a | 0.81 | 0.66 |
| | July only | 21.4 | 17.2 | 19.8 | 9.4 | 0.092 | 0.97 |
| Astrocyte A | all donors | 1.32 | 2.25 | 2.72 | 0.51 | 0.118 | 0.426 |
| | no July | 1.52 | 3.25 | 3.73 | n/a | 0.39 | 0.087 |
| | July only | 0.35 | 0.50 | 0.94 | 0.51 | 0.42 | 0.97 |

At the parent level (all donors): Microglia q 0.0006 (6.04 / 8.68 / 8.42 / 14.9; trend q 0.066), OPC q 0.014, Excitatory_deep q 0.083,
Astrocyte q 0.36 (14.8 / 15.9 / 15.8 / 12.5), Neuron_low_signal q 0.99, Low_signal q 0.89.

- **Microglia B is the one result that holds in every design.** It rises Control -> AD Pure -> AD+LATE (trend q 0.087 without the July
  slides; Microglia as a branch, trend q 0.066 with all donors) and is highest in SORL1, 6.69% against 1.5-3.9% on the same run.
  Exploratory: small groups, and BH treats correlated tests as independent.
- **Astrocytes shift state, not amount.** The branch total is flat; within it A and the stress-state leaf rise and B falls (no July: A
  1.52 -> 3.73, B 6.05 -> 2.19). Only A reaches trend q < 0.10, and only without the July slides. The A/B split depends on the tree
  version: `t` took about half of the old A and B in v4.
- **The SORL1 OPC elevation is mostly the run.** Non-SORL1 donors on the July slides have 9.1-10.3% OPC against 4.5-5.4% elsewhere;
  within the run SORL1 is about 1.5-fold higher and not significant (q 0.27). Microglia A shows the same run shift (6.3-6.7% on July
  slides vs 3.7-5.1%). OPC calls on shallow cells are unreliable (see above), so redo OPC with `opc_confident` before using it.
- **The sink depends on the run.** Non-SORL1 donors have 17-21% in `Neuron_low_signal` on the July slides and 5-7% on the others
  (about 3-fold). The flat group means in the all-donors run (q 0.99) come from averaging the two runs; SORL1 is 9.4%. Do not read
  "no difference by group" as "no depth or run effect".
- **Excitatory and inhibitory neurons are lower in SORL1, but not clearly SORL1-specific.** All donors: SORL1 about 5-fold lower in
  every excitatory leaf (for example `Excitatory_upper` 0.96% vs 2.6-5.7%). July only, the 8 types with q < 0.10 are mostly neuron
  types, which are correlated, and AD Pure is as low as SORL1 for some (`Excitatory_L5_6_NP_like` 0.64 / 0.11 / 0.25 / 0.10).
  The within-slide contrast below shows the neuron deficit is not a run effect.
- **Within-slide contrast (SORL1 vs non-SORL1 AD donors on the same July slide).** The run cannot be repeated, but each of the four
  July slides holds both kinds of donor, so the slide is the block (`within_slide.py`: per slide, mean per-donor % in SORL1 donors
  minus the AD Pure and AD+LATE donors on that slide; 1-2 SORL1 and 1-3 comparator donors per slide). Same sign on all four slides,
  SORL1 against comparator: `Neuron_low_signal` 9.4% vs 18.5% (-6 to -14 points on every slide), Low_signal 18.1% vs 11.4% (+2.8 to
  +17), Microglia B 6.69% vs 3.31% (+2.3 to +7.0), Microglia A 8.19% vs 6.28% (+0.7 to +2.7), and lower in `Excitatory_upper`
  (0.96 vs 2.73), `Excitatory_L4_5_IT_like` (0.50 vs 1.50), the inhibitory leaves, VIP and `Excitatory_CARTPT_GAP43_denovo`. OPC is
  higher on 3 of 4 slides (14.8% vs 9.8%, median +4.9 points; the exception is -1.4). Net: SORL1 donors have roughly 15 points fewer
  neuron-lineage cells and correspondingly more Low_signal, microglia and OPC than the AD donors processed beside them.
- **Depth and donor quality (July slides).** Median counts per cell: five of the six SORL1 donors sit at 191-313 (599 191, 6589 224,
  7068 297, 7796 303, 7181 313; 7137 is 739), against 389-806 for the AD comparators on the same slides (exception: AD Pure 6745 at 140).
  Group medians: SORL1 289, AD Pure 458, AD+LATE 614, Control 654. Restricted to cells with >= 600 counts (104,770 July cells), the
  within-slide pattern persists for Low_signal (31.9% vs 14.2%, slide differences +11 to +25, 4/4), `Neuron_low_signal` (20.1% vs
  32.2%, 4/4), `Excitatory_upper` (1.07 vs 3.28), `Excitatory_L4_5_IT_like` (1.56 vs 3.76) and the inhibitory leaves. Composition alone
  (Low_signal +17.7 points) would take `Excitatory_upper` from 3.28 to about 2.6, so the neuron drop is more than proportional. The
  microglia result shrinks: Microglia B 4.66% vs 1.97% is the same sign on 4/4 slides but only +0.66 and +0.17 on two of them (donors
  599 and 7068 carry it), and Microglia A's median difference is 0.00. OPC falls to 7.5% vs 5.0% (14.8 vs 9.8 with all cells).
- **The pattern is not uniform across SORL1 donors.** Among the >= 600-count cells, Low_signal is 37-48% in four SORL1 donors (599,
  7068, 7796, 6589) and 20.7% and 1.5% in the other two (7137, 7181), which resemble their AD neighbours (6966 13.0%, 7665 5.9%). The four
  high-Low_signal donors are among the shallowest, but depth alone does not predict it: 7181 (median 313) does not follow, and AD Pure
  657 has 29.9% at median 565.
- **Sensitivity: K=21 `a` as its own group (`a_mixed`); a neuron-side leaf (v6) is not pursued.** `a` (31,238 cells) shows astrocyte and
  neuron markers together. The Neuron-routed half (11,693 cells, median 526 counts) has 63-90% of the K=21 neuron level of SNAP25, SYT1
  and STMN2 and 14-20% of the K=21 astrocyte level of SLC1A2, AQP4 and GLUL; the Glia-routed half (17,338 cells, median 327) has 30-49%
  of the neuron level and 30-38% of the astrocyte level. Every marker shifts by about the same factor between the halves, so the first
  split cuts a graded mixture; a leaf under Neuron would name cells that still carry astrocyte signal, and its profile would sit outside
  the K=21 fit. Relabelling `a` as its own group leaves the disease results unchanged: all donors, the same four types have q < 0.10
  (Microglia B 0.0152, OPC 0.0228, `Astrocyte_NEAT1_HSPA1A_denovo` 0.058, Microglia A 0.058; Microglia B group means identical); without
  July, Microglia B and Astrocyte A trend q 0.091 (0.087 before); within-slide, every sign and slide-agreement count is unchanged
  (`Neuron_low_signal` 9.28% vs 18.24%, Low_signal 18.07% vs 11.37%, Microglia B 6.69% vs 3.31%). `a_mixed` is not a SORL1 feature:
  0.81% of cells in SORL1 donors against 3.4-5.2% in the other groups (q 0.14), lower in SORL1 on all four July slides. Whether the
  neuron-like `a` cells are a real hybrid state, two populations merged, or ambient/segmentation mixing is not settled; cell
  coordinates (are they among neurons?) would help. For analyses that need a clean group, report K=21 `a` as "mixed low-signal glia and
  neuron".
- **What cannot be concluded.** A SORL1 effect on neuron or Low_signal shares. Per-type p-values are not possible (with four slides the
  smallest two-sided sign-test p is 0.125, and the agreeing types are correlated), and the contrast cannot separate SORL1 from donor
  RNA quality: the effect sits in four low-depth donors and is absent in two. Next: tissue-quality metadata (post-mortem interval, RNA
  integrity) for 599, 7068, 7796 and 6589, and a donor-level sensitivity analysis that drops donors by a pre-set quality threshold.
  OPC still needs the `opc_confident` flag. Microglia B is the result that holds in every design, but it is carried by two donors here.

### Napari: the v5 columns (not run)

Adds two dropdown columns to the live `_metadata.csv` of the 10 stitched slides, next to the existing ones (`InSituType_K21`, `InSituTree_v4`,
`SEAAD_markers_K21`, ...): `InSituTree_v5` (30 leaves; the L2/3 branch shows as `Excitatory_upper`) and `InSituTree_v5_collapsed`
(`sorl1_hierarchy_v5_collapse.csv`: sink leaves merged, `L4 neuron` -> `Excitatory_L4_5_IT_like`, `r`+`s` -> `Excitatory_L5_IT_like`). It uses the Napari
scripts from PR 92 (`maddie-napari-metadata`: `prepare-sorl1-napari-labels.py`, `validate-napari-label-columns.py`), which are not on this branch, plus
`labels-csv-to-typing-h5.py` and `celltypes-to-napari-metadata.py` from `main`. Do not re-stitch a slide afterwards: the stitcher wipes its output folder,
backups included.

```bash
# 1. on Klone: a labels table with the runs the prepare script needs (tree_v2 is the base InSituTree column and carries the OPC flag)
RUNS="k21:stage4_k21 tree_v2:stage4_insitutree_k21_v2:insitutree_result.h5 tree_v4:stage4_insitutree_k21_v4:insitutree_result.h5 tree_v5:stage4_insitutree_k21_v5:insitutree_result.h5" \
    OPC_FLAG="tree_v2=OPC" OUT_SUBDIR=labels_napari_v5 sbatch pipeline/slurm/54_export_cell_labels.sh
# then, from a Mac terminal (one Duo): scp it down as cell_labels_napari_v5.csv.gz

# 2. on the Mac, in a checkout that has PR 92's scripts: one (cell_id, label) CSV per column
uv run python scripts/prepare-sorl1-napari-labels.py --labels cell_labels_napari_v5.csv.gz \
    --k21-annotations pipeline/reference/denovo_annotations/stage4_k21.csv --tree-run tree_v2 \
    --extra-tree-run tree_v5 --collapse-csv pipeline/reference/sorl1_hierarchy_v5_collapse.csv --out-dir napari_labels_v5
# only insitutree_tree_v5.csv and insitutree_tree_v5_collapsed.csv are used below; the others rebuild the existing columns

# 3. one h5 per new column, then chain the join onto the CURRENT live files (live/<slide>_metadata.csv, downloaded first;
#    flat/<slide>_metadata_file.csv.gz are the export's flat files)
uv run python scripts/labels-csv-to-typing-h5.py --csv napari_labels_v5/insitutree_tree_v5.csv --label-column label --output v5.h5
uv run python scripts/labels-csv-to-typing-h5.py --csv napari_labels_v5/insitutree_tree_v5_collapsed.csv --label-column label --output v5_collapsed.h5
uv run python scripts/celltypes-to-napari-metadata.py --typing-h5 v5.h5 --flatfiles flat --column InSituTree_v5 --merge-into live --out-dir step1
uv run python scripts/celltypes-to-napari-metadata.py --typing-h5 v5_collapsed.h5 --flatfiles flat --column InSituTree_v5_collapsed --merge-into step1 --out-dir step2

# 4. independent check: must report 0 disagreements before anything is uploaded
uv run python scripts/validate-napari-label-columns.py --merged-dir step2 --original-dir live --labels-dir napari_labels_v5 \
    --column InSituTree_v5=insitutree_tree_v5.csv --column InSituTree_v5_collapsed=insitutree_tree_v5_collapsed.csv
```

Then, per slide: copy the live `_metadata.csv` server-side to `_metadata.csv.pre-v5.bak` (check it does not exist first; leave the four existing backups alone),
upload the merged file, download it back and compare md5 of the local files (files over 8 MiB have multipart ETags, so do not compare ETags).
Rollback: copy `_metadata.csv.pre-v5.bak` over `_metadata.csv`. Each upload overwrites a live file that Maddie may have open, so confirm before running it.

## Resolved items

- Slide names `6562 A7 576 …` and `…6802 A6` carry a typo made when the slides were
  entered; the real donors are **526** and **6803**, as in the annotation sheet. The slide
  names are export directory names and stay as they are, so never infer donors from them.
