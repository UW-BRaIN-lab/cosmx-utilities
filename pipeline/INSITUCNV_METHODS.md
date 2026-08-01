# Stage 5 — InSituCNV: Methods

Copy-number inference on CosMx spatial transcriptomics to test whether the transcriptionally
"unresolved" (`Low_signal`) cell fraction is tumor. This documents every design choice,
parameter, and analysis metric for reproducibility and for the manuscript Methods.

---

## 1. Question and rationale

Supervised cell typing (InSituTree, against the CZI GBmap atlas) leaves ~48% of cells as
`Low_signal` — cells with adequate RNA depth and complexity that match no discrete lineage
program on the 6k panel. Because copy-number is a property of the genome, orthogonal to the
expression state, we use CNV to ask whether these cells are malignant (carrying the GBM
signature: chr7 gain, chr10 loss, and additional arm losses) or genuinely normal.

Method: the Moldia **InSituCNV** recipe (`github.com/Moldia/InSituCNV-manuscript`, pinned
commit `e3c7ee3`) — spatial neighborhood smoothing to recover CNV-usable depth, followed by
`infercnvpy` (the Python re-implementation of inferCNV). Adapted from the colleague CLI at
`UW-BRaIN-lab/Glioblastoma/insituCNV-pipeline`, with the deviations documented below.

## 2. Input

- **Cells / labels:** `stage4_insitutree/cosmx_typed.h5ad` — 2,331,655 cells, `obs['cell_type']`
  (InSituTree assignment, incl. `Low_signal`), `obs['Case']` (donor), `obs['Block']` (physical
  tissue block), `obs['Region']` (`Tumor bulk` / `Infiltrating edge` / `Contralateral
  uninvolved`), `obs['slide_id']`.
- **Counts:** raw integer gene counts (`X`), gene probes only (negprobe/falsecode dropped via
  `var['probe_type']`). Panel ≈ 5,900 gene probes → **arm-level CNV resolution, not focal.**
- **Spatial:** per-cell centroids `CenterX_global_px` / `CenterY_global_px` → `obsm['spatial']`.
  These are per-slide (global-px repeats across slides), which motivates §3.

## 3. Per-tissue-section processing (key correctness decision)

The spatial smoothing graph is built **per physical tissue section**, keyed
`tissue_section = "<slide_id>__<Case>__<Block>"`, never cohort- or slide-wide. Rationale:

- Global-px coordinates repeat across slides, so a cohort-level graph bridges patients.
- **Each CosMx slide co-mounts two tissue blocks that are different regions** (e.g. one
  tumor-bulk + one contralateral block, often from the same donor). So even a per-slide (or
  per-donor) graph would smooth counts — and subtract the CNV reference — across physically
  distinct tissues, most damagingly merging a tumor block with its contralateral control.

`(slide, donor, block)` is the finest unit and never merges across slide (coordinate
collision), donor (patient), or block (distinct tissue). inferCNV genomic windows are
gene-position-based and identical across sections, so per-section `X_cnv` matrices concatenate
cleanly for cohort-level comparison. (Note: the upstream colleague wrapper builds a single
cohort-level neighbor graph with no batch key — a cross-tissue/cross-patient confound we
correct here.)

## 4. Gene genomic positions

Genes are annotated with chromosome/start/end from a **full-genome Ensembl BioMart table**
(bundled in the container), keeping **autosomes 1–22 only**. We deliberately do NOT use
InSituCNV's `add_genomic_positions` (which maps against the ~3k-gene `maynard2020_3k` set and
would discard most of the ~5.9k panel, and downloads at runtime). Autosome-only avoids a
sex-mismatch baseline artifact; GBM's chr7/chr10/chr9 events are all autosomal. Result: 5,926
of 6,175 gene-probes mapped (96% of the panel).

## 5. CNV inference recipe (per section)

Following the Moldia manuscript recipe (`run_insitucnv.py`):

1. `layers['counts'] = X` (raw); `sc.pp.normalize_total(target_sum=1e4)` (NOT log-transformed).
2. Spatial neighbor graph `scvelo.pp.neighbors(use_rep='spatial', n_neighbors=20)`, then
   InSituCNV `smooth_data_for_cnv` → `layers['M']` (M = connectivities · normalized counts;
   neighborhood pooling to recover depth).
3. `normalize_total(target_sum=1e4)` + `log1p` on the smoothed layer → `layers['M_log1p']`.
4. `infercnvpy.tl.infercnv(layer='M_log1p', reference=<see §6>, window_size=100, step=10,
   dynamic_threshold=1.5)` → `obsm['X_cnv']` (cells × ~390 genomic windows).

**Smoothing strength (`n_neighbors=20`).** Reduced from InSituCNV's default (100). Heavy
smoothing over-blends in diffusely-infiltrating GBM: it flattens tumor–normal contrast and
contaminates the reference (see §9). Since the primary readout is a group mean over 10⁵–10⁶
cells, heavy per-cell smoothing is unnecessary; 20 balances depth recovery against contrast.
`window_size` reduced 200→100 to match the panel's gene density.

## 6. Reference (diploid baseline) — per-donor matched

inferCNV centers each cell on a diploid reference; the reference choice is the pivotal
scientific decision, and the reference is built explicitly (`build_insitucnv_reference.py`),
not taken from a section's own cells:

- **Cell types (conservative):** confidently non-malignant lineages — oligodendrocyte, neuron,
  all immune (myeloid + lymphoid), all vascular/stroma. **Astrocyte and OPC are excluded**
  because AC-like and OPC-like are malignant Neftel states and the reactive/precursor normals
  are hard to separate cleanly from their tumor look-alikes on this panel.
- **Per-donor matched normal:** the reference is built from each donor's own
  **contralateral-uninvolved** reference-type cells, processed through the identical
  smoothing recipe, averaged per gene into an M_log1p reference profile (passed to inferCNV
  via `reference=`). A cross-donor pooled column ("GLOBAL") is the fallback for donors lacking
  contralateral tissue (all 12 donors here have it). Matched normal is the CNV gold standard —
  it cancels donor germline/batch and, critically, is sourced from tissue with no tumor field
  effect. This was decisive: it un-suppressed the signal ~100× versus per-section
  self-referencing (which contaminates tumor-bulk baselines).

**Reference-choice robustness (immune-free sensitivity).** Re-running with all immune cells
removed from the reference (`insitucnv_reference_types_noimmune.txt`) leaves the per-donor
**median signatures identical to ~3 decimals** — so the result does not depend on immune cells
in the baseline (they are contralateral-sourced and flat). Including immune, if anything,
*raises* the negative-control threshold (see §8), i.e. is the more conservative choice.

## 7. Malignant-signature (primary per-cell metric)

Each cell's genome-wide CNV profile `X_cnv_i` (∈ ℝ^~390 windows; + = gain, − = loss) is scored
by **cosine similarity to the malignant consensus**:

```
signature_i = (X_cnv_i · c) / (‖X_cnv_i‖ · ‖c‖),   c = mean X_cnv over confirmed-malignant cells
```

where the consensus `c` is the mean CNV profile of the 14 confidently-typed malignant states
(9 GBmap Neftel + 5 de-novo). Range −1…+1; it measures whether a cell's **pattern** of
coordinated gains/losses (chr7+, chr10−, chr9−, chr14−…) matches the tumor karyotype,
independent of magnitude.

**Why cosine / direction, not magnitude.** On a targeted panel the inferred CNV *magnitudes*
are small, noisy, and depth-dependent; a non-directional burden (`cnv_score = mean(X_cnv²)`,
our first attempt) failed to separate tumor from normal (malignant 4.5% vs reference ~5%). The
*pattern* — which arms move together — is robust, and cosine to the consensus captures exactly
it, normalizing away magnitude/depth. This is what made the positive control separate.
(`cnv_score` is retained as a secondary output but is not the discriminator.)

**Gains vs losses.** The malignant-signature is a cosine to the FULL genome-wide consensus,
so it uses coordinated **gains AND losses together** and privileges neither direction.
Empirically, in our matched-reference results both are robustly detected in the malignant
states and in tumor-region Low_signal, and **chr7 gain is as strong as or stronger than the
chr10/arm losses** (numbers from the corrected n_neighbors=20 / window=100 run, ~390 windows):
across the core Neftel malignant states chr7 ≈ +0.06 to +0.22 (OPC-like +0.22, MES-like_denovo
+0.20, AC-like +0.17) vs chr10 ≈ −0.04 to −0.12 — lower (chr7 +0.02–0.04) only in the
transcriptionally-flat de-novo stress/hypoxia states, which are CNV-shallow by nature. The
malignant consensus averages chr7 ≈ +0.11 vs chr10 ≈ −0.04 as *whole-chromosome* means,
but at window=100 the structure resolves — **chr7 is a broad both-arm gain** (peaks over the
EGFR/p-arm region and the q-arm past MET, dip near CDK6), while the **chr10 loss is 10q-restricted**
(p-arm near-neutral, q-arm deepening to ≈ −0.08 at the q-terminus through PTEN/MGMT), so the
whole-chromosome chr10 mean understates the actual loss. Low_signal | infiltrating edge
chr7 = +0.113 vs chr10 = −0.047; Low_signal | tumor bulk chr7 = +0.068 vs chr10 = −0.045;
contralateral ≈ 0 for both (chr7 = +0.003, chr10 = +0.001 — control holds).

This is consistent with the detection asymmetry on a sparse targeted panel: the expression
**detection floor** compresses *loss* signal — a gene already near zero counts cannot drop
further, so a hemizygous deletion is confounded with dropout — whereas *gains* have upward
headroom and read cleanly, and high-level focal amplifications (e.g. EGFR) are the most
detectable of all. Pure fold-change math favors losses (2→1 loss, log2 ≈ −1.0, vs 2→3 gain,
≈ +0.58, with gains further buffered by dosage compensation), which is why losses are called
cleanly on *deep whole-transcriptome* scRNA; on a *sparse* panel the floor effect offsets that
advantage, so gains (and amplifications) tend to be the more robust signal. Because the
signature integrates the whole pattern, the classification does not depend on this asymmetry.
We additionally verify panel gene coverage per chromosome against the signal
(`chr_coverage_vs_signal.{csv,png}`: panel gene count vs malignant-consensus CNV per
chromosome), so no signal arm is a gene-density artifact. (See Schmid et al., *Benchmarking
scRNA-seq copy number variation callers*, Nat. Commun. 2025; it benchmarks droplet scRNA-seq,
whereas CosMx is imaging-based/targeted, so the dropout mechanism differs but the low-count /
floor ambiguity is analogous.)

## 8. Thresholds and controls

Signatures are calibrated against controls, not interpreted in absolute terms:

- **Positive control:** the 14 malignant states (must show chr7+/chr10− and high signature).
  Their separation is partly by construction (the consensus is built from them), so it is a
  sanity check, not evidence.
- **Negative controls:** `Low_signal | Contralateral uninvolved` and the diploid reference —
  must stay flat.
- **Malignant-call threshold:** the 95th percentile of the negative controls (**strict**,
  ≈0.75 here). A **sensitive** threshold at the trough between the two modes of the Low_signal
  signature histogram (≈0.45) is reported alongside as a less-conservative estimate. Reported
  fractions are always bracketed strict↔sensitive.

## 9. Field effect and its controls

**Field effect:** in diffusely-infiltrating GBM a diploid cell embedded in tumor has mostly
tumor spatial neighbors, so smoothing pulls ambient tumor RNA into its profile and inferCNV
reports a low-level *apparent* CNV — not a real genomic change. It is region-graded (strong in
dense bulk, weak at the infiltrating edge, absent in contralateral) and lineage-graded
(tumor-resident proliferating/hypoxic/vascular subsets most affected; dispersed quiescent
normals least). Two controls quantify and correct for it:

1. **Within-region contrast** — Low_signal is compared against reference cells in the **same**
   region; both are equally field-exposed, so a Low_signal excess over the region's reference
   (its "field-effect floor" = reference-median signature) is tumor-specific. (`within_region_
   contrast.csv`, `within_region_signature.png` with floor + threshold lines.)
2. **Normal-types-by-region** — the median signature of every reference (non-malignant) type,
   split by region. Real CNV cannot depend on location; a field effect does. Reference types
   are flat in contralateral and elevated only in tumor, proving their apparent CNV is
   contamination, not copy number. (`normal_types_by_region.{csv,png}`.)

Consequence: calls are cleanest at the **infiltrating edge** (Low_signal signature ≫ the ~0.15
edge floor) and are treated as a **lower bound in the dense bulk** (elevated floor, overlapping
distributions). Reducing `n_neighbors` further does not fix the bulk, because in dense tumor
the neighborhood is ~all tumor at 10 or 20 neighbors (the *fraction*, not the count, drives
contamination) while lowering it adds per-cell noise.

## 10. Cohort-level readouts and reporting

- **Group-mean CNV profiles** (per `cell_type`, with `Low_signal` split by region), their
  **cosine-similarity matrix**, and per-chromosome-arm means — the aggregate readout
  appropriate for a noisy targeted panel (per-cell CNV is not clustered).
- **Per-donor / per-region, never pooled-only.** The malignant fraction of Low_signal is
  strongly patient-specific (from ~none to ~90% of a donor's Low_signal). Pooled comparisons
  both hide this and invite Simpson's-paradox artifacts — e.g. a pooled "edge > bulk" gradient
  that holds at every matched RNA-depth bin nonetheless reverses within most donors, so it is
  donor composition, not biology. Per-donor breakdowns (`lowsignal_by_donor*`) and depth-
  stratified diagnostics (`insitucnv_edge_vs_bulk.py`) are reported.

## 11. Outputs

Per section: `<section>_cnv.h5ad` (`obsm['X_cnv']`, `obs['cnv_score']`). Cohort (`compare/`):
`SUMMARY.txt`, `group_mean_cnv.csv`, `cosine_similarity.{csv,png}`, `chr_arm_summary.csv`,
`chr7_chr10.png`, `chr7_chr10_gene_positions.{png,csv}` + `chr7_chr10_windows.csv` (panel gene
positions + malignant-consensus CNV localized along chr7/chr10), `chr_coverage_vs_signal.{csv,png}`,
`chromosome_heatmap.png`, `cnv_score_by_group.csv`, `mal_signature_by_group.csv`,
`within_region_contrast.csv`, `within_region_signature.png`, `normal_types_by_region.{csv,png}`,
`lowsignal_by_region.{csv,png}`, `lowsignal_cnv_by_region[_thr<t>].{csv,png}`,
`lowsignal_by_donor[_thr<t>].{csv,png}`, `mal_signature`/edge-vs-bulk diagnostics.

## 12. Limitations

- Targeted 6k panel → **arm-level, not focal** CNV. Gains and losses are asymmetric on a
  sparse panel (§7): the detection floor compresses single-copy losses while gains (and
  amplifications) read more robustly. In our data both directions are detected (chr7 gain ≈
  chr10 loss or stronger); the signature integrates the whole pattern, so the result does not
  hinge on the asymmetry, and panel coverage of the signal arms is verified.
- **Field effect** in dense tumor limits per-cell specificity there (§9); edge is the clean
  regime, bulk a lower bound.
- Mild **depth-dependence** of the signature; contralateral cells are also shallower (median
  ~187 vs ~400 counts in tumor), though their flatness is corroborated by lineage.
- The malignant consensus is cohort-defined, so positive-control separation is not independent
  evidence; the negative controls and per-donor structure carry the inference.

## 13. Software / reproducibility

CPU-only container `insitucnv.def`: Python 3.11, `infercnvpy==0.6.0`, `scvelo==0.3.3`,
`scanpy==1.11.2`, `anndata==0.11.4` (pinned); Moldia InSituCNV @ `e3c7ee3`; BioMart gene table
baked in. Run order: `86` prep → `86b` build per-donor reference → `87` array (one job per
tissue section) → `88` compare → `89` edge-vs-bulk. Immune-free sensitivity: rerun `86b/87/88`
with `REF_TYPES_BASENAME=insitucnv_reference_types_noimmune.txt` into a separate `STAGE5_DIR`.
