# Stage 4 reference profiles

InSituType (stage 4b) types cells semi-supervised against a genes × cell-types
reference profile matrix, restricted to the CosMx panel genes. One reference per
study; pick it at submit time with `REFERENCE_BASENAME`.

## GBM — CZI GBmap

We use the **CZI GBmap** "core" atlas at **annotation
level 4** (the finest level: 54 cell types — e.g. `AC-like`, `MES-like_hypoxia_MHC`,
`NPC-like_neural`, `TAM`/`Mono` subsets, `Oligodendrocyte`, `Endothelial`, lymphoid
subsets), values = **raw average expression** per cell type.

### Files

- `gbmap_level4_panel.csv` — the GBmap level-4 reference **restricted to the CosMx 6k
  panel genes**. Genes are the row index (HGNC symbols), cell types are columns. This
  is the committed, reproducible artifact and the input to `R/insitutype_typing.R`.

The full ~28k-gene GBmap export is **not** committed (~23 MB; mixes HGNC symbols and
Ensembl IDs). The committed CSV is the panel-restricted subset only.

- `gene_signatures.csv` — gene-signature modules (`module,gene`) for the Stage-5 Phase-2
  continuum-vs-hybrid scoring (`python/insitucnv_hybrid_continuum.py`): the four Neftel GBM
  malignant-state modules (AC/MES/OPC/NPC-like) + Cycling. The two axes (AC↔MES, OPC↔NPC) let the
  script test within-axis co-expression (continuum) vs cross-axis co-expression (two-state hybrid).
  Canonical literature signatures intersected with the CosMx 6k panel (off-panel genes dropped).

### Provenance / how to regenerate

Source (local, not in this repo):
`~/keene-lab/GBM/GBmap/core_GBmap_results/core_GBmap_gene_by_annotation_level_4_raw_avg_expression.csv`

Restrict to the CosMx panel (run once, locally):

```bash
uv run python pipeline/python/prep_insitutype_reference.py \
    --reference-csv ~/keene-lab/GBM/GBmap/core_GBmap_results/core_GBmap_gene_by_annotation_level_4_raw_avg_expression.csv \
    --panel-h5ad <a per-slide stage-1 .h5ad or the stage-3a combined_qc.h5ad> \
    --output pipeline/reference/gbmap_level4_panel.csv
```

Then upload to Kopah so the Slurm typing job can stage it:

```bash
set -a; source pipeline/.env; set +a
AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID" \
AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY" \
S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL" \
    s5cmd cp pipeline/reference/gbmap_level4_panel.csv \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/gbmap_level4_panel.csv"
```

The other GBmap annotation granularities (`..._by_cell_type_...`, `..._level_3_...`)
exist in the same source directory; swap the `--reference-csv` to retype at a coarser level.

## Retina / brain / optic nerve — combined ocular atlas

The retina study is 12 slides, **one donor per slide**, each carrying a cross-section of
**retina + optic nerve + brain**. Its reference is correspondingly a *combined* atlas —
**HRCA v2** (retina) + **Monavarfeshani** (posterior eye / optic nerve) + **Allen Brain
Cell Atlas** (cortex) — built with `pipeline/reference_builder/`
(`examples/retina_brain_ocular.R`), which rescales each source independently and prefixes
type names per source (`HRCA_`, `Mona_`, `Allen_`).

### Files

- `retina_combined_panel.csv` — the combined atlas restricted to the panel and with
  reviewed duplicate merges applied: **5,993 genes × 64 cell types** (97.1% of the 6,175
  panel genes are in the atlas; 68 raw types → 64).
- `retina_type_merges.csv` — the reviewed merge/drop spec: which duplicate source types
  collapse, and which taxonomy QC categories are dropped. See the section below.
- `retina_hierarchy.json` — the InSituTree cell-type hierarchy over those 64 types
  (11 top-level branches, depth 4). Named types only — de novo leaves get added after
  triage, in lockstep with the CosMx-scale profile matrix.
- `cosmx_6k_panel_genes.txt` — the 6,175 gene probes of the **Human RNA 6k Discovery**
  panel, one per line. Derived from an exprMat header with `fov`/`cell_ID` and the
  `Negative*` (20) / `SystemControl*` (324) control probes dropped. Lets
  `prep_insitutype_reference.py` panel-restrict a reference *before* stage 1 exists.
  The GBM and retina studies run the same panel, so this file serves both.

### Provenance / how to regenerate

Source (local, not in this repo):
`~/keene-lab/retina-brain/HRCA-Monavarfeshani-Allen-combined-profile.csv`
(21,040 genes × 68 types).

Two steps — panel-restrict, then apply the reviewed merges:

```bash
uv run python pipeline/python/prep_insitutype_reference.py \
    --reference-csv ~/keene-lab/retina-brain/HRCA-Monavarfeshani-Allen-combined-profile.csv \
    --panel-genes pipeline/reference/cosmx_6k_panel_genes.txt \
    --output pipeline/reference/retina_combined_panel.csv
```

```bash
uv run python pipeline/python/merge_reference_types.py \
    --reference pipeline/reference/retina_combined_panel.csv \
    --spec pipeline/reference/retina_type_merges.csv \
    --output pipeline/reference/retina_combined_panel.csv
```

Marker sanity check on the result (all top-ranked in the expected type): RHO → rod,
NEFL → ganglion, MBP/PLP1 → oligodendrocyte, AQP4 → astrocyte, SLC17A7 → cortical
excitatory, GAD1 → interneuron, P2RY12 → microglia, MLANA → melanocyte,
CLDN5 → endothelium.

### Duplicate source types: what was merged, what competes, and why

Because each source atlas is prefixed rather than merged, one biological cell type recurs
once **per source**. Whether that split is worth keeping is a **signal-to-noise** question,
not a cosine question — two profiles at cosine 0.99 can still differ on a few high-count
genes, and two at 0.95 can differ only on genes that are near-zero in a real cell. So each
duplicate group was scored by *how many expected counts of a 400-count cell sit on genes
that differ ≥2× between copies*:

| group | copies | min cosine | discriminating counts / 400 | decision |
|---|---|---|---|---|
| oligodendrocyte | Mona / Allen | 0.9860 | **10** | **merged** → `Oligodendrocyte` |
| OPC | Mona / Allen | 0.9813 | **7** | **merged** → `OPC` |
| vSMC | Mona / Allen | 0.8777 | 42 | compete |
| microglia | HRCA / Mona / Allen | 0.9609 | 58 | compete |
| fibroblast | Mona / Allen | 0.9130 | 62 | compete |
| endothelial | Mona / Allen | 0.8542 | 64 | compete |
| pericyte | Mona / Allen | 0.8418 | 69 | compete |
| Schwann (myelinating vs not) | Mona | 0.9599 | 100 | compete |
| pigment epithelium | HRCA / Mona | 0.9561 | 108 | compete |
| astrocyte | HRCA / Mona / Allen | 0.9625 | **117** | compete |

Only oligodendrocyte and OPC are genuinely indistinguishable at CosMx depth; keeping those
apart just halves a large population by dataset provenance. The merge spec is
`retina_type_merges.csv`, applied by `python/merge_reference_types.py`.

**The surviving duplicates are confounded, and the reference cannot un-confound them.**
The astrocyte triple has plenty of discriminating signal, but inspect what carries it:

- **NEAT1** is the single largest discriminator (Mona 7.9 / HRCA 4.2 / Allen 1.6 expected
  counts per cell). Per-source means over *all* types are 4.62 / 1.83 / 0.71 — a
  **dataset-wide offset**, not astrocyte biology. NEAT1 is nuclear-retained, so this is
  the single-nucleus-vs-whole-cell assay signature.
- **SLC1A2** is 5.6 in Allen astrocytes and 0.1 in HRCA/Mona. A 56× gap on a core
  pan-astrocyte transporter is hard to read as regional tuning; GJA1 and AQP4 lean the
  same way (Allen 2–4× higher) while GFAP is even across all three.
- **GPC5 / NRXN1** (4.6 / 5.2 in Allen, ~0 elsewhere) plausibly *are* real cortical
  astrocyte identity — both are documented cortical astrocyte genes, and ambient neuronal
  contamination does not explain them (SNAP25 is 0.1 in Allen astrocytes vs 0.9 in Allen
  neurons).

So a real source-level technical offset sits on top of a possibly-real regional signal,
and the reference alone cannot separate them. That has two consequences:

1. **Do not read a flat InSituType astrocyte-copy call as regional identity.** Diagnostic:
   if copy assignment tracks `Region`, the distinction is plausibly real; if it tracks
   **cell area / total counts / nuclear fraction**, it is the NEAT1 assay offset sorting
   cells by segmentation quality. Only the second failure mode looks like success.
2. **The regional question is better answered downstream, on our own tissue.** "Do optic
   nerve astrocytes differ from brain astrocytes?" is a differential-expression test across
   `Region` *within* the astrocyte call — which uses our tissue rather than three labs'
   protocols, has real spatial ground truth, and can discover regional programs instead of
   being limited to what those atlases captured.

### Fine types within a source are also unresolvable at 6k depth

OFF- vs ON-midget ganglion = **0.9985**, cone vs S-cone = 0.9955, H1 vs H2 horizontal =
0.9953, Allen upper- vs deep-layer IT = **0.9969**. Unlike the cross-source duplicates
these are genuine biology from a single dataset, with no inter-atlas confound — they are
simply hard to call on ~400 counts.

### Why the hierarchy, and what it is for

With 64 near-collinear types competing flat, the EM is confident about none of them and
low-count cells fall through into the de novo sink — the mechanism behind the 66–74%
single-cluster collapse in the AtoMx Refit=FALSE calls. `retina_hierarchy.json` addresses
that, but note **what** it buys, which is not "collapsing the astrocytes":

- The **broad call is made first and is robust** — `Astrocyte` vs `Retinal_neuron` vs
  `Vascular` is an easy decision even at 400 counts.
- The **hard call is deferred to its own node, with a readable posterior.** The astrocyte
  triple only ever competes against itself; ON- vs OFF-midget only inside
  `Ganglion → Midget_ganglion`. If those posteriors come out near-uniform, that *is* the
  finding — the distinction is unsupported at this depth, and you collapse to the parent
  without ever having corrupted the broad call. Averaging up front commits to that
  conclusion before there is any evidence for it.
- Cells that fit nothing surface as **low posterior**, not as a de novo cluster.

The shape is 11 top-level branches over 64 leaves, depth 4:

```
Retinal_neuron (27)  Photoreceptor · Bipolar{Rod,Midget,Diffuse,Other} ·
                     Amacrine · Horizontal · Ganglion{Midget,Parasol,Other}
Cortical_neuron (10) Excitatory · Inhibitory · Striatal
Vascular (6)         Endothelial · Mural
Myeloid (4)          microglia ×3 + macrophage
Astrocyte (3) · Oligodendrocyte_lineage (3) · Lymphoid (3) · Pigmented (3)
Peripheral_glia (2) · Fibroblast (2) · Mueller_glia (1)
```

Single-leaf branches (`Mueller_glia`) mirror the GBM hierarchy's proven
`"Neuron": ["Neuron"]`. If `runInSituTree` objects to the depth-4 sub-branches, flatten
`Bipolar` and `Ganglion` to single nodes first — that costs only the isolation of the
hardest pairs.

**The hierarchy currently covers the 64 named types only.** InSituTree does no platform
rescaling, so its `--reference` must be CosMx-scale profiles from one InSituType
`RESCALE=true` run (`python/prep_insitutree_profiles.py`), and the kept de novo leaves must
be added to the hierarchy **in lockstep** with that profile matrix. Until then this file is
not runnable against `85_insitutree.sh`.
