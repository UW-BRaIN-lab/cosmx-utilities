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

- `retina_combined_panel.csv` — the combined atlas restricted to the panel:
  **5,993 genes × 68 cell types** (97.1% of the 6,175 panel genes are in the atlas).
- `cosmx_6k_panel_genes.txt` — the 6,175 gene probes of the **Human RNA 6k Discovery**
  panel, one per line. Derived from an exprMat header with `fov`/`cell_ID` and the
  `Negative*` (20) / `SystemControl*` (324) control probes dropped. Lets
  `prep_insitutype_reference.py` panel-restrict a reference *before* stage 1 exists.
  The GBM and retina studies run the same panel, so this file serves both.

### Provenance / how to regenerate

Source (local, not in this repo):
`~/keene-lab/retina-brain/HRCA-Monavarfeshani-Allen-combined-profile.csv`
(21,040 genes × 68 types).

```bash
uv run python pipeline/python/prep_insitutype_reference.py \
    --reference-csv ~/keene-lab/retina-brain/HRCA-Monavarfeshani-Allen-combined-profile.csv \
    --panel-genes pipeline/reference/cosmx_6k_panel_genes.txt \
    --output pipeline/reference/retina_combined_panel.csv
```

Marker sanity check on the result (all top-ranked in the expected type): RHO → rod,
NEFL → ganglion, MBP/PLP1 → oligodendrocyte, AQP4 → astrocyte, SLC17A7 → cortical
excitatory, GAD1 → interneuron, P2RY12 → microglia, MLANA → melanocyte,
CLDN5 → endothelium.

### ⚠ This reference is strongly collinear — read before interpreting a flat run

Because each source atlas is prefixed rather than merged, the same cell type appears
once **per source**, at near-identical profiles (log1p cosine over the 5,993 panel genes):

| duplicated type | copies | pairwise cosine |
|---|---|---|
| astrocyte | HRCA / Mona / Allen | 0.96–0.98 |
| microglia | HRCA / Mona / Allen | 0.96–0.97 |
| oligodendrocyte + OPC | Mona ×2 / Allen ×2 | up to 0.986 |
| endothelial, fibroblast, pericyte, vSMC | Mona / Allen | 0.84–0.91 |

Fine types *within* a source are also unresolvable at 6k-plex depth: OFF- vs ON-midget
ganglion = **0.9985**, cone vs S-cone = 0.9955, H1 vs H2 horizontal = 0.9953, Allen
upper- vs deep-layer IT = **0.9969**.

Consequences:

1. A **flat** InSituType run splits one biological population across its duplicate
   columns by dataset provenance, not biology — the AtoMx run scattered astrocytes over
   `Mona_astrocyte` (9,199), `HRCA_astrocyte` (1,882) and `Allen_Astrocyte` (1,500) on
   one slide. With 68 near-collinear types the EM is confident about none of them, and
   low-count cells fall through into the de novo sink. This is the mechanism behind the
   66–74% single-de-novo-cluster collapse seen in the Refit=FALSE calls.
2. **InSituTree fixes this properly**, and is the intended endpoint: branch-local
   competition means OFF-midget only competes with its ganglion siblings and the
   astrocyte triple only with each other. See `insitutree_hierarchy.json` for the GBM
   shape; the retina hierarchy is designed after de novo triage.
3. The flat InSituType pass is still **required** — InSituTree does no platform
   rescaling, so its profiles must come from one InSituType `RESCALE=true` run.

So: read a flat retina run as *"which broad population, plus what de novo structure"*,
not as a final fine-type call.
