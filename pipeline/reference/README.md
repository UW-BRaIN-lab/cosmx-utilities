# Stage 4 reference profiles (CZI GBmap)

InSituType (stage 4b) types cells semi-supervised against a genes × cell-types
reference profile matrix. We use the **CZI GBmap "extended"** atlas at **annotation
level 3** (21 cell types — `AC-like`, `Astrocyte`, `B_cell`, `CD4_CD8`, `DC`,
`Endothelial`, `MES-like`, `Mast`, `Mono`, `Mural_cell`, `NK`, `NPC-like`, `Neuron`,
`Neutrophil`, `OPC`, `OPC-like`, `Oligodendrocyte`, `Plasma_B`, `RG`, `TAM-BDM`,
`TAM-MG`), values = **raw average expression** per cell type.

## Files

- `gbmap_extended_level3_panel.csv` — **current reference.** The Extended GBmap level-3
  profiles **restricted to the CosMx 6k panel genes** (5946 genes × 21 cell types).
  Genes are the row index (HGNC symbols), cell types are columns. This is the committed,
  reproducible artifact and the input to `R/insitutype_typing.R` (the default
  `REFERENCE_BASENAME` in `slurm/72_anchor_typing.sh` and `slurm/80_insitutype.sh`).
- `gbmap_level4_panel.csv` — **previous reference**, kept for comparison. The core GBmap
  level-4 profiles restricted to the panel (5962 genes × 54 finer types). Select it by
  setting `REFERENCE_BASENAME=gbmap_level4_panel.csv` at submit time.

The full ~28k-gene GBmap exports are **not** committed (~23 MB each; mix HGNC symbols and
Ensembl IDs). The committed CSVs are the panel-restricted subsets only.

## Provenance / how to regenerate

Source (local, not in this repo):
`~/keene-lab/GBM/GBmap/extended_GBmap_results/extended_GBmap_gene_by_annotation_level_3_raw_avg_expression.csv`

Restrict to the CosMx panel (run once, locally):

```bash
uv run python pipeline/python/prep_insitutype_reference.py \
    --gbmap-csv ~/keene-lab/GBM/GBmap/extended_GBmap_results/extended_GBmap_gene_by_annotation_level_3_raw_avg_expression.csv \
    --panel-h5ad <a per-slide stage-1 .h5ad or the stage-3a combined_qc.h5ad> \
    --output pipeline/reference/gbmap_extended_level3_panel.csv
```

Then upload to Kopah so the Slurm typing job can stage it:

```bash
set -a; source pipeline/.env; set +a
AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID" \
AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY" \
S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL" \
    s5cmd cp pipeline/reference/gbmap_extended_level3_panel.csv \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/gbmap_extended_level3_panel.csv"
```

The other GBmap atlases (`core_GBmap_results/`) and annotation granularities
(`..._by_cell_type_...`, `..._level_4_...`) exist in the same source tree; swap the
`--gbmap-csv` and `--output` to build a reference at a different level.
