# Stage 4 reference profiles (CZI GBmap)

InSituType (stage 4b) types cells semi-supervised against a genes × cell-types
reference profile matrix. We use the **CZI GBmap** "core" atlas at **annotation
level 4** (the finest level: 54 cell types — e.g. `AC-like`, `MES-like_hypoxia_MHC`,
`NPC-like_neural`, `TAM`/`Mono` subsets, `Oligodendrocyte`, `Endothelial`, lymphoid
subsets), values = **raw average expression** per cell type.

## Files

- `gbmap_level4_panel.csv` — the GBmap level-4 reference **restricted to the CosMx 6k
  panel genes**. Genes are the row index (HGNC symbols), cell types are columns. This
  is the committed, reproducible artifact and the input to `R/insitutype_typing.R`.

The full ~28k-gene GBmap export is **not** committed (~23 MB; mixes HGNC symbols and
Ensembl IDs). The committed CSV is the panel-restricted subset only.

## Provenance / how to regenerate

Source (local, not in this repo):
`~/keene-lab/GBM/GBmap/core_GBmap_results/core_GBmap_gene_by_annotation_level_4_raw_avg_expression.csv`

Restrict to the CosMx panel (run once, locally):

```bash
uv run python pipeline/python/prep_insitutype_reference.py \
    --gbmap-csv ~/keene-lab/GBM/GBmap/core_GBmap_results/core_GBmap_gene_by_annotation_level_4_raw_avg_expression.csv \
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
exist in the same source directory; swap the `--gbmap-csv` to retype at a coarser level.
