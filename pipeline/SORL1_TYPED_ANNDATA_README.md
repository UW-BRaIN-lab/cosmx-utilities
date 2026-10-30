# Maddie's SORL1 MTG cohort: typed AnnData (InSituTree v5)

`stage4_insitutree_k21_v5/cosmx_typed.h5ad` under the cohort's Kopah prefix (`cosmx-maddie-sorl1`):
`s3://$KOPAH_BUCKET/cosmx-maddie-sorl1/stage4_insitutree_k21_v5/cosmx_typed.h5ad`. About 1.7 GB (gzip).

793,380 cells that passed QC, 40 donors, 10 slides, with the InSituTree v5 cell types in `obs`. It is the stage-3 clustered AnnData
plus the typing, so it holds expression, Leiden clusters and the UMAP as well.

## Load it

```bash
# needs the lab's Kopah keys (S3_ENDPOINT_URL, AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY)
s5cmd cp "s3://$KOPAH_BUCKET/cosmx-maddie-sorl1/stage4_insitutree_k21_v5/cosmx_typed.h5ad" .
```

```python
import anndata as ad
obs = ad.read_h5ad("cosmx_typed.h5ad", backed="r").obs   # labels only, no expression loaded
adata = ad.read_h5ad("cosmx_typed.h5ad")                 # everything (several GB in memory)
```

## What is in it

- **`X`**: sparse counts, QC-passed cells by all probes. It is the stage-3a counts matrix; check `adata.X.max()` and that the values are
  integers before assuming raw counts.
- **`var["probe_type"]`**: `gene` or `negprobe`. Subset to `gene` for expression analysis; the negative probes are kept because
  InSituType needs them.
- **`obsm`**: `X_pearson_pca` (the batch-aware PCA the clusters were built on) and `X_umap`.
- **`obs`**, the columns to know:

| column | meaning |
|---|---|
| `Region` | donor id. This is the batch. **Slide names are not donors** (each slide holds four donors, and two slide names carry a typo: the donors are 526 and 6803) |
| `slide_id`, `fov` | slide and field of view |
| `case_group`, `case_broad`, `sorl1_mutation` | donor groups from the curated FOV sheet |
| `leiden` | cohort-wide Leiden cluster from stage 3c |
| `total_counts`, `qc_genes_detected`, `qc_area` | depth, genes detected and cell area |
| **`cell_type_v5`** | **the 30 InSituTree v5 leaves with the final names. Use this for cell-level work** |
| **`cell_type_v5_collapsed`** | analysis level: the six neuron sink leaves are `Neuron_low_signal`; `L4 neuron` is `Excitatory_L4_5_IT_like`; `r` and `s` are `Excitatory_L5_IT_like` |
| `cell_type_v5_parent` | branch level (`Astrocyte`, `Microglia`, `Excitatory_deep`, `Vascular`, `OPC`, ...) |
| **`opc_confident`** | called OPC **and** (in Leiden 0 or `total_counts` >= 600). Use it for any OPC result |
| `cell_type`, `insitutype_prob` | the run's own label and its posterior. `cell_type` still says `Excitatory_upper`; prefer `cell_type_v5` |

## Read this before comparing groups

These cause most wrong conclusions in this dataset.

1. **Depth differs a lot by donor.** Median counts per cell range about 6-fold across slides. Four SORL1 donors (599, 7068, 7796, 6589) are
   shallow and put 37-48% of their cells (those with at least 600 counts) in `Low_signal`; two other SORL1 donors (7137, 7181) do not. The
   neuron and `Low_signal` differences between SORL1 and other donors sit in those four donors and cannot be separated from RNA quality.
2. **Run effect.** The July slides (`20260708*`) hold all the SORL1 donors, so SORL1 against the other groups is confounded with the run.
   Compare groups **within** a slide where you can;
   across runs the sink share changes about 3-fold (non-SORL1 donors: 17-21% on the July slides, 5-7% elsewhere).
3. **Use donors, not cells, as the unit.** Per-donor percentages, then group tests; cells within a donor are not independent.
4. **OPC calls on shallow cells are unreliable.** The posterior is near 1 for them. Use `opc_confident`. The OPC elevation in SORL1 is
   largely a run effect.
5. **Some labels are mixed groups, not cell types.**
   - `Excitatory_L2_3_mixed` (about 27k cells): only about 11% is the flat K=21 `L2/3 neuron` call; the rest comes from unresolved clusters. For
     pure L2/3 use the `L2/3 neuron` leaf. Renamed from `Excitatory_upper`.
   - `Astrocyte_low_signal` (about 25k cells) is a sink for cluster `a`, a graded astrocyte/neuron hybrid, not an astrocyte state.
   - The `*_unresolved_*_denovo` leaves and `Low_signal*` are unresolved cells with no biological identity.
   - `L4 neuron` marks layer 5 spatially (Maddie's review), so `Excitatory_L4_5_IT_like` overstates L4. No cluster is a confirmed L4 population.
6. **What held up in the disease comparisons.** Microglia B is higher with disease in every design, but within the July slides two donors (599,
   7068) carry it. Nothing else is firm. The comparisons, tables and caveats are in `pipeline/SORL1_40DONOR_TYPING.md`
   ("Disease comparisons on v5").

## Status of the labels

v5 is the working typing and the deep-excitatory names are confirmed by marker overlap and Maddie's spatial check. It is not formally
declared final: Maddie has reviewed the K=21 cluster identities in Napari but not yet the tree-level labels, and a v6 is possible if that
review finds a spatial pattern. Treat the column names (`*_v5`) as the version stamp. How the typing was built, how it changed from v1 to
v5, and what is still open: `pipeline/SORL1_40DONOR_TYPING.md`.
