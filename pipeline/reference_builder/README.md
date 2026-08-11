# Reference builder: scRNA-seq → InSituType reference profiles

A general, dataset-agnostic toolkit for turning **any** single-cell/nucleus
RNA-seq dataset that has **cell-type labels** into a linear-scale
`genes × cell-types` reference profile matrix for **InSituType semi-supervised
cell typing** of CosMx data. Merge several datasets into one combined reference.

This is the reusable version of the ocular/brain atlas reference (HRCA v2 +
Monavarfeshani + Allen Brain Cell Atlas). It is **not** GBM- or meningioma-
specific — point it at whatever atlas matches your tissue.

> Distinct from `pipeline/python/prep_insitutype_reference.py`, which only
> *restricts an already-built* reference matrix to the CosMx panel genes. This
> toolkit *builds* the matrix from raw counts. A typical flow is: build here →
> panel-restrict with `prep_insitutype_reference.py` → type with
> `pipeline/R/insitutype_typing.R`.

## Files

- `preprocess_h5ad_profiles.py` — large CELLxGENE `.h5ad` → pseudobulk CSV.
  Backed (memory-mapped), chunked, so 3M-cell atlases don't blow up R's memory.
- `build_reference_profile.R` — config-driven driver. Loads each source, derives
  profiles, rescales, intersects genes, merges, runs marker sanity checks, writes
  the final CSV.
- `examples/retina_brain_ocular.R` — a complete worked config (the ocular/brain
  reference) you can copy and adapt.

## Key principles

1. **Linear-scale raw counts only — never log-transformed.**
   `InSituType::getRNAprofiles(x, clust, neg)` expects `x` as cells × genes;
   `neg = 0` for scRNA-seq (no CosMx negative-probe background).
2. **Preprocess big H5ADs in Python.** R's `anndata::read_h5ad()` loads the whole
   matrix into memory and crashes on large files. `preprocess_h5ad_profiles.py`
   reads in `backed='r'` mode and sums counts in 50k-cell chunks.
3. **CELLxGENE stores Ensembl IDs as the index but HUGO symbols in
   `var['feature_name']`.** CosMx uses HUGO, so gene names come from
   `feature_name`; duplicate symbols are summed.
4. **Raw counts may not be in `X`.** The script probes `layers['raw']`,
   `layers['counts']`, `adata.raw`, then `X`, and validates that the chosen slot
   is non-negative integers before trusting it. (This guards the "we thought we
   had raw counts but didn't" failure.)
5. **Merging:** each source is rescaled independently
   (`profiles / quantile(profiles, 0.99) * 1000`), genes are intersected, columns
   are `cbind`-ed. Cell-type names are prefixed per source (`HRCA_`, `Allen_`, …)
   so shared names like `astrocyte` don't collide.
6. **Small GEO datasets (MTX + metadata)** are profiled directly in R via
   `getRNAprofiles` — no Python step needed (see the `mtx` source type).
7. **Sanity-check markers** after merging: known markers should be enriched in the
   expected types (RHO→rods, MBP→oligodendrocytes, SLC17A7→excitatory neurons).

## Two source types

`build_reference_profile.R` reads a config that lists `sources`, each one of:

| type             | input                                               | how profiles are made        |
|------------------|-----------------------------------------------------|------------------------------|
| `pseudobulk_csv` | a CSV from `preprocess_h5ad_profiles.py`            | loaded as-is                 |
| `mtx`            | GEO MTX + gene/cell name files + metadata CSV       | `getRNAprofiles(neg = 0)` in R |

## Usage

**Step 1 — preprocess any large H5AD atlas** (repeat per atlas):

Per-cell labels in `obs` (e.g. HRCA v2):

```bash
uv run python pipeline/reference_builder/preprocess_h5ad_profiles.py \
    --h5ad HRCA_v2_snRNA.h5ad --cell-type-col cell_type \
    --prefix HRCA_ --output HRCA_v2_pseudobulk_profiles.csv
```

One label per file (e.g. Allen superclusters):

```bash
uv run python pipeline/reference_builder/preprocess_h5ad_profiles.py \
    --h5ad Astrocyte.h5ad --label Astrocyte \
    --h5ad Microglia.h5ad --label Microglia \
    --prefix Allen_ --output Allen_pseudobulk_profiles.csv
```

**Step 2 — write a config** (copy `examples/retina_brain_ocular.R`) describing
your sources and `output_path`.

**Step 3 — build the combined reference:**

```bash
Rscript pipeline/reference_builder/build_reference_profile.R my_config.R
```

The output `genes × cell-types` CSV is ready for InSituType semi-supervised
typing (optionally panel-restrict it first with `prep_insitutype_reference.py`).

## References

- InSituType hybrid reference workflow:
  <https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/hybrid-reference-profiles/>
- `InSituType::getRNAprofiles()` — derives profiles from raw scRNA-seq counts.
