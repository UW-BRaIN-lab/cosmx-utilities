################################################################################
# Worked example config for build_reference_profile.R
#
# The ocular/brain reference that this toolkit was generalized from: a combined
# reference for CosMx typing of retina / optic nerve / brain, built from three
# public single-cell datasets. Run:
#
#   Rscript ../build_reference_profile.R retina_brain_ocular.R
#
# Prerequisites (produce the two pseudobulk CSVs first — see README):
#   preprocess_h5ad_profiles.py --h5ad HRCA_v2_snRNA.h5ad --cell-type-col cell_type \
#       --prefix HRCA_ --output HRCA_v2_pseudobulk_profiles.csv
#   preprocess_h5ad_profiles.py --h5ad Astrocyte.h5ad --label Astrocyte \
#       --h5ad Microglia.h5ad --label Microglia ... (one --label per Allen file) \
#       --prefix Allen_ --output Allen_pseudobulk_profiles.csv
#
# Datasets:
#   1. HRCA v2 (Li et al. 2025)          retina snRNA-seq       CELLxGENE H5AD
#   2. Monavarfeshani et al. 2023        posterior eye segment  GEO GSE236566 MTX
#                                        (ON, ONH, sclera, choroid)  + SCP2298 meta
#   3. Allen Brain Cell Atlas v1.0       cortex / white matter  CELLxGENE H5AD
#      (18 superclusters + Vascular)     / vascular / fibroblast
################################################################################

DATA_DIR <- "/mnt/cosmx"   # <-- edit to wherever the inputs live

sources <- list(
  # HRCA v2 — snRNA-seq retina. Too large for R; preprocessed in Python.
  list(
    type   = "pseudobulk_csv",
    path   = file.path(DATA_DIR, "HRCA_v2_pseudobulk_profiles.csv"),
    prefix = "HRCA_"
  ),

  # Monavarfeshani 2023 — posterior eye. Small enough to profile in R directly
  # from the GEO MTX via getRNAprofiles. The gene/cell name CSVs carry a header
  # row ("x"); fread reads it as the first column, which the [[1]] extract uses.
  # SCP_meta.csv has a second "TYPE"/"group" header row -> drop it.
  list(
    type               = "mtx",
    mtx                = file.path(DATA_DIR, "GSE236566_Count_matrix.mtx"),
    genes              = file.path(DATA_DIR, "GSE236566_Count_matrix_genenames.csv"),
    cells              = file.path(DATA_DIR, "GSE236566_Count_matrix_cellnames.csv"),
    meta               = file.path(DATA_DIR, "SCP_meta.csv"),
    meta_id_col        = "NAME",
    meta_drop_type_row = TRUE,
    cell_type_col      = "cell_type__ontology_label",
    prefix             = "Mona_"
  ),

  # Allen Brain Cell Atlas — cortex/white-matter/vascular. Preprocessed in Python
  # (18 superclusters at supercluster level + Vascular at cell_type level).
  list(
    type   = "pseudobulk_csv",
    path   = file.path(DATA_DIR, "Allen_pseudobulk_profiles.csv"),
    prefix = "Allen_"
  )
)

output_path <- file.path(DATA_DIR, "HRCA-Monavarfeshani-Allen-combined-profile.csv")

min_cells_per_type <- 15
rescale_quantile   <- 0.99
rescale_factor     <- 1000

# Optional marker enrichment sanity checks (gene -> expected cell type note).
markers <- list(
  RHO     = "rods",
  OPN1SW  = "S cones",
  RPE65   = "RPE",
  RLBP1   = "Mueller glia / RPE",
  MBP     = "myelinating oligos",
  GFAP    = "astrocytes",
  COL1A1  = "fibroblasts",
  PECAM1  = "endothelium",
  SLC17A7 = "excitatory neurons",
  GAD2    = "inhibitory neurons",
  AQP4    = "astrocytes",
  SLC6A13 = "arachnoid barrier cells",
  DCN     = "leptomeningeal fibroblasts",
  PDGFRA  = "VLMC / OPCs"
)
