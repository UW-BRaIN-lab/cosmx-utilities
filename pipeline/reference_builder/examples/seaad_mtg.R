################################################################################
# Worked example config for build_reference_profile.R: SEA-AD MTG reference.
#
# Single-source reference for CosMx typing of human middle temporal gyrus (MTG),
# built from the SEA-AD "Whole Taxonomy - MTG" atlas (1,378,211 nuclei,
# CELLxGENE collection 1ca90a2d-2943-483d-b678-b809bf464c30). Run:
#
#   Rscript pipeline/reference_builder/build_reference_profile.R \
#       pipeline/reference_builder/examples/seaad_mtg.R
#
# Prerequisite: pipeline/slurm/76_seaad_mtg_reference.sh (downloads the H5AD and
# writes SEAAD_MTG_<LABEL_COL>_pseudobulk_profiles.csv).
#
# SEAAD_DIR is read from the environment so no path is committed.
################################################################################

SEAAD_DIR <- Sys.getenv("SEAAD_DIR")
LABEL_COL <- Sys.getenv("LABEL_COL", unset = "Supertype")
if (!nzchar(SEAAD_DIR)) stop("set SEAAD_DIR to the directory holding the pseudobulk CSV.")

sources <- list(
  list(
    type   = "pseudobulk_csv",
    path   = file.path(SEAAD_DIR, sprintf("SEAAD_MTG_%s_pseudobulk_profiles.csv", LABEL_COL)),
    prefix = "SEAAD_"
  )
)

output_path <- file.path(SEAAD_DIR, sprintf("SEAAD_MTG_%s.profiles.csv", LABEL_COL))

markers <- list(
  SLC17A7 = "excitatory neurons",
  GAD1    = "inhibitory neurons",
  MBP     = "oligodendrocytes",
  GFAP    = "astrocytes",
  CSF1R   = "microglia",
  PDGFRA  = "OPC",
  CLDN5   = "endothelial"
)
