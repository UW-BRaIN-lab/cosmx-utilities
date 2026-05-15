#!/usr/bin/env Rscript
# Extract clinical annotations from one CosMx Seurat .RDS file.
#
# The Case / Block / Region columns in meta.data were added manually in AtoMx
# to map each FOV to its tissue source (each slide carries two tissues). These
# annotations are NOT present in the CosMx flat files, so we extract them once
# into a sidecar CSV that downstream stages left-join into AnnData obs. After
# this runs for all slides, the .RDS files can be archived.
#
# Usage:
#   Rscript pipeline/r/extract_clinical_annotations.R \
#     --input /path/to/seuratObject_<slide>.RDS \
#     --slide-id 7134A77439A6 \
#     --output /path/to/<slide>_clinical.csv

suppressPackageStartupMessages({
  library(optparse)
  library(SeuratObject)
})

option_list <- list(
  make_option("--input",    type = "character", help = "Path to Seurat .RDS file"),
  make_option("--slide-id", type = "character", help = "Slide identifier"),
  make_option("--output",   type = "character", help = "Path to write CSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("input", "slide-id", "output")
missing_args <- setdiff(required_args, names(opt)[!sapply(opt, is.null)])
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

cat("Reading", opt$input, "\n")
obj <- readRDS(opt$input)
md <- obj@meta.data

REQUIRED_COLS <- c("cell_ID", "fov", "Case", "Block", "Region",
                   "cellSegmentationSetName")
missing_cols <- setdiff(REQUIRED_COLS, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing columns in meta.data: ", paste(missing_cols, collapse = ", "))
}

out <- data.frame(
  slide_id              = opt$`slide-id`,
  cell_ID               = as.character(md$cell_ID),
  fov                   = md$fov,
  case                  = md$Case,
  block                 = md$Block,
  region                = md$Region,
  cell_segmentation_set = md$cellSegmentationSetName,
  stringsAsFactors      = FALSE
)

cat("Writing", nrow(out), "cells to", opt$output, "\n")
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
write.csv(out, opt$output, row.names = FALSE)
