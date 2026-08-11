#!/usr/bin/env Rscript
################################################################################
# Build a custom InSituType reference profile from one or more scRNA-seq sources.
#
# Generalized from the ocular/brain atlas project (HRCA v2 + Monavarfeshani +
# Allen Brain Cell Atlas -> one genes x cell-types reference). Any dataset that
# provides raw counts + per-cell type labels can be a source. Two source types:
#
#   type = "pseudobulk_csv"  a genes x cell-types CSV from
#                            preprocess_h5ad_profiles.py (large CELLxGENE H5ADs
#                            preprocessed in Python). Loaded as-is.
#   type = "mtx"             a GEO-style sparse matrix (MTX + gene/cell name
#                            files) + a metadata table with a cell-type column.
#                            Profiles are derived here via
#                            InSituType::getRNAprofiles(neg = 0) — no log
#                            transform, raw counts only.
#
# Each source is rescaled independently to the 99th-percentile convention, genes
# are intersected across sources, and columns are cbind-ed into the final linear-
# scale matrix for InSituType semi-supervised typing. Cell-type names are
# prefixed per source to avoid collisions (e.g. two "astrocyte" columns).
#
# Usage:
#   Rscript build_reference_profile.R <config.R>
# where <config.R> defines `sources`, `output_path`, and optional
# `min_cells_per_type`, `rescale_quantile`, `rescale_factor`, `markers`.
# See examples/retina_brain_ocular.R for a worked configuration.
################################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(InSituType)
})

DEFAULT_MIN_CELLS <- 15
DEFAULT_RESCALE_QUANTILE <- 0.99
DEFAULT_RESCALE_FACTOR <- 1000

# --- config -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript build_reference_profile.R <config.R>")
}
config_env <- new.env()
sys.source(normalizePath(args[[1]]), envir = config_env)

sources <- get0("sources", envir = config_env)
output_path <- get0("output_path", envir = config_env)
if (is.null(sources) || is.null(output_path)) {
  stop("config must define `sources` (a list) and `output_path` (a string).")
}
min_cells_per_type <- get0("min_cells_per_type", envir = config_env,
                           ifnotfound = DEFAULT_MIN_CELLS)
rescale_quantile <- get0("rescale_quantile", envir = config_env,
                         ifnotfound = DEFAULT_RESCALE_QUANTILE)
rescale_factor <- get0("rescale_factor", envir = config_env,
                       ifnotfound = DEFAULT_RESCALE_FACTOR)
markers <- get0("markers", envir = config_env, ifnotfound = list())

# --- helpers ------------------------------------------------------------------

#' Derive genes x cell-types profiles from raw counts + labels via InSituType.
#' @param counts genes x cells (dgCMatrix or matrix), linear-scale raw counts.
#' @param clust  character vector of cell-type labels, length = ncol(counts).
compute_profiles <- function(counts, clust, min_cells, prefix = "") {
  keep <- !is.na(clust) & clust != ""
  counts <- counts[, keep, drop = FALSE]
  clust <- clust[keep]

  type_counts <- table(clust)
  rare <- names(type_counts[type_counts < min_cells])
  if (length(rare) > 0) {
    message(sprintf("  dropping %d types with < %d cells: %s",
                    length(rare), min_cells, paste(rare, collapse = ", ")))
    keep <- !(clust %in% rare)
    counts <- counts[, keep, drop = FALSE]
    clust <- clust[keep]
  }

  # getRNAprofiles wants cells x genes; scRNA-seq has no platform background,
  # so neg = 0 for every cell.
  profiles <- getRNAprofiles(x = Matrix::t(counts), clust = clust,
                             neg = rep(0, length(clust)))
  if (nchar(prefix) > 0) colnames(profiles) <- paste0(prefix, colnames(profiles))
  message(sprintf("  computed %d types x %d genes", ncol(profiles), nrow(profiles)))
  profiles
}

#' Rescale so the 99th-percentile entry maps to `factor`. InSituType ignores
#' column scaling, but a common scale keeps merged sources comparable.
rescale_profiles <- function(profiles, quantile_p, factor) {
  q <- stats::quantile(profiles, quantile_p, na.rm = TRUE)
  if (is.finite(q) && q > 0) profiles <- profiles / as.numeric(q) * factor
  profiles
}

#' Load a source into a rescaled genes x cell-types matrix.
load_source <- function(src) {
  prefix <- if (is.null(src$prefix)) "" else src$prefix
  message(sprintf("Source '%s' (type=%s)", prefix, src$type))

  if (identical(src$type, "pseudobulk_csv")) {
    m <- as.matrix(read.csv(src$path, row.names = 1, check.names = FALSE))
    if (nchar(prefix) > 0 && !all(startsWith(colnames(m), prefix))) {
      colnames(m) <- paste0(prefix, colnames(m))
    }
    message(sprintf("  loaded %d genes x %d types", nrow(m), ncol(m)))

  } else if (identical(src$type, "mtx")) {
    counts <- readMM(src$mtx)
    genes <- fread(src$genes)[[1]]
    cells <- fread(src$cells)[[1]]
    rownames(counts) <- genes
    colnames(counts) <- cells

    meta <- fread(src$meta)
    # Some GEO/SCP exports carry a second header ("TYPE"/"group") row under the
    # column names — drop it if present so labels align to real cells.
    if (!is.null(src$meta_drop_type_row) && isTRUE(src$meta_drop_type_row)) {
      bad <- which(meta[[1]] == "TYPE" |
                     (ncol(meta) > 1 & apply(meta, 1, function(r) any(r == "group"))))
      if (length(bad) > 0) meta <- meta[-bad, ]
    }
    id_col <- if (is.null(src$meta_id_col)) names(meta)[1] else src$meta_id_col
    meta <- meta[match(colnames(counts), meta[[id_col]]), ]
    if (!all(meta[[id_col]] == colnames(counts))) {
      stop(sprintf("  metadata ids do not align to matrix columns for '%s'", prefix))
    }
    m <- compute_profiles(counts, clust = as.character(meta[[src$cell_type_col]]),
                          min_cells = min_cells_per_type, prefix = prefix)
    rm(counts, meta); gc()

  } else {
    stop(sprintf("unknown source type '%s' (use 'pseudobulk_csv' or 'mtx')", src$type))
  }

  rescale_profiles(m, rescale_quantile, rescale_factor)
}

# --- build --------------------------------------------------------------------
message(sprintf("%s, building reference from %d sources", Sys.time(), length(sources)))
profile_list <- lapply(sources, load_source)

shared_genes <- Reduce(intersect, lapply(profile_list, rownames))
message(sprintf("Shared genes across all sources: %d", length(shared_genes)))
if (length(shared_genes) == 0) stop("no genes shared across sources — check gene naming (symbols vs Ensembl).")

combined <- do.call(cbind, lapply(profile_list, function(m) m[shared_genes, , drop = FALSE]))
message(sprintf("Combined reference: %d genes x %d cell types",
                nrow(combined), ncol(combined)))

# --- marker sanity checks -----------------------------------------------------
if (length(markers) > 0) {
  message("Marker sanity checks (top 3 types by expression):")
  for (gene in names(markers)) {
    if (gene %in% rownames(combined)) {
      top <- head(sort(combined[gene, ], decreasing = TRUE), 3)
      message(sprintf("  %-10s (%-28s) %s", gene, markers[[gene]],
                      paste(sprintf("%s=%.1f", names(top), top), collapse = ", ")))
    } else {
      message(sprintf("  %-10s NOT in shared gene set", gene))
    }
  }
}

# --- write --------------------------------------------------------------------
out <- data.table(gene = rownames(combined))
out <- cbind(out, as.data.table(combined))
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
fwrite(out, output_path)
message(sprintf("Wrote %s (%d genes x %d cell types)",
                output_path, nrow(combined), ncol(combined)))
