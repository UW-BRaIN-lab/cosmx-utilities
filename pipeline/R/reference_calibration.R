#!/usr/bin/env Rscript
# Which reference cell types did InSituType actually calibrate, and which profiles came out
# degenerate? Reads the full .rds a fit saves -- a join, not a re-run.
#
# `update_reference_profiles = TRUE` makes updateReferenceProfiles() PIN cells: for each
# reference type it selects the cells that best match that profile (n_anchor_cells,
# min_anchor_cosine, min_anchor_llr) and rescales the profile toward them. A type that
# attracts no cell passing those bars is never calibrated -- the platform correction reaches
# it only through whatever global term the fit shares, and its fitted profile can end up at a
# few expected counts where the data sits in the hundreds. Such a type cannot win a cell, so
# it is invisible in the output: the ONLY place it shows up is here.
#
# Two numbers per reference type:
#   n_pinned          cells pinned to it. 0 = never calibrated.
#   profile_sum       total expected counts of its fitted profile. Compare against the de novo
#                     columns in the same fit, which ARE estimated from the data and so give
#                     the fit's own yardstick for a profile at the right scale.
#
# Terminology: these are InSituType's PINNED CELLS, not the anchor cohort and not the mural
# types -- see CLAUDE.md. The column is called `anchors` in the fit object, which is why it
# appears under that name below.
#
# Usage:
#   reference_calibration.R --result-rds insitutype_result.rds \
#       --reference retina_combined_panel.csv --output-csv reference_calibration.csv

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument: %s", a))
    key <- sub("^--", "", a)
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

# log1p cosine: the same quantity the collinearity table reports, so the two diagnostics can
# be read side by side. Raw cosine between expression profiles is dominated by the few
# highest-expressed genes; log1p spreads the weight without centring.
cosine_log1p <- function(a, b) {
  a <- log1p(a); b <- log1p(b)
  denom <- sqrt(sum(a * a)) * sqrt(sum(b * b))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(a * b) / denom
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--result-rds is required" = !is.null(opt[["result-rds"]]))
stopifnot("--reference is required" = !is.null(opt$reference))
out_csv <- opt[["output-csv"]]

res <- readRDS(opt[["result-rds"]])
cat("objects in the fit:", paste(names(res), collapse = ", "), "\n\n")
stopifnot("the fit has no $profiles" = !is.null(res$profiles))

fitted <- res$profiles
clust <- res$clust

ref_dt <- data.table::fread(opt$reference)
ref <- as.matrix(ref_dt[, -1L])
rownames(ref) <- as.character(ref_dt[[1]])
storage.mode(ref) <- "double"

shared_genes <- intersect(rownames(fitted), rownames(ref))
stopifnot("fit and reference share no genes" = length(shared_genes) > 0)
cat(sprintf("fit: %d genes x %d columns   reference: %d genes x %d types   shared genes: %d\n",
            nrow(fitted), ncol(fitted), nrow(ref), ncol(ref), length(shared_genes)))

reference_types <- colnames(ref)
fitted_types <- colnames(fitted)
named <- intersect(fitted_types, reference_types)
denovo <- setdiff(fitted_types, reference_types)
absent <- setdiff(reference_types, fitted_types)
cat(sprintf("fit carries %d of the reference's %d types, plus %d de novo columns\n",
            length(named), length(reference_types), length(denovo)))
if (length(absent)) {
  cat(sprintf("%d reference type%s ABSENT from the fitted profiles: %s\n",
              length(absent), if (length(absent) == 1) " is" else "s are",
              paste(absent, collapse = ", ")))
}

# --- pinned cells -------------------------------------------------------------
pinned_counts <- integer(0)
if (is.null(res$anchors)) {
  cat("\nNO $anchors IN THE FIT -- it ran with update_reference_profiles=FALSE, so no type\n",
      "was calibrated from pinned cells and the profile scale below is the reference's own.\n",
      sep = "")
} else {
  anchors <- res$anchors[!is.na(res$anchors)]
  pinned_counts <- table(as.character(anchors))
  cat(sprintf("\npinned cells: %d of %d (%.2f%%), spread over %d of %d reference types\n",
              length(anchors), length(clust), 100 * length(anchors) / length(clust),
              length(pinned_counts), length(reference_types)))
}

called_counts <- table(as.character(clust))

per_type <- function(type) {
  in_fit <- type %in% fitted_types
  in_ref <- type %in% reference_types
  sum_fitted <- if (in_fit) sum(fitted[, type]) else NA_real_
  sum_ref <- if (in_ref) sum(ref[shared_genes, type]) else NA_real_
  cos <- if (in_fit && in_ref) {
    cosine_log1p(fitted[shared_genes, type], ref[shared_genes, type])
  } else NA_real_
  list(type = type,
       source = if (in_ref) "reference" else "denovo",
       n_pinned = as.integer(pinned_counts[type] %||% 0L),
       n_called = as.integer(called_counts[type] %||% 0L),
       profile_sum_fitted = sum_fitted,
       profile_sum_reference = sum_ref,
       fitted_over_reference = if (is.na(sum_ref) || sum_ref == 0) NA_real_ else sum_fitted / sum_ref,
       cosine_to_reference = cos)
}

rows <- rbindlist(lapply(union(fitted_types, reference_types), per_type))
rows[is.na(n_pinned), n_pinned := 0L]
rows[is.na(n_called), n_called := 0L]

denovo_median <- rows[source == "denovo" & !is.na(profile_sum_fitted),
                      median(profile_sum_fitted)]
rows[, profile_sum_vs_denovo := profile_sum_fitted / denovo_median]
setorder(rows, source, n_pinned, profile_sum_fitted)

# --- report -------------------------------------------------------------------
cat("\n================ reference types, least calibrated first ================\n")
print(rows[source == "reference",
           .(type, n_pinned, n_called, profile_sum = round(profile_sum_fitted, 1),
             vs_denovo = round(profile_sum_vs_denovo, 3),
             cos_to_ref = round(cosine_to_reference, 4))], nrows = 100)

if (length(denovo)) {
  cat("\n================ de novo columns (the fit's own yardstick) ================\n")
  print(rows[source == "denovo",
             .(type, n_called, profile_sum = round(profile_sum_fitted, 1))], nrows = 100)
  cat(sprintf("\nmedian de novo profile sum: %.1f expected counts\n", denovo_median))
}

n_ref <- rows[source == "reference", .N]
n_calibrated <- rows[source == "reference" & n_pinned > 0, .N]
n_unwinnable <- rows[source == "reference" & n_called == 0, .N]
cat(sprintf("\n%d of %d reference types received pinned cells; %d were never called at all.\n",
            n_calibrated, n_ref, n_unwinnable))

degenerate <- rows[source == "reference" & !is.na(profile_sum_vs_denovo) &
                     profile_sum_vs_denovo < 0.1]
if (nrow(degenerate)) {
  cat(sprintf("\n%d reference profile%s below a TENTH of the de novo median -- %s cannot win\n",
              nrow(degenerate), if (nrow(degenerate) == 1) " sits" else "s sit",
              if (nrow(degenerate) == 1) "it" else "they"))
  cat("a cell against a de novo cluster whatever the biology:\n")
  print(degenerate[, .(type, n_pinned, n_called,
                       profile_sum = round(profile_sum_fitted, 1),
                       vs_denovo = round(profile_sum_vs_denovo, 3))], nrows = 100)
}

if (!is.null(out_csv)) {
  fwrite(rows, out_csv)
  cat(sprintf("\nWrote %s\n", out_csv))
}

cat("\nREAD: a type with n_pinned = 0 was never calibrated. If those types are also the ones\n",
    "with a collapsed profile_sum, the breadth problem is the pinning step, and the lever is\n",
    "either a looser pinning threshold or a reference whose types are far enough apart to\n",
    "attract distinct cells -- NOT a different K.\n", sep = "")
