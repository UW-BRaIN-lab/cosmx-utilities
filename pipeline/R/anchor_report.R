#!/usr/bin/env Rscript
# Stage 4b': anchor report — which reference cell types FAIL TO ANCHOR, and why.
#
# When insitutype() runs with update_reference_profiles = TRUE (our rescale runs), it does
# not rescale the reference wholesale: it first picks ANCHOR CELLS per reference type, then
# rescales each profile from its own anchors. A type that cannot find anchors is dropped
# from the anchor set entirely, so its profile is never adapted to the CosMx platform —
# after which it rarely wins any cell, and those cells fall into the de novo sink instead.
# That is the mechanism behind the retina study's "42 collinear types never anchor" finding
# and, we suspect, a large part of the GBM cohort's Low_signal fraction.
#
# InSituType's anchor filter is a sequence (find_anchor_cells.R upstream):
#   1. cosine to the profile must be >= min_anchor_cosine
#   2. scaled log-likelihood ratio vs the runner-up must be >= min_anchor_llr
#   3. only the top n_anchor_cells candidates per type are kept
#   4. survivors are re-scored against their own centroid; those below min_cosine are cut
#   5. ANY type left with <= insufficient_anchors_thresh anchors loses ALL of them
# Step 5 is the cliff: a type is either anchored or completely unanchored, and nothing in
# the insitutype() output records which happened.
#
# This script reproduces exactly that anchor selection — same functions, and defaults that
# mirror insitutype()'s own (n_anchor_cells = 2000, min_anchor_cosine = 0.3,
# min_anchor_llr = 0.03, insufficient_anchors_thresh = 20, nb_size = 10), NOT the looser
# defaults on find_anchor_cells() itself — and reports per type whether it anchored. It
# does NOT run the EM, so it is far cheaper than a typing run and can be pointed at an
# existing stage-4a input to diagnose a run that has already happened.
#
# Inputs:
#   --input      stage-4a insitutype_input.h5 (prep_insitutype_inputs.py):
#                  /counts/{data,indices,indptr,shape}  gene counts, CSC genes x cells
#                  /genes    gene names    /cell_id  cell ids (col order of /counts)
#                  /neg      per-cell mean negprobe count
#   --reference  panel-restricted reference CSV (genes x cell types; first column "gene").
#   --keep-genes optional kept_genes.txt (73_select_genes.sh) — pass the SAME file the run
#                being diagnosed used, or the anchor picture will not match it.
# Outputs:
#   --output-csv      one row per reference cell type (see build_report()).
#   --output-summary  human-readable summary; also echoed to stdout.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(hdf5r)
  library(InSituType)
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

# Cosine between every pair of reference profiles: the collinearity that makes a type
# unanchorable in the first place.
profile_cosines <- function(ref) {
  norms <- sqrt(colSums(ref^2))
  norms[norms == 0] <- NA_real_
  scaled <- sweep(ref, 2L, norms, "/")
  crossprod(scaled)
}

build_report <- function(ref, anchors, cos_mat, llr_mat, min_cosine,
                         insufficient_anchors_thresh) {
  types <- colnames(ref)
  stopifnot(
    "get_anchor_stats returned columns that do not match the reference types" =
      identical(colnames(cos_mat), types) && identical(colnames(llr_mat), types)
  )
  # Which type each cell looks MOST like on cosine. A type that is never any cell's best
  # match is unanchorable for a different reason than one that is often the best match but
  # always loses the likelihood-ratio tie-break. Cells whose BEST cosine is still under the
  # cut resemble nothing in the reference, so they name no winner: leaving them in would
  # credit them to whichever type max.col's tie-break happens to reach first.
  best_idx <- max.col(cos_mat, ties.method = "first")
  best_cosine <- cos_mat[cbind(seq_len(nrow(cos_mat)), best_idx)]
  best_type <- types[best_idx]
  best_type[best_cosine < min_cosine] <- NA_character_
  anchor_counts <- table(factor(anchors[!is.na(anchors)], levels = types))

  ref_cos <- profile_cosines(ref)
  diag(ref_cos) <- NA_real_

  rows <- lapply(types, function(type) {
    cos_col <- cos_mat[, type]
    above <- cos_col >= min_cosine
    n_above <- sum(above, na.rm = TRUE)
    # Among the cells this type could plausibly claim, which type actually claimed them.
    absorber <- NA_character_
    if (n_above > 0) {
      competing <- table(best_type[above])
      competing <- competing[names(competing) != type]
      if (length(competing) > 0) absorber <- names(which.max(competing))
    }
    nearest <- if (all(is.na(ref_cos[, type]))) NA_integer_ else which.max(ref_cos[, type])
    data.table(
      cell_type = type,
      n_anchors = as.integer(anchor_counts[[type]]),
      anchored = as.integer(anchor_counts[[type]]) > 0L,
      n_cells_above_min_cosine = n_above,
      n_cells_best_on_cosine = sum(best_type == type, na.rm = TRUE),
      max_cosine = suppressWarnings(max(cos_col, na.rm = TRUE)),
      median_llr_above_cosine = if (n_above > 0)
        stats::median(llr_mat[above, type], na.rm = TRUE) else NA_real_,
      empirical_absorber = absorber,
      nearest_reference_type = if (is.na(nearest)) NA_character_ else types[[nearest]],
      nearest_reference_cosine = if (is.na(nearest)) NA_real_ else ref_cos[nearest, type]
    )
  })
  report <- rbindlist(rows)
  report[, insufficient_anchors_thresh := insufficient_anchors_thresh]
  report[order(anchored, -n_anchors)]
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--reference is required" = !is.null(opt$reference))
stopifnot("--output-csv is required" = !is.null(opt[["output-csv"]]))

# Defaults mirror insitutype()'s own anchor arguments, NOT find_anchor_cells()'s looser
# ones (which use n_cells = 500 and min_scaled_llr = 0.01). Matching insitutype() is the
# whole point: this must reproduce what the typing run did.
n_anchor_cells <- as.integer(opt[["n-anchor-cells"]] %||% 2000)
min_anchor_cosine <- as.numeric(opt[["min-anchor-cosine"]] %||% 0.3)
min_anchor_llr <- as.numeric(opt[["min-anchor-llr"]] %||% 0.03)
insufficient_anchors_thresh <- as.integer(opt[["insufficient-anchors-thresh"]] %||% 20)
nb_size <- as.numeric(opt[["nb-size"]] %||% 10)
subsample <- as.integer(opt$subsample %||% 0)
seed <- as.integer(opt$seed %||% 0)
keep_genes_path <- opt[["keep-genes"]]

# --- read the stage-4a input (same layout as insitutype_typing.R) -------------
message(sprintf("%s, reading %s", Sys.time(), opt$input))
f <- H5File$new(opt$input, mode = "r")
shape <- f[["counts/shape"]][]           # c(n_genes, n_cells)
genes <- f[["genes"]][]
cell_id <- f[["cell_id"]][]
counts_gxc <- methods::new(
  "dgCMatrix",
  i = as.integer(f[["counts/indices"]][]),
  p = as.integer(f[["counts/indptr"]][]),
  x = as.double(f[["counts/data"]][]),
  Dim = as.integer(shape),
  Dimnames = list(genes, cell_id)
)
neg <- f[["neg"]][]
names(neg) <- cell_id
f$close_all()

x <- Matrix::t(counts_gxc)               # InSituType wants cells x genes
rm(counts_gxc)
message(sprintf("counts: %d cells x %d genes; neg median %.4f",
                nrow(x), ncol(x), median(neg)))

# Anchor selection is per-type and threshold-based, so a random subsample gives a fast,
# unbiased preview of the anchor picture before committing a full-cohort job.
if (subsample > 0 && subsample < nrow(x)) {
  if (seed != 0) set.seed(seed)
  take <- sort(sample.int(nrow(x), subsample))
  x <- x[take, , drop = FALSE]
  neg <- neg[take]
  message(sprintf("subsampled to %d cells for a preview run", nrow(x)))
}

# --- reference profiles (genes x cell types) ---------------------------------
message(sprintf("%s, reading reference %s", Sys.time(), opt$reference))
ref_dt <- data.table::fread(opt$reference)
ref <- as.matrix(ref_dt[, -1L])
rownames(ref) <- as.character(ref_dt[[1]])
storage.mode(ref) <- "double"
message(sprintf("reference: %d genes x %d cell types", nrow(ref), ncol(ref)))

if (!is.null(keep_genes_path)) {
  keep_genes <- intersect(readLines(keep_genes_path), colnames(x))
  stopifnot("no --keep-genes overlap the panel" = length(keep_genes) > 0)
  message(sprintf("gene pruning: restricting panel %d -> %d genes (--keep-genes %s)",
                  ncol(x), length(keep_genes), keep_genes_path))
  x <- x[, keep_genes, drop = FALSE]
}

shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and reference" = length(shared) > 0)
message(sprintf("shared genes: %d / %d panel, %d / %d reference",
                length(shared), ncol(x), length(shared), nrow(ref)))
x <- x[, shared, drop = FALSE]
ref <- ref[shared, , drop = FALSE]

# --- background, derived exactly as insitutype() derives it from neg ----------
# (utilities.R: s <- rowMeans(counts); bg <- lm(neg ~ s - 1)$fitted). Computing it once
# here lets get_anchor_stats and choose_anchors_from_stats share one pass.
s <- Matrix::rowMeans(x)
bg <- stats::lm(neg ~ s - 1)$fitted
message(sprintf("background: median %.4f (from an intercept-free neg ~ rowMeans fit)",
                median(bg)))

# --- anchor selection --------------------------------------------------------
message(sprintf("%s, get_anchor_stats (size=%g, min_cosine=%g)",
                Sys.time(), nb_size, min_anchor_cosine))
anchorstats <- InSituType::get_anchor_stats(
  counts = x, bg = bg, profiles = ref, size = nb_size, min_cosine = min_anchor_cosine
)

message(sprintf("%s, choose_anchors_from_stats (n_cells=%d, min_scaled_llr=%g, thresh=%d)",
                Sys.time(), n_anchor_cells, min_anchor_llr, insufficient_anchors_thresh))
anchors <- InSituType::choose_anchors_from_stats(
  counts = x, bg = bg, anchorstats = anchorstats,
  n_cells = n_anchor_cells, min_cosine = min_anchor_cosine,
  min_scaled_llr = min_anchor_llr,
  insufficient_anchors_thresh = insufficient_anchors_thresh
)

report <- build_report(ref, anchors, anchorstats$cos, anchorstats$llr,
                       min_anchor_cosine, insufficient_anchors_thresh)
data.table::fwrite(report, opt[["output-csv"]])
message(sprintf("wrote %s", opt[["output-csv"]]))

# --- summary -----------------------------------------------------------------
failed <- report[anchored == FALSE]
lines <- c(
  sprintf("Anchor report: %d reference cell types, %d cells, %d shared genes",
          ncol(ref), nrow(x), length(shared)),
  sprintf("Settings: n_anchor_cells=%d min_anchor_cosine=%g min_anchor_llr=%g thresh=%d",
          n_anchor_cells, min_anchor_cosine, min_anchor_llr, insufficient_anchors_thresh),
  "",
  sprintf("ANCHORED:   %d / %d types (%.1f%%)",
          nrow(report) - nrow(failed), nrow(report),
          100 * (nrow(report) - nrow(failed)) / nrow(report)),
  sprintf("UNANCHORED: %d / %d types (%.1f%%) — profiles never rescaled",
          nrow(failed), nrow(report), 100 * nrow(failed) / nrow(report)),
  sprintf("Total anchor cells: %d of %d (%.2f%%)",
          sum(report$n_anchors), nrow(x), 100 * sum(report$n_anchors) / nrow(x))
)
if (nrow(failed) > 0) {
  lines <- c(lines, "", "Unanchored types (max_cosine, cells above min_cosine, absorber,",
             "                  nearest reference type / its profile cosine):")
  for (i in seq_len(nrow(failed))) {
    r <- failed[i]
    lines <- c(lines, sprintf(
      "  %-38s cos=%.3f  n>=cut=%-8d absorbed_by=%-24s nearest=%s (%.3f)",
      r$cell_type, r$max_cosine, r$n_cells_above_min_cosine,
      r$empirical_absorber, r$nearest_reference_type,
      r$nearest_reference_cosine))
  }
}
summary_text <- paste(lines, collapse = "\n")
cat(summary_text, "\n", sep = "")
if (!is.null(opt[["output-summary"]])) {
  writeLines(lines, opt[["output-summary"]])
  message(sprintf("wrote %s", opt[["output-summary"]]))
}
