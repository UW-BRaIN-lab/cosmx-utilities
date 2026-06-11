#!/usr/bin/env Rscript
# Two-pass cell typing, PASS 2 (per slide): assign one slide's cells to the FIXED
# cohort profile in SUPERVISED mode (InSituType::insitutypeML) — no de novo discovery,
# no reference update. Every slide is typed against the same cohort_profile.csv, so the
# named GBmap types and the de novo a,b,c... labels mean the same thing cohort-wide.
#
# Run as a Slurm array, one task per slide (75_insitutype_supervised.sh).
#
# Inputs:
#   --input      <slide>.h5 (from split_insitutype_inputs.py): /counts (CSC genes x
#                cells), /genes, /cell_id, /neg  — same format as stage 4a.
#   --profile    cohort_profile.csv (from insitutype_profile.R): genes x cell types,
#                first column "gene".
# Output:
#   --output-h5  /cell_id, /cell_type (assignment), /prob (per-cell posterior). Same
#                shape insitutype_typing.R writes, so write_celltypes.py reads it.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(hdf5r)
  library(InSituType)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    out[[key]] <- args[[i + 1]]; i <- i + 2
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input),
          "--profile is required" = !is.null(opt$profile),
          "--output-h5 is required" = !is.null(opt[["output-h5"]]))

# --- read the per-slide input (CSC genes x cells) ----------------------------
message(sprintf("%s, reading %s", Sys.time(), opt$input))
f <- H5File$new(opt$input, mode = "r")
shape <- f[["counts/shape"]][]
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
neg <- f[["neg"]][]; names(neg) <- cell_id
f$close_all()
x <- Matrix::t(counts_gxc); rm(counts_gxc)
message(sprintf("counts: %d cells x %d genes; neg median %.4f",
                nrow(x), ncol(x), median(neg)))

# --- read the fixed cohort profile (genes x cell types) ----------------------
message(sprintf("%s, reading profile %s", Sys.time(), opt$profile))
prof_dt <- data.table::fread(opt$profile)
prof_genes <- as.character(prof_dt[[1]])
profiles <- as.matrix(prof_dt[, -1L]); rownames(profiles) <- prof_genes
storage.mode(profiles) <- "double"
message(sprintf("profile: %d genes x %d cell types", nrow(profiles), ncol(profiles)))

# --- align genes and assign (supervised; profile is FIXED, no rescale) -------
shared <- intersect(colnames(x), rownames(profiles))
stopifnot("no genes shared between slide and profile" = length(shared) > 0)
x <- x[, shared, drop = FALSE]
profiles <- profiles[shared, , drop = FALSE]
message(sprintf("%s, insitutypeML over %d cells vs %d fixed types",
                Sys.time(), nrow(x), ncol(profiles)))

res <- InSituType::insitutypeML(
  x = x,
  neg = neg[rownames(x)],
  reference_profiles = profiles
)

clust <- res$clust[cell_id]
prob <- res$prob[cell_id]
message("cell-type assignment counts:")
print(sort(table(clust), decreasing = TRUE))
message(sprintf("median posterior prob %.3f; %d cells below 0.5",
                median(prob), sum(prob < 0.5)))

message(sprintf("%s, writing %s", Sys.time(), opt[["output-h5"]]))
out <- H5File$new(opt[["output-h5"]], mode = "w")
out[["cell_id"]] <- cell_id
out[["cell_type"]] <- as.character(clust)
out[["prob"]] <- as.numeric(prob)
out$close_all()
message(sprintf("Done. %d cells typed into %d cell types.",
                length(clust), length(unique(clust))))
