#!/usr/bin/env Rscript
# Stage 5f (Phase 2): per-cell TOP-K posterior structure against a FIXED profile set.
#
# The PI's Stage-3 primary discriminator is each cell's top-two profile assignments: adjacent
# malignant states => continuum; a malignant + a neuronal (cross-compartment) pairing =>
# candidate hybrid; scattered / cross-donor top matches => still reference-limited (iterate).
# The compact typing h5 only stores the single top label + prob, and the InSituTree result is
# hierarchical (per-branch competition, not a flat posterior across all leaves). So we recover
# a clean flat posterior here by scoring every cell against the fixed InSituTree profile set
# with InSituType::insitutypeML (ML assignment only — NO de-novo clustering, NO reference
# update): the same profiles InSituTree used, so the posteriors are on one consistent scale.
#
# This is deliberately a re-SCORE, not a re-CLUSTER: it answers "given these profiles, what are
# each cell's top-two?", which is what Stage 3 needs. (The Phase-3 "is 48% a floor?" question —
# derive NEW profiles from the untyped pool and reassign — is a separate, heavier run.)
#
# Inputs:
#   --input      insitutype_input.h5 (stage 4a): /counts (CSC genes x cells), /genes, /cell_id,
#                /neg (per-cell mean negprobe). Same file Stage 4b consumed.
#   --profiles   fixed profile matrix CSV (genes x types; first column "gene"). Default use is
#                pipeline/reference/insitutree_profiles.csv (the 60 InSituTree leaf profiles).
#   --top-k      how many top assignments to emit per cell (default 3).
# Output:
#   --output-csv per-cell table: cell_id, top1_type, top1_prob, top2_type, top2_prob, ...
#                (probabilities are softmax posteriors over the profile set).
#
# Memory: holds cells x genes counts + a cells x types log-likelihood matrix; run on a
# big-memory CPU node (ckpt), like Stage 4b.

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
    a <- args[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument: %s", a))
    key <- sub("^--", "", a)
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]; out[[kv[[1]]]] <- kv[[2]]; i <- i + 1
    } else { out[[key]] <- args[[i + 1]]; i <- i + 2 }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--profiles is required" = !is.null(opt$profiles))
stopifnot("--output-csv is required" = !is.null(opt[["output-csv"]]))
top_k <- as.integer(opt[["top-k"]] %||% 3)

# --- read the stage-4a input (mirrors insitutype_typing.R) -------------------
message(sprintf("%s, reading %s", Sys.time(), opt$input))
f <- H5File$new(opt$input, mode = "r")
shape <- f[["counts/shape"]][]            # c(n_genes, n_cells)
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
message(sprintf("counts: %d cells x %d genes; neg median %.4f", nrow(x), ncol(x), median(neg)))

# --- read the fixed profile set (genes x types) ------------------------------
message(sprintf("%s, reading profiles %s", Sys.time(), opt$profiles))
ref_dt <- data.table::fread(opt$profiles)
ref_genes <- as.character(ref_dt[[1]])
ref <- as.matrix(ref_dt[, -1L]); rownames(ref) <- ref_genes
storage.mode(ref) <- "double"
message(sprintf("profiles: %d genes x %d types", nrow(ref), ncol(ref)))

shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and profiles" = length(shared) > 0)
message(sprintf("shared genes: %d / %d panel, %d / %d profiles",
                length(shared), ncol(x), length(shared), nrow(ref)))
x <- x[, shared, drop = FALSE]
ref <- ref[shared, , drop = FALSE]

# --- flat ML scoring against the fixed profiles ------------------------------
# estimateBackground is not exported in InSituType 2.0; precompute + pass explicitly (same
# runtime shim insitutree_typing.R uses) so the bg = NULL path is never taken.
bg <- tryCatch(InSituType:::estimateBackground(counts = x, neg = neg[rownames(x)]),
               error = function(e) { message("estimateBackground shim failed: ", conditionMessage(e)); NULL })
message(sprintf("%s, insitutypeML: %d cells vs %d fixed profiles", Sys.time(), nrow(x), ncol(ref)))
res <- InSituType::insitutypeML(
  x = x,
  neg = neg[rownames(x)],
  bg = bg,
  reference_profiles = ref
)

ll <- res$logliks
if (is.null(ll)) stop("insitutypeML returned no $logliks; cannot rank top-K posteriors.")
ll <- as.matrix(ll)
if (is.null(colnames(ll))) colnames(ll) <- colnames(ref)
n <- nrow(ll); k <- min(top_k, ncol(ll))
message(sprintf("logliks: %d cells x %d types; extracting top %d", n, ncol(ll), k))

# softmax posteriors over the profile set, using the row max for numerical stability
rmax <- ll[cbind(seq_len(n), max.col(ll, ties.method = "first"))]
lse <- rmax + log(rowSums(exp(sweep(ll, 1, rmax, "-"))))
types <- colnames(ll)

# top-K via repeated max.col with masking (vectorised; avoids a per-row sort over millions)
out <- data.table(cell_id = rownames(ll) %||% cell_id)
work <- ll
for (r in seq_len(k)) {
  idx <- max.col(work, ties.method = "first")
  llr <- work[cbind(seq_len(n), idx)]
  out[[sprintf("top%d_type", r)]] <- types[idx]
  out[[sprintf("top%d_prob", r)]] <- exp(llr - lse)
  work[cbind(seq_len(n), idx)] <- -Inf
}

message(sprintf("%s, writing %s", Sys.time(), opt[["output-csv"]]))
data.table::fwrite(out, opt[["output-csv"]])
message(sprintf("Done. Top-%d posteriors for %d cells against %d profiles.", k, n, ncol(ll)))
