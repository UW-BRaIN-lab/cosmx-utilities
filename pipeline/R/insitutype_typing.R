#!/usr/bin/env Rscript
# Stage 4b: semi-supervised cell typing with InSituType against the CZI GBmap reference.
#
# Cells are matched to the Extended GBmap level-3 reference profiles (named types:
# AC-like, MES-like, TAM-BDM/TAM-MG, Oligodendrocyte, ...) while InSituType ALSO discovers de
# novo clusters (single-letter labels a, b, c, ...) for tumour populations that fit the
# reference poorly. Follows the Bruker CosMx Scratch Space InSituType workflow.
#
# Reference-profile update is gentle on purpose: rescale = TRUE (gene-level platform
# correction for the scRNA-seq -> CosMx shift), refit = FALSE, refinement = FALSE. This
# keeps the GBmap profiles anchored to their biology so off-reference tumour cells stay
# off-reference and surface as the de novo letter clusters, rather than being absorbed
# by an aggressive refit. n_clusts is passed as a RANGE (default 10:20): InSituType then
# auto-selects the optimal number of de novo clusters within it.
#
# Inputs (from stage 4a, prep_insitutype_inputs.py) in one HDF5 file:
#   --input   insitutype_input.h5 with datasets:
#               /counts/{data,indices,indptr,shape}  gene counts, CSC genes x cells
#               /genes    gene names      /cell_id  cell ids (col order of /counts)
#               /neg      per-cell mean negprobe count (InSituType background)
#   --reference  GBmap panel-restricted reference CSV (genes x cell types; first
#                column "gene" is the row index). See pipeline/reference/.
# Outputs:
#   --output-rds  full insitutype() result (clust, prob, profiles, logliks) as .rds —
#                 the canonical R object; load it for a later interactive
#                 InSituType::refineClusters pass to condense over-split de novo clusters.
#   --output-h5   compact result for the Python writeback (stage 4c):
#                 /cell_id, /cell_type (= clust), /prob, /profiles (genes x types),
#                 /profile_types (column names of /profiles).
#
# Memory: holds the cells x genes counts plus a cells x types log-likelihood matrix at
# cohort scale; run on a big-memory CPU node (ckpt), not in a GPU allocation.

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

as_logical <- function(s) tolower(s) %in% c("true", "t", "1", "yes")

# "10:20" -> 10:20; "12" -> 12. The range form triggers InSituType's automatic
# choose-cluster-number step within the range.
parse_n_clusts <- function(s) {
  if (grepl(":", s, fixed = TRUE)) {
    ab <- as.integer(strsplit(s, ":", fixed = TRUE)[[1]])
    stopifnot("--n-clusts range must be 'lo:hi'" = length(ab) == 2L)
    return(seq.int(ab[[1]], ab[[2]]))
  }
  as.integer(s)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--reference is required" = !is.null(opt$reference))
stopifnot("--output-rds is required" = !is.null(opt[["output-rds"]]))
stopifnot("--output-h5 is required" = !is.null(opt[["output-h5"]]))

n_clusts <- parse_n_clusts(opt[["n-clusts"]] %||% "10:20")
update_reference <- as_logical(opt[["update-reference"]] %||% "true")
rescale <- as_logical(opt$rescale %||% "true")
refit <- as_logical(opt$refit %||% "false")
refinement <- as_logical(opt$refinement %||% "false")
n_starts <- as.integer(opt[["n-starts"]] %||% 10)
max_iters <- as.integer(opt[["max-iters"]] %||% 40)
seed <- as.integer(opt$seed %||% 0)

# --- read the stage-4a input -------------------------------------------------
message(sprintf("%s, reading %s", Sys.time(), opt$input))
f <- H5File$new(opt$input, mode = "r")
shape <- f[["counts/shape"]][]           # c(n_genes, n_cells)
genes <- f[["genes"]][]
cell_id <- f[["cell_id"]][]
# CSC genes x cells slots written by stage 4a (sorted, zero-eliminated).
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

# InSituType wants cells x genes.
x <- Matrix::t(counts_gxc)
rm(counts_gxc)
message(sprintf("counts: %d cells x %d genes; neg median %.4f",
                nrow(x), ncol(x), median(neg)))

# --- read the GBmap reference profiles (genes x cell types) -------------------
message(sprintf("%s, reading reference %s", Sys.time(), opt$reference))
ref_dt <- data.table::fread(opt$reference)
ref_genes <- as.character(ref_dt[[1]])           # first column is the gene index
ref <- as.matrix(ref_dt[, -1L])
rownames(ref) <- ref_genes
storage.mode(ref) <- "double"
message(sprintf("reference: %d genes x %d cell types", nrow(ref), ncol(ref)))

# --- align genes (InSituType uses only the shared genes) ----------------------
shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and reference" = length(shared) > 0)
message(sprintf("shared genes: %d / %d panel, %d / %d reference",
                length(shared), ncol(x), length(shared), nrow(ref)))
x <- x[, shared, drop = FALSE]
ref <- ref[shared, , drop = FALSE]

# --- semi-supervised typing ---------------------------------------------------
if (seed != 0) set.seed(seed)
message(sprintf("%s, insitutype: n_clusts=%s update_reference=%s rescale=%s refit=%s refinement=%s",
                Sys.time(),
                if (length(n_clusts) > 1) sprintf("%d:%d", min(n_clusts), max(n_clusts))
                  else as.character(n_clusts),
                update_reference, rescale, refit, refinement))
res <- InSituType::insitutype(
  x = x,
  neg = neg[rownames(x)],
  n_clusts = n_clusts,
  reference_profiles = ref,
  update_reference_profiles = update_reference,
  rescale = rescale,
  refit = refit,
  refinement = refinement,
  n_starts = n_starts,
  max_iters = max_iters
)

message("cell-type assignment counts:")
print(sort(table(res$clust), decreasing = TRUE))

# --- write outputs ------------------------------------------------------------
message(sprintf("%s, writing full result to %s", Sys.time(), opt[["output-rds"]]))
saveRDS(res, opt[["output-rds"]])

# clust / prob are named by cell id; write them in input cell order.
clust <- res$clust[cell_id]
prob <- res$prob[cell_id]

message(sprintf("%s, writing compact result to %s", Sys.time(), opt[["output-h5"]]))
out <- H5File$new(opt[["output-h5"]], mode = "w")
out[["cell_id"]] <- cell_id
out[["cell_type"]] <- as.character(clust)
out[["prob"]] <- as.numeric(prob)
out[["profiles"]] <- matrix(as.numeric(res$profiles), nrow = nrow(res$profiles))
out[["profile_genes"]] <- rownames(res$profiles)
out[["profile_types"]] <- colnames(res$profiles)
out$close_all()

message(sprintf("Done. %d cells typed into %d cell types (%d named + de novo).",
                length(clust), length(unique(clust)),
                sum(colnames(ref) %in% unique(clust))))
