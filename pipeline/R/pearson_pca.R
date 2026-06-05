#!/usr/bin/env Rscript
# Stage 3b: batch-corrected quasipoisson Pearson-residual PCA via scPearsonPCA.
#
# Reproduces the Bruker CosMx vignette's Pearson-residual PCA on the cohort, with
# batch correction at the PATIENT level (the `Case` column) per Patrick Danaher's
# guidance that slide effects are typically minor. Follows the patient-batch example
# in https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/pearsonpca/ :
# pass an HVG subset of the counts together with precomputed all-gene totalcounts and
# per-batch gene frequency, so the SVD (on the genes x genes cross-product) never
# densifies and the R dgCMatrix stays small.
#
# All inputs come from stage 3a (concat_qc_anndata.py) in one HDF5 file:
#   --input   <combined_qc>.pca_input.h5 with datasets:
#               /counts/{data,indices,indptr,shape}  HVG counts, CSC genes x cells
#               /genes        HVG gene names         /cell_id  cell ids (col order)
#               /batch        per-cell batch label   /total_counts  per-cell tc
#               /genefreq/matrix  genes_hvg x batches /genefreq/batch  batch labels
# Output:
#   --output  embedding .h5 with /embedding (cells x npcs), /cell_id, /pc,
#             /loadings (genes x npcs), /gene.
#
# --batch-variable names the obs/batch column (default "Case"). Runs single-threaded
# by default; --ncores parallelises the block-wise embedding projection.
#
# Memory: holds the HVG counts (genes_hvg x cells) plus a working copy as dgCMatrix;
# run on a big-memory CPU node, not in a GPU allocation.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(hdf5r)
  library(scPearsonPCA)
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

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--output is required" = !is.null(opt$output))
batch_variable <- opt[["batch-variable"]] %||% "Case"
npcs <- as.integer(opt$npcs %||% 50)
scale_max <- as.numeric(opt[["scale-max"]] %||% 10)
ncores <- as.integer(opt$ncores %||% 1)

# --- read all stage-3a inputs from the one HDF5 file --------------------------
message(sprintf("%s, reading %s", Sys.time(), opt$input))
f <- H5File$new(opt$input, mode = "r")

shape <- f[["counts/shape"]][]          # c(n_genes_hvg, n_cells)
genes <- f[["genes"]][]
cell_id <- f[["cell_id"]][]
# CSC genes x cells slots written by stage 3a (sorted, zero-eliminated).
x <- methods::new(
  "dgCMatrix",
  i = as.integer(f[["counts/indices"]][]),
  p = as.integer(f[["counts/indptr"]][]),
  x = as.double(f[["counts/data"]][]),
  Dim = as.integer(shape),
  Dimnames = list(genes, cell_id)
)

tc <- f[["total_counts"]][]
names(tc) <- cell_id

batch <- f[["batch"]][]
genefreq_mat <- f[["genefreq/matrix"]][, ]   # genes_hvg x batches
genefreq_batches <- f[["genefreq/batch"]][]
f$close_all()

dimnames(genefreq_mat) <- list(genes, genefreq_batches)

obs <- data.table::data.table(cell_ID = cell_id, batch = batch)
data.table::setnames(obs, "batch", batch_variable)
message(sprintf("counts: %d HVGs x %d cells across %d '%s' batches",
                nrow(x), ncol(x), data.table::uniqueN(batch), batch_variable))

# --- align gene frequency columns to the function's batch ordering ------------
# sparse_quasipoisson_pca_seurat_batch builds its one-hot batch matrix with columns
# in first-appearance order of the batch label across colnames(x), and then asserts
# colnames(grate) == colnames(batch_mat) exactly. data.table::unique preserves first
# appearance, so reorder genefreq columns to match (no-op if 3a already ordered them).
batch_order <- obs[match(colnames(x), cell_ID)][, unique(get(batch_variable))]
stopifnot("genefreq batch labels do not cover the cells' batches" =
            all(batch_order %in% colnames(genefreq_mat)))
grate <- Matrix::Matrix(genefreq_mat[, batch_order, drop = FALSE], sparse = TRUE)

# --- batch-corrected Pearson-residual PCA (HVG x cells; tc + grate over all genes) -
res <- scPearsonPCA::sparse_quasipoisson_pca_seurat_batch(
  x = x,
  obs = obs,
  batch_variable = batch_variable,
  cellid_colname = "cell_ID",
  totalcounts = tc,
  grate = grate,
  scale.max = scale_max,
  do.scale = TRUE,
  do.center = TRUE,
  npcs = npcs,
  ncores = ncores,
  return_seurat_reduction = FALSE,
  verbose = TRUE
)

embedding <- res$cell.embeddings       # cells x npcs, col order of x
loadings <- res$feature.loadings       # genes x npcs
rownames(embedding) <- colnames(x)

# --- write the embedding -------------------------------------------------------
message(sprintf("%s, writing embedding to %s", Sys.time(), opt$output))
out <- H5File$new(opt$output, mode = "w")
out[["embedding"]] <- matrix(as.numeric(embedding), nrow = nrow(embedding))
out[["loadings"]] <- matrix(as.numeric(loadings), nrow = nrow(loadings))
out[["cell_id"]] <- colnames(x)
out[["gene"]] <- rownames(x)
out[["pc"]] <- colnames(loadings)
out$close_all()

message(sprintf("Done. embedding: %d cells x %d PCs.",
                nrow(embedding), ncol(embedding)))
