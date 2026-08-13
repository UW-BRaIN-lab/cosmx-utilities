#!/usr/bin/env Rscript
# Shared reader for the stage-4a InSituType input HDF5 (prep_insitutype_anchor.py /
# prep_insitutype_inputs.py). Both the gene-selection and de-novo-K-sweep helpers read
# the same anchor_input.h5, so the reader lives here instead of being copied into each.
#
# The file holds a CSC genes x cells counts matrix plus per-cell background:
#   /counts/{data,indices,indptr,shape}  gene counts, CSC genes x cells
#   /genes    gene names (row order of /counts)
#   /cell_id  cell ids   (col order of /counts)
#   /neg      per-cell mean negprobe count (InSituType background)
#
# read_anchor_h5() returns InSituType's preferred orientation (cells x genes) so callers
# don't each repeat the transpose. Set max_cells to randomly down-sample columns at read
# time (cheap diagnostics on a modest node); seed makes the draw reproducible.

suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
})

read_anchor_h5 <- function(path, max_cells = NULL, seed = 0L) {
  message(sprintf("%s, reading %s", Sys.time(), path))
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  shape <- f[["counts/shape"]][]              # c(n_genes, n_cells)
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

  # Down-sample cells on the genes x cells CSC (cells are COLUMNS -> cheap column subset)
  # BEFORE transposing, so the big anchor never materializes a full cells x genes copy.
  if (!is.null(max_cells) && max_cells < ncol(counts_gxc)) {
    n_before <- ncol(counts_gxc)
    if (seed != 0L) set.seed(seed)
    keep <- sort(sample.int(n_before, max_cells))
    counts_gxc <- counts_gxc[, keep, drop = FALSE]
    neg <- neg[keep]
    message(sprintf("  down-sampled to %d of %d cells (max_cells=%d)",
                    ncol(counts_gxc), n_before, max_cells))
  }

  # cells x genes for InSituType.
  x <- Matrix::t(counts_gxc)
  rm(counts_gxc)

  message(sprintf("counts: %d cells x %d genes; neg median %.4f",
                  nrow(x), ncol(x), stats::median(neg)))
  list(x = x, neg = neg, genes = colnames(x), cell_id = rownames(x))
}
