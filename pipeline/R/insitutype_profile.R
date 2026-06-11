#!/usr/bin/env Rscript
# Two-pass cell typing, PASS 1 (profile step): turn the anchor's semi-supervised
# InSituType assignment into a single FIXED cohort profile that PASS 2 assigns every
# slide against.
#
# The anchor was typed semi-supervised by insitutype_typing.R (named GBmap types + de
# novo a,b,c... clusters). Here we rebuild clean mean expression profiles from the
# anchor's RAW counts + negprobe background + those labels, via InSituType::getRNAprofiles
# — the same operation Nanostring's CosMx-Cell-Profiles uses to publish reference
# profiles. The result is in CosMx count space, so PASS 2's supervised insitutypeML can
# match each slide's raw counts directly (no further rescale).
#
# Inputs:
#   --anchor-input  anchor_input.h5 (from prep_insitutype_anchor.py): /counts (CSC
#                   genes x cells), /genes, /cell_id, /neg
#   --anchor-typing anchor typing .h5 (from insitutype_typing.R): /cell_id, /cell_type
# Output:
#   --output-csv    cohort_profile.csv: genes x cell types, first column "gene". Same
#                   shape as the GBmap reference CSV, so insitutype_supervised.R reads
#                   it the same way insitutype_typing.R reads --reference.

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
stopifnot("--anchor-input is required" = !is.null(opt[["anchor-input"]]),
          "--anchor-typing is required" = !is.null(opt[["anchor-typing"]]),
          "--output-csv is required" = !is.null(opt[["output-csv"]]))

# --- anchor counts (CSC genes x cells, same layout as stage 4a) --------------
message(sprintf("%s, reading anchor counts %s", Sys.time(), opt[["anchor-input"]]))
f <- H5File$new(opt[["anchor-input"]], mode = "r")
shape <- f[["counts/shape"]][]                  # c(n_genes, n_cells)
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
x <- Matrix::t(counts_gxc)                       # InSituType wants cells x genes
rm(counts_gxc)

# --- anchor labels from PASS-1 semi-supervised typing ------------------------
message(sprintf("%s, reading anchor typing %s", Sys.time(), opt[["anchor-typing"]]))
g <- H5File$new(opt[["anchor-typing"]], mode = "r")
typed_cell_id <- g[["cell_id"]][]
typed_cell_type <- g[["cell_type"]][]
g$close_all()
clust <- setNames(as.character(typed_cell_type), typed_cell_id)[cell_id]
stopifnot("anchor typing does not cover all anchor cells" = !anyNA(clust))
message(sprintf("anchor: %d cells x %d genes; %d cell types",
                nrow(x), ncol(x), length(unique(clust))))

# --- build fixed profile: genes x cell types --------------------------------
message(sprintf("%s, getRNAprofiles over %d types", Sys.time(), length(unique(clust))))
profiles <- InSituType::getRNAprofiles(x = x, neg = neg[rownames(x)], clust = clust)
message(sprintf("profile: %d genes x %d cell types", nrow(profiles), ncol(profiles)))

out <- data.table(gene = rownames(profiles))
out <- cbind(out, as.data.table(profiles))
data.table::fwrite(out, opt[["output-csv"]])
message(sprintf("Wrote %s (%d genes x %d types, incl. de novo).",
                opt[["output-csv"]], nrow(profiles), ncol(profiles)))
