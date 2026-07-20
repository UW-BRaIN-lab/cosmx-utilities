#!/usr/bin/env Rscript
# Stage 4 (alt): SUPERVISED, hierarchical cell typing with InSituTree.
#
# Unlike InSituType (R/insitutype_typing.R), InSituTree does NOT discover de novo
# clusters and does NOT platform-rescale the reference — it forces every cell to a leaf
# of a user-defined cell-type hierarchy, deciding one InSituType::insitutypeML pass per
# internal node so a fine leaf only competes against its siblings within its branch.
# That branch-local competition dissolves the rare-type "attractor" problem that made
# Core GBmap level-4's 54 fine types unusable under a flat refit, so we can use the full
# level-4 type set here. Cells that fit poorly surface as LOW posterior probability
# rather than a de novo cluster (flag low `prob` as unresolved downstream).
#
# Because InSituTree does no scRNA->CosMx rescaling, `full_profiles` must already be in
# CosMx scale AND internally consistent: we build it (python/prep_insitutree_profiles.py)
# from ONE InSituType rescale run's $profiles — Core-L4 named types + our validated
# de-novo tumor states, all in one scaling. See pipeline/reference/insitutree_profiles.csv
# and pipeline/reference/insitutree_hierarchy.json (leaf names must match the profile
# columns).
#
# Inputs (the SAME stage-4a input InSituType uses — no new prep):
#   --input      insitutype_input.h5 with datasets:
#                  /counts/{data,indices,indptr,shape}  gene counts, CSC genes x cells
#                  /genes    gene names      /cell_id  cell ids (col order of /counts)
#                  /neg      per-cell mean negprobe count (background)
#   --reference  InSituTree profile matrix CSV (genes x cell types; first column "gene"
#                is the row index). Every hierarchy leaf must be a column here.
#   --hierarchy  cell-type hierarchy JSON (nested object; leaf groups are string arrays).
# Outputs:
#   --output-rds  full runInSituTree() result (nested per-node list + summaryAnnotation).
#   --output-h5   compact result for the Python writeback (reuses write_celltypes.py):
#                  /cell_id, /cell_type (finest resolved level), /prob.
#   --output-csv  the full multi-level summaryAnnotation (annotLevel_k / probs_annotLevel_k)
#                 for later re-summary at a coarser level without re-running.
#
# Memory: holds cells x genes counts plus one insitutypeML pass per node and a
# cells x (2*levels) summaryAnnotation; run on a big-memory CPU node (ckpt).

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(hdf5r)
  library(jsonlite)
  library(InSituType)
  library(InSituTree)
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

# "GENE1,GENE2" -> c("GENE1","GENE2"); "" / NULL -> character(0).
parse_csv_list <- function(s) {
  if (is.null(s) || !nzchar(s)) return(character(0))
  trimws(strsplit(s, ",", fixed = TRUE)[[1]])
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--reference is required" = !is.null(opt$reference))
stopifnot("--hierarchy is required" = !is.null(opt$hierarchy))
stopifnot("--output-rds is required" = !is.null(opt[["output-rds"]]))
stopifnot("--output-h5 is required" = !is.null(opt[["output-h5"]]))
stopifnot("--output-csv is required" = !is.null(opt[["output-csv"]]))

# InSituTree feature-selection quantiles (0-1); source defaults are 0.5/0.5.
q_abs <- as.numeric(opt[["quantile-absolute"]] %||% "0.5")
q_pct <- as.numeric(opt[["quantile-percent"]] %||% "0.5")
excluded_genes <- parse_csv_list(opt[["excluded-genes"]])
seed <- as.integer(opt$seed %||% 0)

# --- read the stage-4a input (identical layout to InSituType) ----------------
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

# InSituType/InSituTree want cells x genes.
x <- Matrix::t(counts_gxc)
rm(counts_gxc)
message(sprintf("counts: %d cells x %d genes; neg median %.4f",
                nrow(x), ncol(x), median(neg)))

# --- read the reference profile matrix (genes x cell types) -------------------
message(sprintf("%s, reading reference %s", Sys.time(), opt$reference))
ref_dt <- data.table::fread(opt$reference)
ref_genes <- as.character(ref_dt[[1]])           # first column is the gene index
ref <- as.matrix(ref_dt[, -1L])
rownames(ref) <- ref_genes
storage.mode(ref) <- "double"
message(sprintf("reference: %d genes x %d cell types", nrow(ref), ncol(ref)))

# --- align genes (use only the shared genes) ----------------------------------
shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and reference" = length(shared) > 0)
message(sprintf("shared genes: %d / %d panel, %d / %d reference",
                length(shared), ncol(x), length(shared), nrow(ref)))
x <- x[, shared, drop = FALSE]
ref <- ref[shared, , drop = FALSE]

# --- parse and validate the hierarchy ----------------------------------------
# simplifyVector keeps leaf arrays as character vectors; simplifyDataFrame/Matrix FALSE
# keeps internal nodes as nested named lists (the cth shape collapseProfiles expects).
cth <- jsonlite::fromJSON(opt$hierarchy, simplifyVector = TRUE,
                          simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
leaves <- unlist(cth, use.names = FALSE)
missing_leaves <- setdiff(leaves, colnames(ref))
if (length(missing_leaves) > 0) {
  stop(sprintf("hierarchy leaves absent from reference profile columns: %s",
               paste(missing_leaves, collapse = ", ")))
}
dup_leaves <- leaves[duplicated(leaves)]
if (length(dup_leaves) > 0) {
  stop(sprintf("hierarchy leaves listed more than once: %s",
               paste(unique(dup_leaves), collapse = ", ")))
}
message(sprintf("hierarchy: %d leaves over the reference's %d types",
                length(leaves), ncol(ref)))
message("hierarchy structure:")
InSituTree::printTree(cth)

# --- background ---------------------------------------------------------------
# runInSituTree computes background via InSituType::estimateBackground, but that name is
# not exported in the pinned InSituType 2.0 (insitutypeML reaches it internally), so the
# `bg = NULL` path errors. Precompute it with ::: (identical to what insitutypeML uses)
# and pass it explicitly. Drop this shim once a pinned InSituType exports the function.
message(sprintf("%s, estimating background", Sys.time()))
bg <- InSituType:::estimateBackground(counts = x, neg = neg[rownames(x)])

# --- hierarchical supervised typing -------------------------------------------
if (seed != 0) set.seed(seed)
message(sprintf("%s, runInSituTree: quantile_abs=%.2f quantile_pct=%.2f excluded=%d",
                Sys.time(), q_abs, q_pct, length(excluded_genes)))
res <- InSituTree::runInSituTree(
  x = x,
  neg = neg[rownames(x)],
  bg = bg,
  full_profiles = ref,
  cth = cth,
  excluded_genes = excluded_genes,
  quantile_absolute_expression_difference_param = q_abs,
  quantile_percent_expression_difference_param = q_pct,
  return_summary_annotation = TRUE
)

sa <- res$summaryAnnotation
stopifnot("runInSituTree returned no summaryAnnotation" = !is.null(sa))

# The finest resolved label is the highest-numbered annotLevel_* column (summarizeInSituTree
# back-fills NA deeper levels from the parent, so the last column is fully populated for
# every retained cell). Take that column and its paired probability.
level_cols <- grep("^annotLevel_[0-9]+$", colnames(sa), value = TRUE)
stopifnot("no annotLevel_* columns in summaryAnnotation" = length(level_cols) > 0)
last_level <- max(as.integer(sub("^annotLevel_", "", level_cols)))
label_col <- sprintf("annotLevel_%d", last_level)
prob_col <- sprintf("probs_annotLevel_%d", last_level)
message(sprintf("finest resolved level: %d (%s / %s)", last_level, label_col, prob_col))

# Re-index to input cell order (summaryAnnotation rownames are cell ids).
label_by_id <- setNames(as.character(sa[[label_col]]), rownames(sa))
prob_by_id <- setNames(as.numeric(sa[[prob_col]]), rownames(sa))
cell_type <- unname(label_by_id[cell_id])
prob <- unname(prob_by_id[cell_id])

n_na <- sum(is.na(cell_type))
message(sprintf("typed %d cells; %d unresolved at level 1 (NA); %d distinct labels",
                length(cell_type), n_na, length(unique(na.omit(cell_type)))))
message("cell-type assignment counts:")
print(sort(table(cell_type, useNA = "ifany"), decreasing = TRUE))

# --- write outputs ------------------------------------------------------------
message(sprintf("%s, writing full result to %s", Sys.time(), opt[["output-rds"]]))
saveRDS(res, opt[["output-rds"]])

message(sprintf("%s, writing summaryAnnotation to %s", Sys.time(), opt[["output-csv"]]))
sa_out <- data.frame(cell_id = rownames(sa), sa, check.names = FALSE,
                     stringsAsFactors = FALSE)
data.table::fwrite(sa_out, opt[["output-csv"]])

message(sprintf("%s, writing compact result to %s", Sys.time(), opt[["output-h5"]]))
out <- H5File$new(opt[["output-h5"]], mode = "w")
out[["cell_id"]] <- cell_id
# hdf5r cannot write NA character; map unresolved cells to the empty string, which
# write_celltypes.py surfaces as a distinct (empty) category / missing label.
out[["cell_type"]] <- ifelse(is.na(cell_type), "", cell_type)
out[["prob"]] <- as.numeric(prob)   # NA -> NaN in HDF5; read back as float by write_celltypes
out$close_all()

message(sprintf("Done. %d cells typed over a %d-leaf hierarchy.",
                length(cell_type), length(leaves)))
