#!/usr/bin/env Rscript
# Stage 5h (Phase 3): "is 48% a floor?" — recluster the Low_signal pool, derive profiles, and
# reassign, to test whether the untyped fraction is still reference/transfer-limited (would
# collapse on another CosMx-native rescue) or an irreducible core (barely moves).
#
# The 48% Low_signal is already down from ~90% after one CosMx-native rescue: most of the
# initial miss was platform-transfer residue (GBmap is dissociated scRNA in droplet count-
# space). So 48% may still be a way-station. This runs one MORE rescue iteration ON THE
# Low_signal cells only: InSituType semi-supervised against the fixed InSituTree profile set
# (so a cell can snap to an existing named type) PLUS de-novo clustering over an n_clusts range
# (so genuinely off-reference populations surface as new profiles). No rescale, no reference
# update — the profiles are already CosMx-scale (as InSituTree used them).
#
# Caution baked into the interpretation (compare_rescue.py): each rescue iteration also plants
# finer anchors along any continuous structure, which can manufacture boundary ties — so a big
# drop is only meaningful if the rescued cells carry real lineage programs, not just new de-novo
# tiles of the same manifold. The Phase-2 top-two structure is the arbiter.
#
# Inputs:
#   --input      insitutype_input.h5 (stage 4a): /counts (CSC genes x cells), /genes, /cell_id, /neg.
#   --labels-h5  insitutree_result.h5 (stage 4b): /cell_id, /cell_type — selects the Low_signal pool.
#   --profiles   fixed profile matrix CSV (genes x types; default insitutree_profiles.csv).
#   --lowsignal-label  cell_type value marking the pool (default Low_signal).
#   --n-clusts   de-novo cluster range "lo:hi" (default 10:20) — new profiles from the pool.
# Outputs:
#   --output-csv       per Low_signal cell: cell_id, rescue_label, rescue_prob, is_denovo.
#   --output-profiles  the rescued profile matrix (genes x types incl. new de-novo columns) CSV.
#   --output-rds       full insitutype() result for the subset.
#
# Memory: only the Low_signal subset (~1.1M cells) x genes + a cells x types loglik matrix; a
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

parse_n_clusts <- function(s) {
  if (grepl(":", s, fixed = TRUE)) {
    ab <- as.integer(strsplit(s, ":", fixed = TRUE)[[1]])
    stopifnot("--n-clusts range must be 'lo:hi'" = length(ab) == 2L)
    return(seq.int(ab[[1]], ab[[2]]))
  }
  as.integer(s)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("input", "labels-h5", "profiles", "output-csv", "output-profiles")) {
  if (is.null(opt[[k]])) stop(sprintf("--%s is required", k))
}
lowsignal_label <- opt[["lowsignal-label"]] %||% "Low_signal"
n_clusts <- parse_n_clusts(opt[["n-clusts"]] %||% "10:20")
seed <- as.integer(opt$seed %||% 0)

# --- read counts (stage-4a input) --------------------------------------------
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

# --- select the Low_signal pool from the InSituTree labels -------------------
message(sprintf("%s, reading labels %s", Sys.time(), opt[["labels-h5"]]))
lf <- H5File$new(opt[["labels-h5"]], mode = "r")
lab_id <- lf[["cell_id"]][]
lab_ct <- lf[["cell_type"]][]
lf$close_all()
names(lab_ct) <- lab_id
ls_ids <- intersect(cell_id, lab_id[lab_ct == lowsignal_label])
stopifnot("no Low_signal cells found" = length(ls_ids) > 0)
message(sprintf("Low_signal pool: %d cells (of %d typed)", length(ls_ids), length(cell_id)))

x <- Matrix::t(counts_gxc[, ls_ids, drop = FALSE]); rm(counts_gxc)  # cells x genes
neg <- neg[ls_ids]

# --- fixed profile set (genes x types) ---------------------------------------
message(sprintf("%s, reading profiles %s", Sys.time(), opt$profiles))
ref_dt <- data.table::fread(opt$profiles)
ref <- as.matrix(ref_dt[, -1L]); rownames(ref) <- as.character(ref_dt[[1]])
storage.mode(ref) <- "double"
shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and profiles" = length(shared) > 0)
x <- x[, shared, drop = FALSE]; ref <- ref[shared, , drop = FALSE]
message(sprintf("shared genes: %d; profiles: %d types", length(shared), ncol(ref)))

# --- one more rescue iteration on the pool -----------------------------------
# semi-supervised: reference lets a cell snap to an existing type; n_clusts range discovers
# NEW de-novo profiles from the pool. No rescale / no reference update (profiles already
# CosMx-scale, matching how InSituTree consumed them).
if (seed != 0) set.seed(seed)
message(sprintf("%s, insitutype rescue: %d cells, n_clusts=%d:%d, %d reference profiles",
                Sys.time(), nrow(x), min(n_clusts), max(n_clusts), ncol(ref)))
res <- InSituType::insitutype(
  x = x,
  neg = neg[rownames(x)],
  n_clusts = n_clusts,
  reference_profiles = ref,
  update_reference_profiles = FALSE,
  rescale = FALSE,
  refit = FALSE,
  refinement = FALSE
)

message("rescue assignment counts:")
print(sort(table(res$clust), decreasing = TRUE))

clust <- res$clust[ls_ids]
prob <- res$prob[ls_ids]
# de-novo clusters are the single-letter labels InSituType invents (a, b, ...); a named match
# means the cell snapped to an existing InSituTree profile.
is_denovo <- !(clust %in% colnames(ref))

out <- data.table(cell_id = ls_ids, rescue_label = as.character(clust),
                  rescue_prob = as.numeric(prob), is_denovo = as.logical(is_denovo))
data.table::fwrite(out, opt[["output-csv"]])

prof <- as.data.frame(res$profiles)
prof <- cbind(gene = rownames(res$profiles), prof)
data.table::fwrite(prof, opt[["output-profiles"]])
if (!is.null(opt[["output-rds"]])) saveRDS(res, opt[["output-rds"]])

n_named <- sum(!is_denovo)
message(sprintf("Done. %d Low_signal cells reassigned: %d named (%.1f%%), %d de-novo (%.1f%%).",
                length(ls_ids), n_named, 100 * n_named / length(ls_ids),
                sum(is_denovo), 100 * sum(is_denovo) / length(ls_ids)))
