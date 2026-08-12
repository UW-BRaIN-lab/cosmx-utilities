#!/usr/bin/env Rscript
# InSituType gene-pruning test: run cell typing on the FULL 6k panel vs a PRUNED
# informative-gene subset, with everything else identical, and compare how many
# cells fall into de-novo clusters.
#
# Motivation (retina Refit=FALSE collapse): ~66-74% of cells collapsed into one
# flat low-signal de-novo cluster. The InSituType 6k-panel FAQ recommends typing
# on a well-chosen ~3-5k gene subset. This script tests whether pruning reduces
# the collapse, using the SAME custom combined reference (retina + brain + optic
# nerve) and the SAME gentle params (rescale=TRUE, refit=FALSE) as the current run.
#
# Gene selection follows the FAQ's TWO criteria, kept as a UNION so a gene that is
# informative in ONLY ONE tissue (e.g. retinal rod genes) is NOT dropped:
#   (1) solidly above background in the CosMx data, OR
#   (2) moderate-to-high in at least one reference profile.
#
# Run on Hyak (Klone) — this is a full 2.3M-cell InSituType run per panel. Start
# with --slides on a few slides to validate before the full cohort.
#
# Usage:
#   Rscript pipeline/R/insitutype_pruning_test.R \
#     --exprmat-dir /path/to/flatFiles \  # dir of <slide>/<slide>_exprMat_file.csv.gz
#     --reference   /path/to/combined_reference.csv \  # genes x cell types
#     --out-dir     pruning_test \
#     --n-clusts    10:20            # match your AtoMx run
#     [--slides UWA7575eyes,eyes7517]  # optional subset for a quick test
#     [--bg-mult 3] [--ref-quantile 0.5]  # gene-selection knobs (tune to ~3-5k)

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(Matrix); library(InSituType)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
opt <- parse_args(OptionParser(option_list = list(
  make_option("--exprmat-dir", type = "character"),
  make_option("--reference",   type = "character"),
  make_option("--out-dir",     type = "character", default = "pruning_test"),
  make_option("--n-clusts",    type = "character", default = "10:20"),
  make_option("--slides",      type = "character", default = NULL),
  make_option("--bg-mult",     type = "double",    default = 3),   # crit 1 strictness
  make_option("--ref-quantile",type = "double",    default = 0.5), # crit 2 strictness
  make_option("--seed",        type = "integer",   default = 0)
)))
stopifnot(!is.null(opt[["exprmat-dir"]]), !is.null(opt$reference))
dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)
parse_range <- function(s) { p <- as.integer(strsplit(s, ":")[[1]]); if (length(p) == 2) p[1]:p[2] else p }
n_clusts <- parse_range(opt[["n-clusts"]])
is_control <- function(g) grepl("^(Negative|NegPrb|SystemControl|FalseCode)", g, ignore.case = TRUE)

# --- 1. Load exprMat flat files → counts (cells x genes) + per-cell neg --------
files <- Sys.glob(file.path(opt[["exprmat-dir"]], "*", "*_exprMat_file.csv.gz"))
if (!is.null(opt$slides)) {
  keep <- strsplit(opt$slides, ",")[[1]]
  files <- files[basename(dirname(files)) %in% keep]
}
stopifnot("no exprMat files found" = length(files) > 0)
message(sprintf("%s: reading %d exprMat file(s)", Sys.time(), length(files)))

blocks <- list(); negs <- list()
for (f in files) {
  slide <- basename(dirname(f))
  dt <- fread(f)
  idcols <- intersect(c("fov", "cell_ID"), names(dt))
  gcols <- setdiff(names(dt), idcols)
  ctrl <- gcols[is_control(gcols)]
  genes <- setdiff(gcols, ctrl)
  m <- as(as.matrix(dt[, ..genes]), "CsparseMatrix")             # cells x genes
  rownames(m) <- paste(slide, dt$fov, dt$cell_ID, sep = "_")
  blocks[[f]] <- m
  negs[[f]] <- rowMeans(as.matrix(dt[, ..ctrl]))                  # per-cell background
  message(sprintf("  %s: %d cells x %d genes", slide, nrow(m), length(genes)))
}
counts <- do.call(rbind, blocks); neg <- unlist(negs, use.names = TRUE)
names(neg) <- rownames(counts)
message(sprintf("total: %d cells x %d genes; median neg %.4f",
                nrow(counts), ncol(counts), median(neg)))

# --- 2. Reference profiles (genes x cell types) -------------------------------
ref_dt <- fread(opt$reference)
ref <- as.matrix(ref_dt[, -1]); rownames(ref) <- ref_dt[[1]]
shared <- intersect(colnames(counts), rownames(ref))
stopifnot("no shared genes" = length(shared) > 0)
counts <- counts[, shared]; ref <- ref[shared, , drop = FALSE]
message(sprintf("shared genes: %d ; reference cell types: %d", length(shared), ncol(ref)))

# --- 3. Informative-gene selection (FAQ criteria, UNION) ----------------------
gene_mean <- Matrix::colMeans(counts)
bg <- mean(neg)
crit1 <- gene_mean > opt[["bg-mult"]] * bg                       # solidly above bg
ref_max <- apply(ref, 1, max)                                    # max across cell types
crit2 <- ref_max >= quantile(ref_max, opt[["ref-quantile"]])     # moderate-to-high in >=1 profile
keep <- crit1 | crit2                                            # UNION → tissue-specific genes kept
message(sprintf("gene selection: crit1(above-bg)=%d, crit2(ref-informative)=%d, UNION kept=%d of %d",
                sum(crit1), sum(crit2), sum(keep), length(shared)))
message("  → tune --bg-mult / --ref-quantile so 'kept' lands ~3000-5000 (FAQ target)")
writeLines(shared[keep], file.path(opt[["out-dir"]], "kept_genes.txt"))

# --- 4. Run InSituType: FULL panel vs PRUNED panel (identical params) ---------
ref_types <- colnames(ref)
run_it <- function(x, r, tag) {
  message(sprintf("%s: insitutype [%s] %d cells x %d genes, n_clusts=%s",
                  Sys.time(), tag, nrow(x), ncol(x), opt[["n-clusts"]]))
  res <- InSituType::insitutype(
    x = x, neg = neg[rownames(x)], n_clusts = n_clusts,
    reference_profiles = r, update_reference_profiles = TRUE,
    rescale = TRUE, refit = FALSE, refinement = FALSE)
  saveRDS(res, file.path(opt[["out-dir"]], sprintf("insitutype_%s.rds", tag)))
  denovo <- !(res$clust %in% ref_types)
  message(sprintf("  [%s] de-novo fraction: %.1f%% (%d/%d cells)",
                  tag, 100 * mean(denovo), sum(denovo), length(denovo)))
  res
}
full   <- run_it(counts, ref, "full")
pruned <- run_it(counts[, shared[keep]], ref[shared[keep], , drop = FALSE], "pruned")

# --- 5. Compare ---------------------------------------------------------------
denovo_frac <- function(res) mean(!(res$clust %in% ref_types))
tab <- table(full = ifelse(full$clust %in% ref_types, "named", "de-novo"),
             pruned = ifelse(pruned$clust %in% ref_types, "named", "de-novo"))
sink(file.path(opt[["out-dir"]], "summary.txt"))
cat(sprintf("genes: full=%d, pruned=%d\n", length(shared), sum(keep)))
cat(sprintf("de-novo fraction: full=%.1f%%  pruned=%.1f%%\n",
            100 * denovo_frac(full), 100 * denovo_frac(pruned)))
cat("\nnamed/de-novo crosstab (full x pruned):\n"); print(tab)
cat(sprintf("\ncells rescued (de-novo in full → named in pruned): %d (%.1f%%)\n",
            tab["de-novo", "named"], 100 * tab["de-novo", "named"] / length(full$clust)))
sink()
message("wrote ", file.path(opt[["out-dir"]], "summary.txt"))
