#!/usr/bin/env Rscript
# De-novo cluster-number sweep for the anchor reference rebuild.
#
# WHY: the pilot anchor typing ran n_clusts=10:20 and InSituType selected K=20 -- the top
# of the range. InSituType picks K by MINIMUM AIC over the range (chooseClusterNumber:
# AIC = 2*K*n_genes - 2*loglik), evaluated on a small geo-sketched subset (default
# n_chooseclusternumber = 2000). A minimum sitting on the boundary is a CENSORED optimum:
# the true best K is almost certainly higher. This script re-runs that same AIC criterion
# across a WIDER range (default 15:35) and reports the whole curve, so we can see the
# minimum actually turn over instead of guessing.
#
# It calls InSituType::chooseClusterNumber directly -- the exact function insitutype() uses
# internally -- with fixed_profiles = the (un-updated) reference, matching insitutype()'s
# own call. Two knobs make the curve trustworthy at high K where the 2000-cell default is
# noisy: --subset-size (cells per fit) and --n-reps (independent subsamples, averaged).
#
# Pass --keep-genes to sweep on the PRUNED panel (select_informative_genes.R output). Gene
# pruning changes K's AIC penalty directly -- 2*n_genes per added cluster -- so run the
# sweep on both the full and pruned panels and compare where each turns over.
#
# Inputs mirror 72_anchor_typing.sh (anchor_input.h5 + reference CSV).
# Outputs (into --out-dir):
#   k_sweep.csv       n_clusts x {loglik,aic,bic} per rep + mean; best_by_aic flag
#   k_sweep.png       mean AIC and log-likelihood vs K (visual turn-over check)
#   k_sweep_summary.txt  best K by mean AIC/BIC + a CENSORED warning if it hits an edge
#
# Usage:
#   Rscript pipeline/R/choose_denovo_k.R \
#     --input anchor_input.h5 --reference reference.csv --out-dir k_sweep \
#     --n-clusts 15:35 [--keep-genes kept_genes.txt] \
#     [--subset-size 20000] [--n-reps 3] [--max-iters 20] [--max-cells 500000]

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(InSituType)
})
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "read_anchor_h5.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument: %s", a))
    key <- sub("^--", "", a)
    if (grepl("=", key)) { kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]; i <- i + 1
    } else { out[[key]] <- args[[i + 1]]; i <- i + 2 }
  }
  out
}

parse_range <- function(s) {
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
out_dir     <- opt[["out-dir"]] %||% "k_sweep"
n_clusts    <- parse_range(opt[["n-clusts"]] %||% "15:35")
subset_size <- as.integer(opt[["subset-size"]] %||% 20000)
n_reps      <- as.integer(opt[["n-reps"]] %||% 3)
max_iters   <- as.integer(opt[["max-iters"]] %||% 20)
max_cells   <- if (is.null(opt[["max-cells"]])) NULL else as.integer(opt[["max-cells"]])
seed        <- as.integer(opt$seed %||% 0)
keep_path   <- opt[["keep-genes"]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- anchor counts + reference ------------------------------------------------
dat <- read_anchor_h5(opt$input, max_cells = max_cells, seed = seed)
x <- dat$x; neg <- dat$neg

ref_dt <- data.table::fread(opt$reference)
ref <- as.matrix(ref_dt[, -1L]); rownames(ref) <- as.character(ref_dt[[1]])
storage.mode(ref) <- "double"

# --- optional gene pruning (FAQ subset) ---------------------------------------
if (!is.null(keep_path)) {
  keep_genes <- readLines(keep_path)
  keep_genes <- intersect(keep_genes, colnames(x))
  stopifnot("no kept genes overlap the panel" = length(keep_genes) > 0)
  x <- x[, keep_genes, drop = FALSE]
  message(sprintf("pruned panel: %d genes (from --keep-genes %s)", ncol(x), keep_path))
}

shared <- intersect(colnames(x), rownames(ref))
stopifnot("no genes shared between counts and reference" = length(shared) > 0)
message(sprintf("panel %d genes; reference overlap %d genes; AIC penalty per cluster = 2*%d",
                ncol(x), length(shared), ncol(x)))

subset_size <- min(subset_size, nrow(x))

# --- sweep: chooseClusterNumber over the range, n_reps independent subsamples ---
# chooseClusterNumber draws its own random subset_size cells each call, so repeating the
# call gives independent draws; averaging the curves de-noises the AIC minimum at high K.
message(sprintf("%s, K-sweep n_clusts=%d:%d, subset_size=%d, n_reps=%d, max_iters=%d",
                Sys.time(), min(n_clusts), max(n_clusts), subset_size, n_reps, max_iters))
reps <- vector("list", n_reps)
for (r in seq_len(n_reps)) {
  set.seed(seed + r)
  message(sprintf("  rep %d/%d ...", r, n_reps))
  cc <- InSituType::chooseClusterNumber(
    counts = x, neg = neg[rownames(x)],
    assay_type = "rna",
    fixed_profiles = ref,
    n_clusts = n_clusts,
    max_iters = max_iters,
    subset_size = subset_size,
    align_genes = TRUE,
    plotresults = FALSE
  )
  reps[[r]] <- data.table(rep = r, n_clusts = cc$n_clusts,
                          loglik = cc$loglik, aic = cc$aic, bic = cc$bic)
}
long <- data.table::rbindlist(reps)

# --- aggregate + pick best -----------------------------------------------------
agg <- long[, .(loglik = mean(loglik), aic = mean(aic), bic = mean(bic)), by = n_clusts]
data.table::setorder(agg, n_clusts)
best_aic <- agg$n_clusts[which.min(agg$aic)]
best_bic <- agg$n_clusts[which.min(agg$bic)]
agg[, best_by_aic := n_clusts == best_aic]
data.table::fwrite(rbind(
  long[, .(rep = as.character(rep), n_clusts, loglik, aic, bic, best_by_aic = NA)],
  agg[, .(rep = "mean", n_clusts, loglik, aic, bic, best_by_aic)]
), file.path(out_dir, "k_sweep.csv"))

lo <- min(n_clusts); hi <- max(n_clusts)
censored <- best_aic == lo || best_aic == hi
summary_lines <- c(
  sprintf("range swept:        %d:%d  (%d values)", lo, hi, length(n_clusts)),
  sprintf("genes:              %d  (AIC complexity penalty = 2*genes per cluster)", ncol(x)),
  sprintf("subset-size/n-reps: %d / %d", subset_size, n_reps),
  sprintf("best K by mean AIC: %d", best_aic),
  sprintf("best K by mean BIC: %d", best_bic),
  if (censored)
    sprintf("WARNING: best K = %d is at the %s edge of the range -- optimum is CENSORED; extend --n-clusts %s and re-run",
            best_aic, if (best_aic == lo) "LOW" else "HIGH",
            if (best_aic == hi) sprintf("%d:%d", lo, hi + 10) else sprintf("%d:%d", max(1, lo - 10), hi))
  else
    sprintf("OK: AIC turns over at an interior K = %d (not censored)", best_aic)
)
writeLines(summary_lines, file.path(out_dir, "k_sweep_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# --- plot ----------------------------------------------------------------------
png(file.path(out_dir, "k_sweep.png"), width = 1100, height = 500, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 1))
plot(agg$n_clusts, agg$aic, type = "b", pch = 19, xlab = "de-novo clusters (K)",
     ylab = "mean AIC (lower = better)", main = "AIC vs K")
abline(v = best_aic, col = "red", lty = 2); axis(3, at = best_aic, labels = best_aic, col.axis = "red")
plot(agg$n_clusts, agg$loglik, type = "b", pch = 19, xlab = "de-novo clusters (K)",
     ylab = "mean log-likelihood", main = "log-likelihood vs K")
par(op); dev.off()
message("wrote ", file.path(out_dir, "k_sweep.png"))
