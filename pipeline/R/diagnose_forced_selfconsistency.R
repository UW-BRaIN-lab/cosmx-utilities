#!/usr/bin/env Rscript
# Diagnostic: WHY does only 79.3% of already-named cells keep their name under the forced
# supervised run (75)?
#
# The puzzle. 75 scored the anchor cells with insitutypeML against the SAME converged profiles
# the semi-supervised fit (72) produced, minus its 27 de-novo columns. Reading the InSituType
# source, insitutype()'s final phase is itself a call to insitutypeML() on all cells with those
# same profiles, the same estimateBackground() path, one cohort ("all"), and reference_sds NULL
# for RNA. So both runs should be computing the same log-likelihood table — and removing 27
# LOSING columns cannot change an argmax. A cell whose best profile was already a named type
# should keep it. Agreement on the named cells should be 100%, not 79.3%.
#
# So one of those premises is false, and this decomposes which. insitutype() RETURNS $logliks
# (cells x 81), and 72 saved the whole result to anchor_typing.rds — so we can answer from
# stored data, with no re-scoring:
#
#   A. clust  vs  argmax over ALL 81 stored columns
#      Tests "the returned $clust is the argmax of the returned $logliks". If this is not
#      ~100%, insitutype() post-processes its assignment and $logliks is not the whole story.
#   B. clust  vs  argmax over the 54 NAMED stored columns  (named-labelled cells only)
#      The pure column-removal question, entirely inside the anchor fit's own numbers. This
#      MUST be 100% if A holds — it is the same argmax over a subset that still contains the
#      winner. Anything less means the stored logliks are not what produced clust.
#   C. argmax over 54 named stored columns  vs  75's top1_type
#      Isolates OUR re-scoring call. If A and B are clean but C is ~79%, the discrepancy is
#      in how flat_posteriors.R invokes insitutypeML, not in the biology or the reference.
#
# Whichever of A/B/C breaks localizes the cause. The script also reports, for the disagreeing
# cells, the stored top1-vs-top2 log-likelihood margin, so we can see whether they are
# boundary cells (small margin = a near-tie that any numerical difference flips) or confident
# ones (large margin = a real modelling difference).
#
# Inputs:
#   --typing-rds   anchor_typing.rds from 72 (the full insitutype() result: clust, logliks,
#                  profiles). NOT the compact h5 — that drops $logliks.
#   --posteriors   75's supervised_gbmap_posteriors.csv (cell_id, top1_type, top1_prob, ...).
# Outputs:
#   --output-csv   per-semi-supervised-label agreement table for A, B and C.
#   --output-margins-csv  margin distribution for agreeing vs disagreeing cells.
#
# Memory: the stored loglik matrix is ~2.5M x 81 doubles (~1.6GB) plus copies for the two
# max.col passes; ask for well over that.

suppressPackageStartupMessages({
  library(data.table)
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
stopifnot("--typing-rds is required" = !is.null(opt[["typing-rds"]]))
stopifnot("--posteriors is required" = !is.null(opt$posteriors))
stopifnot("--output-csv is required" = !is.null(opt[["output-csv"]]))

# De-novo clusters are InSituType's cluster_name_pool: 1-2 lowercase letters (K=27 overflows
# to "aa"). Named reference types always carry an uppercase letter, digit or separator.
DENOVO_RE <- "^[a-z]{1,2}$"

message(sprintf("%s, reading %s", Sys.time(), opt[["typing-rds"]]))
res <- readRDS(opt[["typing-rds"]])
message("result fields: ", paste(names(res), collapse = ", "))

if (is.null(res$logliks)) {
  stop("anchor_typing.rds has no $logliks — cannot run the diagnostic from stored data. ",
       "Re-score instead, or check that this is the insitutype() result and not a subset.")
}

ll <- as.matrix(res$logliks)
clust <- res$clust
message(sprintf("logliks: %d cells x %d types; clust: %d cells",
                nrow(ll), ncol(ll), length(clust)))

# Align clust to the loglik row order; they should already match but do not assume it.
if (!is.null(rownames(ll)) && !is.null(names(clust))) {
  common <- intersect(rownames(ll), names(clust))
  message(sprintf("cell ids in both logliks and clust: %d", length(common)))
  ll <- ll[common, , drop = FALSE]
  clust <- clust[common]
} else {
  message("WARNING: logliks or clust lack names; assuming identical row order.")
}

types <- colnames(ll)
named_cols <- types[!grepl(DENOVO_RE, types)]
denovo_cols <- types[grepl(DENOVO_RE, types)]
message(sprintf("stored profile columns: %d named + %d de-novo",
                length(named_cols), length(denovo_cols)))
stopifnot("no named columns in the stored logliks" = length(named_cols) > 0)

# --- A: is clust the argmax of the stored logliks? ----------------------------
message(sprintf("%s, A: argmax over all %d stored columns", Sys.time(), ncol(ll)))
argmax_all <- types[max.col(ll, ties.method = "first")]

# --- B: argmax restricted to the named columns (the column-removal question) --
message(sprintf("%s, B: argmax over the %d named columns", Sys.time(), length(named_cols)))
ll_named <- ll[, named_cols, drop = FALSE]
idx_named <- max.col(ll_named, ties.method = "first")
argmax_named <- named_cols[idx_named]

# top1-vs-top2 margin among the named columns: how close was the call?
top1 <- ll_named[cbind(seq_len(nrow(ll_named)), idx_named)]
ll_named[cbind(seq_len(nrow(ll_named)), idx_named)] <- -Inf
top2 <- ll_named[cbind(seq_len(nrow(ll_named)), max.col(ll_named, ties.method = "first"))]
margin <- top1 - top2
rm(ll_named); invisible(gc())

# --- C: our re-scored supervised call ----------------------------------------
message(sprintf("%s, reading %s", Sys.time(), opt$posteriors))
post <- data.table::fread(opt$posteriors, select = c("cell_id", "top1_type"))
data.table::setnames(post, "top1_type", "rescored")

dt <- data.table(cell_id = rownames(ll) %||% names(clust),
                 semisup = as.character(clust),
                 argmax_all = argmax_all,
                 forced_named = argmax_named,
                 margin = margin)
dt <- merge(dt, post, by = "cell_id", all.x = TRUE)
dt[, is_denovo := grepl(DENOVO_RE, semisup)]

n_missing <- dt[is.na(rescored), .N]
if (n_missing > 0) {
  message(sprintf("WARNING: %d cells had no re-scored call; excluded from C.", n_missing))
}

# --- headline numbers ---------------------------------------------------------
pct <- function(x) sprintf("%.2f%%", 100 * mean(x, na.rm = TRUE))
named_cells <- dt[is_denovo == FALSE]

cat("\n================ SELF-CONSISTENCY DECOMPOSITION ================\n")
cat(sprintf("all cells: %d   (named-labelled %d, de-novo-labelled %d)\n",
            nrow(dt), nrow(named_cells), dt[is_denovo == TRUE, .N]))
cat(sprintf("\nA  clust == argmax(all %d stored columns)          : %s\n",
            ncol(ll), pct(dt$semisup == dt$argmax_all)))
cat(  "     expected ~100%. If not, insitutype() post-processes its assignment\n")
cat(sprintf("\nB  clust == argmax(%d named stored columns)         : %s   [named-labelled cells]\n",
            length(named_cols), pct(named_cells$semisup == named_cells$forced_named)))
cat(  "     expected EXACTLY 100%. Removing losing columns cannot change an argmax,\n")
cat(  "     so anything less means the stored logliks did not produce clust.\n")
cat(sprintf("\nC  argmax(named stored) == our re-scored top1        : %s\n",
            pct(dt$forced_named == dt$rescored)))
cat(  "     isolates flat_posteriors.R's insitutypeML call against the fit's own numbers.\n")
cat(sprintf("\nOBSERVED (what 75b reported)  clust == re-scored top1: %s   [named-labelled cells]\n",
            pct(named_cells$semisup == named_cells$rescored)))
cat("===============================================================\n\n")

# --- margins: are the disagreeing cells near-ties? ----------------------------
dt[, agrees_C := forced_named == rescored]
margins <- dt[!is.na(agrees_C), .(
  n = .N,
  margin_p10 = round(quantile(margin, 0.10), 3),
  margin_median = round(median(margin), 3),
  margin_p90 = round(quantile(margin, 0.90), 3)
), by = .(agrees_C)]
cat("Stored top1-vs-top2 loglik margin among the named columns:\n")
print(margins)
cat("  A small margin for the disagreeing cells = near-ties that any numerical\n")
cat("  difference flips (benign). A large margin = a real modelling difference.\n\n")

# --- per-label table ----------------------------------------------------------
per_label <- dt[, .(
  n_cells = .N,
  is_denovo = is_denovo[1],
  A_clust_eq_argmax_all = round(100 * mean(semisup == argmax_all), 2),
  B_clust_eq_forced_named = ifelse(is_denovo[1], NA_real_,
                                   round(100 * mean(semisup == forced_named), 2)),
  C_forced_eq_rescored = round(100 * mean(forced_named == rescored, na.rm = TRUE), 2),
  observed_clust_eq_rescored = ifelse(is_denovo[1], NA_real_,
                                      round(100 * mean(semisup == rescored, na.rm = TRUE), 2)),
  median_margin = round(median(margin), 3)
), by = semisup][order(-n_cells)]

cat("Per semi-supervised label (largest first):\n")
print(head(per_label, 25))

data.table::fwrite(per_label, opt[["output-csv"]])
message(sprintf("\nWrote %s (%d labels)", opt[["output-csv"]], nrow(per_label)))

if (!is.null(opt[["output-margins-csv"]])) {
  data.table::fwrite(margins, opt[["output-margins-csv"]])
  message(sprintf("Wrote %s", opt[["output-margins-csv"]]))
}
message("Done.")
