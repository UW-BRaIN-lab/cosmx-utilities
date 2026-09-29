#!/usr/bin/env Rscript
# Diagnostic: is the fit's returned $clust the argmax of its own returned $logliks?
#
# History. A fully supervised re-score of the anchor cells (insitutypeML against the converged
# profiles minus the 27 de-novo columns) kept the name of only 79.3% of already-named cells.
# Reading the InSituType source, that should be 100%: insitutype()'s final phase is itself an
# insitutypeML() call on all cells with those same profiles, and removing 27 LOSING columns
# cannot change an argmax. That re-score was dropped (it is on the tag pr79-archive-09f2397);
# this script answers the question from stored data instead. insitutype() RETURNS $logliks
# (cells x 81) and 72 saved the whole result to anchor_typing.rds, so no re-scoring is needed:
#
#   A. clust  vs  argmax over ALL 81 stored columns
#      Tests "the returned $clust is the argmax of the returned $logliks". If this is not
#      ~100%, insitutype() post-processes its assignment and $logliks is not the whole story.
#      (It is not: 75o shows the gap is InSituType's pinned cells, whose labels are overwritten.)
#   B. clust  vs  argmax over the 54 NAMED stored columns  (named-labelled cells only)
#      The pure column-removal question, entirely inside the anchor fit's own numbers. This
#      MUST equal A's agreement on named-labelled cells: the same argmax over a subset that
#      still contains the winner.
#
# The script also reports, for the mismatching cells, the stored top1-vs-top2 log-likelihood
# margin among the named columns, so we can see whether they are boundary cells (small margin)
# or confident ones (large margin).
#
# Inputs:
#   --typing-rds   anchor_typing.rds from 72 (the full insitutype() result: clust, logliks,
#                  profiles). NOT the compact h5 — that drops $logliks.
# Outputs:
#   --output-csv   per-semi-supervised-label agreement table for A and B.
#   --output-margins-csv  margin distribution for matching vs mismatching cells.
#   --output-mismatch-csv one row per cell whose assigned label is NOT its best-scoring
#                       profile, with the loglik deficit, worst first.
#   --output-forced-csv   THE PER-CELL FORCED-NAMED CALL: (cell_id, top1_type, top1_prob,
#                  top2_type, top2_prob) plus `margin`, read straight off the anchor fit's own
#                  stored logliks over the named columns — "what the fit would call each cell
#                  with the de-novo option removed". Feed it to 75b, 75f and the 75d-75m views.
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
n_cells <- nrow(ll_named)
top1 <- ll_named[cbind(seq_len(n_cells), idx_named)]
# Softmax over the named columns, row max for numerical stability — the same posterior
# flat_posteriors.R reports, so the emitted CSV keeps the (cell_id, top1_type, ...) schema.
lse <- top1 + log(rowSums(exp(ll_named - top1)))
prob1 <- exp(top1 - lse)
ll_named[cbind(seq_len(n_cells), idx_named)] <- -Inf
idx2 <- max.col(ll_named, ties.method = "first")
top2 <- ll_named[cbind(seq_len(n_cells), idx2)]
argmax_named_2 <- named_cols[idx2]
prob2 <- exp(top2 - lse)
margin <- top1 - top2
rm(ll_named); invisible(gc())

# Per-cell loglik of the label the cell WAS given, and of the one that scored best. Both are
# read by row position, so dt must stay in ll's row order (do not merge/sort it by cell_id).
assigned_idx <- match(as.character(clust), types)
if (anyNA(assigned_idx)) {
  message(sprintf("WARNING: %d cells carry a label with no loglik column; their deficit is NA.",
                  sum(is.na(assigned_idx))))
}
row_seq <- seq_len(nrow(ll))
assigned_loglik <- ifelse(is.na(assigned_idx), NA_real_,
                          ll[cbind(row_seq, ifelse(is.na(assigned_idx), 1L, assigned_idx))])
argmax_loglik <- ll[cbind(row_seq, max.col(ll, ties.method = "first"))]

dt <- data.table(cell_id = rownames(ll) %||% names(clust),
                 semisup = as.character(clust),
                 argmax_all = argmax_all,
                 assigned_loglik = assigned_loglik,
                 argmax_loglik = argmax_loglik,
                 forced_named = argmax_named,
                 top1_prob = prob1,
                 top2_type = argmax_named_2,
                 top2_prob = prob2,
                 margin = margin)
dt[, is_denovo := grepl(DENOVO_RE, semisup)]

# --- headline numbers ---------------------------------------------------------
pct <- function(x) sprintf("%.2f%%", 100 * mean(x, na.rm = TRUE))
named_cells <- dt[is_denovo == FALSE]

cat("\n================ SELF-CONSISTENCY DECOMPOSITION ================\n")
cat(sprintf("all cells: %d   (named-labelled %d, de-novo-labelled %d)\n",
            nrow(dt), nrow(named_cells), dt[is_denovo == TRUE, .N]))
cat(sprintf("\nA  clust == argmax(all %d stored columns)          : %s\n",
            ncol(ll), pct(dt$semisup == dt$argmax_all)))
cat(  "     would be ~100% if clust were the argmax. The shortfall is insitutype()'s\n")
cat(  "     post-processing: the pinned cells' labels are overwritten (see 75o).\n")
cat(sprintf("\nB  clust == argmax(%d named stored columns)         : %s   [named-labelled cells]\n",
            length(named_cols), pct(named_cells$semisup == named_cells$forced_named)))
cat(  "     equals A on these cells: removing losing columns cannot change an argmax,\n")
cat(  "     so any difference from A would mean the stored logliks did not produce clust.\n")
cat("===============================================================\n\n")

# --- margins: are the mismatching cells near-ties? ----------------------------
dt[, mismatch := semisup != argmax_all]
margins <- dt[, .(
  n = .N,
  margin_p10 = round(quantile(margin, 0.10), 3),
  margin_median = round(median(margin), 3),
  margin_p90 = round(quantile(margin, 0.90), 3)
), by = .(mismatch)]
cat("Stored top1-vs-top2 loglik margin among the named columns:\n")
print(margins)
cat("  A small margin for the mismatching cells = near-ties that any numerical\n")
cat("  difference flips (benign). A large margin = the label contradicts the evidence.\n\n")

# --- per-label table ----------------------------------------------------------
per_label <- dt[, .(
  n_cells = .N,
  is_denovo = is_denovo[1],
  A_clust_eq_argmax_all = round(100 * mean(semisup == argmax_all), 2),
  B_clust_eq_forced_named = ifelse(is_denovo[1], NA_real_,
                                   round(100 * mean(semisup == forced_named), 2)),
  median_margin = round(median(margin), 3)
), by = semisup][order(-n_cells)]

cat("Per semi-supervised label (largest first):\n")
print(head(per_label, 25))

data.table::fwrite(per_label, opt[["output-csv"]])
message(sprintf("\nWrote %s (%d labels)", opt[["output-csv"]], nrow(per_label)))

if (!is.null(opt[["output-forced-csv"]])) {
  forced <- data.table(cell_id = dt$cell_id, top1_type = dt$forced_named,
                       top1_prob = dt$top1_prob, top2_type = dt$top2_type,
                       top2_prob = dt$top2_prob, margin = dt$margin)
  data.table::fwrite(forced, opt[["output-forced-csv"]])
  message(sprintf("Wrote %s (%d cells) — the CORRECTED forced-named call, from the fit's own logliks",
                  opt[["output-forced-csv"]], nrow(forced)))
}

if (!is.null(opt[["output-margins-csv"]])) {
  data.table::fwrite(margins, opt[["output-margins-csv"]])
  message(sprintf("Wrote %s", opt[["output-margins-csv"]]))
}
if (!is.null(opt[["output-mismatch-csv"]])) {
  # A above reports the RATE at which the assigned label differs from the best-scoring profile.
  # This is the list itself, with the DEFICIT — how many log-likelihood units better the winning
  # profile scored — because that is what separates a near-tie, which any numerical wobble would
  # flip, from a cell whose label the fit's own evidence plainly contradicts.
  dt[, deficit := argmax_loglik - assigned_loglik]
  mism <- dt[semisup != argmax_all][order(-deficit)]
  cat(sprintf("\n================ LABEL vs ARGMAX MISMATCHES ================\n"))
  cat(sprintf("%d of %d cells (%.2f%%) carry a label that is NOT their best-scoring profile.\n",
              nrow(mism), nrow(dt), 100 * nrow(mism) / nrow(dt)))
  if (nrow(mism) > 0) {
    cat(sprintf("deficit (loglik units the argmax beats the assigned label by): "))
    cat(sprintf("median %.1f, p90 %.1f, max %.1f\n",
                median(mism$deficit, na.rm = TRUE),
                quantile(mism$deficit, 0.9, na.rm = TRUE),
                max(mism$deficit, na.rm = TRUE)))
    pairs <- mism[, .(n = .N, median_deficit = round(median(deficit, na.rm = TRUE), 1)),
                  by = .(assigned = semisup, argmax = argmax_all)][order(-n)]
    cat("\nWhere the mismatched cells would have gone (top 20 label -> argmax pairs):\n")
    print(head(pairs, 20))
    by_label <- dt[, .(n_cells = .N, n_mismatch = sum(semisup != argmax_all),
                       pct_mismatch = round(100 * mean(semisup != argmax_all), 1)),
                   by = .(assigned = semisup)][order(-pct_mismatch)][n_cells >= 100]
    cat("\nBy assigned label, worst first (labels with >=100 cells):\n")
    print(head(by_label, 20))
  }
  cat("===========================================================\n\n")
  data.table::fwrite(
    mism[, .(cell_id, assigned = semisup, argmax = argmax_all,
             assigned_loglik = round(assigned_loglik, 3),
             argmax_loglik = round(argmax_loglik, 3),
             deficit = round(deficit, 3))],
    opt[["output-mismatch-csv"]])
  message(sprintf("Wrote %s (%d cells)", opt[["output-mismatch-csv"]], nrow(mism)))
}

message("Done.")
