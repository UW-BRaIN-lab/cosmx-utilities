#!/usr/bin/env Rscript
# Would the InSituType FAQ's 80% posterior threshold catch the flat cells if there were no flat leaf?
#
# The b decision for the InSituTree rebuild is: carry b's flat profile as a Low_signal_denovo
# leaf, or drop it. InSituTree has no Low_signal outcome of its own, so without the leaf every
# b cell is forced onto a named leaf, and "Low_signal" becomes whatever we flag afterwards. The
# FAQ's recipe is to call cells below a posterior of 0.8 "unclassified". But b's cells are
# near-ties (median top1-vs-top2 margin 4.1 against 25.1 for confident cells), and with
# thousands of genes a margin of 4 can still read as a posterior near 0.98. If most b cells
# clear 0.8 anyway, the threshold alone will not flag them and a margin cut is needed as well.
#
# This answers it from the fit's own stored 2.54M x 81 log-likelihood table, exactly as 75c
# derives the forced call: drop the de-novo columns, take the posterior over the named ones.
#
# THE ONE THING THAT CAN BE WRONG is the posterior formula. InSituType's `prob` may be a plain
# softmax of the log-likelihoods, or it may add a log prior of cluster frequencies. So both are
# tested against the stored `prob`, on cells whose label IS their argmax (pinned cells carry an
# overridden label, so their `prob` is not comparable), and the one that reproduces it is used.
# If neither does, the outputs are still written but every row says formula_validated = FALSE.
#
# Inputs:
#   --typing-rds  anchor_typing.rds from 72 (clust, prob, logliks). NOT the h5, which drops
#                 $logliks.
# Options:
#   --thresholds  posterior cut-offs to report (default 0.8,0.95,0.99)
#   --margins     top1-vs-top2 log-likelihood cut-offs to report (default 2,5,10,25)
#   --focus       letter to cross-tabulate posterior against margin (default b)
#   --chunk       rows per block, to bound memory (default 250000)
#   --tolerance   mean |formula - stored prob| that counts as reproduced (default 1e-3)
# Outputs (--output-dir):
#   posterior_formula_check.csv      how well each candidate formula reproduces stored `prob`
#   named_posterior_by_label.csv     per label: share below each posterior and margin cut
#   <focus>_posterior_by_margin.csv  the focus letter: posterior bin x margin bin
#   named_posterior_per_cell.csv.gz  per cell: label, named-only top1/top2, posterior, margin

suppressPackageStartupMessages({ library(data.table) })

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
stopifnot("--output-dir is required" = !is.null(opt[["output-dir"]]))
outdir <- opt[["output-dir"]]; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

THRESHOLDS <- as.numeric(strsplit(opt$thresholds %||% "0.8,0.95,0.99", ",")[[1]])
MARGIN_CUTS <- as.numeric(strsplit(opt$margins %||% "2,5,10,25", ",")[[1]])
FOCUS <- opt$focus %||% "b"
CHUNK <- as.integer(opt$chunk %||% 250000)
TOL <- as.numeric(opt$tolerance %||% 1e-3)
# Same convention as 75c: de-novo clusters are 1-2 lowercase letters (K=27 overflows to "aa").
DENOVO_RE <- "^[a-z]{1,2}$"
POSTERIOR_BREAKS <- c(0, 0.5, 0.8, 0.95, 0.99, 1)

message(sprintf("%s, reading %s", Sys.time(), opt[["typing-rds"]]))
res <- readRDS(opt[["typing-rds"]])
if (is.null(res$logliks)) stop("anchor_typing.rds has no $logliks -- this must be the full ",
                               "insitutype() result, not the compact h5.")
if (is.null(res$prob)) stop("anchor_typing.rds has no $prob, so the posterior formula cannot ",
                            "be validated against the fit.")
ll <- as.matrix(res$logliks)
clust <- res$clust
stored_prob <- res$prob
if (!is.null(names(clust))) {
  common <- intersect(rownames(ll), names(clust))
  ll <- ll[common, , drop = FALSE]; clust <- clust[common]
  if (!is.null(names(stored_prob))) stored_prob <- stored_prob[common]
}
rm(res); invisible(gc())
types <- colnames(ll)
named_cols <- types[!grepl(DENOVO_RE, types)]
message(sprintf("logliks: %d cells x %d types (%d named, %d de novo)",
                nrow(ll), ncol(ll), length(named_cols), ncol(ll) - length(named_cols)))
stopifnot("no named columns found" = length(named_cols) >= 2)

# Record the package's own posterior function if the container ships it, so a formula mismatch
# is diagnosable from the log rather than a guess. Not exported, hence the namespace lookup.
f <- tryCatch(get("logliks2probs", envir = asNamespace("InSituType")), error = function(e) NULL)
if (!is.null(f)) { message("---- InSituType:::logliks2probs ----"); print(f); message("----") }

# Posterior of the winner, the runner-up, and the top1-vs-top2 log-likelihood margin, computed
# block by block. Subtracting the row maximum before exp() is what keeps this finite: raw
# log-likelihoods run to the thousands and exp() of them underflows to zero.
posterior_block <- function(block, logprior) {
  if (!is.null(logprior)) block <- sweep(block, 2, logprior, `+`)
  n <- nrow(block)
  i1 <- max.col(block, ties.method = "first")
  top1 <- block[cbind(seq_len(n), i1)]
  masked <- block; masked[cbind(seq_len(n), i1)] <- -Inf
  i2 <- max.col(masked, ties.method = "first")
  top2 <- block[cbind(seq_len(n), i2)]
  lse <- top1 + log(rowSums(exp(block - top1)))
  list(i1 = i1, i2 = i2, prob1 = exp(top1 - lse), prob2 = exp(top2 - lse),
       margin = top1 - top2)
}

run_blocks <- function(cols, logprior) {
  n <- nrow(ll)
  out <- list(i1 = integer(n), i2 = integer(n), prob1 = numeric(n), prob2 = numeric(n),
              margin = numeric(n))
  for (start in seq(1, n, by = CHUNK)) {
    rows <- start:min(start + CHUNK - 1, n)
    r <- posterior_block(ll[rows, cols, drop = FALSE], logprior)
    for (k in names(out)) out[[k]][rows] <- r[[k]]
  }
  out
}

# ---- 1. which formula reproduces the stored prob? ---------------------------
freq <- table(factor(clust, levels = types)) / length(clust)
candidates <- list(
  softmax = NULL,
  softmax_with_frequency_prior = log(pmax(as.numeric(freq), .Machine$double.xmin)))

full_top <- types[max.col(ll, ties.method = "first")]
is_argmax <- full_top == clust
message(sprintf("label is the stored argmax for %.1f%% of cells; validating on those only ",
                100 * mean(is_argmax)),
        "(pinned cells carry an overridden label, so their prob is not comparable)")
check_rows <- which(is_argmax & is.finite(stored_prob))
if (length(check_rows) > 200000) { set.seed(1); check_rows <- sort(sample(check_rows, 200000)) }

check <- rbindlist(lapply(names(candidates), function(nm) {
  r <- posterior_block(ll[check_rows, , drop = FALSE], candidates[[nm]])
  err <- abs(r$prob1 - stored_prob[check_rows])
  data.table(formula = nm, n_cells = length(check_rows),
             mean_abs_error = mean(err), median_abs_error = median(err),
             p99_abs_error = quantile(err, 0.99), max_abs_error = max(err))
}))
# Judge on the MEAN, not the median. The two formulas only disagree on near-tie cells, which are
# a minority, so both medians sit at ~0 whichever formula the fit used; the mean is what those
# cells move.
check[, reproduces_stored_prob := mean_abs_error <= TOL]
fwrite(check, file.path(outdir, "posterior_formula_check.csv"))
print(check)

best <- check[which.min(mean_abs_error)]
VALIDATED <- best$reproduces_stored_prob
FORMULA <- best$formula
if (!VALIDATED) {
  message("\n*** Neither candidate reproduces the stored prob (best mean |error| ",
          sprintf("%.3g", best$mean_abs_error), " for ", FORMULA, "). Continuing with it, but ",
          "every output row is marked formula_validated = FALSE. ***\n")
} else if (all(check$reproduces_stored_prob)) {
  message(sprintf(paste0("stored prob reproduced by BOTH formulas (too few near-ties to tell ",
                         "them apart); using %s, the closer (mean |error| %.2g)"),
                  FORMULA, best$mean_abs_error))
} else {
  message(sprintf("stored prob reproduced by %s only (mean |error| %.2g)",
                  FORMULA, best$mean_abs_error))
}

# ---- 2. the counterfactual: posterior over the named columns only -----------
# The margin is always the LIKELIHOOD-only top1-vs-top2 gap, the quantity the 4.1-vs-25.1
# comparison is about. A prior, if the fit uses one, enters the posterior only.
message(sprintf("%s, named-only posterior over %d columns", Sys.time(), length(named_cols)))
lik <- run_blocks(named_cols, NULL)
nm <- lik
if (FORMULA != "softmax") {
  # The fit's own frequencies are the wrong prior here: a named type that rarely won against the
  # letters would carry a crushing prior, although b's cells would redistribute onto it once the
  # letters are gone. So re-estimate from the named-only calls themselves, with a pseudocount so
  # a type no cell picks is penalised but not excluded.
  counts <- tabulate(lik$i1, nbins = length(named_cols))
  fr <- (counts + 1) / (sum(counts) + length(named_cols))
  message("prior re-estimated from the named-only calls")
  nm <- run_blocks(named_cols, log(fr))
}

cells <- data.table(
  cell_id = rownames(ll), label = clust,
  label_kind = ifelse(grepl(DENOVO_RE, clust), "de novo", "named"),
  stored_prob = as.numeric(stored_prob),
  named_top1 = named_cols[nm$i1], named_top1_prob = nm$prob1,
  named_top2 = named_cols[nm$i2], named_top2_prob = nm$prob2,
  named_margin = lik$margin)
rm(nm, lik); invisible(gc())
fwrite(cells, file.path(outdir, "named_posterior_per_cell.csv.gz"))

# ---- 3. per label: what share would each cut-off flag? ----------------------
by_label <- cells[, c(
  .(label_kind = label_kind[1], n = .N,
    median_named_posterior = median(named_top1_prob),
    median_named_margin = median(named_margin),
    median_stored_prob = median(stored_prob, na.rm = TRUE),
    pct_stored_prob_below_0.8 = round(100 * mean(stored_prob < 0.8, na.rm = TRUE), 2)),
  setNames(lapply(THRESHOLDS, function(t) round(100 * mean(named_top1_prob < t), 2)),
           sprintf("pct_named_posterior_below_%g", THRESHOLDS)),
  setNames(lapply(MARGIN_CUTS, function(m) round(100 * mean(named_margin < m), 2)),
           sprintf("pct_named_margin_below_%g", MARGIN_CUTS))),
  by = label]
by_label[, `:=`(formula = FORMULA, formula_validated = VALIDATED)]
setorder(by_label, -n)
fwrite(by_label, file.path(outdir, "named_posterior_by_label.csv"))

message("\nPer label, largest first (named-only posterior and margin):")
show <- c("label", "label_kind", "n", "median_named_posterior", "median_named_margin",
          sprintf("pct_named_posterior_below_%g", THRESHOLDS[1]),
          sprintf("pct_named_margin_below_%g", MARGIN_CUTS))
print(by_label[, ..show], nrows = 100)

# ---- 4. the focus letter: posterior against margin --------------------------
if (FOCUS %in% cells$label) {
  f_cells <- cells[label == FOCUS]
  f_cells[, posterior_bin := cut(named_top1_prob, POSTERIOR_BREAKS, include.lowest = TRUE,
                                 right = FALSE)]
  f_cells[, margin_bin := cut(named_margin, c(0, MARGIN_CUTS, Inf), include.lowest = TRUE,
                              right = FALSE)]
  grid <- f_cells[, .(n = .N), by = .(posterior_bin, margin_bin)]
  grid[, pct_of_letter := round(100 * n / sum(n), 2)]
  setorder(grid, posterior_bin, margin_bin)
  fwrite(grid, file.path(outdir, sprintf("%s_posterior_by_margin.csv", FOCUS)))
  message(sprintf("\n%s (%d cells): named-only posterior bin x margin bin", FOCUS,
                  nrow(f_cells)))
  print(dcast(grid, posterior_bin ~ margin_bin, value.var = "pct_of_letter", fill = 0))
  confident_but_tied <- f_cells[, mean(named_top1_prob >= THRESHOLDS[1] &
                                         named_margin < MARGIN_CUTS[2])]
  message(sprintf(paste0("\nREAD: %.1f%% of %s clears the %g posterior threshold while its top ",
                         "two named types are under %g log-units apart. Those are the cells a ",
                         "posterior cut alone would leave as confident named calls."),
                  100 * confident_but_tied, FOCUS, THRESHOLDS[1], MARGIN_CUTS[2]))
} else {
  message(sprintf("focus letter %s not among the labels; skipping the grid", FOCUS))
}
message(sprintf("\n%s, wrote outputs to %s", Sys.time(), outdir))
