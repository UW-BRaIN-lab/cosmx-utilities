#!/usr/bin/env Rscript
# Tests for named_posterior_threshold.R -- the named-only posterior and the formula check.
#
# Synthetic loglik tables are run through the real script, planting the situations it exists to
# tell apart:
#
#   FORMULA         the stored prob is built WITH a frequency prior, so the check must pick
#                   that candidate over the plain softmax, and flag neither when prob is noise.
#   PINNED CELLS    cells whose label is not their argmax must not enter the validation.
#   NEAR-TIE LETTER a letter whose top two named types sit ~4 log-units apart: its named
#                   posterior is ~0.98, so it clears 0.8 while failing a margin cut of 5 --
#                   the exact case the job exists to measure.
#   CLEAR LETTER    a letter with a wide named margin clears both.
#
# Run:  Rscript pipeline/R/tests/test_named_posterior_threshold.R

suppressPackageStartupMessages({ library(data.table) })

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT <- file.path(dirname(dirname(normalizePath(this_file))), "named_posterior_threshold.R")
stopifnot("cannot locate the script under test" = file.exists(SCRIPT))

NAMED <- c("OPC-like", "AC-like", "Neuron", "Oligodendrocyte")
LETTERS_DN <- c("b", "n")
N_EACH <- 3000

softmax_prob1 <- function(ll, logprior = NULL) {
  if (!is.null(logprior)) ll <- sweep(ll, 2, logprior, `+`)
  top <- apply(ll, 1, max)
  exp(top - (top + log(rowSums(exp(ll - top)))))
}

# b: near-tie between OPC-like and AC-like (gap ~4), and the letter itself wins by a mile.
# n: named Neuron wins by ~40 over the rest, letter n wins overall.
# Named cells: Oligodendrocyte cells whose own column wins cleanly.
make_fit <- function(path, prob_mode = c("prior", "noise"), n_pinned = 0) {
  prob_mode <- match.arg(prob_mode)
  set.seed(42)
  cols <- c(NAMED, LETTERS_DN)
  block <- function(label, fill) {
    m <- matrix(-3000 + rnorm(N_EACH * length(cols), 0, 1), nrow = N_EACH,
                dimnames = list(sprintf("%s_%d", label, seq_len(N_EACH)), cols))
    fill(m)
  }
  b <- block("b", function(m) {
    m[, "OPC-like"] <- -1000 + rnorm(N_EACH, 0, 0.3)
    m[, "AC-like"] <- m[, "OPC-like"] - 4
    m[, "b"] <- -900; m })
  n <- block("n", function(m) {
    m[, "Neuron"] <- -1000 + rnorm(N_EACH, 0, 0.3)
    m[, "AC-like"] <- m[, "Neuron"] - 40
    m[, "n"] <- -900; m })
  o <- block("oligo", function(m) {
    m[, "Oligodendrocyte"] <- -1000 + rnorm(N_EACH, 0, 0.3)
    m[, "Neuron"] <- m[, "Oligodendrocyte"] - 30; m })
  # Near-ties between two full-menu types of very different frequency: the only cells on which
  # a frequency prior changes the posterior, so the only ones that can tell the formulas apart.
  tie <- block("tie", function(m) {
    m <- m[seq_len(600), , drop = FALSE]
    m[, "Oligodendrocyte"] <- -1000 + rnorm(nrow(m), 0, 0.3)
    m[, "AC-like"] <- m[, "Oligodendrocyte"] - 0.5; m })
  ll <- rbind(b, n, o, tie)
  clust <- setNames(cols[max.col(ll, ties.method = "first")], rownames(ll))
  freq <- table(factor(clust, levels = cols)) / length(clust)
  prob <- if (prob_mode == "prior") softmax_prob1(ll, log(pmax(as.numeric(freq), 1e-300)))
          else runif(nrow(ll))
  names(prob) <- rownames(ll)
  # Pinned cells: overwrite the label to something that is NOT their argmax, and give them a
  # prob no formula would reproduce. They must be excluded from the validation.
  if (n_pinned > 0) {
    pinned <- rownames(ll)[seq_len(n_pinned)]
    clust[pinned] <- "Neuron"; prob[pinned] <- 0.01
  }
  saveRDS(list(logliks = ll, clust = clust, prob = prob), path)
}

run <- function(prob_mode, n_pinned = 0) {
  dir <- tempfile("npt_"); dir.create(dir)
  rds <- file.path(dir, "fit.rds")
  make_fit(rds, prob_mode, n_pinned)
  out <- file.path(dir, "out")
  status <- system2("Rscript", c(SCRIPT, "--typing-rds", rds, "--output-dir", out,
                                 "--chunk", "1000"), stdout = TRUE, stderr = TRUE)
  cat(paste(status, collapse = "\n"), "\n")
  list(check = fread(file.path(out, "posterior_formula_check.csv")),
       by_label = fread(file.path(out, "named_posterior_by_label.csv")),
       grid = fread(file.path(out, "b_posterior_by_margin.csv")),
       cells = fread(file.path(out, "named_posterior_per_cell.csv.gz")))
}

failures <- 0
expect <- function(ok, what) {
  cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", what))
  if (!isTRUE(ok)) failures <<- failures + 1
}

cat("\n===== stored prob built with a frequency prior, 200 pinned cells =====\n")
r <- run("prior", n_pinned = 200)
print(r$check)
chosen <- r$check[which.min(mean_abs_error)]
expect(chosen$formula == "softmax_with_frequency_prior",
       "formula check picks the frequency-prior candidate")
expect(isTRUE(chosen$reproduces_stored_prob),
       "chosen formula reproduces the stored prob despite 200 pinned cells")
expect(all(r$by_label$formula_validated), "every output row is marked validated")

b_row <- r$by_label[label == "b"]
print(b_row)
expect(b_row$median_named_margin > 3.5 && b_row$median_named_margin < 4.5,
       sprintf("b's named margin is the planted ~4 (got %.2f)", b_row$median_named_margin))
expect(b_row$pct_named_posterior_below_0.8 < 1,
       sprintf("b clears the 0.8 posterior cut (%.2f%% below)", b_row$pct_named_posterior_below_0.8))
expect(b_row$pct_named_margin_below_5 > 99,
       sprintf("...but fails a margin cut of 5 (%.2f%% below)", b_row$pct_named_margin_below_5))

n_row <- r$by_label[label == "n"]
expect(n_row$pct_named_margin_below_10 < 1 && n_row$pct_named_posterior_below_0.99 < 1,
       "the clear letter n clears both the posterior and the margin cuts")
expect(all(r$cells[label == "n", named_top1] == "Neuron"),
       "n's named-only call is Neuron, the planted runner-up behind the letter")
expect(all(r$cells[label == "b", named_top1] %in% c("OPC-like", "AC-like")),
       "b's named-only call is one of its two planted near-tie types")

grid <- r$grid
expect(abs(sum(grid$pct_of_letter) - 100) < 0.1, "b's posterior x margin grid sums to 100%")

cat("\n===== stored prob is noise: nothing should validate =====\n")
r2 <- run("noise")
print(r2$check)
expect(!any(r2$check$reproduces_stored_prob), "neither formula validates against noise")
expect(!any(r2$by_label$formula_validated), "outputs are marked formula_validated = FALSE")

cat(sprintf("\n%d failure(s)\n", failures))
quit(status = if (failures) 1 else 0)
