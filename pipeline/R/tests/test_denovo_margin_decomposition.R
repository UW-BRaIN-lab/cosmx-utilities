#!/usr/bin/env Rscript
# Tests for denovo_margin_decomposition.R Part 1 -- the model-free verdict (no S3 / network).
#
# The pair that matters is the two worlds Part 1 exists to tell apart, planted as synthetic
# loglik tables and run through the real script:
#
#   ONE CONTINUUM  cells drawn from a single population, so the argmax boundary at zero cuts
#                  straight through the peak -> ratio_at_zero must be HIGH.
#   TWO POPULATIONS a genuine gap either side of zero -> ratio_at_zero must be LOW.
#
# Without both, a statistic that always said "arbitrary cut" would be useless, and that is
# exactly the conclusion this script is being used to argue for.
#
# Run:  Rscript pipeline/R/tests/test_denovo_margin_decomposition.R

suppressPackageStartupMessages({ library(data.table) })

# Rscript exposes the test's own path as --file=; the script under test is its parent dir.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT <- file.path(dirname(dirname(normalizePath(this_file))),
                    "denovo_margin_decomposition.R")
stopifnot("cannot locate the script under test" = file.exists(SCRIPT))

NAMED <- c("Endo_capilar", "Endo_arterial", "Pericyte", "Astrocyte")
N <- 4000

# Build a fake insitutype() result whose margin (letter minus destination) follows `margins`.
# Every other named column is pushed well below both so the argmax is always one of the two.
make_rds <- function(path, margins) {
  n <- length(margins)
  ids <- sprintf("SL01_F%d_C%d", seq_len(n) %% 100 + 1, seq_len(n))
  base <- -500 + rnorm(n, 0, 5)
  ll <- matrix(-5000, nrow = n, ncol = length(NAMED) + 1,
               dimnames = list(ids, c(NAMED, "l")))
  ll[, "Endo_capilar"] <- base
  ll[, "l"] <- base + margins
  # clust is the argmax over all columns, which is the letter exactly when the margin is > 0.
  clust <- ifelse(margins > 0, "l", "Endo_capilar")
  names(clust) <- ids
  saveRDS(list(logliks = ll, clust = clust), path)
  invisible(path)
}

run <- function(margins, tag) {
  tmp <- file.path(tempdir(), tag); dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  rds <- file.path(tmp, "typing.rds"); make_rds(rds, margins)
  out <- file.path(tmp, "out")
  res <- system2("Rscript", c(SCRIPT, "--typing-rds", rds, "--letter", "l",
                              "--destination", "Endo_capilar", "--output-dir", out),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(out, "density_at_zero.csv"))) {
    cat(paste(res, collapse = "\n"), "\n"); stop(tag, ": no verdict written")
  }
  fread(file.path(out, "density_at_zero.csv"))[scale == "margin"]
}

set.seed(11)

# --- world 1: one population, boundary slices its peak ------------------------
one <- run(rnorm(N, mean = 0, sd = 40), "continuum")
cat(sprintf("one continuum : ratio_at_zero %.3f  trough %s\n",
            one$ratio_at_zero, one$trough_at_zero))
stopifnot("a single population must leave heavy density at the boundary" =
            one$ratio_at_zero > 0.8)
stopifnot("a single population must not report a trough" = !isTRUE(one$trough_at_zero))

# --- world 2: two populations with a real gap ---------------------------------
gap <- c(rnorm(N / 2, mean = -300, sd = 40), rnorm(N / 2, mean = 300, sd = 40))
two <- run(gap, "twopop")
cat(sprintf("two populations: ratio_at_zero %.3f  trough %s\n",
            two$ratio_at_zero, two$trough_at_zero))
stopifnot("a real gap must empty the density at the boundary" = two$ratio_at_zero < 0.2)

# --- the two must be far apart, not marginally different ----------------------
stopifnot("the statistic must separate the two worlds decisively" =
            one$ratio_at_zero > 4 * two$ratio_at_zero)

# --- a lopsided split (what the real data looks like) still reports honestly ---
# 91% of cells go to the letter, from ONE population -- the boundary is far out in the tail but
# still inside a single continuum, so the density there must not be mistaken for a gap.
lop <- run(rnorm(N, mean = 55, sd = 40), "lopsided")
cat(sprintf("lopsided 1-pop : ratio_at_zero %.3f  trough %s\n",
            lop$ratio_at_zero, lop$trough_at_zero))
stopifnot("a lopsided single population must still show real density at the boundary" =
            lop$ratio_at_zero > 0.3)
stopifnot("a lopsided single population is not a trough" = !isTRUE(lop$trough_at_zero))

cat("all tests passed\n")
