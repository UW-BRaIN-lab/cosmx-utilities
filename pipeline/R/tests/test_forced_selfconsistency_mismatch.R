#!/usr/bin/env Rscript
# Tests for 75c's label-vs-argmax mismatch list.
#
# The one that matters is ALIGNMENT. dt is merged with the re-scored posteriors by cell_id,
# and data.table's merge re-sorts by that key — so a per-cell loglik looked up after the merge
# pairs each cell with a DIFFERENT cell's likelihoods, silently. The fixture therefore gives
# cells ids whose sorted order is not their row order, and checks the deficit of named cells,
# which only comes out right if the lookup happened before the merge.
#
#   Rscript pipeline/R/tests/test_forced_selfconsistency_mismatch.R

suppressMessages(library(data.table))
# Rscript exposes the test's own path as --file=; the script under test is its parent dir.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT <- file.path(dirname(dirname(normalizePath(this_file))),
                    "diagnose_forced_selfconsistency.R")
stopifnot("cannot locate the script under test" = file.exists(SCRIPT))

set.seed(23)
N <- 300
NAMED <- c("Endo_capilar", "Endo_arterial", "Pericyte")
TYPES <- c(NAMED, "l", "c")

tmp <- file.path(tempdir(), "mismatch"); dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

# Ids that sort into a different order than the rows they sit on.
ids <- sprintf("CELL_%s", sample(sprintf("%04d", seq_len(N))))
ll <- matrix(rnorm(N * length(TYPES), -600, 30), nrow = N,
             dimnames = list(ids, TYPES))

# Everyone is assigned their argmax...
clust <- TYPES[max.col(ll, ties.method = "first")]
names(clust) <- ids

# ...except 40 cells, deliberately relabelled to a type that scores WORSE by a known amount.
mismatched <- sort(sample(seq_len(N), 40))
wrong <- character(length(mismatched))
planted <- numeric(length(mismatched))
for (k in seq_along(mismatched)) {
  i <- mismatched[[k]]
  best <- which.max(ll[i, ])
  alt <- setdiff(seq_along(TYPES), best)[1]
  ll[i, alt] <- ll[i, best] - (10 * k)   # a deficit unique to this cell
  wrong[[k]] <- TYPES[[alt]]
  planted[[k]] <- 10 * k
}
clust[mismatched] <- wrong

saveRDS(list(logliks = ll, clust = clust), file.path(tmp, "typing.rds"))
# The re-scored posteriors the script joins on; deliberately written in yet another order.
fwrite(data.table(cell_id = sample(ids), top1_type = NAMED[[1]]),
       file.path(tmp, "posteriors.csv"))

out_csv <- file.path(tmp, "mismatch.csv")
res <- system2("Rscript", c(SCRIPT,
                            "--typing-rds", file.path(tmp, "typing.rds"),
                            "--posteriors", file.path(tmp, "posteriors.csv"),
                            "--output-csv", file.path(tmp, "by_label.csv"),
                            "--output-mismatch-csv", out_csv),
               stdout = TRUE, stderr = TRUE)
if (!file.exists(out_csv)) { cat(paste(res, collapse = "\n"), "\n"); stop("no mismatch csv") }

got <- fread(out_csv)
expected <- data.table(cell_id = ids[mismatched], assigned = wrong, deficit = planted)

stopifnot("every planted mismatch must be listed, and nothing else" =
            setequal(got$cell_id, expected$cell_id) && nrow(got) == length(mismatched))
cat(sprintf("listed %d of %d planted mismatches\n", nrow(got), length(mismatched)))

m <- merge(got, expected, by = "cell_id", suffixes = c("", "_want"))
stopifnot("the assigned label must be the one the cell actually carries" =
            all(m$assigned == m$assigned_want))
# THE ALIGNMENT CHECK: a deficit read off the wrong row would not reproduce the planted value.
stopifnot("the deficit must match the planted gap" =
            max(abs(m$deficit - m$deficit_want)) < 1e-6)
cat(sprintf("deficits reproduce exactly (max error %.2e)\n",
            max(abs(m$deficit - m$deficit_want))))

stopifnot("argmax must never equal the assigned label" = all(got$argmax != got$assigned))
stopifnot("deficit must be positive by construction" = all(got$deficit > 0))
stopifnot("the list must be sorted worst-first" = !is.unsorted(rev(got$deficit)))
cat("all tests passed\n")
