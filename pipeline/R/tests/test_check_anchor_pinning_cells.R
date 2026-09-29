#!/usr/bin/env Rscript
# Test for check_anchor_pinning.R's optional pinned-cells export (third argument).
#
# The export is what denovo_vs_native_pseudobulk.py --pinned-csv reads, so it must list EXACTLY
# the cells whose $anchors entry is non-NA, with the label they carry, and nothing else.
#
#   Rscript pipeline/R/tests/test_check_anchor_pinning_cells.R

suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT <- file.path(dirname(dirname(normalizePath(this_file))), "check_anchor_pinning.R")
stopifnot("cannot locate the script under test" = file.exists(SCRIPT))

set.seed(5)
N <- 400
TYPES <- c("Endo_arterial", "Pericyte", "l")
ids <- sprintf("CELL_%s", sample(sprintf("%04d", seq_len(N))))   # sorted order != row order
ll <- matrix(rnorm(N * 3, -500, 20), N, dimnames = list(ids, TYPES))
clust <- TYPES[max.col(ll, ties.method = "first")]; names(clust) <- ids

pinned_rows <- sort(sample(seq_len(N), 30))
anchors <- rep(NA_character_, N); names(anchors) <- ids
anchors[pinned_rows] <- "Endo_arterial"
clust[pinned_rows] <- "Endo_arterial"       # the overwrite: label ignores the likelihood

tmp <- file.path(tempdir(), "pin"); dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(clust = clust, logliks = ll, anchors = anchors), file.path(tmp, "fit.rds"))
res <- system2("Rscript", c(SCRIPT, file.path(tmp, "fit.rds"), file.path(tmp, "by_label.csv"),
                            file.path(tmp, "pinned.csv")), stdout = TRUE, stderr = TRUE)
if (!file.exists(file.path(tmp, "pinned.csv"))) { cat(paste(res, collapse = "\n")); stop("no export") }

got <- fread(file.path(tmp, "pinned.csv"))
stopifnot("columns" = identical(names(got), c("cell_id", "pinned_type", "label")))
stopifnot("exactly the pinned cells" = setequal(got$cell_id, ids[pinned_rows]) && nrow(got) == 30)
stopifnot("pinned type is carried" = all(got$pinned_type == "Endo_arterial"))
stopifnot("label is the assigned label" = all(got$label == clust[got$cell_id]))
cat(sprintf("exported %d pinned cells, all matching\n", nrow(got)))

# Without the third argument nothing is written and the script still succeeds.
res2 <- system2("Rscript", c(SCRIPT, file.path(tmp, "fit.rds"), file.path(tmp, "by_label2.csv")),
                stdout = TRUE, stderr = TRUE)
stopifnot("still works without the export" = file.exists(file.path(tmp, "by_label2.csv")))
cat("all tests passed\n")
