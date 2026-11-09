#!/usr/bin/env Rscript
# Test for reference_calibration.R.
#
# The fixture encodes the failure the script exists to find: a reference type that attracted
# no pinned cell and whose fitted profile collapsed to a fraction of the de novo scale. Such a
# type is INVISIBLE in the typing output -- it wins nothing, so it appears nowhere -- and the
# table must name it anyway.
#
#   Rscript pipeline/R/tests/test_reference_calibration.R

suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT <- file.path(dirname(dirname(normalizePath(this_file))), "reference_calibration.R")
stopifnot("cannot locate the script under test" = file.exists(SCRIPT))

set.seed(11)
N_CELLS <- 500
N_GENES <- 40
GENES <- sprintf("GENE%02d", seq_len(N_GENES))
# calibrated: pinned cells, profile at the de novo scale. starved: no pinned cell, profile
# collapsed. absent: in the reference, never fitted at all.
REF_TYPES <- c("Calibrated_glia", "Starved_neuron", "Absent_type")
DENOVO <- c("a", "b")

reference <- matrix(abs(rnorm(N_GENES * length(REF_TYPES), 8, 3)),
                    N_GENES, dimnames = list(GENES, REF_TYPES))

# The fit carries two of the three reference types plus the de novo columns.
fitted_types <- c("Calibrated_glia", "Starved_neuron", DENOVO)
fitted <- cbind(
  Calibrated_glia = reference[, "Calibrated_glia"] * 1.05,   # rescaled, same scale
  Starved_neuron  = reference[, "Starved_neuron"] * 0.01,    # collapsed: 1% of reference
  a               = abs(rnorm(N_GENES, 8, 3)),
  b               = abs(rnorm(N_GENES, 8, 3)))
rownames(fitted) <- GENES
stopifnot(identical(colnames(fitted), fitted_types))

ids <- sprintf("CELL_%04d", seq_len(N_CELLS))
clust <- sample(c("Calibrated_glia", DENOVO), N_CELLS, replace = TRUE,
                prob = c(0.3, 0.4, 0.3))
names(clust) <- ids
anchors <- rep(NA_character_, N_CELLS); names(anchors) <- ids
pinned_rows <- sample(seq_len(N_CELLS), 40)
anchors[pinned_rows] <- "Calibrated_glia"
clust[pinned_rows] <- "Calibrated_glia"

tmp <- file.path(tempdir(), "refcal")
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
ref_csv <- file.path(tmp, "reference.csv")
fwrite(data.table(gene = GENES, as.data.table(reference)), ref_csv)
rds <- file.path(tmp, "fit.rds")
saveRDS(list(clust = clust, profiles = fitted, anchors = anchors), rds)

out_csv <- file.path(tmp, "reference_calibration.csv")
log <- system2("Rscript", c(SCRIPT, "--result-rds", rds, "--reference", ref_csv,
                            "--output-csv", out_csv), stdout = TRUE, stderr = TRUE)
if (!file.exists(out_csv)) { cat(paste(log, collapse = "\n"), "\n"); stop("no output csv") }
res <- fread(out_csv)
setkey(res, type)

cat(paste(log, collapse = "\n"), "\n\n--- assertions ---\n")

stopifnot("every reference type and de novo column must appear" =
            setequal(res$type, union(REF_TYPES, fitted_types)))
stopifnot("the pinned type must carry its pinned count" =
            res["Calibrated_glia", n_pinned] == 40L)
stopifnot("a type with no pinned cell must read 0, not NA" =
            res["Starved_neuron", n_pinned] == 0L)
stopifnot("an unfitted reference type must read 0 pinned and 0 called" =
            res["Absent_type", n_pinned] == 0L && res["Absent_type", n_called] == 0L)
stopifnot("de novo columns must be labelled as such" =
            all(res[DENOVO, source] == "denovo"))
stopifnot("reference columns must be labelled as such" =
            all(res[REF_TYPES, source] == "reference"))

# The collapsed profile is the finding: ~1% of its reference scale, far under the de novo
# median, while the calibrated one sits beside the de novo columns.
starved <- res["Starved_neuron", profile_sum_vs_denovo]
calibrated <- res["Calibrated_glia", profile_sum_vs_denovo]
cat(sprintf("starved vs de novo: %.4f   calibrated vs de novo: %.3f\n", starved, calibrated))
stopifnot("the collapsed profile must fall below a tenth of the de novo median" = starved < 0.1)
stopifnot("the calibrated profile must not" = calibrated > 0.5)
stopifnot("fitted/reference must recover the 1% collapse" =
            abs(res["Starved_neuron", fitted_over_reference] - 0.01) < 1e-6)
# Scaling a profile does not rotate it, so the log1p cosine stays high even when collapsed --
# which is exactly why the SUM is the diagnostic and the cosine is not.
stopifnot("cosine to the reference must be reported for fitted reference types" =
            all(!is.na(res[c("Calibrated_glia", "Starved_neuron"), cosine_to_reference])))
stopifnot("the script must name the starved type in its report" =
            any(grepl("Starved_neuron", log)))
stopifnot("the headline must count the calibrated types" =
            any(grepl("1 of 3 reference types received pinned cells", log)))

cat("\nall reference_calibration tests passed\n")
