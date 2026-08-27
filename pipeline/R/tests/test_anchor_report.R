#!/usr/bin/env Rscript
# Unit tests for anchor_report.R's reporting logic.
#
# The InSituType calls in anchor_report.R need the pinned package (only present in
# insitutype.sif), but the part most likely to be wrong is the bookkeeping AROUND them:
# which types count as anchored, which type absorbed an unanchored type's cells, and the
# profile-collinearity lookup. Those are pure functions of the cos/llr matrices and the
# anchor vector, so they are tested here against a synthetic fixture with a known answer —
# no InSituType, no HDF5, no cluster.
#
# Run:  Rscript pipeline/R/tests/test_anchor_report.R

suppressPackageStartupMessages(library(data.table))

SCRIPT <- file.path(dirname(dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))),
  "anchor_report.R")

# Load ONLY the function definitions: evaluating the whole script would run its argument
# parsing and demand InSituType. An assignment counts as a definition when its right-hand
# side is literally a `function(...)` expression, so nothing else can execute.
load_functions <- function(path) {
  env <- new.env(parent = globalenv())
  for (expr in parse(path)) {
    if (!is.call(expr)) next
    if (!as.character(expr[[1]])[[1]] %in% c("<-", "=")) next
    rhs <- expr[[3]]
    if (is.call(rhs) && identical(as.character(rhs[[1]])[[1]], "function")) {
      eval(expr, env)
    }
  }
  env
}

env <- load_functions(SCRIPT)
for (needed in c("profile_cosines", "build_report")) {
  if (!exists(needed, envir = env, inherits = FALSE)) {
    stop(sprintf("FAIL: %s not found in %s", needed, SCRIPT))
  }
}

failures <- 0L
check <- function(label, actual, expected, tolerance = 0) {
  ok <- if (is.numeric(expected) && is.numeric(actual) && tolerance > 0) {
    !is.na(actual) && abs(actual - expected) <= tolerance
  } else {
    identical(actual, expected)
  }
  if (!ok) {
    failures <<- failures + 1L
    cat(sprintf("FAIL  %-52s got %s, expected %s\n", label,
                paste(format(actual), collapse = ","),
                paste(format(expected), collapse = ",")))
  } else {
    cat(sprintf("ok    %s\n", label))
  }
}

# --- fixture -----------------------------------------------------------------
# Astro and Astro_reactive are deliberately near-collinear: the situation that makes a
# type unanchorable. Astro always edges it out on cosine, so Astro_reactive should come
# back unanchored with Astro named as its empirical absorber.
GENES <- paste0("G", 1:6)
TYPES <- c("Astro", "Astro_reactive", "Neuron", "Micro", "Rare")
ref <- matrix(
  c(10, 9, 1, 1, 1, 1,          # Astro
    10.2, 8.8, 1.1, 1, 1, 1,    # Astro_reactive — nearly the same direction
    1, 1, 10, 9, 1, 1,          # Neuron
    1, 1, 1, 1, 10, 9,          # Micro
    9, 1, 1, 10, 1, 1),         # Rare
  nrow = length(GENES), dimnames = list(GENES, TYPES)
)

N_CELLS <- 100L
MIN_COSINE <- 0.3
cos_mat <- matrix(0.05, nrow = N_CELLS, ncol = length(TYPES),
                  dimnames = list(paste0("c", seq_len(N_CELLS)), TYPES))
cos_mat[1:70, "Astro"] <- 0.90            # Astro wins every cell it shares with reactive
cos_mat[1:70, "Astro_reactive"] <- 0.85   # plausible on 70 cells, best on none
cos_mat[71:90, "Neuron"] <- 0.80
cos_mat[91:95, "Micro"] <- 0.70
# "Rare" is left at 0.05 everywhere: no cell ever clears min_cosine for it.

llr_mat <- matrix(0.5, nrow = N_CELLS, ncol = length(TYPES),
                  dimnames = dimnames(cos_mat))
llr_mat[1:70, "Astro_reactive"] <- 0.001  # loses the likelihood-ratio tie-break too

# Anchors as choose_anchors_from_stats would return them: types at or below the
# insufficient-anchors threshold have already had every anchor stripped to NA.
anchors <- rep(NA_character_, N_CELLS)
names(anchors) <- rownames(cos_mat)
anchors[1:40] <- "Astro"
anchors[71:90] <- "Neuron"

report <- env$build_report(ref, anchors, cos_mat, llr_mat, MIN_COSINE, 20L)
row_for <- function(type) report[cell_type == type]

# --- anchored vs unanchored --------------------------------------------------
check("row per reference type", nrow(report), length(TYPES))
check("Astro anchored", row_for("Astro")$anchored, TRUE)
check("Astro anchor count", row_for("Astro")$n_anchors, 40L)
check("Neuron anchored", row_for("Neuron")$anchored, TRUE)
check("Neuron anchor count", row_for("Neuron")$n_anchors, 20L)
check("Astro_reactive unanchored", row_for("Astro_reactive")$anchored, FALSE)
check("Micro unanchored", row_for("Micro")$anchored, FALSE)
check("Rare unanchored", row_for("Rare")$anchored, FALSE)
check("unanchored types counted", nrow(report[anchored == FALSE]), 3L)

# --- the "why": plausible cells, who won them, and collinearity --------------
check("Astro_reactive cells above cut", row_for("Astro_reactive")$n_cells_above_min_cosine, 70L)
check("Astro_reactive best on cosine for none", row_for("Astro_reactive")$n_cells_best_on_cosine, 0L)
check("Astro_reactive absorbed by Astro", row_for("Astro_reactive")$empirical_absorber, "Astro")
check("Astro best on cosine for its 70", row_for("Astro")$n_cells_best_on_cosine, 70L)
check("Astro_reactive max cosine", row_for("Astro_reactive")$max_cosine, 0.85, tolerance = 1e-12)
check("Astro_reactive median llr", row_for("Astro_reactive")$median_llr_above_cosine,
      0.001, tolerance = 1e-12)

# Cells 96-100 clear min_cosine for nothing, so no type may claim them as its best match
# (an earlier version credited them to whichever column max.col's tie-break reached first).
check("only cells clearing the cut name a best type",
      sum(report$n_cells_best_on_cosine), 95L)

# A type no cell resembles has no absorber and no llr to report, rather than a wrong one.
check("Rare has no cells above cut", row_for("Rare")$n_cells_above_min_cosine, 0L)
check("Rare has no absorber", row_for("Rare")$empirical_absorber, NA_character_)
check("Rare has no median llr", row_for("Rare")$median_llr_above_cosine, NA_real_)

# Nearest reference type must be the collinear partner, not the type itself.
check("Astro nearest is Astro_reactive", row_for("Astro")$nearest_reference_type, "Astro_reactive")
check("Astro_reactive nearest is Astro", row_for("Astro_reactive")$nearest_reference_type, "Astro")
check("collinear pair cosine > 0.99", row_for("Astro")$nearest_reference_cosine > 0.99, TRUE)
check("Neuron nearest is not itself", row_for("Neuron")$nearest_reference_type != "Neuron", TRUE)

# --- profile_cosines directly ------------------------------------------------
pc <- env$profile_cosines(ref)
check("profile cosine matrix is types x types", dim(pc), c(5L, 5L))
check("self-cosine is 1", pc["Neuron", "Neuron"], 1, tolerance = 1e-12)
check("cosine is symmetric", pc["Astro", "Neuron"], pc["Neuron", "Astro"], tolerance = 1e-12)
check("collinear pair near 1", pc["Astro", "Astro_reactive"] > 0.999, TRUE)

# A zero profile must not produce NaN cosines that silently read as "not similar".
ref_zero <- cbind(ref, Empty = rep(0, length(GENES)))
pc_zero <- env$profile_cosines(ref_zero)
check("zero profile yields NA, not NaN", is.na(pc_zero["Astro", "Empty"]), TRUE)

cat(sprintf("\n%s: %d checks failed\n", if (failures == 0L) "PASS" else "FAIL", failures))
quit(status = if (failures == 0L) 0L else 1L)
