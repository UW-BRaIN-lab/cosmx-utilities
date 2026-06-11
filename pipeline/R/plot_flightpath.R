#!/usr/bin/env Rscript
# Stage 4 QC: InSituType "flightpath" plot from the saved typing result.
#
# The flightpath lays every cell out by its posterior log-likelihoods across clusters:
# confident cells sit at their cluster's vertex, ambiguous cells pile toward the middle,
# and each cluster is labelled with its mean assignment confidence. It's the canonical
# InSituType QC for "are these calls trustworthy and are the clusters distinct?" — see
# https://github.com/Nanostring-Biostats/InSituType/blob/main/FAQs.md
#
# Reads the full result saved by insitutype_typing.R (stage 4b):
#   --input   insitutype_result.rds  (list with $logliks cells x clusters + $profiles)
# Writes:
#   --output  a PNG of the flightpath.
#
# flightpath_plot() only consumes $logliks and $profiles, so this works off the .rds
# alone. At cohort scale the full cell count is too dense to render, so a random subset
# is plotted (--n-cells, default 1e5); the layout is unchanged, just fewer points drawn.
#
# Run in the InSituType container (APPTAINER_INSITUTYPE), e.g.:
#   apptainer exec "$APPTAINER_INSITUTYPE" Rscript pipeline/R/plot_flightpath.R \
#       --input insitutype_result.rds --output flightpath.png

suppressPackageStartupMessages({
  library(InSituType)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument: %s", a))
    key <- sub("^--", "", a)
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--output is required" = !is.null(opt$output))
n_cells <- as.integer(opt[["n-cells"]] %||% 100000)
seed <- as.integer(opt$seed %||% 1)

message(sprintf("%s, reading %s", Sys.time(), opt$input))
res <- readRDS(opt$input)
stopifnot("result is missing $logliks" = !is.null(res$logliks))
n <- nrow(res$logliks)

# Subsample cells for a renderable plot; the cluster layout is driven by $profiles, so
# the vertices are unchanged — we're just drawing fewer points.
set.seed(seed)
idx <- if (n > n_cells) sort(sample.int(n, n_cells)) else seq_len(n)
sub <- list(logliks = res$logliks[idx, , drop = FALSE], profiles = res$profiles)
message(sprintf("flightpath on %d / %d cells across %d clusters",
                length(idx), n, ncol(res$logliks)))

p <- InSituType::flightpath_plot(insitutype_result = sub, showclusterconfidence = TRUE)

ggplot2::ggsave(opt$output, p, width = 12, height = 10, dpi = 150)
message(sprintf("%s, wrote %s", Sys.time(), opt$output))
