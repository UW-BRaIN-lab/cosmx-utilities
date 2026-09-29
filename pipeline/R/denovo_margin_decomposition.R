#!/usr/bin/env Rscript
# Why did a cell choose a de-novo letter over a named GBmap type — and is the boundary between
# them a real gap or an arbitrary cut through one population?
#
# 75d compares the two groups on markers and its header records the catch: a cell landed in the
# letter BECAUSE it fit the letter's profile better, so "they differ" is definitional. 75h then
# showed the split is not donor-driven and not batch-driven: ~91% of vascular cells go to `l` in
# every one of the 38 donors. So the remaining question is about the decision itself.
#
# The assignment is argmax over the 81 stored log-likelihoods, so for a letter L and a
# destination D the quantity that decided every cell is
#
#     margin_i = loglik(i, L) - loglik(i, D)
#
# and A cells (clust == L) have margin > 0 while B cells (clust == D) have margin <= 0 BY
# CONSTRUCTION. That is why "do A and B separate" is not a question worth asking here: they are
# separated at zero by definition. What is NOT definitional is where zero falls in the DENSITY.
#
#   * If the pooled margin density is unimodal with substantial mass AT zero, the boundary
#     slices through the bulk of one continuous population -- B is just the tail that happened
#     to fall inside GBmap's basin, and `l` is not a distinct cell type.
#   * If the density has a trough at zero with a mode either side, the boundary sits in a
#     natural gap and the two are genuinely different populations.
#
# This answers that from the stored logliks alone -- no model, no counts, nothing to get wrong.
# (A per-gene decomposition of the margin, which needs InSituType's likelihood recomputed, was
# tried and never reproduced the stored margin; it is not shipped. See the tag
# pr79-archive-09f2397 for the attempts.)
#
# Inputs:
#   --typing-rds  anchor_typing.rds from 72 (clust, logliks, profiles). NOT the h5, which drops
#                 $logliks.
#   --counts-h5   anchor_input.h5 (stage 4a), OPTIONAL: adds the per-count margin (depth-
#                 normalised) beside the raw one. Omit and depth is NA.
#   --letter      de-novo letter, e.g. l
#   --destination named GBmap type to compare against, e.g. Endo_capilar
# Outputs (--output-dir):
#   margin_summary.csv     per-group margin stats, raw and per-count
#   margin_histogram.csv   binned pooled density, for plotting on the Mac
#   density_at_zero.csv    the verdict statistic: density at the boundary vs the modes

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(hdf5r)
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
stopifnot("--letter is required" = !is.null(opt$letter))
stopifnot("--destination is required" = !is.null(opt$destination))
stopifnot("--output-dir is required" = !is.null(opt[["output-dir"]]))
outdir <- opt[["output-dir"]]; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

LETTER <- opt$letter
DEST   <- opt$destination
# Same convention as 75c: de-novo clusters are 1-2 lowercase letters (K=27 overflows to "aa").
DENOVO_RE <- "^[a-z]{1,2}$"

message(sprintf("%s, reading %s", Sys.time(), opt[["typing-rds"]]))
res <- readRDS(opt[["typing-rds"]])
if (is.null(res$logliks)) stop("anchor_typing.rds has no $logliks -- this must be the full ",
                               "insitutype() result, not the compact h5.")
ll <- as.matrix(res$logliks)
clust <- res$clust
if (!is.null(rownames(ll)) && !is.null(names(clust))) {
  common <- intersect(rownames(ll), names(clust))
  ll <- ll[common, , drop = FALSE]; clust <- clust[common]
}
types <- colnames(ll)
stopifnot("letter not among the stored profile columns" = LETTER %in% types)
stopifnot("destination not among the stored profile columns" = DEST %in% types)
message(sprintf("logliks: %d cells x %d types", nrow(ll), ncol(ll)))

# The forced call, read off the fit's own stored table exactly as 75c does -- never re-scored.
named_cols <- types[!grepl(DENOVO_RE, types)]
forced <- named_cols[max.col(ll[, named_cols, drop = FALSE], ties.method = "first")]

in_a <- clust == LETTER & forced == DEST     # the letter's cells GBmap would call DEST
in_b <- clust == DEST                        # cells that took DEST without needing a letter
message(sprintf("A (%s->%s): %d cells    B (%s native): %d cells",
                LETTER, DEST, sum(in_a), DEST, sum(in_b)))
stopifnot("group A is empty" = sum(in_a) > 0, "group B is empty" = sum(in_b) > 0)

# ---- the margin, and where zero sits in its density --------------------------
keep <- in_a | in_b
margin <- ll[keep, LETTER] - ll[keep, DEST]
group  <- ifelse(in_a[keep], sprintf("A: %s->%s", LETTER, DEST), sprintf("B: %s native", DEST))

# The raw margin scales with sequencing depth -- a deep cell accumulates a larger difference at
# the same per-transcript preference -- so the per-count margin is reported beside it.
depth <- if (!is.null(res$profiles) && !is.null(attr(res, "totalcounts"))) {
  attr(res, "totalcounts")[keep]
} else rep(NA_real_, sum(keep))
if (!is.null(opt[["counts-h5"]])) {
  fh <- H5File$new(opt[["counts-h5"]], mode = "r")
  cid <- fh[["cell_id"]][]
  shp <- fh[["counts/shape"]][]
  cnt <- methods::new("dgCMatrix", i = as.integer(fh[["counts/indices"]][]),
                      p = as.integer(fh[["counts/indptr"]][]),
                      x = as.double(fh[["counts/data"]][]),
                      Dim = as.integer(shp), Dimnames = list(fh[["genes"]][], cid))
  neg_all <- fh[["neg"]][]; names(neg_all) <- cid
  fh$close_all()
  tot <- Matrix::colSums(cnt)
  depth <- tot[rownames(ll)[keep]]
  message(sprintf("counts: %d genes x %d cells; median depth %.0f",
                  nrow(cnt), ncol(cnt), median(depth, na.rm = TRUE)))
}

# Is clust actually the argmax of the stored logliks for these cells? 75c measured 98.9%
# COHORT-WIDE, but that is dominated by large clusters; a small named type competing against a
# large de-novo one is exactly where the exceptions would hide. If clust is not the argmax then
# the margin does not explain the assignment for that group, and the boundary is not at zero.
top_all <- types[max.col(ll[keep, , drop = FALSE], ties.method = "first")]
is_argmax <- top_all == clust[keep]
argmax_tbl <- data.table(group = group, is_argmax = is_argmax)[
  , .(n = .N, pct_clust_is_argmax = round(100 * mean(is_argmax), 1)), by = group]
fwrite(argmax_tbl, file.path(outdir, "argmax_consistency.csv"))
message("\nIs the stored argmax the assigned label?")
print(argmax_tbl)
message("Anything well below 100% means insitutype did not assign that group by the stored\n",
        "logliks, so for those cells the margin is not the reason they went where they did.\n")

dt <- data.table(cell_id = rownames(ll)[keep], group = group, margin = margin,
                 depth = as.numeric(depth))
dt[, margin_per_count := margin / pmax(depth, 1)]

# Enough quantiles to draw a box or a violin from the summary alone, not just a median and a
# 10/90 pair -- the PI's lead figure is the per-group DISTRIBUTION, and a 3-number summary
# cannot show whether a group is tight against the boundary or spread across it.
qs <- c(.01, .05, .1, .25, .5, .75, .9, .95, .99)
summ <- dt[, c(.(n = .N, median_margin = median(margin), q10 = quantile(margin, .1),
                 q90 = quantile(margin, .9), median_depth = median(depth, na.rm = TRUE),
                 median_margin_per_count = median(margin_per_count, na.rm = TRUE),
                 pct_above_zero = round(100 * mean(margin > 0), 1)),
               setNames(as.list(quantile(margin, qs)),
                        sprintf("margin_q%02d", round(qs * 100)))), by = group]
fwrite(summ, file.path(outdir, "margin_summary.csv"))
print(summ)

# Density of the pooled margin. The boundary is at zero by construction; the question is whether
# zero lands in a trough (two populations) or in the bulk (one continuum, arbitrarily cut).
verdicts <- list()
group_curves <- list()
for (scale_name in c("margin", "margin_per_count")) {
  sub <- dt[is.finite(get(scale_name))]
  v <- sub[[scale_name]]
  if (length(v) < 100) next
  # Trim the extreme tails so the bandwidth is not set by outliers.
  lim <- quantile(v, c(0.001, 0.999))
  inside <- v >= lim[1] & v <= lim[2]
  d <- density(v[inside], n = 2048)
  at0 <- approx(d$x, d$y, xout = 0)$y
  peak <- max(d$y)
  # A trough at the boundary means a genuine gap. `local_min_x` locates it, but do NOT decide
  # on "is the argmin near zero": when the gap is CLEAN the density there is flat at ~0, so
  # which.min picks an arbitrary point inside the flat region and the test fails on exactly the
  # clearest cases (a planted two-mode fixture with a 140-unit gap came out FALSE). Ask the
  # structural question instead -- is there a real mode on EACH side of zero, and is the
  # density at zero far below both of them?
  inner <- which(d$x > quantile(v, .05) & d$x < quantile(v, .95))
  local_min <- if (length(inner) > 10) d$x[inner][which.min(d$y[inner])] else NA_real_
  left <- d$y[d$x < 0]; right <- d$y[d$x > 0]
  flanked <- length(left) > 10 && length(right) > 10 &&
    max(left) > 0.1 * peak && max(right) > 0.1 * peak
  verdicts[[scale_name]] <- data.table(
    scale = scale_name, density_at_zero = at0, peak_density = peak,
    ratio_at_zero = at0 / peak, mode_x = d$x[which.max(d$y)],
    local_min_x = local_min,
    trough_at_zero = isTRUE(flanked && at0 / min(max(left), max(right)) < 0.5))
  fwrite(data.table(scale = scale_name, x = d$x, density = d$y),
         file.path(outdir, sprintf("margin_histogram_%s.csv", scale_name)))

  # The pooled curve above answers "where does zero fall". It cannot answer the PI's question,
  # which is how the FORCED and NATIVE groups sit relative to each other -- so also emit one
  # curve per group, on ONE SHARED BANDWIDTH AND GRID. density() chooses its own bandwidth per
  # vector, so a 3,456-cell group and a 61,721-cell one would be smoothed differently and the
  # reader would be comparing the smoothing rather than the cells.
  for (g in unique(sub$group)) {
    vg <- v[inside & sub$group == g]
    if (length(vg) < 20) next
    dg <- density(vg, bw = d$bw, from = lim[[1]], to = lim[[2]], n = 512)
    group_curves[[paste(scale_name, g)]] <- data.table(
      scale = scale_name, group = g, n = length(vg), x = dg$x, density = dg$y)
  }
}
if (length(group_curves)) {
  fwrite(rbindlist(group_curves), file.path(outdir, "margin_density_by_group.csv"))
}
verdict <- rbindlist(verdicts)
fwrite(verdict, file.path(outdir, "density_at_zero.csv"))
print(verdict)
message("\nREAD: ratio_at_zero near 1 = the boundary cuts through the PEAK of one population,\n",
        "so the letter and the named type are one continuum split arbitrarily. A ratio near 0\n",
        "with a trough at zero = a genuine gap, i.e. two populations.\n")
