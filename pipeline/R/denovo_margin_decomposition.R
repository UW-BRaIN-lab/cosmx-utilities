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
# PART 1 answers that from the stored logliks alone -- no model, no counts, nothing to get wrong.
# It always runs.
#
# PART 2 asks WHICH GENES drove the decision. The likelihood is a sum over genes, so the margin
# decomposes exactly, per gene. That requires recomputing InSituType's likelihood, which is the
# one thing here that can be wrong -- so it is GATED ON A VALIDATION: the per-gene contributions
# must re-sum to the stored margin. If they do not, the script reports the discrepancy and exits
# non-zero rather than emitting a decomposition built on the wrong model. Part 1's outputs are
# already written by then.
#
# Inputs:
#   --typing-rds  anchor_typing.rds from 72 (clust, logliks, profiles). NOT the h5, which drops
#                 $logliks.
#   --counts-h5   anchor_input.h5 (stage 4a) -- PART 2 ONLY. Omit to run Part 1 alone.
#   --letter      de-novo letter, e.g. l
#   --destination named GBmap type to compare against, e.g. Endo_capilar
# Outputs (--output-dir):
#   margin_summary.csv     per-group margin stats, raw and per-count
#   margin_histogram.csv   binned pooled density, for plotting on the Mac
#   density_at_zero.csv    the verdict statistic: density at the boundary vs the modes
#   gene_contributions.csv PART 2, only if validation passes

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
SUBSAMPLE <- as.integer(opt$subsample %||% 20000)
NB_SIZE   <- as.numeric(opt[["nb-size"]] %||% 10)
TOL       <- as.numeric(opt$tolerance %||% 1e-4)
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

# ---- PART 1: the margin, and where zero sits in its density ------------------
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

summ <- dt[, .(n = .N, median_margin = median(margin), q10 = quantile(margin, .1),
               q90 = quantile(margin, .9), median_depth = median(depth, na.rm = TRUE),
               median_margin_per_count = median(margin_per_count, na.rm = TRUE)), by = group]
fwrite(summ, file.path(outdir, "margin_summary.csv"))
print(summ)

# Density of the pooled margin. The boundary is at zero by construction; the question is whether
# zero lands in a trough (two populations) or in the bulk (one continuum, arbitrarily cut).
verdicts <- list()
for (scale_name in c("margin", "margin_per_count")) {
  v <- dt[[scale_name]]
  v <- v[is.finite(v)]
  if (length(v) < 100) next
  # Trim the extreme tails so the bandwidth is not set by outliers.
  lim <- quantile(v, c(0.001, 0.999))
  d <- density(v[v >= lim[1] & v <= lim[2]], n = 2048)
  at0 <- approx(d$x, d$y, xout = 0)$y
  peak <- max(d$y)
  # A trough at the boundary means a genuine gap. Look for a local minimum within the middle
  # half of the range, and report how deep it is relative to the surrounding modes.
  inner <- which(d$x > quantile(v, .05) & d$x < quantile(v, .95))
  local_min <- if (length(inner) > 10) d$x[inner][which.min(d$y[inner])] else NA_real_
  verdicts[[scale_name]] <- data.table(
    scale = scale_name, density_at_zero = at0, peak_density = peak,
    ratio_at_zero = at0 / peak, mode_x = d$x[which.max(d$y)],
    local_min_x = local_min,
    trough_at_zero = isTRUE(abs(local_min) < diff(range(d$x)) * 0.02) && at0 / peak < 0.5)
  fwrite(data.table(scale = scale_name, x = d$x, density = d$y),
         file.path(outdir, sprintf("margin_histogram_%s.csv", scale_name)))
}
verdict <- rbindlist(verdicts)
fwrite(verdict, file.path(outdir, "density_at_zero.csv"))
print(verdict)
message("\nREAD: ratio_at_zero near 1 = the boundary cuts through the PEAK of one population,\n",
        "so the letter and the named type are one continuum split arbitrarily. A ratio near 0\n",
        "with a trough at zero = a genuine gap, i.e. two populations.\n")

# ---- PART 2: which genes decided it ------------------------------------------
if (is.null(opt[["counts-h5"]])) {
  message("No --counts-h5; skipping the per-gene decomposition (Part 1 outputs are written).")
  quit(save = "no", status = 0)
}
if (is.null(res$profiles)) stop("no $profiles in the rds; cannot decompose per gene.")

message(sprintf("%s, PART 2: per-gene decomposition", Sys.time()))
# Record the exact model InSituType uses, so a validation failure is diagnosable rather than a
# guessing game. lldist is not exported; reaching it with ::: is the same shim 72 uses for
# estimateBackground.
for (fn in c("lldist", "lls_rna")) {
  f <- tryCatch(get(fn, envir = asNamespace("InSituType")), error = function(e) NULL)
  if (is.null(f)) { message("WARNING: could not reach InSituType:::", fn); next }
  message("---- InSituType:::", fn, " (", paste(names(formals(f)), collapse = ", "), ") ----")
  print(body(f))
  message("---- end ", fn, " ----")
}
# lls_rna is the branch lldist actually takes when bg is a per-cell vector, which is our case, so
# it -- not the else branch printed from lldist -- holds the scaling that has to be reproduced.

prof <- as.matrix(res$profiles)
stopifnot("letter missing from $profiles" = LETTER %in% colnames(prof),
          "destination missing from $profiles" = DEST %in% colnames(prof))

set.seed(1)
idx_a <- sample(which(in_a), min(SUBSAMPLE, sum(in_a)))
idx_b <- sample(which(in_b), min(SUBSAMPLE, sum(in_b)))
cells <- rownames(ll)[c(idx_a, idx_b)]
shared <- intersect(rownames(prof), rownames(cnt))
message(sprintf("decomposing %d cells over %d shared genes", length(cells), length(shared)))

x <- Matrix::t(cnt[shared, cells, drop = FALSE])          # cells x genes
bg <- InSituType:::estimateBackground(counts = x, neg = neg_all[cells])
message(sprintf("bg per cell: median %.4f, range %.4f-%.4f",
                median(bg), min(bg), max(bg)))
pl <- prof[shared, LETTER]; pd <- prof[shared, DEST]

# The model, read off the lldist source printed above rather than assumed:
#
#     bgsub_i = sum_g max(count_ig - bg_i, 0)     per-cell, background-subtracted total
#     s_i^k   = bgsub_i / sum(profile_k)          scaling is PER CLUSTER, not per cell alone
#     yhat_ig = s_i^k * profile_gk + bg_i
#     ll_ig   = dnbinom(count_ig, mu = yhat_ig, size, log = TRUE)
#
# The first version of this script used s_i = rowSums(counts) — raw depth, with no background
# subtraction and no division by the profile total — and the gate caught it at a median relative
# error of 296. The per-cluster divisor is the part that is easy to miss: sum(profile_k) differs
# between clusters, so s is not a property of the cell alone.
# lls_rna is compiled, so its scaling cannot be read -- but it takes `bgsub` as an ARGUMENT
# rather than deriving it, which is what makes an exact decomposition possible. Reimplementing
# the likelihood was the wrong approach twice (median relative error 296, then 0.81); instead
# call the package's own code, holding bgsub fixed at its all-genes value while subsetting the
# genes. Then the per-gene terms are the package's, not an approximation of them.
#
# bgsub itself IS visible, in the lldist source: background is subtracted from each nonzero
# entry, floored at zero, and summed per cell.
bgsub_of <- function(counts_mat, bgv) {
  b <- counts_mat
  b@x <- pmax(b@x - bgv[b@i + 1], 0)
  Matrix::rowSums(b)
}
bgsub <- bgsub_of(x, bg)
pair <- prof[shared, c(LETTER, DEST), drop = FALSE]
stored <- ll[cells, LETTER] - ll[cells, DEST]

# Step 1: can the package's own lldist reproduce the stored logliks at all? If not, the stored
# table came from different inputs (gene set, background, profiles) and no decomposition of ours
# will match it -- that is a different problem and must not be papered over.
sweep <- data.table(size = numeric(), err = numeric())
best <- list(size = NA_real_, err = Inf, ll = NULL)
for (sz in unique(c(NB_SIZE, 10, 1, 0.5, 100))) {
  got <- tryCatch(
    InSituType:::lldist(x = pair, mat = x, bg = bg, size = sz, digits = 12,
                        assay_type = "rna"),
    error = function(e) { message("  lldist(size=", sz, ") failed: ",
                                  conditionMessage(e)); NULL })
  if (is.null(got)) next
  m <- got[, LETTER] - got[, DEST]
  e <- median(abs(m - stored) / pmax(abs(stored), 1))
  sweep <- rbind(sweep, data.table(size = sz, err = e))
  message(sprintf("  lldist size %-6s median relative error %.4g", format(sz), e))
  if (e < best$err) best <- list(size = sz, err = e, ll = got)
}
fwrite(sweep, file.path(outdir, "lldist_size_sweep.csv"))
if (!is.finite(best$err) || best$err > TOL) {
  stop(sprintf(paste0("InSituType's OWN lldist does not reproduce the stored margin ",
                      "(best median relative error %.4g at size %s). The stored logliks were ",
                      "produced from different inputs -- most likely a different gene set, a ",
                      "different background, or profiles updated after the logliks were saved. ",
                      "Compare the fit's gene panel and bg against what this job staged before ",
                      "decomposing anything. Part 1's outputs are valid and already written."),
               best$err, format(best$size)))
}
SIZE <- best$size
message(sprintf("lldist reproduces the stored margin at size %s (median relative error %.3g)",
                format(SIZE), best$err))

# Step 2: the same call gene by gene, with bgsub pinned to its all-genes value.
message(sprintf("%s, decomposing %d genes", Sys.time(), length(shared)))
contrib <- matrix(NA_real_, nrow = length(cells), ncol = length(shared),
                  dimnames = list(cells, shared))
for (gi in seq_along(shared)) {
  g <- shared[gi]
  r <- InSituType:::lls_rna(mat = x[, g, drop = FALSE], bgsub = bgsub,
                            x = pair[g, , drop = FALSE], bg = bg, size_dnb = SIZE)
  contrib[, gi] <- r[, LETTER] - r[, DEST]
  if (gi %% 500 == 0) message(sprintf("  %d/%d genes", gi, length(shared)))
}

# THE GATE, now asking only whether the per-gene calls sum to the whole-gene call -- a question
# about separability, not about whether we guessed the model right.
recomputed <- rowSums(contrib)
err <- abs(recomputed - stored) / pmax(abs(stored), 1)
message(sprintf("validation: median relative error %.3g, 90th pct %.3g, max %.3g",
                median(err), quantile(err, .9), max(err)))
if (median(err) > TOL) {
  fwrite(data.table(cell_id = cells, stored = stored, recomputed = recomputed, rel_err = err),
         file.path(outdir, "validation_failure.csv"))
  stop(sprintf(paste0("the per-gene terms do not sum to the whole (median relative error %.3g). ",
                      "lls_rna is then not separable per gene at fixed bgsub -- decompose by ",
                      "leave-one-gene-out against the full call instead. Part 1 is unaffected."),
               median(err)))
}

per_gene <- data.table(gene = shared,
                       mean_contrib_A = colMeans(contrib[seq_along(idx_a), , drop = FALSE]),
                       mean_contrib_B = colMeans(contrib[-seq_along(idx_a), , drop = FALSE]))
per_gene[, difference := mean_contrib_A - mean_contrib_B]
setorder(per_gene, -mean_contrib_A)
fwrite(per_gene, file.path(outdir, "gene_contributions.csv"))
message("\nTop 15 genes pushing A cells toward the letter:")
print(head(per_gene, 15))
message("\nTop 15 pushing them toward the named type:")
print(head(per_gene[order(mean_contrib_A)], 15))
message(sprintf("\nWrote %s", file.path(outdir, "gene_contributions.csv")))
