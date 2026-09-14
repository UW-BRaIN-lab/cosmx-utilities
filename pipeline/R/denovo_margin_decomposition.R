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
ld <- tryCatch(InSituType:::lldist, error = function(e) NULL)
if (!is.null(ld)) {
  message("InSituType:::lldist formals: ", paste(names(formals(ld)), collapse = ", "))
  message("---- lldist source ----"); print(body(ld)); message("---- end ----")
} else {
  message("WARNING: could not reach InSituType:::lldist; proceeding with the documented model.")
}

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
s  <- Matrix::rowSums(x)                                  # per-cell scaling
pl <- prof[shared, LETTER]; pd <- prof[shared, DEST]

# Per-gene log-likelihood contribution under each profile. Negative binomial with the package's
# size parameter; expected counts are the cell's scaling times the profile, plus background.
xm <- as.matrix(x)
mu_l <- outer(s, pl) + bg
mu_d <- outer(s, pd) + bg
contrib <- dnbinom(xm, mu = mu_l, size = NB_SIZE, log = TRUE) -
           dnbinom(xm, mu = mu_d, size = NB_SIZE, log = TRUE)

# THE GATE. The per-gene terms must re-sum to the margin the fit actually stored.
recomputed <- rowSums(contrib)
stored <- ll[cells, LETTER] - ll[cells, DEST]
err <- abs(recomputed - stored) / pmax(abs(stored), 1)
message(sprintf("validation: median relative error %.3g, 90th pct %.3g, max %.3g",
                median(err), quantile(err, .9), max(err)))
if (median(err) > TOL) {
  fwrite(data.table(cell_id = cells, stored = stored, recomputed = recomputed, rel_err = err),
         file.path(outdir, "validation_failure.csv"))
  stop(sprintf(paste0("per-gene decomposition does NOT reproduce the stored margin ",
                      "(median relative error %.3g > %.3g). The likelihood model above is wrong ",
                      "-- read the lldist source printed earlier and correct mu/size/scaling. ",
                      "Part 1 outputs are valid and already written; validation_failure.csv has ",
                      "the per-cell discrepancy."), median(err), TOL))
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
