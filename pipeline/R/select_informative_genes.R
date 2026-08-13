#!/usr/bin/env Rscript
# InSituType 6k-panel gene selection (FAQ "Which genes to use").
#
# NanoString's 6000-plex FAQ recommends typing on a well-chosen ~3000-5000 gene subset
# rather than the full panel. A gene is RETAINED if EITHER criterion holds (UNION, so a
# gene informative in only one biology is not dropped):
#   (1) solidly above background in the CosMx data  (mean count > bg-mult x mean neg), OR
#   (2) moderate-to-high in at least one reference profile  (max over types >= a quantile).
#
# This is the CHEAP tuning step: it only computes per-gene means over the anchor and the
# reference, so run it first and sweep --bg-mult / --ref-quantile until the UNION lands in
# the 3000-5000 band. The expensive de-novo K-sweep (choose_denovo_k.R) and the anchor
# re-fit (insitutype_typing.R --keep-genes) then consume kept_genes.txt.
#
# Inputs mirror 72_anchor_typing.sh: the stage-4a anchor_input.h5 and the panel-restricted
# reference CSV (genes x cell types, first column "gene").
#
# Outputs (into --out-dir):
#   kept_genes.txt          one retained gene per line  -> insitutype_typing.R --keep-genes
#   excluded_genes.csv      comma-joined dropped genes  -> 85b EXCLUDED_GENES (full-cohort typing)
#   gene_selection.csv      per-gene stats + both criteria + kept flag (re-threshold offline)
#   gene_selection_summary.txt  crit counts + whether the union is in the target band
#
# Usage:
#   Rscript pipeline/R/select_informative_genes.R \
#     --input anchor_input.h5 --reference reference.csv --out-dir gene_selection \
#     [--bg-mult 3] [--ref-quantile 0.5] [--target-min 3000] [--target-max 5000] \
#     [--max-cells 500000]   # subsample for a faster mean estimate; means are stable

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "read_anchor_h5.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument: %s", a))
    key <- sub("^--", "", a)
    if (grepl("=", key)) { kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]; i <- i + 1
    } else { out[[key]] <- args[[i + 1]]; i <- i + 2 }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot("--input is required" = !is.null(opt$input))
stopifnot("--reference is required" = !is.null(opt$reference))
out_dir     <- opt[["out-dir"]] %||% "gene_selection"
bg_mult     <- as.numeric(opt[["bg-mult"]] %||% 3)
ref_quantile<- as.numeric(opt[["ref-quantile"]] %||% 0.5)
target_min  <- as.integer(opt[["target-min"]] %||% 3000)
target_max  <- as.integer(opt[["target-max"]] %||% 5000)
max_cells   <- if (is.null(opt[["max-cells"]])) NULL else as.integer(opt[["max-cells"]])
seed        <- as.integer(opt$seed %||% 0)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- anchor counts (cells x genes) + per-cell background ----------------------
dat <- read_anchor_h5(opt$input, max_cells = max_cells, seed = seed)
gene_mean <- Matrix::colMeans(dat$x)               # mean count per gene across cells
bg <- mean(dat$neg)                                # cohort mean background rate
message(sprintf("panel genes: %d ; cohort mean background (neg): %.4f", length(gene_mean), bg))

# --- reference profiles (genes x cell types) ----------------------------------
ref_dt <- data.table::fread(opt$reference)
ref_genes <- as.character(ref_dt[[1]])
ref <- as.matrix(ref_dt[, -1L]); rownames(ref) <- ref_genes
storage.mode(ref) <- "double"
ref_max_by_gene <- apply(ref, 1, max)              # per gene: max expression across types
message(sprintf("reference: %d genes x %d cell types", nrow(ref), ncol(ref)))

# --- FAQ criteria, evaluated over the PANEL genes -----------------------------
panel <- names(gene_mean)
ref_max <- ref_max_by_gene[panel]                  # NA for panel genes absent from reference
ref_max[is.na(ref_max)] <- 0
# crit2 threshold uses the reference's own distribution (genes present in the reference),
# so "moderate-to-high in >=1 profile" is judged against real profile magnitudes.
ref_cut <- as.numeric(stats::quantile(ref_max_by_gene, ref_quantile))

crit1 <- gene_mean > bg_mult * bg                  # solidly above background
crit2 <- ref_max >= ref_cut                        # moderate-to-high in >=1 profile
kept  <- crit1 | crit2                             # UNION

report <- data.table(
  gene = panel,
  mean_count = as.numeric(gene_mean),
  above_bg = crit1,
  ref_max = as.numeric(ref_max),
  in_reference = panel %in% ref_genes,
  ref_informative = crit2,
  kept = kept
)
data.table::fwrite(report, file.path(out_dir, "gene_selection.csv"))
writeLines(panel[kept], file.path(out_dir, "kept_genes.txt"))
writeLines(paste(panel[!kept], collapse = ","), file.path(out_dir, "excluded_genes.csv"))

n_keep <- sum(kept)
in_band <- n_keep >= target_min && n_keep <= target_max
summary_lines <- c(
  sprintf("panel genes:            %d", length(panel)),
  sprintf("bg-mult:                %.2f  (background rate = %.4f)", bg_mult, bg),
  sprintf("ref-quantile:           %.2f  (ref_max cutoff = %.4f)", ref_quantile, ref_cut),
  sprintf("crit1 above-background: %d", sum(crit1)),
  sprintf("crit2 ref-informative:  %d", sum(crit2)),
  sprintf("UNION kept:             %d  (excluded %d)", n_keep, length(panel) - n_keep),
  sprintf("target band [%d, %d]: %s", target_min, target_max,
          if (in_band) "IN BAND" else "OUT OF BAND -- retune --bg-mult / --ref-quantile"),
  if (!in_band && n_keep > target_max)
    "  too many kept: raise --bg-mult and/or --ref-quantile"
  else if (!in_band && n_keep < target_min)
    "  too few kept: lower --bg-mult and/or --ref-quantile" else NULL
)
writeLines(summary_lines, file.path(out_dir, "gene_selection_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
message("wrote ", file.path(out_dir, "kept_genes.txt"))
