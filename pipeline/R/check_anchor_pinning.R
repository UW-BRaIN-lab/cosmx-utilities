#!/usr/bin/env Rscript
# Are the label-vs-argmax mismatches simply InSituType's anchor cells?
#
# The 2.0 source settles the mechanism on paper. insitutype() is called here with anchors=NULL,
# but when update_reference_profiles=TRUE it DERIVES anchors — updateReferenceProfiles() picks
# the cells that best match each reference profile (n_anchor_cells, min_anchor_cosine,
# min_anchor_llr), uses them to rescale the reference, and returns them — after which
#     anchors <- update_result$anchors
# makes the argument non-NULL, and the final line of insitutype() is
#     out$clust[!is.na(anchors)] <- anchors[names(out$clust[!is.na(anchors)])]
# a hard overwrite that ignores the likelihoods entirely.
#
# The fit saves $anchors, so this is a join rather than a re-run. If the mismatched cells ARE
# the anchored ones, the "native groups were not assigned by the likelihoods" reading is wrong:
# those cells are the ones InSituType deliberately selected as the BEST match to each reference
# profile, which is close to the opposite of a weakly-supported residue.

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
rds <- args[[1]]
out_csv <- if (length(args) > 1) args[[2]] else NULL

res <- readRDS(rds)
cat("objects in the fit:", paste(names(res), collapse = ", "), "\n\n")

if (is.null(res$anchors)) {
  cat("NO $anchors IN THE FIT — the pinning hypothesis is REJECTED for this run.\n")
  quit(status = 0)
}

anchors <- res$anchors
clust <- res$clust
ll <- res$logliks
types <- colnames(ll)

common <- intersect(rownames(ll), names(clust))
ll <- ll[common, , drop = FALSE]; clust <- clust[common]
anchors <- anchors[common]          # named over the same cells; NA where not an anchor

is_anchor <- !is.na(anchors)
argmax <- types[max.col(ll, ties.method = "first")]
is_mismatch <- clust != argmax

cat(sprintf("cells: %d   anchors: %d (%.2f%%)   label!=argmax: %d (%.2f%%)\n\n",
            length(clust), sum(is_anchor), 100 * mean(is_anchor),
            sum(is_mismatch), 100 * mean(is_mismatch)))

tab <- table(anchor = is_anchor, mismatch = is_mismatch)
cat("anchor x mismatch:\n"); print(tab); cat("\n")

if (sum(is_anchor) > 0) {
  cat(sprintf("of ANCHORED cells,     %.1f%% have label != argmax\n",
              100 * mean(is_mismatch[is_anchor])))
  cat(sprintf("of NON-anchored cells, %.1f%% have label != argmax\n",
              100 * mean(is_mismatch[!is_anchor])))
  cat(sprintf("\nof MISMATCHED cells,   %.1f%% are anchors\n",
              100 * mean(is_anchor[is_mismatch])))
  cat(sprintf("anchored cells whose label equals their anchor type: %.1f%%\n\n",
              100 * mean(clust[is_anchor] == as.character(anchors[is_anchor]))))
}

per <- data.table(label = clust, anchor = is_anchor, mismatch = is_mismatch)[
  , .(n = .N, pct_anchored = round(100 * mean(anchor), 1),
      pct_mismatch = round(100 * mean(mismatch), 1)), by = label][
  n >= 100][order(-pct_mismatch)]
cat("by assigned label (>=100 cells), worst mismatch first:\n")
print(head(per, 20))

if (!is.null(out_csv)) {
  fwrite(per, out_csv)
  cat(sprintf("\nWrote %s\n", out_csv))
}
cat("\nREAD: if pct_anchored tracks pct_mismatch, the mismatches are the pinning step and\n",
    "nothing more — and the affected cells are the reference's own best exemplars.\n", sep = "")
