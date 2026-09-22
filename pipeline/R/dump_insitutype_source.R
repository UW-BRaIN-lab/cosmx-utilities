#!/usr/bin/env Rscript
# Print the installed InSituType's assignment path. See 75n_insitutype_source.sh for why.
#
# Deliberately prints WHOLE function bodies rather than grepped lines: the question is about
# ORDER — whether the matrix that clust is taken from is the one that gets returned — and a
# grep drops exactly the surrounding lines that settle it.

suppressPackageStartupMessages(library(InSituType))

cat("InSituType version:", as.character(packageVersion("InSituType")), "\n")
cat("library path:", dirname(system.file(package = "InSituType")), "\n\n")

show_fn <- function(name) {
  cat(strrep("=", 78), "\n", name, "\n", strrep("=", 78), "\n", sep = "")
  f <- tryCatch(get(name, envir = asNamespace("InSituType")), error = function(e) NULL)
  if (is.null(f)) { cat("  (not found in the namespace)\n\n"); return(invisible()) }
  cat(paste(deparse(f), collapse = "\n"), "\n\n")
}

for (fn in c("insitutypeML", "update_logliks_with_cohort_freqs", "insitutype")) show_fn(fn)

cat(strrep("=", 78), "\n", "exported + internal object names\n", strrep("=", 78), "\n", sep = "")
print(sort(ls(asNamespace("InSituType"))))
