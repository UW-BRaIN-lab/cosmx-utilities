#!/usr/bin/env Rscript
# Print the installed InSituType's assignment path. See 75n_insitutype_source.sh for why.
#
# Deliberately prints WHOLE bodies rather than grepped lines: the question is about ORDER —
# whether the matrix `clust` is taken from is the one that gets returned — and a grep drops
# exactly the surrounding lines that settle it.
#
# In InSituType 2.0 both insitutype() and insitutypeML() are S4 GENERICS, so deparsing the
# object by name yields only the dispatch stub. The implementation lives in the methods, which
# is what this walks.

suppressPackageStartupMessages(library(InSituType))

NS <- asNamespace("InSituType")
cat("InSituType version:", as.character(packageVersion("InSituType")), "\n")
cat("library path:", dirname(system.file(package = "InSituType")), "\n\n")

rule <- function(ch = "=") cat(strrep(ch, 78), "\n", sep = "")

show_body <- function(label, f) {
  rule("-")
  cat(label, "\n")
  rule("-")
  cat(paste(deparse(f), collapse = "\n"), "\n\n")
}

show_fn <- function(name) {
  rule(); cat(name, "\n"); rule()
  obj <- tryCatch(get(name, envir = NS), error = function(e) NULL)
  if (is.null(obj)) { cat("  (not found in the namespace)\n\n"); return(invisible()) }

  if (isGeneric(name, where = NS) || methods::is(obj, "standardGeneric")) {
    sigs <- tryCatch(findMethods(name, where = NS), error = function(e) NULL)
    if (is.null(sigs) || length(sigs) == 0) {
      sigs <- tryCatch(findMethods(name), error = function(e) NULL)
    }
    cat("S4 generic with", length(sigs), "method(s):",
        paste(names(sigs), collapse = ", "), "\n\n")
    for (sig in names(sigs)) show_body(sprintf("%s  method: %s", name, sig), sigs[[sig]])
  } else {
    show_body(name, obj)
  }
}

# updateReferenceProfiles is the one that SELECTS the pinned cells — everything about which
# cells get their label overwritten is decided in there — so it is dumped alongside the two
# generics and the refinement passes.
for (fn in c("insitutypeML", "insitutype", "updateReferenceProfiles",
             "update_logliks_with_cohort_freqs", "refineClusters", "refineAnchors",
             "chooseAnchors")) show_fn(fn)

rule(); cat("exported + internal object names\n"); rule()
print(sort(ls(NS)))
