# The export's guard against being read as raw counts (docs/PLAN.md C3f, C5's K11; milestone
# M13): Seurat 5.5.1's Read10X(), the usual way R reads a 10x directory, on every directory of a
# `scrust cluster --export` output and of a `scrust sctransform --corrected-counts` output. It must
# fail everywhere but the export's counts/, which it must read as the input counts: the same
# dimensions, entries and values as the input directory read the same way. Reads plain files
# only (no .rds); run through sbatch in the production R env (validation/slurm/m13_guard.sbatch).
#
#   Rscript validation/R/guard_read10x.R <out.tsv> <export dir> <input counts dir> [<sct dir> ...]
#
# Writes one row per directory tried: the directory (relative to its output), whether Read10X
# read it, the error or the dimensions, and whether that is what the guard requires. Exit status
# 1 if any directory behaves otherwise.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: guard_read10x.R <out.tsv> <export dir> <input counts dir> [<sct dir> ...]")
suppressPackageStartupMessages({
  library("Matrix")
  library("Seurat")
})
out <- args[1]
export <- normalizePath(args[2])
input <- normalizePath(args[3])
others <- if (length(args) > 3) normalizePath(args[-(1:3)]) else character()

dirs_under <- function(root) {
  d <- list.dirs(root, full.names = TRUE, recursive = TRUE)
  d[!grepl("/stages(/|$)", d)]
}

try_read <- function(d) {
  tryCatch(list(ok = TRUE, m = suppressWarnings(Read10X(data.dir = d)), error = ""),
           error = function(e) list(ok = FALSE, m = NULL, error = conditionMessage(e)))
}

rows <- list()
pass <- TRUE
want_counts <- Read10X(data.dir = input)
for (root in c(export, others)) {
  for (d in dirs_under(root)) {
    rel <- if (d == root) "." else substring(d, nchar(root) + 2)
    is_counts <- root == export && rel == "counts"
    r <- try_read(d)
    detail <- r$error
    ok <- if (is_counts) r$ok else !r$ok
    if (r$ok) {
      m <- r$m
      detail <- sprintf("read: %d x %d, %d entries", nrow(m), ncol(m), length(m@x))
      if (is_counts) {
        same <- identical(dim(m), dim(want_counts)) && identical(dimnames(m), dimnames(want_counts)) &&
          identical(m@i, want_counts@i) && identical(m@p, want_counts@p) && identical(m@x, want_counts@x)
        detail <- paste0(detail, if (same) "; identical to the input read by Read10X" else "; DIFFERS from the input")
        ok <- ok && same
      }
    }
    pass <- pass && ok
    rows[[length(rows) + 1]] <- data.frame(output = basename(root), dir = rel, expected = if (is_counts) "read" else "fail",
                                           read = r$ok, as_required = ok, detail = gsub("[\t\n]", " ", detail),
                                           stringsAsFactors = FALSE)
    cat(sprintf("%-28s %-12s %-6s %s\n", basename(root), rel, if (ok) "ok" else "WRONG", detail))
  }
}
tab <- do.call(rbind, rows)
write.table(tab, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Seurat %s Read10X: %d directories, %d as required\n", as.character(packageVersion("Seurat")), nrow(tab), sum(tab$as_required)))
if (!pass) quit(status = 1)
