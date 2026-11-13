# Reads the singlets/ directories of Rust runs (milestone M6b) back with R's readers, as a
# downstream R step would, and compares each with the saved workflow object's RNA counts as
# exported (data/doubletfinder/<object>/saved/singlets/, export_doubletfinder_singlets.R):
# - lib/io.R's exact reader (mtx_read_10x): must be identical(), slot for slot, explicit zeros
#   included;
# - Seurat::Read10X(), the usual way R reads a 10x directory: the same matrix, values and
#   dimnames, with whatever its readMM() conversion does to explicit zeros recorded.
# Reads only plain files.
#
#   Rscript validation/R/read10x_singlets.R <label> <profile> <object> [<object> ...]

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
setwd(normalizePath(file.path(here, "..", "..")))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: read10x_singlets.R <label> <profile> <object> [<object> ...]")
suppressPackageStartupMessages({
  library("Matrix")
  library("Seurat")
})
label <- args[1]
profile <- args[2]
ok <- TRUE
for (obj in args[-(1:2)]) {
  d <- file.path("data", "doubletfinder", obj, "rust", label, profile, "run", "singlets")
  ref <- mtx_read_10x(file.path("data", "doubletfinder", obj, "saved", "singlets"))
  exact <- mtx_read_10x(d)
  same_exact <- identical(exact, ref, num.eq = FALSE)
  rm(exact)
  r10x <- Read10X(d)
  same_values <- identical(dim(r10x), dim(ref)) && identical(dimnames(r10x), dimnames(ref)) &&
    identical(drop0(r10x), drop0(ref), num.eq = FALSE)
  same_slots <- identical(r10x, ref, num.eq = FALSE)
  cat(sprintf("%s %s: %d x %d, %d entries (%d explicit zeros); exact reader identical: %s; Read10X: same matrix %s, every stored entry kept %s (%d entries)\n",
              obj, profile, nrow(ref), ncol(ref), length(ref@x), sum(ref@x == 0), same_exact, same_values,
              same_slots, length(r10x@x)))
  ok <- ok && same_exact && same_values
  rm(ref, r10x)
  invisible(gc())
}
if (!ok) quit(status = 1)
