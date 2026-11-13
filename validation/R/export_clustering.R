# One-time export of a clustered object's chain input and results (docs/PLAN.md C4; milestone
# M8). Reads .rds files, so run it only through validation/slurm/export_clustering.sbatch, never on
# a login node.
#
#   Rscript validation/R/export_clustering.R --object <name> --rds <saved .rds> --markers <csv> \
#     --out <data/clustering/<object>> --level full|standard|light|labels --counts-from <dir>|-
#
# <out>/input/: the chain's input, the object's RNA counts layer as saved (C1):
#   matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz  10x Matrix Market, entries exactly as stored
#       (explicit zeros included), cells and features in object order (lib/io.R mtx_write_10x).
#       With --counts-from <dir>, an earlier export of the same counts (M2's input/ of a
#       DoubletFinder input, M6b's saved/singlets/ of a DoubletFinder output) is reused instead:
#       the object's counts must equal that directory's, read back exactly, slot for slot, and the
#       three files become relative symbolic links to it.
#   assay.txt     v3 or v5: the class of the object's RNA assay (Assay or Assay5), which decides
#                 SCTransform's dispatch (C1); the reference runs build the same
#   export.json   provenance, dimensions, value checks, the round trip, and what was reused
# <out>/saved/: what the chain left in the object, in lib/clustering.R's layout (cl_write_object,
#   at --level), so a reference run compares with it file by file; plus
#   object.json       facts and checks (lib/clustering.R), the RNA counts' and object's provenance,
#                     and the parameters the object's command log recorded for the chain's calls
#   markergenes.csv   a copy of the workflow's combined marker table for this object
# <out>/export_run.json: timings, peak RSS, SLURM job and environment fingerprint.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))

args <- df_args(list(object = NULL, rds = NULL, markers = NULL, out = NULL, level = NULL, `counts-from` = NULL))
cl_setup()
t0 <- proc.time()[["elapsed"]]
elapsed <- function() proc.time()[["elapsed"]] - t0
in_dir <- file.path(args$out, "input")
saved_dir <- file.path(args$out, "saved")
if (file.exists(file.path(saved_dir, "object.json"))) stop(saved_dir, " exists; not overwriting")
dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(saved_dir, recursive = TRUE, showWarnings = FALSE)

cat("reading", args$rds, "\n")
o <- readRDS(args$rds)
read_s <- elapsed()
version <- as.character(SeuratObject::Version(o))
cl_check_names(o)
counts <- LayerData(o, assay = "RNA", layer = "counts")
if (!inherits(counts, "dgCMatrix")) stop("RNA counts layer is a ", class(counts)[1], ", not a dgCMatrix")
invisible(validObject(counts))
if (!identical(colnames(counts), colnames(o))) stop("RNA counts columns are not the object's cells")
if (anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts))) stop("duplicated cells or features")
vst_attr <- intersect(CL_VST_ATTR, colnames(o@meta.data))
if (length(vst_attr)) stop("meta.data has columns vst() would use as cell attributes: ", paste(vst_attr, collapse = ", "))
rna_class <- class(o[["RNA"]])[1]
assay <- switch(rna_class, Assay = "v3", Assay5 = "v5", stop("RNA assay class ", rna_class))
x <- counts@x
info <- list(
  object = args$object,
  source = normalizePath(args$rds),
  source_sha256 = sha256_file(args$rds),
  source_bytes = file.size(args$rds),
  seurat_object_version = version,
  rna_assay_class = rna_class,
  assay = assay,
  rna_layers = Layers(o[["RNA"]]),
  n_features = nrow(counts),
  n_cells = ncol(counts),
  nnz = length(x),
  explicit_zeros = sum(x == 0),
  non_integer_values = sum(x != round(x)),
  min_value = min(x),
  max_value = max(x)
)

# Input: the RNA counts, written or reused -----------------------------------------------------
files <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
if (identical(args$`counts-from`, "-")) {
  info$mtx_field <- mtx_write_10x(counts, in_dir, software = "scRust validation/R/export_clustering.R")
  info$counts_from <- NULL
} else {
  src <- normalizePath(args$`counts-from`)
  rel_dir <- system2("realpath", c("--relative-to", shQuote(normalizePath(in_dir)), shQuote(src)), stdout = TRUE)
  for (f in files) {
    target <- file.path(in_dir, f)
    if (file.exists(target) || nzchar(Sys.readlink(target))) unlink(target)
    if (!file.symlink(file.path(rel_dir, f), target)) stop("cannot link ", target)
    if (!identical(normalizePath(target), normalizePath(file.path(src, f)))) stop("link ", target, " does not resolve to ", src)
  }
  info$counts_from <- sub(paste0("^", normalizePath(file.path(here, "..", "..")), "/"), "", src)
  info$counts_from_sha256 <- as.list(vapply(files, function(f) sha256_file(file.path(src, f)), character(1)))
}
writeLines(assay, file.path(in_dir, "assay.txt"))
back <- mtx_read_10x(in_dir)
same <- identical(back, counts, num.eq = FALSE)
if (!same) {
  if (!identical(back@i, counts@i) || !identical(back@p, counts@p) || !identical(dimnames(back), dimnames(counts))) {
    stop("round trip changed the matrix structure or dimnames")
  }
  stop("round trip changed ", sum(back@x != counts@x), " values")
}
info$round_trip_identical <- same
rm(back, counts, x)
invisible(gc())
input_s <- elapsed()

# Saved: what the chain left -------------------------------------------------------------------
facts <- cl_write_object(o, saved_dir, args$level, software = "scRust validation/R/export_clustering.R")
cmd_params <- lapply(o@commands, function(cmd) {
  p <- slot(cmd, "params")
  p[vapply(p, function(v) is.null(v) || is.atomic(v), logical(1))]
})
if (!file.copy(args$markers, file.path(saved_dir, "markergenes.csv"), overwrite = FALSE, copy.date = TRUE)) stop("cannot copy the marker table")
json_write(c(facts, list(object = args$object, source = info$source, source_sha256 = info$source_sha256,
                         markers_source = normalizePath(args$markers),
                         markers_sha256 = sha256_file(args$markers),
                         commands = names(o@commands), command_params = cmd_params)),
           file.path(saved_dir, "object.json"))
json_write(info, file.path(in_dir, "export.json"))
json_write(list(timing_s = list(read_rds = read_s, input = input_s - read_s, saved = elapsed() - input_s,
                                total = elapsed()),
                peak_rss_kib = peak_rss_kib(), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                fingerprint = fingerprint(packages = c("Seurat", "SeuratObject", "sctransform", "glmGamPoi",
                                                       "irlba", "RcppAnnoy", "bluster", "BiocNeighbors",
                                                       "Matrix", "Rcpp"))),
           file.path(args$out, "export_run.json"))
cat("exported", args$object, ":", info$n_features, "features x", info$n_cells, "cells (", assay, "),",
    facts$n_clusters, "clusters,", facts$sct_n_features, "SCT features, level", args$level,
    if (!is.null(info$counts_from)) paste("; counts reused from", info$counts_from) else "",
    "; peak RSS", peak_rss_kib(), "KiB,", round(elapsed()), "s\n")
