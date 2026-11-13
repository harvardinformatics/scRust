# One-time export of the singlets' counts in the workflow's saved DoubletFinder outputs
# (milestone M6b; docs/PLAN.md A3 and B4): the RNA assay's "counts" layer of
# seurat_doubletfinder_<object>.rds. The workflow's doubletfinder.R keeps the singlets with
# subset(); its post-hoc SCTransform and clustering add assays and leave the RNA counts as they
# are. So this layer is DoubletFinder's count-level result, what `scrust doubletfinder` writes to
# singlets/. Reads .rds files, so run it only through
# validation/slurm/export_doubletfinder_singlets.sbatch, never on a login node.
#
#   Rscript validation/R/export_doubletfinder_singlets.R --object <name> \
#     --saved <seurat_doubletfinder_<object>.rds> --export <data/doubletfinder/<object>> --out <dir>
#
# <out>/singlets/: matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz (lib/io.R mtx_write_10x:
#   entries exactly as stored, explicit zeros included, cells and features in object order; the
#   feature lines are the input export's), and export.json (provenance and checks).
# <out>/export_singlets_run.json: timings, peak RSS, SLURM job and environment fingerprint.
#
# Fatal checks:
#   - the .rds is the one milestone M2 exported (sha256 = saved/saved.json's source_sha256);
#   - its cells are saved/cells.tsv's kept cells, in order, and its features the input's, in order;
#   - the written directory reads back identical (mtx_read_10x; identical(), slot for slot).
# Recorded, not fatal: whether the layer equals the input export's counts restricted to the kept
# cells, slot for slot (stored entries, explicit zeros and their order included). That is what
# `scrust doubletfinder` writes; the comparison shows whether R's subset() keeps the stored
# entries as they are.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))

args <- df_args(list(object = NULL, saved = NULL, export = NULL, out = NULL))
suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
})
t0 <- proc.time()[["elapsed"]]
elapsed <- function() proc.time()[["elapsed"]] - t0
out_dir <- file.path(args$out, "singlets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Provenance: the object M2 exported ----------------------------------------------------------
saved_info <- jsonlite::fromJSON(file.path(args$export, "saved", "saved.json"))
sha <- sha256_file(args$saved)
if (!identical(sha, saved_info$source_sha256)) {
  stop(args$saved, ": sha256 ", sha, " is not that of the object M2 exported, ", saved_info$source_sha256)
}

# The saved object's RNA counts ----------------------------------------------------------------
cat("reading", args$saved, "\n")
o <- readRDS(args$saved)
read_s <- elapsed()
version <- as.character(SeuratObject::Version(o))
if (package_version(version) < "5.0") stop("SeuratObject version ", version, " < 5.0")
counts <- LayerData(o, assay = "RNA", layer = "counts")
if (!inherits(counts, "dgCMatrix")) stop("counts layer is a ", class(counts)[1], ", not a dgCMatrix")
invisible(validObject(counts))
cells <- rownames(o@meta.data)
rna_layers <- Layers(o[["RNA"]])
default_assay <- DefaultAssay(o)
rm(o)
invisible(gc())
if (!identical(colnames(counts), cells)) stop("counts columns are not in meta.data order")
kept <- read.delim(file.path(args$export, "saved", "cells.tsv"), colClasses = "character", quote = "")$barcode
if (!identical(cells, kept)) stop("the saved object's cells are not saved/cells.tsv's kept cells, in order")

# The input export, re-read exactly (lib/io.R), restricted to the kept cells slot for slot ------
input <- mtx_read_10x(file.path(args$export, "input"))
if (!identical(rownames(counts), rownames(input))) stop("features differ from the input's (names or order)")
feature_ids <- read.delim(file.path(args$export, "input", "features.tsv.gz"), header = FALSE,
                          colClasses = "character", quote = "", comment.char = "")[[1]]
j <- match(kept, colnames(input))
if (anyNA(j)) stop("kept cells missing from the input")
nper <- diff(input@p)[j]
idx <- sequence(nper, from = input@p[j] + 1L)
expected <- new("dgCMatrix", i = input@i[idx], p = as.integer(c(0L, cumsum(nper))), x = input@x[idx],
                Dim = c(nrow(input), length(j)), Dimnames = list(rownames(input), kept))
rm(input, idx)
invisible(gc())
equals_input_subset <- identical(counts, expected, num.eq = FALSE)
subset_detail <- list(
  nnz = length(counts@x),
  nnz_input_subset = length(expected@x),
  explicit_zeros = sum(counts@x == 0),
  explicit_zeros_input_subset = sum(expected@x == 0),
  same_structure = identical(counts@i, expected@i) && identical(counts@p, expected@p),
  same_after_drop0 = if (equals_input_subset) TRUE else identical(drop0(counts), drop0(expected), num.eq = FALSE)
)
rm(expected)
invisible(gc())
compare_s <- elapsed()

# Write, then read back ------------------------------------------------------------------------
field <- mtx_write_10x(counts, out_dir, feature_ids = feature_ids,
                       software = "scRust validation/R/export_doubletfinder_singlets.R")
write_s <- elapsed()
back <- mtx_read_10x(out_dir)
if (!identical(back, counts, num.eq = FALSE)) stop("round trip changed the matrix")
rm(back)
verify_s <- elapsed()

x <- counts@x
info <- list(
  object = args$object,
  source = normalizePath(args$saved),
  source_sha256 = sha,
  source_bytes = file.size(args$saved),
  seurat_object_version = version,
  rna_layers = rna_layers,
  default_assay = default_assay,
  n_features = nrow(counts),
  n_cells = ncol(counts),
  nnz = length(x),
  explicit_zeros = sum(x == 0),
  non_integer_values = sum(x != round(x)),
  min_value = if (length(x)) min(x) else NA,
  max_value = if (length(x)) max(x) else NA,
  mtx_field = field,
  round_trip_identical = TRUE,
  equals_input_subset = equals_input_subset,
  input_subset = subset_detail
)
json_write(info, file.path(out_dir, "export.json"))
json_write(list(timing_s = list(read_rds = read_s, compare = compare_s - read_s, write = write_s - compare_s,
                                verify = verify_s - write_s, total = elapsed()),
                peak_rss_kib = peak_rss_kib(), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                fingerprint = fingerprint()),
           file.path(args$out, "export_singlets_run.json"))
cat("exported", args$object, "singlets:", info$n_features, "features x", info$n_cells, "cells,", info$nnz,
    "entries (", info$explicit_zeros, "explicit zeros ); equals the input's counts restricted to the kept cells:",
    equals_input_subset, "; peak RSS", peak_rss_kib(), "KiB,", round(elapsed()), "s\n")
