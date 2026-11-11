# One-time export of DoubletFinder's inputs and of the compared fields of the workflow's saved
# outputs (docs/PLAN.md A4 and B4). Reads .rds files, so run it only through
# validation/slurm/export_doubletfinder.sbatch, never on a login node.
#
#   Rscript validation/R/export_doubletfinder.R --object <name> --input <seurat_<object>.rds> \
#     --saved <seurat_doubletfinder_<object>.rds> --out <dir>
#
# <out>/input/: what DoubletFinder 2.0.6 reads with sct = TRUE (paramSweep.R, doubletFinder.R):
#   the RNA assay's "counts" layer (LayerData(seu, assay = "RNA", layer = "counts"), SeuratObject
#   >= 5.0) and, for the workflow's homotypic adjustment, seurat_clusters.
#   matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz  10x (Cell Ranger v3) Matrix Market, entries
#       exactly as stored in the dgCMatrix (explicit zeros included: SoupX objects have them, and
#       sctransform's min_variance, median(umi@x), counts them), cells and features in object
#       order (cells = rownames(meta.data), the order DoubletFinder uses)
#   seurat_clusters.tsv         barcode, seurat_clusters
#   seurat_clusters_levels.txt  the factor's levels in order (modelHomotypic sums in this order)
#   export.json                 provenance, dimensions, value checks and the round-trip result
# <out>/export_run.json: timings, peak RSS, SLURM job and environment fingerprint of the export
# <out>/saved/: the fields of the saved workflow output that are compared (plan B4)
#   cells.tsv     barcode, pANN of every kept cell (the output keeps singlets only), output order
#   doublets.txt  input barcodes minus output barcodes, in input order: the doublet set
#   saved.json    the pANN and DF.classifications column names, and pN, pK and nExp parsed from
#                 them (pK both as R's label and as the number as.numeric() gives)
#
# Values must round-trip exactly: the export re-reads the matrix with the same reader the
# reference runs use (lib/io.R mtx_read_10x) and checks identical(..., num.eq = FALSE), i.e.
# bitwise, against the object's own matrix; it also re-reads the TSVs as text.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))

args <- df_args(list(object = NULL, input = NULL, saved = NULL, out = NULL))
suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
})
t0 <- proc.time()[["elapsed"]]
elapsed <- function() proc.time()[["elapsed"]] - t0
in_dir <- file.path(args$out, "input")
saved_dir <- file.path(args$out, "saved")
dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(saved_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs ----------------------------------------------------------------------------------------
cat("reading", args$input, "\n")
seu <- readRDS(args$input)
read_s <- elapsed()
version <- as.character(SeuratObject::Version(seu))
if (package_version(version) < "5.0") stop("SeuratObject version ", version, " < 5.0: DoubletFinder would use GetAssayData")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
if (!inherits(counts, "dgCMatrix")) stop("counts layer is a ", class(counts)[1], ", not a dgCMatrix")
invisible(validObject(counts))
barcodes <- rownames(seu@meta.data)
if (!identical(colnames(counts), barcodes)) {
  stop("counts columns are not in meta.data order; DoubletFinder would reorder them")
}
if (anyDuplicated(barcodes) || anyDuplicated(rownames(counts))) stop("duplicated barcodes or features")
clusters <- seu@meta.data$seurat_clusters
if (!is.factor(clusters)) stop("seurat_clusters is not a factor")
feature_meta <- tryCatch(seu[["RNA"]][[]], error = function(e) NULL)
x <- counts@x
info_input <- list(
  object = args$object,
  source = normalizePath(args$input),
  source_sha256 = sha256_file(args$input),
  source_bytes = file.size(args$input),
  seurat_object_version = version,
  rna_assay_class = class(seu[["RNA"]])[1],
  rna_layers = Layers(seu[["RNA"]]),
  default_assay = DefaultAssay(seu),
  n_features = nrow(counts),
  n_cells = ncol(counts),
  nnz = length(x),
  explicit_zeros = sum(x == 0),
  non_integer_values = sum(x != round(x)),
  min_value = min(x),
  max_value = max(x),
  feature_metadata_columns = if (is.null(feature_meta)) character() else colnames(feature_meta),
  cluster_levels = levels(clusters),
  cluster_sizes = as.integer(table(clusters))
)
rm(seu)
invisible(gc())

field <- mtx_write_10x(counts, in_dir)
writeLines(levels(clusters), file.path(in_dir, "seurat_clusters_levels.txt"))
tsv_write(data.frame(barcode = barcodes, seurat_clusters = as.character(clusters)),
          file.path(in_dir, "seurat_clusters.tsv"))
write_s <- elapsed()

# Round trip: the exact reader must give back the identical matrix, bit for bit.
back <- mtx_read_10x(in_dir)
same <- identical(back, counts, num.eq = FALSE)
if (!same) {
  # identical() failing on values alone would point at R's decimal parser (non-integer values,
  # see lib/io.R); the structure must match in any case.
  if (!identical(back@i, counts@i) || !identical(back@p, counts@p) || !identical(dimnames(back), dimnames(counts))) {
    stop("round trip changed the matrix structure or dimnames")
  }
  stop("round trip changed ", sum(back@x != counts@x), " values")
}
rm(back)
cl_back <- read.delim(file.path(in_dir, "seurat_clusters.tsv"), colClasses = "character", quote = "")
lev_back <- readLines(file.path(in_dir, "seurat_clusters_levels.txt"))
if (!identical(factor(cl_back$seurat_clusters, levels = lev_back), clusters) ||
    !identical(cl_back$barcode, barcodes)) stop("seurat_clusters round trip failed")
verify_s <- elapsed()
info_input$mtx_field <- field
info_input$round_trip_identical <- same
timing_s <- list(read_rds = read_s, write = write_s - read_s, verify = verify_s - write_s)
input_barcodes <- barcodes
rm(counts, x)
invisible(gc())

# Saved workflow output ------------------------------------------------------------------------
cat("reading", args$saved, "\n")
o <- readRDS(args$saved)
meta <- o@meta.data
rm(o)
invisible(gc())
pann_col <- grep("^pANN_", colnames(meta), value = TRUE)
class_col <- grep("^DF\\.classifications_", colnames(meta), value = TRUE)
if (length(pann_col) != 1 || length(class_col) != 1) stop("expected exactly one pANN and one DF.classifications column")
parts <- strsplit(sub("^pANN_", "", pann_col), "_", fixed = TRUE)[[1]]
if (length(parts) != 3) stop("cannot parse ", pann_col)
if (!identical(sub("^DF\\.classifications_", "", class_col), sub("^pANN_", "", pann_col))) {
  stop("pANN and classification columns disagree: ", pann_col, " vs ", class_col)
}
kept <- rownames(meta)
if (!all(kept %in% input_barcodes)) stop("output has cells that are not in the input")
if (!identical(kept, input_barcodes[input_barcodes %in% kept])) stop("output changed the cell order")
if (!all(meta[[class_col]] == "Singlet")) stop("saved output keeps non-singlets")
doublets <- input_barcodes[!input_barcodes %in% kept]
pann <- meta[[pann_col]]
if (!is.double(pann)) stop("pANN column is not double")
tsv_write(data.frame(barcode = kept, pANN = pann), file.path(saved_dir, "cells.tsv"))
writeLines(doublets, file.path(saved_dir, "doublets.txt"))
back <- read.delim(file.path(saved_dir, "cells.tsv"), colClasses = "character", quote = "")
if (!identical(back$barcode, kept) || !identical(back$pANN, fmt17(pann))) stop("saved cells.tsv round trip failed")
info_saved <- list(
  object = args$object,
  source = normalizePath(args$saved),
  source_sha256 = sha256_file(args$saved),
  source_bytes = file.size(args$saved),
  pANN_column = pann_col,
  classification_column = class_col,
  pN_label = parts[1],
  pK_label = parts[2],
  nExp_label = parts[3],
  pK = as.numeric(parts[2]),
  nExp = as.numeric(parts[3]),
  n_input = length(input_barcodes),
  n_kept = length(kept),
  n_doublets = length(doublets)
)
if (info_saved$n_doublets != info_saved$nExp) stop("doublet count ", info_saved$n_doublets, " != nExp ", info_saved$nExp)

json_write(info_input, file.path(in_dir, "export.json"))
json_write(list(timing_s = c(timing_s, list(total = elapsed())), peak_rss_kib = peak_rss_kib(),
                slurm_job_id = Sys.getenv("SLURM_JOB_ID"), fingerprint = fingerprint()),
           file.path(args$out, "export_run.json"))
json_write(info_saved, file.path(saved_dir, "saved.json"))
cat("exported", args$object, ":", info_input$n_features, "features x", info_input$n_cells, "cells,",
    info_input$nnz, "entries (", info_input$explicit_zeros, "explicit zeros ), kept", info_saved$n_kept,
    ", column", pann_col, ", peak RSS", peak_rss_kib(), "KiB,", round(elapsed()), "s\n")
