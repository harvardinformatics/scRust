#!/usr/bin/env Rscript
# Milestone M4 benchmark, R side (docs/PLAN.md B5, B6): Seurat 5.5.1's SCTransform() on one
# DoubletFinder context's SCT input, built from the object's counts (the export in input/) and
# the context's draws in an instrumented run's stages/ (S1), as DoubletFinder 2.0.6 builds it:
#   data <- counts[, real.cells]; doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
#   colnames(doublets) <- paste("X", 1:n_doublets, sep = ""); cbind(data, doublets)
#   seu <- CreateSeuratObject(counts = data_wdoublets); seu <- SCTransform(seu)
# SCTransform() alone is timed (proc.time), and its results are written in the layout of
# `scrust sctransform --from-stages ... --context ...` (lib/sct_outputs.R), so the two compare byte
# for byte under the same OpenBLAS kernel. Also writes timing.json (seconds, peak RSS).
#
# Reads no .rds, but the largest contexts need over 100 GB: run it through
# validation/slurm/m4_benchmark.sbatch.
#   Rscript validation/R/sctransform_benchmark.R <input dir> <stages dir> <context> <out dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: sctransform_benchmark.R <input dir> <stages dir> <context> <out dir>")
stopifnot(R.version$major == "4", R.version$minor == "5.3")
t_start <- proc.time()[["elapsed"]]
suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
    library(Seurat)
})
stopifnot(packageVersion("Seurat") == "5.5.1", packageVersion("sctransform") == "0.4.3",
          packageVersion("glmGamPoi") == "1.22.0")
source(file.path("validation", "R", "lib", "io.R"))
source(file.path("validation", "R", "lib", "sct_outputs.R"))
input <- args[[1]]; stages <- args[[2]]; ctx <- args[[3]]; out <- args[[4]]
if (file.exists(file.path(out, "timing.json"))) stop(out, " holds a previous run")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- mtx_read_10x(input)
features <- read.delim(file.path(input, "features.tsv.gz"), header = FALSE, colClasses = "character",
                       quote = "", comment.char = "")
id_of <- setNames(features[[1]], features[[2]])
index <- read.delim(file.path(stages, "index.tsv"), colClasses = "character")
block <- function(file, context) {
    row <- index[index$file == file & index$ctx == context, ]
    if (nrow(row) != 1) stop("no block ", file, " / ", context)
    a <- npy_read(file.path(stages, file))
    a <- if (is.null(dim(a))) matrix(a, ncol = 1) else a
    a[as.integer(row$row_start) + seq_len(as.integer(row$nrow)), , drop = FALSE]
}
# the sweep's real cells are its subsample when DoubletFinder drew one (over 10,000 cells)
has_subsample <- any(index$file == "S1_subsample.npy" & index$ctx == "sweep")
real_idx <- if (ctx != "final" && has_subsample) block("S1_subsample.npy", "sweep")[, 1] + 1L else seq_len(ncol(counts))
pairs <- block("S1_doublets.npy", ctx)
real.cells <- colnames(counts)[real_idx]
data <- counts[, real.cells]
real.cells1 <- real.cells[pairs[, 1] + 1L]
real.cells2 <- real.cells[pairs[, 2] + 1L]
doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
colnames(doublets) <- paste("X", seq_len(nrow(pairs)), sep = "")
data_wdoublets <- cbind(data, doublets)
rm(counts, data, doublets)
invisible(gc())
t_input <- proc.time()[["elapsed"]]
rss_input <- peak_rss_kib()

captured <- sct_capture()
seu <- CreateSeuratObject(counts = data_wdoublets)
rm(data_wdoublets)
t0 <- proc.time()
seu <- SCTransform(seu)
t1 <- proc.time()
v <- captured$vst_out
stopifnot(!is.null(v), identical(v$variable_features, VariableFeatures(seu)))
rm(seu)
sct_write_outputs(v, out, id_of)
json_write(list(context = ctx, n_cells = nrow(v$cell_attr), n_genes = nrow(v$gene_attr),
                input_s = t_input - t_start, sctransform_elapsed_s = (t1 - t0)[["elapsed"]],
                sctransform_cpu_s = (t1 - t0)[["user.self"]] + (t1 - t0)[["sys.self"]],
                total_s = proc.time()[["elapsed"]] - t_start, peak_rss_after_input_kib = rss_input,
                peak_rss_kib = peak_rss_kib(), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                openblas_coretype = Sys.getenv("OPENBLAS_CORETYPE")),
           file.path(out, "timing.json"))
cat("SCTransform:", round((t1 - t0)[["elapsed"]], 1), "s; written to", normalizePath(out), "\n")
