#!/usr/bin/env Rscript
# Seurat 5.5.1's RunPCA() on an SCTransform output directory in `scrust sctransform`'s layout
# (docs/PLAN.md B3b, milestone M5), written in the layout of `scrust pca`, so the two compare byte
# for byte under the same OpenBLAS kernel: the reference for `scrust pca` (tests/fixtures/pca/
# cli_<object>/) and the R side of the M5 benchmark (validation/slurm/m5_benchmark.sbatch).
#
# The input: scale_data.npy (variable features x cells), variable_features.tsv (id, name; the
# rows, in VariableFeatures order) and cell_attr.tsv (barcode first). They become a Seurat object
# whose default assay "SCT" is a v3 Assay holding that scale.data and those VariableFeatures, as
# DoubletFinder's SCTransform() leaves it, so RunPCA() takes DoubletFinder's path:
# RunPCA.Seurat -> RunPCA.Assay -> PrepDR -> RunPCA.default -> irlba (an SCTAssay is an Assay,
# and Seurat has no RunPCA method of its own for it). RunPCA() runs with its defaults and
# npcs (DoubletFinder: 10), and is timed alone (proc.time).
#
# Outputs: embeddings.tsv, loadings.tsv, stdev.tsv, random_seed.txt (as `scrust pca` writes
# them; sprintf("%.17g") text), timing.json (seconds, peak RSS, total.variance, irlba's d, iter
# and mprod).
#
# Reads no .rds; the largest inputs (2.5 GB) need about 10 GB: run them through sbatch. From the
# repository root, with the production R env, the hermetic flags and the comparison's kernel:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 \
#     "$R_ENV/bin/Rscript" validation/R/runpca_reference.R <sct dir> <out dir> [npcs]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:3) stop("usage: runpca_reference.R <sct dir> <out dir> [npcs]")
stopifnot(R.version$major == "4", R.version$minor == "5.3")
t_start <- proc.time()[["elapsed"]]
suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
    library(Seurat)
})
stopifnot(packageVersion("Seurat") == "5.5.1", packageVersion("SeuratObject") == "5.4.0",
          packageVersion("irlba") == "2.4.1")
source(file.path("validation", "R", "lib", "io.R"))
sct <- args[[1]]; out <- args[[2]]
npcs <- if (length(args) == 3) as.integer(args[[3]]) else 50L
if (file.exists(file.path(out, "timing.json"))) stop(out, " holds a previous run")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

vf <- read.delim(file.path(sct, "variable_features.tsv"), colClasses = "character", quote = "", comment.char = "")
cells <- read.delim(file.path(sct, "cell_attr.tsv"), colClasses = "character", quote = "", comment.char = "")$barcode
sd <- npy_read(file.path(sct, "scale_data.npy"))
stopifnot(identical(dim(sd), c(nrow(vf), length(cells))))
dimnames(sd) <- list(vf$name, cells)
id_of <- setNames(vf$id, vf$name)
# a v3 assay: empty counts of the right shape, then scale.data and the variable features
empty <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = dim(sd), dimnames = dimnames(sd))
assay <- suppressWarnings(CreateAssayObject(counts = empty))
assay <- SetAssayData(assay, layer = "scale.data", new.data = sd)
VariableFeatures(assay) <- vf$name
rm(sd, empty)
seu <- CreateSeuratObject(counts = assay, assay = "SCT")
rm(assay)
stopifnot(inherits(seu[["SCT"]], "Assay"), !inherits(seu[["SCT"]], "Assay5"), DefaultAssay(seu) == "SCT")
invisible(gc())
t_input <- proc.time()[["elapsed"]]
rss_input <- peak_rss_kib()

captured <- new.env()
suppressMessages(trace("irlba", where = asNamespace("irlba"), print = FALSE,
                       exit = bquote(assign("out", returnValue(), envir = .(captured)))))
t0 <- proc.time()
seu <- RunPCA(seu, npcs = npcs)
t1 <- proc.time()
rng_after <- .Random.seed
untrace("irlba", where = asNamespace("irlba"))

pca <- seu[["pca"]]
emb <- Embeddings(pca)
lod <- Loadings(pca)
tsv <- function(name, header, cols) {
    lines <- do.call(paste, c(cols, sep = "\t"))
    writeLines(c(paste(header, collapse = "\t"), lines), file.path(out, name), useBytes = TRUE)
}
tsv("embeddings.tsv", c("barcode", colnames(emb)), c(list(rownames(emb)), lapply(seq_len(ncol(emb)), function(k) fmt17(emb[, k]))))
tsv("loadings.tsv", c("id", "name", colnames(lod)),
    c(list(unname(id_of[rownames(lod)]), rownames(lod)), lapply(seq_len(ncol(lod)), function(k) fmt17(lod[, k]))))
tsv("stdev.tsv", c("pc", "stdev"), list(colnames(emb), fmt17(Stdev(pca))))
writeLines(as.character(rng_after), file.path(out, "random_seed.txt"))
ir <- captured$out
json_write(list(n_cells = nrow(emb), n_features = nrow(vf), n_pca_features = nrow(lod), npcs = ncol(emb),
                total_variance = Misc(pca, "total.variance"), irlba_d = ir$d, irlba_iter = ir$iter,
                irlba_mprod = ir$mprod, input_s = t_input - t_start,
                runpca_elapsed_s = (t1 - t0)[["elapsed"]],
                runpca_cpu_s = (t1 - t0)[["user.self"]] + (t1 - t0)[["sys.self"]],
                total_s = proc.time()[["elapsed"]] - t_start, peak_rss_after_input_kib = rss_input,
                peak_rss_kib = peak_rss_kib(), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                openblas_coretype = Sys.getenv("OPENBLAS_CORETYPE")),
           file.path(out, "timing.json"))
cat("RunPCA:", round((t1 - t0)[["elapsed"]], 1), "s; written to", normalizePath(out), "\n")
