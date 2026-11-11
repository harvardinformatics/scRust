#!/usr/bin/env Rscript
# Reference vectors for the sctransform crate and DoubletFinder's SCT input (docs/PLAN.md B3a,
# milestone M4): the installed packages' own code (Matrix 1.7-6, SeuratObject 5.4.0, Seurat
# 5.5.1, sctransform 0.4.3, glmGamPoi 1.22.0) on small synthetic inputs, so that each rule the
# port copies is tested directly, not only through the values it happened to produce on the
# reference data.
#
# Reads no .rds, so it may run on a login node. Run from the repository root with the production
# R env, the hermetic flags and the portable OpenBLAS profile (some glmGamPoi functions call
# BLAS through Armadillo):
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 \
#     "$R_ENV/bin/Rscript" validation/R/sctransform_vectors.R tests/fixtures/sctransform
#
# Output: tab-separated files; lines starting with "#" describe the columns. Doubles are written as
# their IEEE-754 bit patterns (16 hex digits), lists comma-separated, "-" for an empty list.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path("tests", "fixtures", "sctransform")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

stopifnot(R.version$major == "4", R.version$minor == "5.3")
suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
    library(Seurat)
})
stopifnot(packageVersion("Matrix") == "1.7.6", packageVersion("SeuratObject") == "5.4.0",
          packageVersion("Seurat") == "5.5.1", packageVersion("sctransform") == "0.4.3",
          packageVersion("glmGamPoi") == "1.22.0")

hex <- function(x) {
    x <- as.double(x)
    if (!length(x)) return(character())
    b <- writeBin(x, raw(), size = 8, endian = "big")
    apply(matrix(as.character(b), nrow = 8), 2, paste, collapse = "")
}
lst <- function(v) if (length(v)) paste(v, collapse = ",") else "-"
hexl <- function(x) lst(hex(x))

write_tsv <- function(name, header, rows) {
    path <- file.path(out, name)
    writeLines(c(paste0("# ", header), unlist(rows)), path)
    cat(sprintf("%-22s %6d lines\n", name, length(unlist(rows))))
}

# ---------------------------------------------------------------------------------------------
# DoubletFinder's SCT input and the matrix vst() works on, built by the installed code from a
# counts matrix with explicit zeros:
#   data_wdoublets <- cbind(data, (data[, rc1] + data[, rc2]) / 2)        (DoubletFinder 2.0.6)
#   CreateSeuratObject(counts = data_wdoublets)                          (parallel_paramSweep)
#   CreateAssayObject(LayerData(seu, layer = "counts"))                   (SCTransform.StdAssay)
#   as.sparse(GetAssayData(assay, layer = "counts"))                      (SCTransform.Assay/.default)
#   umi[genes, ] with genes <- rownames(umi)[rowSums(umi >= 0.01) >= min_cells]   (vst)
# For each case: the input's slots, the doublet pairs, then i, p and x of the SCT input and of the
# filtered matrix, the kept rows, and median(umi@x), row_gmean and row_var of the filtered matrix.

sparse_case <- function(seed, nr, nc, nd, density, zero_frac, min_cells) {
    set.seed(seed)
    m <- rsparsematrix(nr, nc, density = density, rand.x = function(n) as.double(rpois(n, 3) + 1))
    m <- as(m, "CsparseMatrix")
    # turn some stored entries into explicit zeros (as SoupX's rounding leaves them)
    z <- runif(length(m@x)) < zero_frac
    m@x[z] <- 0
    dimnames(m) <- list(paste0("g", seq_len(nr)), paste0("c", seq_len(nc)))
    real.cells <- colnames(m)
    rc1 <- sample(real.cells, nd, replace = TRUE)
    rc2 <- sample(real.cells, nd, replace = TRUE)
    doublets <- (m[, rc1] + m[, rc2]) / 2
    colnames(doublets) <- paste("X", 1:nd, sep = "")
    dw <- cbind(m, doublets)
    seu <- suppressWarnings(CreateSeuratObject(counts = dw))
    layer <- LayerData(seu, assay = "RNA", layer = "counts")
    a3 <- suppressWarnings(CreateAssayObject(layer))
    umi <- as.sparse(GetAssayData(a3, layer = "counts"))
    stopifnot(identical(umi, dw))
    genes <- rownames(umi)[rowSums(umi >= 0.01) >= min_cells]
    umi2 <- umi[genes, ]
    paste(nr, nc, lst(m@i), lst(m@p), hexl(m@x),
          lst(match(rc1, real.cells) - 1L), lst(match(rc2, real.cells) - 1L),
          lst(dw@i), lst(dw@p), hexl(dw@x),
          lst(match(genes, rownames(umi)) - 1L), lst(umi2@i), lst(umi2@p), hexl(umi2@x),
          hex(median(umi2@x)), hexl(sctransform:::row_gmean(umi2, eps = 1)),
          hexl(sctransform:::row_var(umi2)), sep = "\t")
}

write_tsv("sct_input.tsv",
          paste("nrow\tncol\ti\tp\tx\tpairs a\tpairs b\tSCT input i\tp\tx\tkept rows\tfiltered i\tp\tx",
                "median(umi@x)\trow_gmean\trow_var", sep = "\t"),
          list(sparse_case(1, 40, 25, 9, 0.3, 0.2, 3),
               sparse_case(2, 60, 30, 10, 0.15, 0.5, 5),
               sparse_case(3, 25, 12, 4, 0.6, 0.05, 2),
               sparse_case(4, 80, 40, 14, 0.05, 0.3, 5),
               sparse_case(5, 30, 20, 30, 0.4, 0.0, 1)))

cat("written to", normalizePath(out), "\n")
