#!/usr/bin/env Rscript
# Reference outputs for `scrust sctransform` (docs/PLAN.md B3a, milestone M4): Seurat 5.5.1's
# SCTransform() with its defaults on a 10x counts directory, called as DoubletFinder calls it
# (CreateSeuratObject(counts), then SCTransform(seu)). What SCTransform.default returns is
# captured by a passive trace() and written in the layout of scrust sctransform's outputs
# (sprintf("%.17g") text, npy_write), so the two compare byte for byte: cell_attr.tsv,
# gene_attr.tsv, model_pars.tsv, model_pars_fit.tsv, variable_features.tsv, scale_data.npy (and
# scale_data.npy.sha256, which tests/fixtures/sctransform/cli_<object>/ keeps instead of the
# file).
#
# Reads no .rds, so it may run on a login node for small inputs. Run from the repository root
# with the production R env, the hermetic flags and the OpenBLAS kernel of the comparison (the
# residuals' products are BLAS calls):
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 \
#     "$R_ENV/bin/Rscript" validation/R/sctransform_cli_reference.R <10x dir> <out dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: sctransform_cli_reference.R <10x dir> <out dir>")
stopifnot(R.version$major == "4", R.version$minor == "5.3")
suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
    library(Seurat)
})
stopifnot(packageVersion("Seurat") == "5.5.1", packageVersion("sctransform") == "0.4.3",
          packageVersion("glmGamPoi") == "1.22.0")
source(file.path("validation", "R", "lib", "io.R"))
out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X(args[[1]])
features <- read.delim(gzfile(list.files(args[[1]], "^features\\.tsv(\\.gz)?$", full.names = TRUE)),
                       header = FALSE, colClasses = "character")
stopifnot(identical(make.unique(features[[2]]), rownames(counts)))
id_of <- setNames(features[[1]], rownames(counts))

captured <- new.env()
suppressMessages(trace("SCTransform.default", where = asNamespace("Seurat"), print = FALSE,
                       exit = bquote(assign("vst_out", returnValue(), envir = .(captured)))))
seu <- CreateSeuratObject(counts = counts)
seu <- SCTransform(seu)
v <- captured$vst_out
stopifnot(!is.null(v), identical(v$variable_features, VariableFeatures(seu)))

tsv <- function(name, header, cols) {
    lines <- do.call(paste, c(cols, sep = "\t"))
    writeLines(c(paste(header, collapse = "\t"), lines), file.path(out, name), useBytes = TRUE)
}
ca <- v$cell_attr
tsv("cell_attr.tsv", c("barcode", "umi", "log_umi"), list(rownames(ca), fmt17(ca$umi), fmt17(ca$log_umi)))
ga <- v$gene_attr
tsv("gene_attr.tsv", c("id", "name", "detection_rate", "gmean", "amean", "variance", "residual_mean", "residual_variance"),
    list(id_of[rownames(ga)], rownames(ga), fmt17(ga$detection_rate), fmt17(ga$gmean), fmt17(ga$amean),
         fmt17(ga$variance), fmt17(ga$residual_mean), fmt17(ga$residual_variance)))
pars <- function(name, m) {
    tsv(name, c("id", "name", "theta", "(Intercept)", "log_umi"),
        list(id_of[rownames(m)], rownames(m), fmt17(m[, "theta"]), fmt17(m[, "(Intercept)"]), fmt17(m[, "log_umi"])))
}
pars("model_pars.tsv", v$model_pars)
pars("model_pars_fit.tsv", v$model_pars_fit)
vf <- v$variable_features
tsv("variable_features.tsv", c("id", "name"), list(id_of[vf], vf))
stopifnot(identical(rownames(v$y), vf))
npy_write(file.path(out, "scale_data.npy"), unname(v$y))
# the fixture keeps the digest of the (large) scale.data file rather than the file
writeLines(unname(tools::sha256sum(file.path(out, "scale_data.npy"))), file.path(out, "scale_data.npy.sha256"))
cat("written to", normalizePath(out), "\n")
