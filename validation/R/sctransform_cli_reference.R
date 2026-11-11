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
source(file.path("validation", "R", "lib", "sct_outputs.R"))
out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X(args[[1]])
features <- read.delim(gzfile(list.files(args[[1]], "^features\\.tsv(\\.gz)?$", full.names = TRUE)),
                       header = FALSE, colClasses = "character")
stopifnot(identical(make.unique(features[[2]]), rownames(counts)))
id_of <- setNames(features[[1]], rownames(counts))

captured <- sct_capture()
seu <- CreateSeuratObject(counts = counts)
seu <- SCTransform(seu)
v <- captured$vst_out
stopifnot(!is.null(v), identical(v$variable_features, VariableFeatures(seu)))
sct_write_outputs(v, out, id_of)
cat("written to", normalizePath(out), "\n")
