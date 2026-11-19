# Milestone M15b: percent.mt as the production run saved it (docs/PLAN.md C5's M15b row). Reads one
# saved workflow object (.rds) and writes its cells' percent.mt, with nCount_RNA and nFeature_RNA,
# as %.17g text, for comparison with scrust's percent.mt computed by the workflow's pattern. Every
# workflow script computes it before SCTransform:
#   obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "(?i)^mt-")
# The script also computes it again on the saved object's RNA counts (the same call) and records
# whether that gives the saved column bit for bit, which features the pattern matches (Seurat's row
# names), and the counts' digest against an export of the same counts when one is given. Reading a
# saved object is not a new R run of the workflow. Run through validation/slurm/m15b_saved.sbatch
# (it reads .rds files):
#   Rscript validation/R/export_percent_mt.R --rds <saved .rds> --out <dir> [--counts-from <10x dir>]
# Output: <out>/percent_mt.tsv (barcode, nCount_RNA, nFeature_RNA, percent.mt) and
# <out>/percent_mt.json.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args

args <- df_args(list(rds = NULL, out = NULL, `counts-from` = ""))
suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
})
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
o <- readRDS(args$rds)
md <- o@meta.data
if (!"percent.mt" %in% colnames(md)) stop(args$rds, ": no percent.mt column")
if (!identical(rownames(md), colnames(o))) stop("meta.data rows differ from the cells")
tsv_write(data.frame(barcode = colnames(o), nCount_RNA = md$nCount_RNA, nFeature_RNA = md$nFeature_RNA,
                     percent.mt = md$percent.mt, stringsAsFactors = FALSE, check.names = FALSE),
          file.path(args$out, "percent_mt.tsv"))
again <- unname(PercentageFeatureSet(o, pattern = "(?i)^mt-", assay = "RNA"))
features <- grep("(?i)^mt-", rownames(o[["RNA"]]), value = TRUE)
facts <- list(
  rds = args$rds,
  n_cells = ncol(o),
  rna_assay_class = class(o[["RNA"]])[1],
  default_assay = DefaultAssay(o),
  percent_mt_class = class(md$percent.mt)[1],
  recomputed_identical = identical(again, unname(md$percent.mt)),
  pattern = "(?i)^mt-",
  n_features = length(features),
  features = I(features),
  cells_zero = sum(md$percent.mt == 0),
  cells_nan = sum(is.nan(md$percent.mt)),
  max = if (ncol(o)) max(md$percent.mt) else NA
)
if (nzchar(args$`counts-from`)) {
  cnt <- LayerData(o, assay = "RNA", layer = "counts")
  ex <- mtx_read_10x(args$`counts-from`)
  facts$counts_from <- args$`counts-from`
  facts$counts_identical_to_export <- identical(unname(cnt), unname(ex)) && identical(colnames(cnt), colnames(ex))
}
json_write(facts, file.path(args$out, "percent_mt.json"))
cat(args$rds, ":", ncol(o), "cells,", length(features), "features match the pattern; recomputed identical:",
    facts$recomputed_identical, "\n")
