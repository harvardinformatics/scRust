# Milestone M15b: how doubletfinder.R made its saved percent.mt. The script subsets its input object
# (a clustered object whose default assay is SCT) to the singlets and calls
# PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-") with no assay, so the default
# assay's: the SCT assay's counts layer (the corrected counts) over nCount_SCT, not the RNA counts
# over nCount_RNA. This script computes both on the saved input object, takes the singlets (the
# saved DoubletFinder object's cells, from validation/R/export_percent_mt.R's percent_mt.tsv) and
# compares each with the saved column as %.17g text (so bit for bit). With --rust, also with
# scrust's percent.mt of the same singlets (validation/slurm/m15b_saved_rust.sbatch's re-clustering,
# its metadata.tsv), so the RNA-based value R computes is checked against scrust's. Reads .rds
# files; not a new R run of the workflow. Through validation/slurm/m15b_percent_mt_sct.sbatch:
#   Rscript validation/R/percent_mt_sct.R --rds <saved input .rds> --saved <percent_mt.tsv> --out <json>
#       [--rust <scrust's metadata.tsv>]

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args

args <- df_args(list(rds = NULL, saved = NULL, out = NULL, rust = NULL))
suppressPackageStartupMessages(library("Seurat"))
o <- readRDS(args$rds)
saved <- read.delim(args$saved, colClasses = "character")
cells <- saved$barcode
if (!all(cells %in% colnames(o))) stop("the saved object's cells are not all the input's")
res <- list(rds = args$rds, saved = args$saved, n_cells = length(cells), default_assay = DefaultAssay(o))
rust <- NULL
if (!is.null(args$rust)) {
  rust <- read.delim(args$rust, colClasses = "character")
  if (!identical(rust$barcode, cells)) stop("scrust's cells are not the saved object's, in order")
  res$rust <- args$rust
}
for (assay in c("SCT", "RNA")) {
  pm <- PercentageFeatureSet(o, pattern = "(?i)^mt-", assay = assay)
  txt <- fmt17(unname(pm[cells]))
  res[[assay]] <- list(n_features = length(grep("(?i)^mt-", rownames(o[[assay]]))),
                       identical_to_saved = identical(txt, saved$percent.mt),
                       differing = sum(txt != saved$percent.mt))
  if (!is.null(rust)) {
    res[[assay]]$identical_to_scrust <- identical(txt, rust$percent.mt)
    res[[assay]]$differing_from_scrust <- sum(txt != rust$percent.mt)
  }
}
json_write(res, args$out)
cat(args$rds, ": default assay", res$default_assay, "; SCT identical", res$SCT$identical_to_saved, "(", res$SCT$differing,
    "differ ); RNA identical", res$RNA$identical_to_saved, "(", res$RNA$differing, "differ )",
    if (!is.null(rust)) paste("; scrust: RNA identical", res$RNA$identical_to_scrust, "(", res$RNA$differing_from_scrust, "differ ), SCT identical",
                              res$SCT$identical_to_scrust, "(", res$SCT$differing_from_scrust, "differ )"), "\n")
