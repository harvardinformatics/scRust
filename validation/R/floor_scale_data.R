# The portable profile's SCT scale.data against the saved object's, element by element, for M8's
# R-vs-R floor (docs/validation/clustering_m8.md), where the reference runs did not write
# scale.data in full (dump level light: the three largest tier-P objects). Reads .rds files, so
# run it only through validation/slurm/floor_scale_data.sbatch.
#
#   Rscript validation/R/floor_scale_data.R --object <name> --rds <saved .rds> \
#     --input <export input dir> --out <json>
#
# It runs the chain's first call, SCTransform(seed.use = 12345, verbose = FALSE)
# (lib/clustering.R), on the export's counts under the job's OpenBLAS kernel, then reads the saved
# object's SCT assay and compares: scale.data element by element (values, differing, largest
# absolute difference, largest difference relative to the saved value); the corrected counts
# (identical or not); and the SCTModel's feature attributes column by column. It also records the
# sha256 of the new scale.data's doubles, which must equal the digest the profile's reference run
# recorded (summary.json: scale_data_sha256), so the comparison is the run's own scale.data. On an
# object whose runs kept scale.data in full (soupx_tenx_SRR18256555), its numbers equal the
# element-by-element comparison of the saved and run scale_data.npy (checked in M8).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))

args <- df_args(list(object = NULL, rds = NULL, input = NULL, out = NULL))
cl_setup()
t0 <- proc.time()[["elapsed"]]
assay <- readLines(file.path(args$input, "assay.txt"))
seu <- cl_load_input(args$input, assay)$seurat
seu <- SCTransform(seu, seed.use = WORKFLOW_SEED, verbose = FALSE)
sct <- seu[["SCT"]]
sd_run <- GetAssayData(sct, layer = "scale.data")
cnt_run <- GetAssayData(sct, layer = "counts")
fa_run <- slot(sct, "SCTModel.list")[[1]]@feature.attributes
rm(seu, sct)
invisible(gc())
sct_s <- local({
  o <- readRDS(args$rds)
  s <- o[["SCT"]]
  rm(o)
  s
})
invisible(gc())
sd_saved <- GetAssayData(sct_s, layer = "scale.data")
cnt_saved <- GetAssayData(sct_s, layer = "counts")
fa_saved <- slot(sct_s, "SCTModel.list")[[1]]@feature.attributes
rm(sct_s)
invisible(gc())
if (!identical(dimnames(sd_run), dimnames(sd_saved))) stop("scale.data's rows or columns differ")

# element by element, in column blocks to bound memory
differing <- 0
max_abs <- 0
max_rel <- 0
step <- 5000
for (s in seq(1, ncol(sd_saved), by = step)) {
  cols <- s:min(ncol(sd_saved), s + step - 1)
  a <- sd_saved[, cols, drop = FALSE]
  b <- sd_run[, cols, drop = FALSE]
  d <- a != b
  if (any(d)) {
    differing <- differing + sum(d)
    diff <- abs(a[d] - b[d])
    max_abs <- max(max_abs, diff)
    nz <- a[d] != 0
    if (any(nz)) max_rel <- max(max_rel, diff[nz] / abs(a[d][nz]))
  }
}
# The feature attributes: NA (the step-1 columns of genes outside step 1's subset) equals NA; NA
# against a number counts as a difference; the sizes are over the values that are numbers in both.
if (!identical(rownames(fa_saved), rownames(fa_run))) stop("the feature attributes' genes differ")
fa <- list()
for (col in union(colnames(fa_saved), colnames(fa_run))) {
  x <- fa_saved[[col]]
  y <- fa_run[[col]]
  if (is.null(x) || is.null(y) || !is.numeric(x) || !is.numeric(y)) {
    fa[[col]] <- list(identical = identical(x, y))
    next
  }
  both <- !is.na(x) & !is.na(y)
  dv <- both & x != y
  d <- dv | xor(is.na(x), is.na(y))
  nz <- dv & x != 0
  fa[[col]] <- list(identical = identical(x, y), differing = sum(d), na_mismatches = sum(xor(is.na(x), is.na(y))),
                    max_abs = if (any(dv)) max(abs(x[dv] - y[dv])) else 0,
                    max_rel = if (any(nz)) max(abs(x[nz] - y[nz]) / abs(x[nz])) else 0)
}
res <- list(object = args$object, profile = Sys.getenv("SCRUST_PROFILE"),
            scale_data = list(dim = dim(sd_saved), values = length(sd_saved), differing = differing,
                              max_abs = max_abs, max_rel = max_rel,
                              run_sha256 = raw_sha256(sd_run), saved_sha256 = raw_sha256(sd_saved)),
            corrected_counts_identical = identical(cnt_run, cnt_saved),
            feature_attributes = fa,
            seconds = proc.time()[["elapsed"]] - t0, peak_rss_kib = peak_rss_kib(),
            slurm_job_id = Sys.getenv("SLURM_JOB_ID"), fingerprint = fingerprint(packages = CL_PACKAGES))
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
json_write(res, args$out)
cat(args$object, ": scale.data", differing, "of", length(sd_saved), "values differ, max abs", max_abs,
    "max rel", max_rel, "; corrected counts identical", res$corrected_counts_identical, "\n")
