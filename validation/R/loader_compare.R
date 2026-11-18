# K11 (docs/PLAN.md C3g, C4 "Loader tests", C5; milestone M14): the R loaders' objects against
# the saved object. Builds a Seurat object (scrustio::load_seurat) and a SingleCellExperiment
# (scrustio::load_sce) from one scrust output and compares each with the saved workflow object
# (.rds), component by component: the counts, the SCT layers and model, the chain's metadata
# columns, the identities, the PCA and the graphs. Also M17's comparison (C4 "Full panel"):
# with --loaders seurat --files no --round-trip no it compares one export with its saved .rds
# and writes the tables, after which the export's large files can go. Reads .rds files: run it
# only through sbatch (validation/slurm/m14_loaders.sbatch, m14_m17.sbatch).
#
#   Rscript validation/R/loader_compare.R --object <name> --export <scrust output> --rds <saved .rds> \
#     --out <dir> [--loaders seurat,sce] [--files yes|no] [--round-trip yes|no] \
#     [--kind cluster|doubletfinder] [--tmp <dir>]
#
# Grading (as M11-M13 grade, docs/validation/clustering.md "Kernel-gated grading"): exact when the
# export's recorded OpenBLAS kernel (summary.json run.fingerprint.openblas_corename) is the saved
# objects' (SapphireRapids, every workflow job's); otherwise ("floor") the BLAS-free components
# stay exact and the BLAS-dependent ones may differ within [clustering.floor.kernel] of
# validation/tolerances.toml (scale.data and the model's residual statistics, S2d; the PCA, S3,
# compared after aligning each component's sign; the silhouette widths, K8; the purities, K9);
# everything the floor does not list (the counts, the corrected counts and their log1p, the
# labels, levels and identities, the graphs, the total variance) stays exact.
#
# --files yes: the loaders load exactly what the files hold. Every value the Seurat loader read
# from a text file is formatted again with %.17g and compared with the file's text (%.17g is
# one-to-one on doubles, so equal text is the same double), entry by entry for the Matrix Market
# files, column by column for the tables, and .npy files byte for byte; the summary's numbers by
# their shortest decimal form.
# --round-trip yes: scrustio::save_rds() (the one-line command) writes each object as .rds, and
# readRDS() must give an object identical to the loader's.
#
# --kind doubletfinder: the output of scrust doubletfinder against the saved DoubletFinder
# object: the singlets' raw counts and DoubletFinder's two metadata columns.
#
# Output in <dir>: k11.tsv (one row per loader and component: status identical, within floor,
# beyond floor, differing or not rebuilt, with the statistics), files.tsv, roundtrip.tsv and
# compare.json (mode, kernels, counts by status, timings and peak RSS per step).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "loaders.R"))
root <- normalizePath(file.path(here, "..", ".."))
args <- df_args(list(object = NULL, export = NULL, rds = NULL, out = NULL, loaders = "seurat,sce",
                     files = "yes", `round-trip` = "yes", kind = "cluster", tmp = tempdir(),
                     tolerances = file.path(root, "validation", "tolerances.toml"),
                     `saved-kernel` = "SapphireRapids"))
ld_setup(root)
loaders <- strsplit(args$loaders, ",", fixed = TRUE)[[1]]
if ("sce" %in% loaders) suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("Seurat"))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

t0 <- proc.time()[["elapsed"]]
steps <- list()
mark <- function(name) steps[[name]] <<- list(elapsed_s = proc.time()[["elapsed"]] - t0, peak_rss_kib = ld_vmhwm_kib())

m <- read_manifest(args$export)
export_kernel <- m$summary$run$fingerprint$openblas_corename
if (is.null(export_kernel)) export_kernel <- NA_character_
mode <- if (identical(export_kernel, args$`saved-kernel`)) "exact" else "floor"
floor <- ld_floor(args$tolerances)
mark("manifest")

rows <- list()
add <- function(loader, component, ident, stats = NULL, floor_name = "", blas = TRUE, detail = "", status = NULL) {
  fe <- if (nzchar(floor_name)) floor[[floor_name]] else NULL
  if (nzchar(floor_name) && is.null(fe)) stop("no floor entry ", floor_name)
  st <- if (!is.null(status)) status else ld_grade(ident, stats, if (blas) fe else NULL, mode)
  g <- function(k) if (is.null(stats[[k]])) NA else stats[[k]]
  rows[[length(rows) + 1]] <<- data.frame(
    object = args$object, loader = loader, component = component,
    class = if (blas) "blas_dependent" else "blas_free", floor = floor_name, status = st,
    n = g("n"), differing = g("differing"), max_abs = g("max_abs"), max_scaled = g("max_scaled"),
    detail = gsub("[\t\n]", " ", detail), stringsAsFactors = FALSE)
}

# A whole object (the SCT assay, its model, the PCA): identical, or, in floor mode, graded as its
# parts were (rows of this loader whose component starts with one of `prefixes`).
add_whole <- function(loader, component, ident, prefixes) {
  if (isTRUE(ident) || mode != "floor") return(add(loader, component, ident))
  parts <- Filter(function(r) r$loader == loader && any(startsWith(r$component, prefixes)), rows)
  st <- vapply(parts, `[[`, "", "status")
  status <- if (all(st %in% c("identical", "within floor"))) "within floor" else "differing"
  add(loader, component, ident, status = status, detail = paste("graded by its", length(parts), "parts"))
}

CHAIN <- c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "seurat_clusters",
           "silhouette_width", "neighborhood_purity")
COL_FLOOR <- c(silhouette_width = "K8_silhouette_width", neighborhood_purity = "K9_purity")
COL_BLAS_FREE <- c("nCount_RNA", "nFeature_RNA")

# The metadata columns, the identities -----------------------------------------------------------
cmp_column <- function(L, col, a, b) {
  ident <- identical(a, b)
  if (is.factor(a) || is.factor(b)) {
    d <- if (is.factor(a) && is.factor(b) && identical(levels(a), levels(b))) "" else "levels differ"
    add(L, paste("metadata", col), ident, list(n = length(b), differing = sum(as.character(a) != as.character(b))),
        detail = d, blas = !col %in% COL_BLAS_FREE)
  } else {
    fl <- if (col %in% names(COL_FLOOR)) COL_FLOOR[[col]] else ""
    st <- if (length(a) == length(b) && is.numeric(a) && is.numeric(b)) ld_num_stats(a, b) else NULL
    d <- if (!identical(typeof(a), typeof(b))) paste("types", typeof(a), typeof(b)) else ""
    add(L, paste("metadata", col), ident, st, floor_name = fl, blas = !col %in% COL_BLAS_FREE, detail = d)
  }
}

# The SCTModel's slots, from a list (SCE) or the S4 object (Seurat) ------------------------------
cmp_model <- function(L, fa, ca, other, m2) {
  fb <- m2@feature.attributes
  ra <- rownames(fa); rb <- rownames(fb)
  add(L, "SCTModel feature.attributes rows and columns", identical(ra, rb) && identical(colnames(fa), colnames(fb)),
      blas = FALSE, detail = paste(ncol(fa), "columns"))
  # the floor's S2d_residual_stats is one array of both columns, cbind(residual_mean, residual_variance)
  rcols <- c("residual_mean", "residual_variance")
  rs <- if (all(rcols %in% colnames(fa)) && nrow(fa) == nrow(fb)) {
    ld_num_stats(as.matrix(fa[, rcols]), as.matrix(fb[, rcols]))
  } else NULL
  for (col in colnames(fb)) {
    a <- fa[[col]]; b <- fb[[col]]
    fl <- if (col %in% rcols) "S2d_residual_stats" else ""
    blas <- !col %in% c("detection_rate", "gmean", "variance")
    st <- if (col %in% rcols && !is.null(rs)) rs else if (is.numeric(a) && is.numeric(b) && length(a) == length(b)) ld_num_stats(a, b) else list(n = length(b))
    add(L, paste("SCTModel feature.attributes", col), identical(a, b), st, floor_name = fl, blas = blas,
        detail = if (col %in% rcols) "statistics over residual_mean and residual_variance together, as the floor's S2d_residual_stats" else "")
  }
  cb <- m2@cell.attributes
  add(L, "SCTModel cell.attributes", identical(as.data.frame(ca), cb), list(n = nrow(cb)), blas = FALSE,
      detail = paste(colnames(cb), collapse = ","))
  # cells_step1 on its own (also within cell.attributes): the loaders draw it again (M14)
  if ("cells_step1" %in% colnames(cb)) {
    a <- ca[["cells_step1"]]; b <- cb[["cells_step1"]]
    add(L, "SCTModel cells_step1", identical(a, b),
        list(n = length(b), differing = if (length(a) == length(b)) sum(a != b) else NA), blas = FALSE,
        detail = paste0(sum(b), " of ", length(b), " cells in the saved model: ",
                        if (all(b)) "every cell (at most ncells, no draw)" else "drawn again from the manifest's seed and ncells"))
  }
  for (s in c("clips", "umi.assay", "model", "arguments", "median_umi")) {
    add(L, paste("SCTModel", s), identical(other[[s]], methods::slot(m2, s)), blas = s %in% c("clips", "arguments"))
  }
}

# Seurat ------------------------------------------------------------------------------------------
cmp_seurat <- function(obj, o) {
  L <- "seurat"
  a <- LayerData(obj, assay = "RNA", layer = "counts"); b <- LayerData(o, assay = "RNA", layer = "counts")
  add(L, "RNA counts", identical(a, b), ld_sparse_stats(a, b), blas = FALSE)
  if (inherits(o[["RNA"]], "Assay5")) {
    add(L, "RNA assay (Assay5, as a whole)", identical(obj[["RNA"]], o[["RNA"]]), blas = FALSE)
  } else {
    add(L, "RNA assay (Assay5, as a whole)", NA, blas = FALSE, status = "not rebuilt",
        detail = paste0("the saved RNA assay is a v3 ", class(o[["RNA"]])[1], " (emptyDrops path: as.Seurat()), with ",
                        "its own data layer and meta.features; the loader makes the Assay5 CreateSeuratObject() makes"))
  }
  add(L, "metadata orig.ident", identical(obj$orig.ident, o$orig.ident), blas = FALSE)
  for (col in CHAIN) cmp_column(L, col, obj@meta.data[[col]], o@meta.data[[col]])
  add(L, "identities (active.ident)", identical(obj@active.ident, o@active.ident))
  add(L, "default assay", identical(DefaultAssay(obj), DefaultAssay(o)), blas = FALSE)
  s1 <- obj[["SCT"]]; s2 <- o[["SCT"]]
  add(L, "SCT counts", identical(s1@counts, s2@counts), ld_sparse_stats(s1@counts, s2@counts))
  add(L, "SCT data", identical(s1@data, s2@data), ld_sparse_stats(s1@data, s2@data))
  sd_names <- identical(dimnames(s1@scale.data), dimnames(s2@scale.data))
  add(L, "SCT scale.data rows and columns", sd_names, blas = FALSE,
      detail = "rows in the SCT assay's gene order, as Seurat keeps them")
  add(L, "SCT scale.data", identical(s1@scale.data, s2@scale.data),
      if (sd_names) ld_num_stats(s1@scale.data, s2@scale.data) else NULL, floor_name = "S2d_scale_data")
  add(L, "SCT variable features", identical(s1@var.features, s2@var.features), blas = FALSE)
  same_slots <- identical(class(s1), class(s2)) && identical(s1@key, s2@key) && identical(s1@assay.orig, s2@assay.orig) &&
    identical(s1@meta.features, s2@meta.features) && identical(s1@misc, s2@misc)
  add(L, "SCT assay class, key, assay.orig, meta.features, misc", same_slots, blas = FALSE)
  add(L, "SCT models (one, model1)", identical(names(s1@SCTModel.list), names(s2@SCTModel.list)), blas = FALSE,
      detail = paste(length(s1@SCTModel.list), "loaded,", length(s2@SCTModel.list), "saved"))
  m1 <- s1@SCTModel.list[[1]]; m2 <- s2@SCTModel.list[[1]]
  cmp_model(L, m1@feature.attributes, m1@cell.attributes,
            list(clips = m1@clips, umi.assay = m1@umi.assay, model = m1@model, arguments = m1@arguments,
                 median_umi = m1@median_umi), m2)
  add_whole(L, "SCTModel (as a whole)", identical(m1, m2), "SCTModel ")
  add_whole(L, "SCT assay (as a whole)", identical(s1, s2), c("SCT ", "SCTModel "))
  p1 <- obj[["pca"]]; p2 <- o[["pca"]]
  e1 <- Embeddings(p1); e2 <- Embeddings(p2)
  add(L, "pca embeddings rows and columns", identical(dimnames(e1), dimnames(e2)), blas = FALSE)
  nd <- ncol(e2); k30 <- seq_len(min(30, nd))
  add(L, "pca embeddings, dims 1:30", identical(e1[, k30], e2[, k30]), ld_num_stats(e1[, k30], e2[, k30], align = TRUE),
      floor_name = "S3_embeddings_dims_1_30")
  if (nd > 30) {
    k <- 31:nd
    add(L, "pca embeddings, dims 31:50", identical(e1[, k], e2[, k]), ld_num_stats(e1[, k], e2[, k], align = TRUE),
        floor_name = "S3_embeddings_dims_31_50")
  }
  l1 <- Loadings(p1); l2 <- Loadings(p2)
  add(L, "pca loadings rows and columns", identical(dimnames(l1), dimnames(l2)), blas = FALSE)
  add(L, "pca loadings", identical(l1, l2), ld_num_stats(l1, l2, align = TRUE), floor_name = "S3_loadings")
  add(L, "pca stdev", identical(p1@stdev, p2@stdev), ld_num_stats(p1@stdev, p2@stdev), floor_name = "S3_sdev")
  add(L, "pca total.variance (misc)", identical(p1@misc, p2@misc),
      ld_num_stats(p1@misc$total.variance, p2@misc$total.variance))
  same <- identical(p1@key, p2@key) && identical(p1@assay.used, p2@assay.used) && identical(p1@global, p2@global) &&
    identical(p1@jackstraw, p2@jackstraw) && identical(p1@feature.loadings.projected, p2@feature.loadings.projected)
  add(L, "pca key, assay.used, global, jackstraw, projected", same, blas = FALSE)
  add_whole(L, "pca (as a whole)", identical(p1, p2), "pca ")
  for (g in c("SCT_nn", "SCT_snn")) {
    add(L, paste("graph", g), identical(obj[[g]], o[[g]]), ld_sparse_stats(obj[[g]], o[[g]]),
        detail = paste("assay.used", obj[[g]]@assay.used))
  }
  add(L, "project.name, SeuratObject version", identical(obj@project.name, o@project.name) &&
        identical(obj@version, o@version), blas = FALSE)
  # what the loader does not rebuild
  if ("umap" %in% names(o@reductions)) add(L, "umap reduction", NA, status = "not rebuilt", blas = FALSE,
                                             detail = "RunUMAP is not part of the chain scrust reproduces (A7b; it changes nothing downstream, M8)")
  if (length(o@commands)) add(L, "commands", NA, status = "not rebuilt", blas = FALSE,
                              detail = paste0("SeuratCommand records of the R calls (", paste(names(o@commands), collapse = ", "), "); scrust records its run in summary.json"))
  extra <- setdiff(colnames(o@meta.data), colnames(obj@meta.data))
  if (length(extra)) add(L, "other metadata columns", NA, status = "not rebuilt", blas = FALSE,
                         detail = paste0("columns of other workflow steps, not of the chain: ", paste(extra, collapse = ", ")))
  added <- setdiff(colnames(obj@meta.data), colnames(o@meta.data))
  if (length(added)) add(L, "metadata columns the export adds", NA, status = "added", blas = FALSE,
                         detail = paste0("bluster's other outputs, which the workflow drops: ", paste(added, collapse = ", ")))
  for (s in c("tools", "misc", "neighbors", "images")) {
    if (length(methods::slot(o, s))) add(L, paste("slot", s), NA, status = "not rebuilt", blas = FALSE,
                                         detail = paste(names(methods::slot(o, s)), collapse = ", "))
  }
}

# SingleCellExperiment ------------------------------------------------------------------------------
cmp_sce <- function(sce, o) {
  L <- "sce"
  a <- assay(sce, "counts"); b <- LayerData(o, assay = "RNA", layer = "counts")
  add(L, "counts", identical(a, b), ld_sparse_stats(a, b), blas = FALSE)
  cd <- as.data.frame(colData(sce))
  add(L, "colData rows (cells)", identical(rownames(cd), colnames(o)), blas = FALSE)
  for (col in CHAIN) cmp_column(L, col, cd[[col]], o@meta.data[[col]])
  sct <- altExp(sce, "SCT"); s2 <- o[["SCT"]]
  a <- assay(sct, "corrected_counts"); add(L, "altExp SCT corrected_counts", identical(a, s2@counts), ld_sparse_stats(a, s2@counts))
  a <- assay(sct, "logcounts"); add(L, "altExp SCT logcounts", identical(a, s2@data), ld_sparse_stats(a, s2@data))
  sd <- assay(altExp(sce, "SCT_scale_data"), "scaledata")
  vf <- metadata(sct)$variable_features
  add(L, "altExp SCT variable features (metadata)", identical(vf, VariableFeatures(s2)), blas = FALSE)
  same_rows <- setequal(rownames(sd), rownames(s2@scale.data)) && identical(colnames(sd), colnames(s2@scale.data))
  add(L, "altExp SCT_scale_data rows (VariableFeatures order) and columns", same_rows && identical(rownames(sd), VariableFeatures(s2)),
      blas = FALSE)
  b <- s2@scale.data[rownames(sd), , drop = FALSE]
  add(L, "altExp SCT_scale_data scaledata", identical(sd, b), if (same_rows) ld_num_stats(sd, b) else NULL, floor_name = "S2d_scale_data")
  md <- metadata(sct)$SCTModel
  fa <- as.data.frame(rowData(sct), optional = TRUE)
  ca <- as.data.frame(colData(sct), optional = TRUE)
  cmp_model(L, fa, ca, md, s2@SCTModel.list[[1]])
  p2 <- o[["pca"]]
  e <- reducedDim(sce, "PCA")
  e1 <- e; attributes(e1) <- attributes(e)[c("dim", "dimnames")]
  e2 <- Embeddings(p2)
  add(L, "reducedDim PCA rows and columns", identical(dimnames(e1), dimnames(e2)), blas = FALSE)
  k30 <- seq_len(min(30, ncol(e2)))
  add(L, "reducedDim PCA, dims 1:30", identical(e1[, k30], e2[, k30]), ld_num_stats(e1[, k30], e2[, k30], align = TRUE),
      floor_name = "S3_embeddings_dims_1_30")
  if (ncol(e2) > 30) {
    k <- 31:ncol(e2)
    add(L, "reducedDim PCA, dims 31:50", identical(e1[, k], e2[, k]), ld_num_stats(e1[, k], e2[, k], align = TRUE),
        floor_name = "S3_embeddings_dims_31_50")
  }
  rot <- attr(e, "rotation")
  add(L, "reducedDim PCA rotation", identical(rot, Loadings(p2)), ld_num_stats(rot, Loadings(p2), align = TRUE), floor_name = "S3_loadings")
  add(L, "reducedDim PCA stdev", identical(attr(e, "stdev"), p2@stdev), ld_num_stats(attr(e, "stdev"), p2@stdev), floor_name = "S3_sdev")
  add(L, "reducedDim PCA total.variance", identical(attr(e, "total.variance"), p2@misc$total.variance))
  for (g in c("SCT_nn", "SCT_snn")) {
    a <- colPair(sce, g, asSparse = TRUE)
    b <- methods::as(o[[g]], "dgCMatrix"); dimnames(b) <- NULL
    a <- methods::as(a, "CsparseMatrix")
    ident <- identical(dim(a), dim(b)) && identical(a@i, b@i) && identical(a@p, b@p) && identical(as.double(a@x), b@x)
    add(L, paste("colPair", g), ident, ld_sparse_stats(a, b))
  }
}

# DoubletFinder --------------------------------------------------------------------------------------
cmp_doubletfinder <- function(obj, o, L) {
  if (L == "seurat") {
    a <- LayerData(obj, assay = "RNA", layer = "counts"); md <- obj@meta.data
  } else {
    a <- assay(obj, "counts"); md <- as.data.frame(colData(obj), optional = TRUE)
  }
  b <- LayerData(o, assay = "RNA", layer = "counts")
  add(L, "singlets' raw counts", identical(a, b), ld_sparse_stats(a, b), blas = FALSE)
  res <- m$summary$result
  for (col in c(res$pANN_column, res$classification_column)) {
    if (!col %in% colnames(o@meta.data)) {
      add(L, paste("metadata", col), FALSE, detail = "the saved object has no such column")
      next
    }
    x <- md[[col]]; y <- o@meta.data[[col]]
    st <- if (is.numeric(x)) ld_num_stats(x, y) else list(n = length(y), differing = sum(x != y))
    add(L, paste("metadata", col), identical(x, y), st)
  }
  if (L == "seurat") {
    for (col in c("nCount_RNA", "nFeature_RNA")) cmp_column(L, col, md[[col]], o@meta.data[[col]])
  }
}

# The files (--files yes): text identity of every value the Seurat loader read ----------------------
frows <- list()
fadd <- function(file, component, n, differing, detail = "") {
  frows[[length(frows) + 1]] <<- data.frame(object = args$object, file = file, component = component, values = n,
                                            differing = differing, status = if (differing == 0) "identical" else "differing",
                                            detail = detail, stringsAsFactors = FALSE)
}

# A Matrix Market file against a loaded dgCMatrix, entry by entry in stored order, in pieces.
files_mm <- function(role, mat, label) {
  e <- manifest_entries(m, role)[[1]]
  path <- file.path(m$dir, e$path)
  con <- gzfile(path, "r"); on.exit(close(con))
  line <- readLines(con, n = 1)
  while (startsWith(line, "%")) line <- readLines(con, n = 1)
  nnz <- as.numeric(strsplit(trimws(line), " +")[[1]])[3]
  mat <- methods::as(mat, "dgCMatrix")
  jj <- rep.int(seq_len(ncol(mat)), diff(mat@p))
  done <- 0; bad <- 0
  while (done < nnz) {
    k <- min(5e6, nnz - done)
    tr <- scan(con, what = list("", "", ""), nmax = k, quiet = TRUE)
    idx <- (done + 1):(done + k)
    bad <- bad + sum(tr[[1]] != as.character(mat@i[idx] + 1L) | tr[[2]] != as.character(jj[idx]) | tr[[3]] != fmt17(mat@x[idx]))
    done <- done + k
  }
  fadd(e$path, label, nnz, bad, "row, column and value of every entry, in stored order")
}

# Columns of a TSV file (read as text) against loaded numeric vectors, formatted with %.17g.
files_tsv <- function(role, cols, label) {
  e <- manifest_entries(m, role)[[1]]
  df <- utils::read.delim(file.path(m$dir, e$path), colClasses = "character", quote = "", comment.char = "",
                          check.names = FALSE, na.strings = character())
  n <- 0; bad <- 0
  for (nm in names(cols)) {
    n <- n + length(cols[[nm]])
    bad <- bad + sum(df[[nm]] != fmt17(cols[[nm]]))
  }
  fadd(e$path, label, n, bad, paste(names(cols), collapse = ", "))
}

# A summary.json number (one key per line, as scrust writes it) against the loaded double, which
# jsonlite parsed with C's strtod: R's own parser must read the same double from the text.
files_json_number <- function(key, x, label) {
  txt <- readLines(file.path(m$dir, "summary.json"), warn = FALSE)
  hits <- grep(paste0("\"", key, "\": *-?[0-9]"), txt, value = TRUE)
  vals <- as.numeric(sub(".*: *(-?[0-9.eE+-]+),?[[:space:]]*$", "\\1", hits))
  bad <- if (length(vals)) sum(vals != x) else 1
  fadd("summary.json", label, max(1, length(vals)), bad, paste("result's", key, "parsed by R and by jsonlite"))
}

files_check <- function(obj) {
  if (args$kind == "doubletfinder") {
    files_mm("raw_counts", LayerData(obj, assay = "RNA", layer = "counts"), "singlets' raw counts")
    calls <- read.delim(file.path(m$dir, manifest_entries(m, "doublet_calls")[[1]]$path), colClasses = "character")
    keep <- match(colnames(obj), calls$barcode)
    col <- m$summary$result$pANN_column
    fadd(manifest_entries(m, "doublet_calls")[[1]]$path, "pANN", length(keep), sum(calls$pANN[keep] != fmt17(obj@meta.data[[col]])))
    return(invisible())
  }
  files_mm("raw_counts", LayerData(obj, assay = "RNA", layer = "counts"), "RNA counts")
  s <- obj[["SCT"]]
  files_mm("sct_counts", s@counts, "SCT counts")
  files_mm("sct_data", s@data, "SCT data")
  files_mm("graph_nn", obj[["SCT_nn"]], "graph SCT_nn")
  files_mm("graph_snn", obj[["SCT_snn"]], "graph SCT_snn")
  # scale_data.npy: rows in VariableFeatures order in the file
  e <- manifest_entries(m, "sct_scale_data")[[1]]
  sd <- s@scale.data[s@var.features, , drop = FALSE]
  con <- file(file.path(m$dir, e$path), "rb")
  readBin(con, "raw", 8); hl <- readBin(con, "integer", size = 2, signed = FALSE, endian = "little"); readBin(con, "raw", hl)
  n <- length(sd); done <- 0; bad <- 0
  while (done < n) {
    k <- min(2^26, n - done)
    v <- readBin(con, "double", k, size = 8, endian = "little")
    bad <- bad + sum(v != sd[(done + 1):(done + k)] | is.na(v) != is.na(sd[(done + 1):(done + k)]))
    done <- done + k
  }
  close(con)
  fadd(e$path, "SCT scale.data (VariableFeatures order)", n, bad, "every double, as bytes")
  mdl <- s@SCTModel.list[[1]]
  fa <- mdl@feature.attributes
  files_tsv("sct_feature_attributes", as.list(fa[, c("detection_rate", "gmean", "variance", "residual_mean", "residual_variance")]),
            "SCTModel feature.attributes")
  files_tsv("sct_model_pars_fit", list(theta = fa$theta, `(Intercept)` = fa$`(Intercept)`, log_umi = fa$log_umi),
            "SCTModel theta, (Intercept), log_umi")
  e1 <- manifest_entries(m, "sct_model_pars")[[1]]
  st <- utils::read.delim(file.path(m$dir, e1$path), colClasses = "character", check.names = FALSE, quote = "")
  rn <- seurat_feature_names(st$name, warn = FALSE)
  fadd(e1$path, "SCTModel step1_theta, step1_(Intercept), step1_log_umi", 3 * nrow(st),
       sum(st$theta != fmt17(fa[rn, "step1_theta"])) + sum(st$`(Intercept)` != fmt17(fa[rn, "step1_(Intercept)"])) +
         sum(st$log_umi != fmt17(fa[rn, "step1_log_umi"])))
  ca <- mdl@cell.attributes
  files_tsv("sct_cell_attributes", list(umi = ca$umi, log_umi = ca$log_umi), "SCTModel cell.attributes")
  md <- obj@meta.data
  num <- intersect(c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "silhouette_width", "neighborhood_purity"), colnames(md))
  files_tsv("cell_metadata", as.list(md[, num]), "metadata (numeric columns)")
  p <- obj[["pca"]]
  emb <- Embeddings(p)
  files_tsv("pca_embeddings", setNames(lapply(seq_len(ncol(emb)), function(k) emb[, k]), colnames(emb)), "pca embeddings")
  ld <- Loadings(p)
  files_tsv("pca_loadings", setNames(lapply(seq_len(ncol(ld)), function(k) ld[, k]), colnames(ld)), "pca loadings")
  files_tsv("pca_stdev", list(stdev = p@stdev), "pca stdev")
  files_json_number("total_variance", p@misc$total.variance, "pca total.variance")
  files_json_number("min_variance", mdl@arguments$set_min_var, "SCTModel set_min_var")
}

# Run ----------------------------------------------------------------------------------------------------
built <- list()
if ("seurat" %in% loaders) {
  built$seurat <- load_seurat(m)
  mark("load_seurat")
}
if ("sce" %in% loaders) {
  built$sce <- load_sce(m)
  mark("load_sce")
}
if (identical(args$files, "yes")) {
  if (is.null(built$seurat)) stop("--files yes needs the Seurat loader")
  files_check(built$seurat)
  mark("files")
}
o <- readRDS(args$rds)
mark("read_rds")
if (args$kind == "doubletfinder") {
  for (L in names(built)) cmp_doubletfinder(built[[L]], o, L)
} else {
  if (!is.null(built$seurat)) cmp_seurat(built$seurat, o)
  if (!is.null(built$sce)) cmp_sce(built$sce, o)
}
rm(o); invisible(gc())
mark("compare")
rt <- list()
if (identical(args$`round-trip`, "yes")) {
  dir.create(args$tmp, recursive = TRUE, showWarnings = FALSE)
  for (L in names(built)) {
    f <- file.path(args$tmp, paste0(args$object, "_", L, ".rds"))
    x <- save_rds(m$dir, f, object = L)
    back <- readRDS(f)
    rt[[L]] <- data.frame(object = args$object, loader = L, rds_bytes = file.size(f),
                          identical_to_loader = identical(back, built[[L]]),
                          identical_to_rebuilt = identical(back, x), stringsAsFactors = FALSE)
    unlink(f)
    rm(x, back); invisible(gc())
  }
  mark("round_trip")
}

tab <- do.call(rbind, rows)
write.table(tab, file.path(args$out, "k11.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (length(frows)) write.table(do.call(rbind, frows), file.path(args$out, "files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(rt)) write.table(do.call(rbind, rt), file.path(args$out, "roundtrip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
counts <- as.list(table(factor(tab$status, levels = c("identical", "within floor", "beyond floor", "differing", "not rebuilt", "added"))))
by_loader <- lapply(split(tab, tab$loader), function(t) as.list(table(factor(t$status, levels = c("identical", "within floor", "beyond floor", "differing", "not rebuilt", "added")))))
failing <- tab$status %in% c("beyond floor", "differing")
files_ok <- if (length(frows)) all(vapply(frows, function(r) r$differing == 0, logical(1))) else NA
rt_ok <- if (length(rt)) all(vapply(rt, function(r) r$identical_to_loader && r$identical_to_rebuilt, logical(1))) else NA
json_write(list(object = args$object, kind = args$kind, export = normalizePath(args$export), rds = normalizePath(args$rds),
                export_kernel = export_kernel, saved_kernel = args$`saved-kernel`, mode = mode,
                pass = !any(failing) && !isFALSE(files_ok) && !isFALSE(rt_ok),
                components = nrow(tab), status = counts, by_loader = by_loader,
                failing = I(paste(tab$loader[failing], tab$component[failing], sep = ": ")),
                files_identical = files_ok, round_trip_identical = rt_ok,
                steps = steps, peak_rss_kib = ld_vmhwm_kib(), elapsed_s = proc.time()[["elapsed"]] - t0,
                slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                versions = list(scrustio = as.character(packageVersion("scrustio")), Seurat = as.character(packageVersion("Seurat")),
                                SeuratObject = as.character(packageVersion("SeuratObject")),
                                SingleCellExperiment = if ("sce" %in% loaders) as.character(packageVersion("SingleCellExperiment")) else NULL)),
           file.path(args$out, "compare.json"))
cat(sprintf("%s (%s, %s mode): %d components: %s; files %s; round trip %s; %.0f s, peak RSS %.2f GiB\n",
            args$object, args$kind, mode, nrow(tab),
            paste(names(counts), unlist(counts), sep = " ", collapse = ", "),
            if (is.na(files_ok)) "not checked" else if (files_ok) "identical" else "DIFFERING",
            if (is.na(rt_ok)) "not run" else if (rt_ok) "identical" else "DIFFERING",
            proc.time()[["elapsed"]] - t0, ld_vmhwm_kib() / 2^20))
if (any(failing)) {
  print(tab[failing, c("loader", "component", "status", "n", "differing", "max_abs", "max_scaled", "detail")], row.names = FALSE)
  quit(status = 1)
}
