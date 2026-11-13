# The workflow's clustering chain (docs/PLAN.md C1) and the plain-file layout of its results,
# shared by export_clustering.R (fields of the saved objects), reference_clustering.R (direct
# runs) and instrumented_clustering.R. Requires lib/io.R.
#
# Copied from the pipeline repo (scRNAseq-preprocessing; read at commit ea053fa, where these
# scripts last changed in 8d34de1 on 2026-10-30, before the test-reference run of 2026-11-06 and
# the production run of 2026-11-08/09):
# - workflow/scripts/silhouette_utils.R: set_workflow_seed(), require_min_cells_for_pca() and
#   add_silhouette_to_metadata(), verbatim except two timing marks in the last (marked below);
# - the chain's calls, identical in the nine scripts that run it (tenx2seuratrds.R, emptydrops.R,
#   soupx.R, decontx.R, cellbender2seurat.R, doubletfinder.R, scdblfinder.R,
#   posthocfilter_threshold.R; posthocfilter_mad.R is not configured), and their
#   write_cluster_metadata(), verbatim;
# - the packages tenx2seuratrds.R loads (Seurat, glmGamPoi) and its future option.
#
# The chain's input is the object's RNA counts layer as saved (C1). The scripts' other work
# (reading 10x files, emptyDrops, SoupX, decontX, DoubletFinder, scDblFinder, subset(),
# percent.mt) is not part of the chain and is not run: SCTransform reads only the counts and the
# cell names (Seurat passes meta.data to vst() as cell_attr, which uses only the columns umi,
# gene, log_umi, log_gene, umi_per_gene and log_umi_per_gene if present; export_clustering.R
# checks that no saved object has them), and every later step reads only what the chain made.
# R's random stream is reseeded inside every step that draws from it (C2), so the work the scripts
# do before the chain does not affect it.

WORKFLOW_SEED <- 12345L

# From workflow/scripts/silhouette_utils.R (verbatim) -------------------------------------------

set_workflow_seed <- function(seed = Sys.getenv("SCRNASEQ_PREPROCESS_SEED", "12345")) {
  seed <- suppressWarnings(as.integer(seed))
  if (length(seed) != 1 || is.na(seed)) {
    stop("SCRNASEQ_PREPROCESS_SEED must be an integer")
  }
  set.seed(seed)
  seed
}

require_min_cells_for_pca <- function(seurat_obj, context = "", npcs = 50L) {
  n_cells <- ncol(seurat_obj)
  if (n_cells <= npcs) {
    prefix <- if (nzchar(context)) paste0(context, ": ") else ""
    stop(sprintf(
      paste0(
        "%s[LOW_QUALITY_SAMPLE] only %d cell(s) remain - too few to compute %d principal ",
        "components (RunPCA and downstream UMAP/clustering require more cells than PCs). ",
        "This usually means upstream QC/filtering removed nearly all cells for this ",
        "sample; consider dataset-specific thresholds or excluding this sample."
      ),
      prefix, n_cells, npcs
    ), call. = FALSE)
  }
  invisible(n_cells)
}

# Timing marks: the run scripts replace this with their own; the copy below calls it twice.
cl_mark <- function(name) invisible(NULL)

add_silhouette_to_metadata <- function(
    seurat_obj,
    cluster_col = "seurat_clusters",
    reduction = "pca",
    dims = 1:30,
    sil_col_name = "silhouette_width",
    purity_col_name = "neighborhood_purity",
    purity_k = 50,
    purity_weighted = TRUE,
    purity_BNPARAM = NULL,
    purity_BPPARAM = NULL
) {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }

  emb <- Embeddings(seurat_obj, reduction = reduction)
  dims <- dims[dims <= ncol(emb)]
  if (length(dims) == 0) {
    stop(paste("No", reduction, "dimensions available for silhouette calculation"))
  }

  clust <- seurat_obj@meta.data[[cluster_col]]
  valid_cells <- !is.na(clust)
  sil_values <- rep(NA_real_, length(clust))
  purity_values <- rep(NA_real_, length(clust))

  if (sum(valid_cells) < 2 || length(unique(clust[valid_cells])) < 2) {
    seurat_obj@meta.data[[sil_col_name]] <- sil_values
    seurat_obj@meta.data[[purity_col_name]] <- purity_values
    return(seurat_obj)
  }

  message(sprintf(
    "Calculating approximate silhouette widths for %d cells with bluster::approxSilhouette",
    sum(valid_cells)
  ))
  sil <- bluster::approxSilhouette(
    emb[valid_cells, dims, drop = FALSE],
    clust[valid_cells]
  )
  cl_mark("silhouette")  # scRust: timing mark, not in the workflow
  sil_values[valid_cells] <- as.numeric(sil$width)
  seurat_obj@meta.data[[sil_col_name]] <- sil_values

  message(sprintf(
    "Calculating neighborhood purities for %d cells with bluster::neighborPurity",
    sum(valid_cells)
  ))
  purity <- bluster::neighborPurity(
    emb[valid_cells, dims, drop = FALSE],
    clust[valid_cells],
    k = purity_k,
    weighted = purity_weighted,
    BNPARAM = purity_BNPARAM,
    BPPARAM = purity_BPPARAM
  )
  cl_mark("purity")  # scRust: timing mark, not in the workflow
  purity_values[valid_cells] <- as.numeric(purity$purity)
  seurat_obj@meta.data[[purity_col_name]] <- purity_values

  seurat_obj
}

# From workflow/scripts/tenx2seuratrds.R (verbatim) ----------------------------------------------

write_cluster_metadata <- function(seurat_obj, nclusters_output, cluster_ids_output) {
  cluster_ids <- levels(Idents(seurat_obj))
  if (is.null(cluster_ids) || length(cluster_ids) == 0) {
    cluster_ids <- sort(unique(as.character(Idents(seurat_obj))))
  }
  cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
  writeLines(cluster_ids, con = cluster_ids_output)
  writeLines(as.character(length(cluster_ids)), con = nclusters_output)
}

# Setup, input and the chain ---------------------------------------------------------------------

cl_setup <- function() {
  suppressPackageStartupMessages({
    library("Seurat")
    library("glmGamPoi")
  })
  options(future.globals.maxSize = 16 * 1024^3)
}

# The object the chain starts from, built from an export's counts (export_clustering.R, or an M2,
# M6b or M7 export it reuses): the identical dgCMatrix (lib/io.R's exact reader) wrapped with
# CreateSeuratObject(), as the workflow's scripts build their objects. assay "v5" makes an Assay5
# (every object made by CreateSeuratObject or subset from one); "v3" makes a v3 Assay, as the
# emptyDrops objects hold (emptydrops.R: as.Seurat() of a SingleCellExperiment), so SCTransform
# dispatches to SCTransform.Assay directly. A path ending in .rds reads a saved object instead and
# keeps only its RNA counts (for diagnosis).
cl_load_input <- function(input, assay = c("v5", "v3")) {
  assay <- match.arg(assay)
  if (grepl("\\.rds$", input)) {
    o <- readRDS(input)
    counts <- LayerData(o, assay = "RNA", layer = "counts")
    rm(o)
  } else {
    counts <- mtx_read_10x(input)
  }
  if (assay == "v3") {
    old <- options(Seurat.object.assay.version = "v3")
    on.exit(options(old), add = TRUE)
  }
  seu <- CreateSeuratObject(counts = counts)
  want <- if (assay == "v3") "Assay" else "Assay5"
  if (!identical(class(seu[["RNA"]])[1], want)) stop("RNA assay is a ", class(seu[["RNA"]])[1], ", not ", want)
  if (!identical(LayerData(seu, assay = "RNA", layer = "counts"), counts)) {
    stop("CreateSeuratObject() changed the counts matrix")
  }
  if (!identical(colnames(seu), colnames(counts))) stop("cells reordered")
  list(seurat = seu, source = if (grepl("\\.rds$", input)) "rds" else "export")
}

cl_rng_digest <- function() {
  s <- rng_state()
  if (is.null(s)) return("none")
  digest::digest(paste(s, collapse = ","), algo = "sha256", serialize = FALSE)
}

# The workflow's calls (C1), in order. `umap = FALSE` leaves out RunUMAP (M8's check that it
# changes nothing downstream). mark(name) is called after each step. FindClusters' printed
# output (the modularity optimizer's numbers) and its messages (GroupSingletons) are captured
# and still printed. Returns the object and what the steps printed, with a digest of R's
# .Random.seed after each step.
cl_chain <- function(seu, umap = TRUE, mark = function(name) invisible(NULL)) {
  rng <- list()
  seu <- SCTransform(seu, seed.use = WORKFLOW_SEED, verbose = FALSE)
  rng$SCTransform <- cl_rng_digest()
  mark("SCTransform")
  require_min_cells_for_pca(seu, context = "reference_clustering")
  seu <- RunPCA(seu, seed.use = WORKFLOW_SEED, verbose = FALSE)
  rng$RunPCA <- cl_rng_digest()
  mark("RunPCA")
  if (umap) {
    seu <- RunUMAP(seu, dims = 1:30, seed.use = WORKFLOW_SEED)
    rng$RunUMAP <- cl_rng_digest()
    mark("RunUMAP")
  }
  seu <- FindNeighbors(seu, dims = 1:30)
  rng$FindNeighbors <- cl_rng_digest()
  mark("FindNeighbors")
  msgs <- character()
  printed <- utils::capture.output(
    withCallingHandlers(
      seu <- FindClusters(seu, random.seed = WORKFLOW_SEED),
      message = function(m) msgs <<- c(msgs, conditionMessage(m))),
    split = TRUE)
  rng$FindClusters <- cl_rng_digest()
  mark("FindClusters")
  seu <- add_silhouette_to_metadata(seu)
  rng$add_silhouette_to_metadata <- cl_rng_digest()
  list(seurat = seu, printed = printed, messages = msgs, rng = rng)
}

# The optimizer's numbers as Seurat prints them, and GroupSingletons' message.
cl_parse_printed <- function(printed, messages) {
  one <- function(re) {
    hit <- grep(re, printed, value = TRUE)
    if (length(hit) != 1) stop("expected one line matching ", re, ", found ", length(hit))
    sub(re, "\\1", hit)
  }
  sg <- grep("^[0-9]+ singletons identified\\. [0-9]+ final clusters\\.", trimws(messages), value = TRUE)
  comm <- as.integer(one("^Number of communities: ([0-9]+)$"))
  list(
    nodes = as.integer(one("^Number of nodes: ([0-9]+)$")),
    edges = as.integer(one("^Number of edges: ([0-9]+)$")),
    max_modularity = one("^Maximum modularity in 10 random starts: ([0-9.]+)$"),
    communities = comm,
    singletons = if (length(sg)) as.integer(sub("^([0-9]+) singletons.*", "\\1", sg[1])) else 0L,
    final_clusters = if (length(sg)) as.integer(sub(".*identified\\. ([0-9]+) final clusters\\..*", "\\1", sg[1])) else comm
  )
}

# Plain files --------------------------------------------------------------------------------------

# The chain's metadata columns (C1), in this order; the saved objects' other columns are not part
# of the chain.
CL_META <- c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8",
             "seurat_clusters", "silhouette_width", "neighborhood_purity")

# Names of the cell attributes vst() would take from meta.data instead of computing them.
CL_VST_ATTR <- c("umi", "gene", "log_umi", "log_gene", "umi_per_gene", "log_umi_per_gene")

gz_lines <- function(lines, path) {
  con <- gzfile(path, "w")
  writeLines(lines, con, useBytes = TRUE)
  close(con)
}

# A sparse matrix (dgCMatrix, or a Graph) as Matrix Market, entries exactly as stored (column by
# column, rows in stored order), values through fmt17, written in chunks of columns. The header
# says integer when every value is integer-valued. `comments` go after the banner. Row and column
# names are written separately.
mtx_write_plain <- function(m, path, comments = character(), chunk_nnz = 5e6) {
  x <- m@x
  if (any(!is.finite(x))) stop("mtx_write_plain: non-finite values in ", path)
  if (any(x == 0 & 1 / x < 0)) stop("mtx_write_plain: negative zeros in ", path)
  field <- if (all(x == round(x))) "integer" else "real"
  con <- gzfile(path, "w")
  writeLines(c(sprintf("%%%%MatrixMarket matrix coordinate %s general", field),
               paste0("% ", comments),
               sprintf("%d %d %d", nrow(m), ncol(m), length(x))), con)
  p <- m@p
  nc <- ncol(m)
  j <- 1L
  while (j <= nc) {
    k <- j
    while (k < nc && p[k + 2] - p[j] <= chunk_nnz) k <- k + 1L
    lo <- p[j] + 1
    hi <- p[k + 1]
    if (hi >= lo) {
      cols <- rep.int(j:k, diff(p[j:(k + 1)]))
      writeLines(paste(m@i[lo:hi] + 1L, cols, fmt17(x[lo:hi])), con)
    }
    j <- k + 1L
  }
  close(con)
  invisible(field)
}

# sha256 of a numeric vector's doubles as little-endian bytes (what a .npy holds after its header),
# for matrices too large to write in full; written in chunks.
raw_sha256 <- function(x) {
  x <- as.double(x)
  n <- length(x)
  if (n == 0) return(digest::digest(raw(0), algo = "sha256", serialize = FALSE))
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f))
  con <- file(f, "wb")
  step <- 5e7
  for (s in seq(1, n, by = step)) writeBin(x[s:min(n, s + step - 1)], con, size = 8, endian = "little")
  close(con)
  digest::digest(file = f, algo = "sha256")
}

cl_check_names <- function(seu) {
  if (!identical(rownames(seu@meta.data), colnames(seu))) stop("meta.data rows differ from the cells")
}

# Writes what the chain leaves (C1, "What the step leaves") from a Seurat object into `dir`, in one
# layout for the saved objects (export_clustering.R) and the reference runs, so the two compare
# file by file. Levels:
#   labels    metadata.tsv and the levels files only
#   light     + sct/ (corrected counts in full; the log1p data and scale.data as digests), pca/,
#             graphs/
#   standard  + sct/scale_data.npy
#   full      + sct/log1p_corrected_counts.mtx.gz
# Files:
#   metadata.tsv                  barcode and CL_META, in cell order (doubles exact, %.17g)
#   seurat_clusters_levels.txt, SCT_snn_res.0.8_levels.txt, idents_levels.txt
#   sct/corrected_counts.mtx.gz   the SCT assay's counts layer (integer Matrix Market, entries as
#                                 stored; rows sct/sct_features.tsv.gz, columns sct/sct_cells.tsv.gz)
#   sct/log1p_corrected_counts.mtx.gz   its data layer (real; only at level full)
#   sct/scale_data.npy            its scale.data layer (rows in stored order,
#                                 sct/scale_data_features.tsv; columns the cells)
#   sct/variable_features.tsv     VariableFeatures, in order
#   sct/model/feature_attributes.tsv, sct/model/cell_attributes.tsv, sct/model/model.json
#                                 the SCTModel's slots
#   pca/embeddings.npy, pca/loadings.npy, pca/loadings_features.tsv, pca/stdev.npy, pca/pca.json
#   graphs/nn.mtx.gz, graphs/snn.mtx.gz  SCT_nn and SCT_snn as stored (real Matrix Market, values
#                                 %.17g; rows and columns are the cells), graphs/graphs.json
# Returns the facts for object.json: dimensions, checks and digests.
cl_write_object <- function(seu, dir, level = c("full", "standard", "light", "labels"),
                            software = "scRust validation/R/lib/clustering.R") {
  level <- match.arg(level)
  cl_check_names(seu)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cells <- colnames(seu)
  md <- seu@meta.data
  facts <- list(level = level, n_cells = length(cells))
  missing <- setdiff(CL_META, colnames(md))
  if (length(missing)) stop("metadata columns missing: ", paste(missing, collapse = ", "))
  meta <- data.frame(barcode = cells, stringsAsFactors = FALSE)
  for (col in CL_META) meta[[col]] <- md[[col]]
  tsv_write(meta, file.path(dir, "metadata.tsv"))
  for (col in c("seurat_clusters", "SCT_snn_res.0.8")) {
    if (!is.factor(md[[col]])) stop(col, " is not a factor")
    writeLines(levels(md[[col]]), file.path(dir, paste0(col, "_levels.txt")))
  }
  writeLines(levels(Idents(seu)), file.path(dir, "idents_levels.txt"))
  facts$idents_equal_seurat_clusters <- identical(unname(Idents(seu)), unname(md$seurat_clusters)) ||
    identical(as.character(Idents(seu)), as.character(md$seurat_clusters))
  facts$n_clusters <- nlevels(md$seurat_clusters)
  facts$cluster_sizes <- as.integer(table(md$seurat_clusters))
  facts$default_assay <- DefaultAssay(seu)
  facts$assays <- as.list(vapply(Assays(seu), function(a) class(seu[[a]])[1], character(1)))
  facts$metadata_vst_attr_names <- I(intersect(CL_VST_ATTR, colnames(md)))
  if (level == "labels") return(facts)

  # SCT assay
  sct <- seu[["SCT"]]
  sdir <- file.path(dir, "sct")
  dir.create(file.path(sdir, "model"), recursive = TRUE, showWarnings = FALSE)
  cnt <- GetAssayData(sct, layer = "counts")
  dat <- GetAssayData(sct, layer = "data")
  sd <- GetAssayData(sct, layer = "scale.data")
  if (!inherits(cnt, "dgCMatrix") || !inherits(dat, "dgCMatrix")) stop("SCT counts or data is not a dgCMatrix")
  if (!identical(colnames(cnt), cells) || !identical(colnames(sd), cells)) stop("SCT layers' cells differ from the object's")
  facts$sct_n_features <- nrow(cnt)
  facts$sct_nnz <- length(cnt@x)
  facts$sct_counts_integer <- all(cnt@x == round(cnt@x))
  facts$sct_counts_min <- min(cnt@x)
  facts$sct_counts_explicit_zeros <- sum(cnt@x == 0)
  facts$sct_data_is_log1p_counts <- identical(dat@i, cnt@i) && identical(dat@p, cnt@p) &&
    identical(dimnames(dat), dimnames(cnt)) && identical(dat@x, log1p(cnt@x))
  facts$nCount_SCT_is_colSums <- identical(md$nCount_SCT, unname(Matrix::colSums(cnt)))
  facts$nFeature_SCT_is_nonzero <- identical(md$nFeature_SCT, unname(Matrix::colSums(cnt > 0)))
  comments <- c(sprintf("metadata_json: {\"software_version\": \"%s\", \"format_version\": 2}", software),
                "Seurat SCT assay, counts layer: the corrected counts of sctransform::correct(), NOT raw counts.",
                "Rows: sct_features.tsv.gz; columns: sct_cells.tsv.gz.")
  mtx_write_plain(cnt, file.path(sdir, "corrected_counts.mtx.gz"), comments)
  gz_lines(rownames(cnt), file.path(sdir, "sct_features.tsv.gz"))
  gz_lines(cells, file.path(sdir, "sct_cells.tsv.gz"))
  facts$sct_counts_x_sha256 <- raw_sha256(cnt@x)
  facts$sct_data_x_sha256 <- raw_sha256(dat@x)
  if (level == "full") {
    mtx_write_plain(dat, file.path(sdir, "log1p_corrected_counts.mtx.gz"),
                    c(comments[1], "Seurat SCT assay, data layer: log1p of the corrected counts, NOT raw counts.",
                      comments[3]))
  }
  facts$scale_data_dim <- dim(sd)
  facts$scale_data_sha256 <- raw_sha256(sd)
  writeLines(rownames(sd), file.path(sdir, "scale_data_features.tsv"))
  if (level %in% c("full", "standard")) npy_write(file.path(sdir, "scale_data.npy"), unname(sd))
  vf <- VariableFeatures(sct)
  writeLines(vf, file.path(sdir, "variable_features.tsv"))
  facts$n_variable_features <- length(vf)
  facts$scale_data_rows_are_variable_features <- setequal(rownames(sd), vf)
  facts$scale_data_rows_in_variable_feature_order <- identical(rownames(sd), vf)
  models <- slot(sct, "SCTModel.list")
  facts$sct_models <- length(models)
  m <- models[[1]]
  fa <- m@feature.attributes
  tsv_write(cbind(data.frame(name = rownames(fa), stringsAsFactors = FALSE), fa),
            file.path(sdir, "model", "feature_attributes.tsv"))
  ca <- m@cell.attributes
  tsv_write(cbind(data.frame(barcode = rownames(ca), stringsAsFactors = FALSE), ca),
            file.path(sdir, "model", "cell_attributes.tsv"))
  facts$model_feature_rows_are_sct_features <- identical(rownames(fa), rownames(cnt))
  facts$model_cell_rows_are_cells <- identical(rownames(ca), cells)
  json_write(list(clips = m@clips, umi.assay = m@umi.assay, model = m@model, median_umi = m@median_umi,
                  arguments = m@arguments, feature_attribute_columns = colnames(fa),
                  cell_attribute_columns = colnames(ca)),
             file.path(sdir, "model", "model.json"))

  # PCA
  pdir <- file.path(dir, "pca")
  dir.create(pdir, showWarnings = FALSE)
  p <- seu[["pca"]]
  emb <- Embeddings(p)
  if (!identical(rownames(emb), cells)) stop("PCA embeddings' cells differ")
  npy_write(file.path(pdir, "embeddings.npy"), unname(emb))
  ld <- Loadings(p)
  npy_write(file.path(pdir, "loadings.npy"), unname(ld))
  writeLines(rownames(ld), file.path(pdir, "loadings_features.tsv"))
  npy_write(file.path(pdir, "stdev.npy"), Stdev(p))
  json_write(list(key = Key(p), assay_used = DefaultAssay(p), npcs = ncol(emb),
                  total_variance = p@misc$total.variance, global = p@global,
                  loadings_projected_dim = dim(p@feature.loadings.projected),
                  loadings_features_are_variable_features = identical(rownames(ld), vf)),
             file.path(pdir, "pca.json"))

  # Graphs
  gdir <- file.path(dir, "graphs")
  dir.create(gdir, showWarnings = FALSE)
  gfacts <- list()
  for (g in c("SCT_nn", "SCT_snn")) {
    G <- seu@graphs[[g]]
    if (is.null(G)) stop("graph ", g, " missing")
    if (!identical(rownames(G), cells) || !identical(colnames(G), cells)) stop("graph ", g, "'s names differ from the cells")
    name <- if (g == "SCT_nn") "nn" else "snn"
    mtx_write_plain(G, file.path(gdir, paste0(name, ".mtx.gz")),
                    c(sprintf("metadata_json: {\"software_version\": \"%s\", \"format_version\": 2}", software),
                      sprintf("Seurat graph %s as stored; rows and columns are the cells (metadata.tsv order).", g)))
    colidx <- rep.int(seq_len(ncol(G)), diff(G@p))
    same_col <- colidx[-1] == colidx[-length(colidx)]
    gfacts[[name]] <- list(class = class(G)[1], assay_used = DefaultAssay(G), dim = dim(G), nnz = length(G@x),
                           symmetric = Matrix::isSymmetric(as(G, "dgCMatrix")),
                           rows_sorted_within_columns = all(diff(G@i)[same_col] > 0),
                           diagonal_entries = sum(Matrix::diag(as(G, "dgCMatrix")) != 0),
                           distinct_values = length(unique(G@x)))
  }
  json_write(gfacts, file.path(gdir, "graphs.json"))
  facts
}

# The run driver ---------------------------------------------------------------------------------

cl_vmrss_kib <- function() {
  st <- readLines("/proc/self/status")
  as.numeric(sub("^VmRSS:\\s+([0-9]+) kB$", "\\1", grep("^VmRSS:", st, value = TRUE)))
}

CL_PACKAGES <- c("Seurat", "SeuratObject", "sctransform", "glmGamPoi", "irlba", "uwot", "RcppAnnoy",
                 "RcppEigen", "bluster", "BiocNeighbors", "assorthead", "Matrix", "future",
                 "future.apply", "Rcpp")

# A reference run (reference_clustering.R, instrumented_clustering.R): loads the export's counts,
# runs the chain, and writes the outputs. args: input, out, level, umap (yes/no), assay
# (v3/v5/auto: auto reads input/assay.txt). `hooks` (instrumented runs) is a list of functions
# install(seurat), finish(seurat, umap) called around the chain.
# Outputs in args$out: cl_write_object()'s layout at args$level; cluster_ids.txt and
# nclusters.txt (write_cluster_metadata); summary.json with a deterministic "result" section
# (equal between runs that agree: the object's facts, the optimizer's printed numbers and
# GroupSingletons' message, a digest of .Random.seed after each step) and a "run" section
# (timings, peak RSS after each step, current RSS, fingerprint).
cl_main <- function(args, kind, hooks = NULL) {
  umap <- switch(args$umap, yes = TRUE, no = FALSE, stop("--umap must be yes or no"))
  assay <- if (identical(args$assay, "auto")) readLines(file.path(args$input, "assay.txt")) else args$assay
  dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
  cl_setup()
  phases <- list()
  t0 <- proc.time()[["elapsed"]]
  mark <- function(name) {
    phases[[name]] <<- list(elapsed_s = proc.time()[["elapsed"]] - t0, peak_rss_kib = peak_rss_kib(),
                            rss_kib = cl_vmrss_kib())
  }
  assign("cl_mark", mark, envir = globalenv())
  inp <- cl_load_input(args$input, assay)
  mark("load")
  seu <- inp$seurat
  run_info <- list(kind = kind, input = args$input, input_source = inp$source, outdir = args$out,
                   level = args$level, umap = umap, assay = assay,
                   start_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                   slurm_job_id = Sys.getenv("SLURM_JOB_ID"), profile = Sys.getenv("SCRUST_PROFILE"),
                   fingerprint = fingerprint(packages = CL_PACKAGES))
  rm(inp)
  if (!is.null(hooks)) hooks$install(seu)
  res <- cl_chain(seu, umap = umap, mark = mark)
  rm(seu)
  seu <- res$seurat
  mark("chain_end")
  if (!is.null(hooks)) hooks$finish(seu, umap)
  write_cluster_metadata(seu, file.path(args$out, "nclusters.txt"), file.path(args$out, "cluster_ids.txt"))
  facts <- cl_write_object(seu, args$out, args$level, software = paste0("scRust ", kind))
  mark("write")
  result <- c(facts, list(
    rna_assay_class = class(seu[["RNA"]])[1],
    seed = WORKFLOW_SEED,
    optimizer = cl_parse_printed(res$printed, res$messages),
    rng_after = res$rng,
    rng_end_sha256 = cl_rng_digest()))
  json_write(list(result = result,
                  run = c(run_info, list(phases = phases, peak_rss_kib = peak_rss_kib(),
                                         printed = I(res$printed), messages = I(res$messages)))),
             file.path(args$out, "summary.json"))
  cat(kind, ":", facts$n_cells, "cells,", facts$n_clusters, "clusters, optimizer", result$optimizer$nodes,
      "nodes", result$optimizer$edges, "edges, Q", result$optimizer$max_modularity, ",",
      result$optimizer$singletons, "singletons; peak RSS", peak_rss_kib(), "KiB,",
      round(proc.time()[["elapsed"]] - t0), "s\n")
  invisible(result)
}
