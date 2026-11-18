# Seurat objects (Seurat 5.5.1, SeuratObject 5.4.0), as docs/PLAN.md C3g lays them out:
# the raw counts in the RNA assay's counts layer; an SCTAssay named SCT with the corrected counts
# in counts, their log1p in data, the residuals in scale.data and exactly one SCTModel;
# meta.data with factors in their levels' order; the pca reduction; the SCT_nn and SCT_snn
# graphs; identities from seurat_clusters; SCT the default assay. Each piece is made with the
# SeuratObject and Seurat calls the workflow's chain makes it with (SCTransform, RunPCA,
# FindNeighbors, FindClusters), so the object holds what theirs holds.

.need <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) stop("scrustio: this needs the R package ", p, call. = FALSE)
  }
}

# The versions whose classes the loaders build and are tested with; others are used with a warning.
.TESTED <- c(Seurat = "5.5.1", SeuratObject = "5.4.0")

.check_versions <- function(tested = .TESTED) {
  for (p in names(tested)) {
    v <- as.character(utils::packageVersion(p))
    if (!identical(v, tested[[p]])) {
      warning("scrustio builds ", p, " ", tested[[p]], "'s classes and is tested with it; ", p, " ", v,
              " is installed, untested: check the object", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.cells_of <- function(m) {
  for (role in c("raw_counts", "sct_counts", "cell_metadata", "rna_cell_metadata")) {
    e <- .entry(m, role, required = FALSE)
    if (is.null(e)) next
    ids <- if (role %in% c("raw_counts", "sct_counts")) e$columns else e$rows
    return(.identifiers(m, ids))
  }
  stop("scrustio: this ", m$command, " output names no cells", call. = FALSE)
}

#' Build a Seurat object from a scrust output
#'
#' From `scrust cluster --export` (the raw counts, the SCT assay with its model, the metadata,
#' the PCA, the graphs and the identities), `scrust sctransform --corrected-counts` (the SCT
#' assay and the counts' per-cell columns; that output holds no raw counts) or
#' `scrust doubletfinder` (the singlets' raw counts and DoubletFinder's calls as metadata
#' columns). Files are found through the manifest, by role, and each matrix goes only into a
#' slot meant for its scale (slot_table()).
#'
#' @param dir the output directory.
#' @param check_sha256 check every listed file's SHA-256 first.
#' @return a Seurat object.
load_seurat <- function(dir, check_sha256 = FALSE) {
  .need(c("SeuratObject", "Seurat", "Matrix"))
  .check_versions()
  m <- if (inherits(dir, "scrust_manifest")) dir else read_manifest(dir, check_sha256 = check_sha256)
  old <- options(Seurat.object.assay.version = "v5")
  on.exit(options(old), add = TRUE)
  feat <- .counts_features(m)
  cells <- .cells_of(m)
  obj <- NULL
  if (.has_role(m, "raw_counts")) {
    counts <- read_for_slot(m, "raw_counts", "Seurat RNA counts")
    if (!identical(colnames(counts), cells)) stop("scrustio: the raw counts' cells differ from the output's", call. = FALSE)
    rownames(counts) <- feat$seurat
    obj <- SeuratObject::CreateSeuratObject(counts = counts)
  }
  if (.has_role(m, "sct_counts")) {
    sct <- .sct_assay(m, feat, cells)
    if (is.null(obj)) {
      obj <- SeuratObject::CreateSeuratObject(counts = sct, assay = "SCT")
    } else {
      obj[["SCT"]] <- sct
    }
    SeuratObject::DefaultAssay(obj) <- "SCT"
  }
  if (is.null(obj)) {
    stop("scrustio: this ", m$command, " output holds neither raw counts nor SCTransform's layers ",
         "(scrust cluster --export holds both); add its per-cell columns to an object with add_scrust()",
         call. = FALSE)
  }
  obj <- .add_metadata_seurat(obj, m)
  if (.has_role(m, "pca_embeddings")) obj[["pca"]] <- .pca_reduction(m, feat, cells)
  for (g in list(c("graph_nn", "SCT_nn"), c("graph_snn", "SCT_snn"))) {
    if (!.has_role(m, g[1])) next
    G <- read_for_slot(m, g[1], paste("Seurat graph", g[2]))
    if (!identical(rownames(G), cells) || !identical(colnames(G), cells)) stop("scrustio: graph ", g[2], "'s cells differ", call. = FALSE)
    G <- SeuratObject::as.Graph(G)
    SeuratObject::DefaultAssay(G) <- "SCT"
    obj[[g[2]]] <- G
  }
  if ("seurat_clusters" %in% colnames(obj[[]])) {
    SeuratObject::Idents(obj) <- "seurat_clusters"
  }
  obj
}

# The SCT assay: CreateAssayObject() on the corrected counts, VariableFeatures, the data and
# scale.data layers set as SCTransform.Assay sets them (SetAssayData() keeps scale.data's rows in
# the assay's gene order), and the one SCTModel; the class and assay.orig as SCTransform.Seurat
# leaves them.
.sct_assay <- function(m, feat, cells) {
  tabs <- .sct_tables(m, feat)
  sct_names_written <- .identifiers(m, .entry(m, "sct_counts")$rows)
  if (!identical(sct_names_written, tabs$gene_attr$name)) {
    stop("scrustio: the SCT features (sct_features) differ from the model's genes (gene_attr)", call. = FALSE)
  }
  genes <- tabs$gene_names
  cc <- read_for_slot(m, "sct_counts", "Seurat SCT counts")
  dd <- read_for_slot(m, "sct_data", "Seurat SCT data")
  if (!identical(colnames(cc), cells) || !identical(colnames(dd), cells)) stop("scrustio: the SCT layers' cells differ", call. = FALSE)
  if (!identical(cc@i, dd@i) || !identical(cc@p, dd@p)) stop("scrustio: the SCT data layer's entries differ from the counts layer's", call. = FALSE)
  if (!identical(dd@x, log1p(cc@x))) {
    stop("scrustio: the SCT data layer is not log1p of the counts layer in ", sum(dd@x != log1p(cc@x)),
         " values", call. = FALSE)
  }
  dimnames(cc) <- list(genes, cells)
  dimnames(dd) <- list(genes, cells)
  vf <- .read_tsv(.full_path(m, .entry(m, "sct_variable_features")))
  vf_names <- .map_features(feat, vf$name, vf$id, "variable features")
  sd <- read_for_slot(m, "sct_scale_data", "Seurat SCT scale.data")
  if (!identical(rownames(sd), vf$name) || !identical(colnames(sd), cells)) {
    stop("scrustio: scale.data's rows or columns differ from the variable features or the cells", call. = FALSE)
  }
  dimnames(sd) <- list(vf_names, cells)
  assay <- SeuratObject::CreateAssayObject(counts = cc)
  SeuratObject::VariableFeatures(assay) <- vf_names
  assay <- SeuratObject::SetAssayData(assay, layer = "data", new.data = dd)
  assay <- SeuratObject::SetAssayData(assay, layer = "scale.data", new.data = sd)
  model <- .sct_model(m, tabs, cells)
  methods::new(methods::getClass("SCTAssay", where = asNamespace("Seurat")), assay,
               SCTModel.list = list(model1 = model), assay.orig = "RNA")
}

# The pca reduction, as RunPCA.default makes it: CreateDimReducObject() with the embeddings, the
# loadings (rows: the PCA's features, in VariableFeatures order), the standard deviations, the
# key PC_, the assay SCT and misc$total.variance.
.pca_reduction <- function(m, feat, cells) {
  for (s in c("Seurat pca cell.embeddings", "Seurat pca feature.loadings", "Seurat pca stdev")) {
    e <- .entry(m, .SLOTS[[s]]$role)
    .check_slot(s, e$role, e$scale, e$path, e$format)
  }
  p <- .pca_parts(m, feat, cells, to = "seurat")
  misc <- if (is.null(p$total_variance)) list() else list(total.variance = p$total_variance)
  SeuratObject::CreateDimReducObject(embeddings = p$embeddings, loadings = p$loadings, assay = "SCT",
                                     stdev = p$stdev, key = "PC_", misc = misc)
}

# The PCA's files: the embeddings (cells x components), the loadings (features x components),
# the standard deviations (unnamed) and the total variance (from the summary's result).
.pca_parts <- function(m, feat, cells, to) {
  e_emb <- .entry(m, "pca_embeddings")
  e_ld <- .entry(m, "pca_loadings")
  e_sd <- .entry(m, "pca_stdev")
  for (e in list(e_emb, e_ld, e_sd)) .check_shape(m, e)
  emb_hdr <- strsplit(.read_lines_n(.full_path(m, e_emb), 1), "\t", fixed = TRUE)[[1]]
  pcs <- emb_hdr[-1]
  emb_df <- .read_tsv(.full_path(m, e_emb), numeric = pcs)
  if (!identical(emb_df$barcode, cells)) stop("scrustio: the PCA embeddings' cells differ", call. = FALSE)
  emb <- as.matrix(emb_df[, pcs, drop = FALSE])
  dimnames(emb) <- list(cells, pcs)
  ld_df <- .read_tsv(.full_path(m, e_ld), numeric = pcs)
  ld <- as.matrix(ld_df[, pcs, drop = FALSE])
  dimnames(ld) <- list(.map_features(feat, ld_df$name, ld_df$id, "PCA features", to), pcs)
  sd_df <- .read_tsv(.full_path(m, e_sd), numeric = "stdev")
  if (!identical(sd_df$pc, pcs)) stop("scrustio: the standard deviations' components differ from the embeddings'", call. = FALSE)
  res <- m$summary$result
  tv <- if (!is.null(res$pca$total_variance)) res$pca$total_variance else res$total_variance
  list(embeddings = emb, loadings = ld, stdev = as.numeric(sd_df$stdev),
       total_variance = if (is.null(tv)) NULL else as.numeric(tv))
}

# meta.data: the output's per-cell columns (and DoubletFinder's calls), each with Seurat's type,
# factors in their levels' order. Columns Seurat computes itself when the assays are added
# (nCount_RNA, nFeature_RNA, nCount_SCT, nFeature_SCT) must equal the output's.
.add_metadata_seurat <- function(obj, m) {
  md <- .read_metadata(m)
  calls <- .read_doublet_calls(m)
  cells <- colnames(obj)
  if (!is.null(calls)) {
    missing <- setdiff(cells, rownames(calls))
    if (length(missing)) stop("scrustio: ", length(missing), " cells have no DoubletFinder call", call. = FALSE)
    calls <- calls[cells, , drop = FALSE]
    md <- if (is.null(md)) calls else cbind(md, calls)
  }
  if (is.null(md)) return(obj)
  if (!identical(rownames(md), cells)) stop("scrustio: the metadata's cells differ from the object's", call. = FALSE)
  have <- obj[[]]
  for (col in intersect(colnames(md), colnames(have))) {
    if (!col %in% c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT")) {
      stop("scrustio: the output's column ", col, " would replace the object's", call. = FALSE)
    }
    if (!identical(as.double(md[[col]]), as.double(have[[col]]))) {
      stop("scrustio: the output's ", col, " differs from what Seurat computes from its matrices", call. = FALSE)
    }
  }
  new_cols <- setdiff(colnames(md), colnames(have))
  if (length(new_cols)) obj <- SeuratObject::AddMetaData(obj, metadata = md[, new_cols, drop = FALSE])
  obj
}
