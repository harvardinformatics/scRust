# SingleCellExperiment objects, as docs/PLAN.md C3g lays them out: the raw counts as the counts
# assay; altExp(sce, "SCT") with the assays corrected_counts and logcounts (the model genes only),
# the SCTModel's feature attributes as its rowData, its cell attributes as its colData and its
# other slots in its metadata; altExp(sce, "SCT_scale_data") with the assay scaledata (the
# variable features, in VariableFeatures order); colData; reducedDim(sce, "PCA") with the loadings
# as its "rotation" attribute; the colPairs SCT_nn and SCT_snn. Features are named by the features
# file's names, made unique; rowData holds their ids, names and types.

#' Build a SingleCellExperiment from a scrust output
#'
#' From `scrust cluster --export`, `scrust sctransform --corrected-counts` (no raw counts: the
#' main experiment then has no rows, and the SCT data are in the altExps) or
#' `scrust doubletfinder` (the singlets' raw counts and DoubletFinder's calls as colData). Each
#' matrix goes only into a slot meant for its scale (slot_table()).
#'
#' @param dir the output directory.
#' @param check_sha256 check every listed file's SHA-256 first.
#' @return a SingleCellExperiment.
load_sce <- function(dir, check_sha256 = FALSE) {
  .need(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "Matrix"))
  m <- if (inherits(dir, "scrust_manifest")) dir else read_manifest(dir, check_sha256 = check_sha256)
  feat <- .counts_features(m)
  cells <- .cells_of(m)
  md <- .read_metadata(m)
  calls <- .read_doublet_calls(m)
  if (!is.null(calls)) {
    missing <- setdiff(cells, rownames(calls))
    if (length(missing)) stop("scrustio: ", length(missing), " cells have no DoubletFinder call", call. = FALSE)
    calls <- calls[cells, , drop = FALSE]
    md <- if (is.null(md)) calls else cbind(md, calls)
  }
  if (is.null(md)) md <- data.frame(row.names = cells)
  if (!identical(rownames(md), cells)) stop("scrustio: the metadata's cells differ from the output's", call. = FALSE)
  cd <- S4Vectors::DataFrame(md, check.names = FALSE)
  if (.has_role(m, "raw_counts")) {
    counts <- read_for_slot(m, "raw_counts", "SingleCellExperiment counts")
    if (!identical(colnames(counts), cells)) stop("scrustio: the raw counts' cells differ from the output's", call. = FALSE)
    rownames(counts) <- feat$sce
    rd <- S4Vectors::DataFrame(ID = feat$id, Symbol = feat$name, Type = feat$type, row.names = feat$sce)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = cd, rowData = rd)
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(colData = cd)
  }
  if (.has_role(m, "sct_counts")) {
    tabs <- .sct_tables(m, feat, to = "sce")
    if (!identical(.identifiers(m, .entry(m, "sct_counts")$rows), tabs$gene_attr$name)) {
      stop("scrustio: the SCT features (sct_features) differ from the model's genes (gene_attr)", call. = FALSE)
    }
    cc <- read_for_slot(m, "sct_counts", "SingleCellExperiment altExp SCT corrected_counts")
    dd <- read_for_slot(m, "sct_data", "SingleCellExperiment altExp SCT logcounts")
    if (!identical(cc@i, dd@i) || !identical(cc@p, dd@p) || !identical(dd@x, log1p(cc@x))) {
      stop("scrustio: the SCT data layer is not log1p of the counts layer", call. = FALSE)
    }
    dimnames(cc) <- list(tabs$gene_names, cells)
    dimnames(dd) <- list(tabs$gene_names, cells)
    parts <- .sct_model_parts(m, tabs, cells)
    vf <- .read_tsv(.full_path(m, .entry(m, "sct_variable_features")))
    vf_names <- .map_features(feat, vf$name, vf$id, "variable features", "sce")
    sct <- SingleCellExperiment::SingleCellExperiment(
      assays = list(corrected_counts = cc, logcounts = dd),
      rowData = S4Vectors::DataFrame(parts$feature.attributes, check.names = FALSE),
      colData = S4Vectors::DataFrame(parts$cell.attributes, check.names = FALSE),
      metadata = list(SCTModel = parts[c("clips", "umi.assay", "model", "arguments", "median_umi")],
                      variable_features = vf_names))
    SingleCellExperiment::altExp(sce, "SCT") <- sct
    sd <- read_for_slot(m, "sct_scale_data", "SingleCellExperiment altExp SCT_scale_data scaledata")
    if (!identical(rownames(sd), vf$name) || !identical(colnames(sd), cells)) {
      stop("scrustio: scale.data's rows or columns differ from the variable features or the cells", call. = FALSE)
    }
    dimnames(sd) <- list(vf_names, cells)
    SingleCellExperiment::altExp(sce, "SCT_scale_data") <- SingleCellExperiment::SingleCellExperiment(
      assays = list(scaledata = sd))
  }
  if (.has_role(m, "pca_embeddings")) {
    for (s in c("SingleCellExperiment reducedDim PCA", "SingleCellExperiment reducedDim PCA rotation",
                "SingleCellExperiment reducedDim PCA stdev")) {
      e <- .entry(m, .SLOTS[[s]]$role)
      .check_slot(s, e$role, e$scale, e$path, e$format)
    }
    pca <- .pca_parts(m, feat, cells, to = "sce")
    emb <- pca$embeddings
    attr(emb, "rotation") <- pca$loadings
    attr(emb, "stdev") <- pca$stdev
    if (!is.null(pca$total_variance)) attr(emb, "total.variance") <- pca$total_variance
    SingleCellExperiment::reducedDim(sce, "PCA") <- emb
  }
  for (g in list(c("graph_nn", "SCT_nn"), c("graph_snn", "SCT_snn"))) {
    if (!.has_role(m, g[1])) next
    G <- read_for_slot(m, g[1], paste("SingleCellExperiment colPair", g[2]))
    if (!identical(rownames(G), cells) || !identical(colnames(G), cells)) stop("scrustio: graph ", g[2], "'s cells differ", call. = FALSE)
    if (any(G@x == 0)) stop("scrustio: graph ", g[2], " stores zeros, which a colPair cannot hold", call. = FALSE)
    SingleCellExperiment::colPair(sce, g[2]) <- G
  }
  S4Vectors::metadata(sce)$scrust <- list(command = m$command,
                                          scrust_version = m$provenance$scrust_version,
                                          command_line = unlist(m$provenance$command_line))
  sce
}
