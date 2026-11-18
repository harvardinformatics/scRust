# Feature names as Seurat holds them. scrust names features by the 10x features file's column 2
# as written; Seurat's Read10X() makes them unique (make.unique(), unique.features = TRUE) and
# SeuratObject's CreateAssayObject() and CreateSeuratObject() replace "_" and "|" with "-"
# (CheckFeaturesNames()). SCTransform then runs on those names, so every list of features of an
# output (the SCT features, the variable features, the model's genes, the PCA's loadings) is
# mapped onto the counts' features by position and renamed the same way.

#' Seurat's names for a features file's names
#'
#' @param names the features' names (column 2 of features.tsv), in row order.
#' @param warn warn, as Seurat does, when underscores or pipes are replaced.
#' @return the names Seurat gives these rows.
seurat_feature_names <- function(names, warn = TRUE) {
  out <- if (anyDuplicated(names)) make.unique(names) else names
  if (any(grepl("_", out, fixed = TRUE))) {
    if (warn) warning("Feature names cannot have underscores ('_'), replacing with dashes ('-')", call. = FALSE)
    out <- gsub("_", "-", out, fixed = TRUE)
  }
  if (any(grepl("|", out, fixed = TRUE))) {
    if (warn) warning("Feature names cannot have pipe characters ('|'), replacing with dashes ('-')", call. = FALSE)
    out <- gsub("|", "-", out, fixed = TRUE)
  }
  out
}

# The counts' features (ids, names as written, Seurat's names), from the raw counts' features
# file (role raw_counts_features), or NULL if the output holds no raw counts.
.counts_features <- function(m) {
  e <- .entry(m, "raw_counts_features", required = FALSE)
  if (is.null(e)) return(NULL)
  lines <- .read_lines(.full_path(m, e))
  f <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(f, `[`, "", 1)
  nm <- vapply(f, function(x) if (length(x) >= 2) x[2] else x[1], "")
  type <- vapply(f, function(x) if (length(x) >= 3) x[3] else NA_character_, "")
  list(id = ids, name = nm, type = type, seurat = seurat_feature_names(nm, warn = FALSE),
       sce = if (anyDuplicated(nm)) make.unique(nm) else nm)
}

# Seurat's names (to = "seurat") or the SingleCellExperiment's (to = "sce": the names as
# written, made unique) for a list of features given by name (and, where names repeat, by id),
# as rows of the counts. With no raw counts in the output (scrust sctransform's), the names are
# used as written, renamed by the same rules.
.map_features <- function(feat, names, ids = NULL, what = "features", to = "seurat") {
  if (is.null(feat)) {
    return(if (to == "seurat") seurat_feature_names(names, warn = FALSE) else if (anyDuplicated(names)) make.unique(names) else names)
  }
  if (!anyDuplicated(feat$name)) {
    k <- match(names, feat$name)
  } else if (!is.null(ids) && !anyDuplicated(feat$id)) {
    k <- match(ids, feat$id)
  } else {
    stop("scrustio: the counts' feature names repeat and the ", what, " cannot be matched to them by id",
         call. = FALSE)
  }
  if (anyNA(k)) {
    stop("scrustio: ", sum(is.na(k)), " of the ", what, " are not among the counts' features (first: ",
         names[which(is.na(k))[1]], ")", call. = FALSE)
  }
  feat[[to]][k]
}
