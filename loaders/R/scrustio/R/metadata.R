# Per-cell tables: scrust cluster's metadata.tsv (role cell_metadata), scrust sctransform's
# cell_metadata.tsv (rna_cell_metadata), scrust doubletfinder's cells.tsv (doublet_calls), with
# each column's type as Seurat holds it and each factor in its levels' order.

# Columns of metadata.tsv with a known type (docs/output_manifest.md, role cell_metadata);
# percent.mt (M15b) when the output computed it.
.META_DOUBLE <- c("nCount_RNA", "percent.mt", "nCount_SCT", "silhouette_width", "neighborhood_purity")
.META_INTEGER <- c("nFeature_RNA", "nFeature_SCT")

# Factor levels, from the files of role cluster_levels: each holds one column's levels in order,
# one per line, and is named after its column (`<column>_levels.txt`, docs/output_manifest.md);
# the manifest gives no other link from a levels file to its column.
.cluster_levels <- function(m) {
  out <- list()
  for (e in manifest_entries(m, "cluster_levels")) {
    col <- sub("_levels\\.txt$", "", basename(e$path))
    if (identical(col, basename(e$path))) stop("scrustio: levels file ", e$path, " is not named <column>_levels.txt", call. = FALSE)
    out[[col]] <- .read_lines(.full_path(m, e))
  }
  out
}

# A label column as a factor in its levels' order; an error if a label is not a level.
.as_factor <- function(x, levels, col) {
  bad <- setdiff(unique(x), levels)
  if (length(bad)) stop("scrustio: column ", col, " holds labels missing from its levels: ", paste(head(bad, 5), collapse = ", "), call. = FALSE)
  factor(x, levels = levels)
}

# The metadata table of a scrust cluster or scrust sctransform output, typed, with the cells as
# row names. Label columns with a levels file become factors in that order; silhouette_other and
# purity_maximum, which name clusters, take seurat_clusters' levels.
.read_metadata <- function(m) {
  e <- .entry(m, "cell_metadata", required = FALSE)
  if (is.null(e)) e <- .entry(m, "rna_cell_metadata", required = FALSE)
  if (is.null(e)) return(NULL)
  .check_shape(m, e)
  path <- .full_path(m, e)
  hdr <- strsplit(.read_lines_n(path, 1), "\t", fixed = TRUE)[[1]]
  levels <- .cluster_levels(m)
  numeric <- intersect(hdr, .META_DOUBLE)
  integer <- intersect(hdr, .META_INTEGER)
  df <- .read_tsv(path, numeric = numeric, integer = integer)
  rownames(df) <- df$barcode
  df$barcode <- NULL
  for (col in names(df)) {
    if (col %in% names(levels)) {
      df[[col]] <- .as_factor(df[[col]], levels[[col]], col)
    } else if (col %in% c("silhouette_other", "purity_maximum") && "seurat_clusters" %in% names(levels)) {
      x <- df[[col]]
      bad <- setdiff(unique(x[!is.na(x)]), levels[["seurat_clusters"]])
      if (length(bad)) stop("scrustio: column ", col, " names clusters missing from seurat_clusters' levels", call. = FALSE)
      df[[col]] <- factor(x, levels = levels[["seurat_clusters"]])
    } else if (is.character(df[[col]]) && !col %in% c(numeric, integer)) {
      df[[col]] <- utils::type.convert(df[[col]], as.is = TRUE, na.strings = "NA")
    }
  }
  df
}

# scrust doubletfinder's calls (role doublet_calls: barcode, pANN, classification), named as
# DoubletFinder 2.0.6 names its metadata columns (pANN_<pN>_<pK>_<nExp> and
# DF.classifications_<pN>_<pK>_<nExp>, which the summary gives).
.read_doublet_calls <- function(m) {
  e <- .entry(m, "doublet_calls", required = FALSE)
  if (is.null(e)) return(NULL)
  .check_shape(m, e)
  df <- .read_tsv(.full_path(m, e), numeric = "pANN")
  res <- m$summary$result
  pann <- if (is.character(res$pANN_column)) res$pANN_column else "pANN"
  cls <- if (is.character(res$classification_column)) res$classification_column else "DF.classifications"
  out <- data.frame(df$pANN, df$classification, stringsAsFactors = FALSE)
  names(out) <- c(pann, cls)
  rownames(out) <- df$barcode
  out
}
