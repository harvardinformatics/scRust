#' Write a scrust output as an .rds file
#'
#' Builds the object with load_seurat() or load_sce() and saves it with saveRDS(). As one line:
#'
#'     Rscript -e 'scrustio::save_rds("<scrust output>", "<object.rds>")'
#'     Rscript -e 'scrustio::save_rds("<scrust output>", "<sce.rds>", object = "sce")'
#'
#' @param dir the scrust output directory.
#' @param file the .rds file to write.
#' @param object "seurat" (a Seurat object) or "sce" (a SingleCellExperiment).
#' @param check_sha256 check every listed file's SHA-256 first.
#' @return the object, invisibly.
save_rds <- function(dir, file, object = c("seurat", "sce"), check_sha256 = FALSE) {
  object <- match.arg(object)
  x <- switch(object,
              seurat = load_seurat(dir, check_sha256 = check_sha256),
              sce = load_sce(dir, check_sha256 = check_sha256))
  saveRDS(x, file)
  invisible(x)
}

#' Add a scrust output's per-cell columns to an existing object
#'
#' Adds the per-cell columns of a scrust output to a Seurat object's meta.data or a
#' SingleCellExperiment's colData, matching cells by barcode: scrust doubletfinder's calls (named
#' as DoubletFinder names them), or scrust cluster's metadata (cluster labels as factors in their
#' levels' order, the silhouette widths and purities). Cells the output does not cover get NA;
#' columns the object already has are left as they are, with a message.
#'
#' @param object a Seurat object or a SingleCellExperiment.
#' @param dir the scrust output directory.
#' @return the object with the new columns.
add_scrust <- function(object, dir) {
  m <- if (inherits(dir, "scrust_manifest")) dir else read_manifest(dir)
  md <- .read_metadata(m)
  calls <- .read_doublet_calls(m)
  tabs <- Filter(Negate(is.null), list(md, calls))
  if (!length(tabs)) stop("scrustio: this ", m$command, " output has no per-cell table", call. = FALSE)
  cells <- colnames(object)
  for (tab in tabs) {
    hit <- match(cells, rownames(tab))
    if (all(is.na(hit))) stop("scrustio: none of the object's cells are in the output", call. = FALSE)
    if (inherits(object, "Seurat")) {
      have <- colnames(object[[]])
    } else if (methods::is(object, "SingleCellExperiment")) {
      have <- colnames(SummarizedExperiment::colData(object))
    } else {
      stop("scrustio: add_scrust() takes a Seurat object or a SingleCellExperiment", call. = FALSE)
    }
    skip <- intersect(colnames(tab), have)
    if (length(skip)) message("scrustio: the object already has ", paste(skip, collapse = ", "), "; left as they are")
    add <- tab[hit, setdiff(colnames(tab), have), drop = FALSE]
    rownames(add) <- cells
    if (!ncol(add)) next
    if (inherits(object, "Seurat")) {
      object <- SeuratObject::AddMetaData(object, metadata = add)
    } else {
      for (col in colnames(add)) SummarizedExperiment::colData(object)[[col]] <- add[[col]]
    }
  }
  object
}
