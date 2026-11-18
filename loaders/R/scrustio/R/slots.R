# The rule (docs/PLAN.md C3g): normalised data goes only into slots meant for normalised data.
# The raw counts always come from the raw counts (counts/ of an export, singlets/ of scrust
# doubletfinder), and the corrected counts never go into a raw-count slot (Seurat's RNA counts,
# the SingleCellExperiment's counts).
#
# Every slot a loader fills from a matrix file takes exactly one role and one scale. A file is
# put into a slot only through read_for_slot(), which refuses any other role or scale, and also
# refuses, for a raw-count slot, a Matrix Market file whose comment lines say "NOT raw counts"
# (scrust marks every normalised matrix so), and, for a slot of SCTransform's layers, one whose
# comment lines do not say it.

.SLOTS <- list(
  `Seurat RNA counts`                      = list(role = "raw_counts", format = "matrix_market", scale = "raw_counts", raw = TRUE),
  `Seurat SCT counts`                      = list(role = "sct_counts", format = "matrix_market", scale = "corrected_counts", normalised = TRUE),
  `Seurat SCT data`                        = list(role = "sct_data", format = "matrix_market", scale = "log1p_corrected_counts", normalised = TRUE),
  `Seurat SCT scale.data`                  = list(role = "sct_scale_data", format = "npy", scale = "residuals"),
  `Seurat pca cell.embeddings`             = list(role = "pca_embeddings", format = "tsv", scale = "pca_embeddings"),
  `Seurat pca feature.loadings`            = list(role = "pca_loadings", format = "tsv", scale = "pca_loadings"),
  `Seurat pca stdev`                       = list(role = "pca_stdev", format = "tsv", scale = "pca_stdev"),
  `Seurat graph SCT_nn`                    = list(role = "graph_nn", format = "matrix_market", scale = "graph_nn"),
  `Seurat graph SCT_snn`                   = list(role = "graph_snn", format = "matrix_market", scale = "graph_snn"),
  `SingleCellExperiment counts`            = list(role = "raw_counts", format = "matrix_market", scale = "raw_counts", raw = TRUE),
  `SingleCellExperiment altExp SCT corrected_counts` = list(role = "sct_counts", format = "matrix_market", scale = "corrected_counts", normalised = TRUE),
  `SingleCellExperiment altExp SCT logcounts` = list(role = "sct_data", format = "matrix_market", scale = "log1p_corrected_counts", normalised = TRUE),
  `SingleCellExperiment altExp SCT_scale_data scaledata` = list(role = "sct_scale_data", format = "npy", scale = "residuals"),
  `SingleCellExperiment reducedDim PCA`    = list(role = "pca_embeddings", format = "tsv", scale = "pca_embeddings"),
  `SingleCellExperiment reducedDim PCA rotation` = list(role = "pca_loadings", format = "tsv", scale = "pca_loadings"),
  `SingleCellExperiment reducedDim PCA stdev` = list(role = "pca_stdev", format = "tsv", scale = "pca_stdev"),
  `SingleCellExperiment colPair SCT_nn`    = list(role = "graph_nn", format = "matrix_market", scale = "graph_nn"),
  `SingleCellExperiment colPair SCT_snn`   = list(role = "graph_snn", format = "matrix_market", scale = "graph_snn")
)

#' The slots the R loaders fill from matrix files, with the one role and scale each takes
#'
#' @return a data frame: slot, role, scale, format, and whether the slot holds raw counts.
slot_table <- function() {
  data.frame(slot = names(.SLOTS),
             role = vapply(.SLOTS, `[[`, "", "role"),
             scale = vapply(.SLOTS, `[[`, "", "scale"),
             format = vapply(.SLOTS, `[[`, "", "format"),
             raw_counts = vapply(.SLOTS, function(s) isTRUE(s$raw), logical(1)),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Whether `slot` takes a file of this role, scale and format; an error naming them if not.
.check_slot <- function(slot, role, scale, path = "a file", format = NULL) {
  s <- .SLOTS[[slot]]
  if (is.null(s)) stop("scrustio: no slot named ", slot, call. = FALSE)
  if (!is.null(format) && !identical(format, s$format)) {
    stop("scrustio: refusing to put ", path, " (format ", format, ") into ", slot, ", which takes only ",
         s$format, " files", call. = FALSE)
  }
  if (!identical(role, s$role) || !identical(scale, s$scale)) {
    stop("scrustio: refusing to put ", path, " (role ", role, ", scale ", scale, ") into ", slot,
         ", which takes only role ", s$role, " with scale ", s$scale,
         if (isTRUE(s$raw)) ": raw counts come only from the raw counts (counts/ or singlets/)" else "",
         call. = FALSE)
  }
  invisible(TRUE)
}

# The content rule, on a Matrix Market file's comment lines.
.check_comments <- function(slot, comments, path) {
  s <- .SLOTS[[slot]]
  marked <- any(grepl("NOT raw counts", comments, fixed = TRUE))
  if (isTRUE(s$raw) && marked) {
    stop("scrustio: refusing to put ", path, " into ", slot, ": its comment lines say it is NOT raw counts",
         call. = FALSE)
  }
  if (isTRUE(s$normalised) && !marked) {
    stop("scrustio: refusing to put ", path, " into ", slot, ": scrust marks SCTransform's layers ",
         "\"NOT raw counts\" in their comment lines, and this file is not so marked", call. = FALSE)
  }
  invisible(TRUE)
}

#' Read a listed file for a slot, refusing it if the slot is not meant for it
#'
#' @param m a manifest (read_manifest()).
#' @param entry the entry to read (or a role held by exactly one file).
#' @param slot a slot of slot_table().
#' @return what read_entry() returns, for matrix files with their row and column names.
read_for_slot <- function(m, entry, slot) {
  e <- if (is.character(entry)) .entry(m, entry) else entry
  .check_slot(slot, e$role, e$scale, e$path, e$format)
  if (identical(e$format, "matrix_market")) {
    r <- .read_mm(.full_path(m, e))
    .check_comments(slot, r$comments, e$path)
    .check_shape(m, e)
    mm <- r$matrix
    if (!is.null(e$shape) && !identical(as.numeric(dim(mm)), e$shape[1:2])) {
      stop("scrustio: ", e$path, " is ", nrow(mm), " x ", ncol(mm), "; the manifest's shape is ",
           paste(e$shape, collapse = " x "), call. = FALSE)
    }
    if (!is.null(e$nnz) && length(mm@x) != e$nnz) {
      stop("scrustio: ", e$path, " stores ", length(mm@x), " entries; the manifest lists ", e$nnz, call. = FALSE)
    }
    dimnames(mm) <- list(.identifiers(m, e$rows), .identifiers(m, e$columns))
    return(mm)
  }
  read_entry(m, e)
}
