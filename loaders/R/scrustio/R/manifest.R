# The output manifest: summary.json's `files` section, schema scrust-files version 1
# (docs/output_manifest.md in the scRust repository). Every file is found through it, by its
# role, never by its name.

.SCHEMA <- "scrust-files"
.SCHEMA_VERSION <- 1L

# Roles and scales of schema version 1 (docs/output_manifest.md, "Roles" and "Scales"). A role or
# scale outside these is refused when a loader is asked to put that file into a slot; reading a
# manifest that lists one is allowed (fields and files may be added within a version).
.ROLES <- c("summary", "readme", "raw_counts", "raw_counts_cells", "raw_counts_features",
            "rna_cell_metadata", "cell_metadata", "cluster_labels", "cluster_levels",
            "cluster_columns", "workflow_cluster_ids", "workflow_nclusters", "optimizer_ids",
            "sct_cell_attributes", "sct_feature_attributes", "sct_model_pars", "sct_model_pars_fit",
            "sct_variable_features", "sct_scale_data", "sct_counts", "sct_data", "sct_features",
            "sct_cells", "pca_embeddings", "pca_loadings", "pca_stdev", "rng_state", "knn_index",
            "knn_distance", "graph_nn", "graph_snn", "annoy_index", "doublet_calls",
            "doublet_sweep_stats", "doublet_bcmvn", "stage_dump", "stage_index", "stage_scalars",
            "start_partitions")
.SCALES <- c("raw_counts", "corrected_counts", "log1p_corrected_counts", "residuals",
             "pca_embeddings", "pca_loadings", "pca_stdev", "graph_nn", "graph_snn",
             "neighbour_indices", "neighbour_distances", "none")
.FORMATS <- c("matrix_market", "tsv", "npy", "json", "text", "annoy")

#' Read and check the manifest of a scrust output directory
#'
#' Reads `summary.json`, refuses a manifest whose `schema` is not "scrust-files" or whose
#' `version` is not 1, and checks that every listed file exists with its listed size (and, with
#' `check_sha256 = TRUE`, its SHA-256). Files are then found by role (manifest_entries()).
#'
#' @param dir the output directory (it holds summary.json).
#' @param check_sha256 also check every listed file's SHA-256 (reads every file; a large export
#'   holds gigabytes).
#' @return a `scrust_manifest`: the directory, the command, the provenance, the entries and the
#'   whole summary.json.
read_manifest <- function(dir, check_sha256 = FALSE) {
  if (!dir.exists(dir)) stop("scrustio: no directory ", dir, call. = FALSE)
  dir <- normalizePath(dir, mustWork = TRUE)
  path <- file.path(dir, "summary.json")
  if (!file.exists(path)) stop("scrustio: ", dir, " has no summary.json: not a scrust output", call. = FALSE)
  js <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  f <- js$files
  if (is.null(f)) {
    stop("scrustio: ", path, " has no `files` section (the manifest): written by a scrust older ",
         "than 0.1.0's M13, or not by scrust", call. = FALSE)
  }
  if (!identical(f$schema, .SCHEMA)) {
    stop("scrustio: ", path, ": manifest schema ", deparse(f$schema), ", not \"", .SCHEMA,
         "\"; refusing it", call. = FALSE)
  }
  v <- f$version
  if (!is.numeric(v) || length(v) != 1 || !identical(as.integer(v), .SCHEMA_VERSION) || v != .SCHEMA_VERSION) {
    stop("scrustio: ", path, ": manifest schema version ", deparse(v), "; this loader reads version ",
         .SCHEMA_VERSION, " only. Refusing it", call. = FALSE)
  }
  entries <- f$entries
  if (!is.list(entries) || !length(entries)) stop("scrustio: ", path, ": the manifest lists no files", call. = FALSE)
  seen <- character()
  for (k in seq_along(entries)) {
    e <- entries[[k]]
    for (field in c("path", "role", "scale", "format")) {
      if (!is.character(e[[field]]) || length(e[[field]]) != 1 || !nzchar(e[[field]])) {
        stop("scrustio: ", path, ": entry ", k, " has no valid `", field, "`", call. = FALSE)
      }
    }
    p <- e$path
    if (startsWith(p, "/") || grepl("(^|/)\\.\\.(/|$)", p) || grepl("\\\\", p)) {
      stop("scrustio: ", path, ": entry path ", p, " is not a relative path inside the output", call. = FALSE)
    }
    if (p %in% seen) stop("scrustio: ", path, ": ", p, " is listed twice", call. = FALSE)
    seen <- c(seen, p)
    e$shape <- if (is.null(e$shape)) NULL else as.numeric(unlist(e$shape))
    entries[[k]] <- e
    full <- file.path(dir, p)
    if (!file.exists(full)) stop("scrustio: ", p, " is listed in the manifest but missing from ", dir, call. = FALSE)
    if (!is.null(e$bytes) && file.size(full) != e$bytes) {
      stop("scrustio: ", p, " has ", file.size(full), " bytes; the manifest lists ", e$bytes, call. = FALSE)
    }
  }
  m <- structure(list(dir = dir, command = f$command, provenance = f$provenance, entries = entries,
                      summary = js, cache = new.env(parent = emptyenv())),
                 class = "scrust_manifest")
  if (isTRUE(check_sha256)) .check_sha256(m)
  m
}

.check_sha256 <- function(m, entries = m$entries) {
  for (e in entries) {
    if (is.null(e$sha256)) next
    got <- unname(tools::sha256sum(file.path(m$dir, e$path)))
    if (!identical(got, e$sha256)) {
      stop("scrustio: ", e$path, ": SHA-256 ", got, ", the manifest lists ", e$sha256, call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @export
print.scrust_manifest <- function(x, ...) {
  roles <- vapply(x$entries, `[[`, "", "role")
  cat("scrust output (", x$command, "), ", x$dir, "\n", sep = "")
  cat("manifest: schema ", .SCHEMA, " version ", .SCHEMA_VERSION, ", ", length(roles), " files\n", sep = "")
  tab <- sort(table(roles[!roles %in% c("stage_dump", "stage_index", "stage_scalars")]), decreasing = TRUE)
  cat("roles:", paste0(names(tab), ifelse(tab > 1, paste0(" (", tab, ")"), ""), collapse = ", "), "\n")
  invisible(x)
}

#' The manifest's entries, all of them or those of one role
#'
#' @param m a manifest (read_manifest()).
#' @param role a role (docs/output_manifest.md), or NULL for every entry.
#' @return a list of entries (each a list: path, role, scale, format, gzip, header, shape, nnz,
#'   rows, columns, description, bytes, sha256).
manifest_entries <- function(m, role = NULL) {
  if (is.null(role)) return(m$entries)
  Filter(function(e) identical(e$role, role), m$entries)
}

# The one entry with this role: an error if there is none (when required) or several.
.entry <- function(m, role, required = TRUE) {
  hits <- manifest_entries(m, role)
  if (length(hits) > 1) {
    stop("scrustio: ", length(hits), " files have role ", role, " (",
         paste(vapply(hits, `[[`, "", "path"), collapse = ", "), "); expected one", call. = FALSE)
  }
  if (!length(hits)) {
    if (required) stop("scrustio: this ", m$command, " output has no file with role ", role, call. = FALSE)
    return(NULL)
  }
  hits[[1]]
}

.has_role <- function(m, role) length(manifest_entries(m, role)) > 0

.entry_at <- function(m, path) {
  for (e in m$entries) if (identical(e$path, path)) return(e)
  stop("scrustio: the manifest names ", path, " as an identifier file but does not list it", call. = FALSE)
}

.full_path <- function(m, e) file.path(m$dir, e$path)

# The identifiers a `rows` or `columns` field points to: the listed file, its column by header
# name (header true) or by 1-based number. Cached per file and column.
.identifiers <- function(m, ids) {
  if (is.null(ids)) return(NULL)
  key <- paste0(ids$file, "\r", ids$column)
  if (!is.null(m$cache[[key]])) return(m$cache[[key]])
  e <- .entry_at(m, ids$file)
  out <- .read_id_column(.full_path(m, e), isTRUE(e$header), ids$column)
  assign(key, out, envir = m$cache)
  out
}

#' Check a manifest's identifiers against its shapes (and, if asked, the SHA-256)
#'
#' For every entry with `rows` or `columns`, reads the identifiers and checks that there are as
#' many as the entry's `shape` says. read_entry() and the loaders make the same check for each
#' file they read.
check_manifest <- function(m, sha256 = FALSE) {
  for (e in m$entries) .check_shape(m, e)
  if (isTRUE(sha256)) .check_sha256(m)
  invisible(TRUE)
}

.check_shape <- function(m, e) {
  if (is.null(e$shape)) return(invisible(TRUE))
  if (!is.null(e$rows)) {
    n <- length(.identifiers(m, e$rows))
    if (n != e$shape[1]) {
      stop("scrustio: ", e$path, ": the manifest's shape gives ", e$shape[1], " rows, but its row ",
           "identifiers (", e$rows$file, ", column ", e$rows$column, ") number ", n, call. = FALSE)
    }
  }
  if (!is.null(e$columns)) {
    if (length(e$shape) < 2) stop("scrustio: ", e$path, ": column identifiers but a one-dimensional shape", call. = FALSE)
    n <- length(.identifiers(m, e$columns))
    if (n != e$shape[2]) {
      stop("scrustio: ", e$path, ": the manifest's shape gives ", e$shape[2], " columns, but its column ",
           "identifiers (", e$columns$file, ", column ", e$columns$column, ") number ", n, call. = FALSE)
    }
  }
  invisible(TRUE)
}
