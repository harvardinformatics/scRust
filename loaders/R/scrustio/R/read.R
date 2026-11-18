# Readers for the formats of docs/output_manifest.md, in base R and Matrix.
#
# Exactness: scrust writes doubles in text with 17 significant digits (C's %.17g), so a correctly
# rounded reader recovers every double. R's own parser (R_strtod, used by scan() and read.delim())
# accumulates the digits in an 80-bit long double and divides by an exact power of ten; for at
# most 17 significant digits of a normal double that leaves the result within the half-ulp the
# 17 digits guarantee, so it returns the double that was written. The loaders' comparisons with
# the saved objects (milestone M14's K11) check this on every value of every reference object.
# .npy files are read as bytes.

# Lines of a text file, gzipped or not.
.read_lines <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

# One column of a TSV file, as text: by header name (header = TRUE) or 1-based number.
.read_id_column <- function(path, header, column) {
  if (!header && (identical(column, 1) || identical(column, 1L))) {
    first <- .read_lines(path)
    return(sub("\t.*$", "", first))
  }
  df <- utils::read.delim(path, header = header, sep = "\t", quote = "", comment.char = "",
                          colClasses = "character", check.names = FALSE, na.strings = character(),
                          stringsAsFactors = FALSE)
  if (is.character(column)) {
    if (!column %in% names(df)) stop("scrustio: ", path, " has no column ", column, call. = FALSE)
    return(df[[column]])
  }
  if (column > ncol(df)) stop("scrustio: ", path, " has no column ", column, call. = FALSE)
  df[[column]]
}

# A TSV file with a header line. `numeric` names the columns read as doubles (all others are
# text); `integer` those read as integers. "NA" is read as NA; the identifier columns (barcode,
# id, name, pc) must have none.
.read_tsv <- function(path, numeric = character(), integer = character()) {
  hdr <- strsplit(.read_lines_n(path, 1), "\t", fixed = TRUE)[[1]]
  classes <- rep("character", length(hdr))
  classes[hdr %in% numeric] <- "numeric"
  classes[hdr %in% integer] <- "integer"
  df <- utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                          colClasses = classes, check.names = FALSE, na.strings = "NA",
                          stringsAsFactors = FALSE)
  if (!identical(names(df), hdr)) stop("scrustio: ", path, ": unexpected header", call. = FALSE)
  for (col in intersect(hdr, c("barcode", "id", "name", "pc"))) {
    if (anyNA(df[[col]])) stop("scrustio: ", path, ": NA in identifier column ", col, call. = FALSE)
  }
  df
}

.read_lines_n <- function(path, n) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  readLines(con, n = n, warn = FALSE)
}

# Matrix Market (coordinate, general): the header, its comment lines and the entries. The
# entries become a dgCMatrix with the entries in the order stored (column-major, rows increasing
# within a column, as scrust and Cell Ranger write them), explicit zeros kept, as Seurat holds
# them; any other order is sorted (no position may repeat). Values are parsed as doubles
# (integer fields exactly).
.read_mm <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  banner <- readLines(con, n = 1, warn = FALSE)
  w <- strsplit(trimws(banner), "[[:space:]]+")[[1]]
  if (length(w) < 5 || !identical(w[1], "%%MatrixMarket") || tolower(w[2]) != "matrix" ||
      tolower(w[3]) != "coordinate" || tolower(w[5]) != "general") {
    stop("scrustio: ", path, ": unsupported Matrix Market header: ", banner, call. = FALSE)
  }
  field <- tolower(w[4])
  if (!field %in% c("integer", "real", "pattern")) stop("scrustio: ", path, ": unsupported field ", field, call. = FALSE)
  comments <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("scrustio: ", path, ": no size line", call. = FALSE)
    if (!startsWith(line, "%")) break
    comments <- c(comments, line)
  }
  dims <- as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]])
  if (length(dims) != 3 || anyNA(dims)) stop("scrustio: ", path, ": bad size line: ", line, call. = FALSE)
  nr <- dims[1]; nc <- dims[2]; nnz <- dims[3]
  what <- if (field == "pattern") list(0L, 0L) else list(0L, 0L, 0)
  tr <- scan(con, what = what, nmax = nnz, quiet = TRUE)
  if (length(tr[[1]]) != nnz) stop("scrustio: ", path, ": ", nnz, " entries announced, ", length(tr[[1]]), " read", call. = FALSE)
  i <- tr[[1]]; j <- tr[[2]]
  x <- if (field == "pattern") rep(1, nnz) else tr[[3]]
  rm(tr)
  if (nnz && (min(i) < 1 || max(i) > nr || min(j) < 1 || max(j) > nc)) {
    stop("scrustio: ", path, ": an entry lies outside the ", nr, " x ", nc, " matrix", call. = FALSE)
  }
  sorted <- !is.unsorted(j)
  if (sorted && nnz > 1) {
    same <- j[-1] == j[-nnz]
    sorted <- !any(i[-1][same] <= i[-nnz][same])
  }
  if (!sorted) {
    o <- order(j, i)
    i <- i[o]; j <- j[o]; x <- x[o]
    if (nnz > 1 && any(i[-1] == i[-nnz] & j[-1] == j[-nnz])) {
      stop("scrustio: ", path, ": a position is stored twice", call. = FALSE)
    }
  }
  p <- c(0L, cumsum(tabulate(j, nbins = nc)))
  m <- methods::new("dgCMatrix", i = i - 1L, p = as.integer(p), x = as.double(x),
                    Dim = as.integer(c(nr, nc)))
  list(matrix = m, field = field, comments = comments, stored_order = sorted)
}

# A .npy file (format 1.0 to 3.0, little-endian <f8, <f4 or <i4, any order): a vector, or a
# matrix in R's column-major layout. Read in pieces, so an array larger than 2^31 values loads.
.read_npy <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 6)
  if (!identical(magic, as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59)))) stop("scrustio: ", path, ": not a .npy file", call. = FALSE)
  ver <- as.integer(readBin(con, "raw", 2))
  hlen <- if (ver[1] == 1L) {
    readBin(con, "integer", size = 2, signed = FALSE, endian = "little")
  } else {
    readBin(con, "integer", size = 4, endian = "little")
  }
  header <- rawToChar(readBin(con, "raw", hlen))
  descr <- sub(".*'descr': *'([^']+)'.*", "\\1", header)
  fortran <- grepl("'fortran_order': *True", header)
  shape_txt <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
  shape <- as.numeric(strsplit(gsub(" ", "", shape_txt), ",")[[1]])
  shape <- shape[!is.na(shape)]
  n <- if (length(shape)) prod(shape) else 1
  spec <- switch(descr, "<f8" = list("double", 8L), "<f4" = list("double", 4L), "<i4" = list("integer", 4L),
                 stop("scrustio: ", path, ": unsupported dtype ", descr, call. = FALSE))
  x <- if (spec[[1]] == "double") double(n) else integer(n)
  piece <- 2^27
  done <- 0
  while (done < n) {
    k <- min(piece, n - done)
    v <- readBin(con, spec[[1]], k, size = spec[[2]], endian = "little")
    if (length(v) != k) stop("scrustio: ", path, ": the file ends early", call. = FALSE)
    x[(done + 1):(done + k)] <- v
    done <- done + k
  }
  if (length(shape) > 1) {
    if (fortran) {
      dim(x) <- shape
    } else {
      x <- aperm(array(x, rev(shape)))
    }
  }
  x
}

#' Read one listed file
#'
#' Reads a file the manifest lists, by its format: a dgCMatrix (matrix_market), a data frame
#' (tsv), an array (npy), a list (json) or lines (text), with the row and column names the
#' manifest's `rows` and `columns` point to, after checking their number against `shape`. No
#' slot is involved: read_for_slot() puts a matrix into a slot.
#'
#' @param m a manifest (read_manifest()).
#' @param entry an entry of the manifest, or a role held by exactly one file.
read_entry <- function(m, entry) {
  e <- if (is.character(entry)) .entry(m, entry) else entry
  .check_shape(m, e)
  path <- .full_path(m, e)
  out <- switch(e$format,
    matrix_market = {
      r <- .read_mm(path)
      mm <- r$matrix
      if (!is.null(e$shape) && !identical(as.numeric(dim(mm)), e$shape[1:2])) {
        stop("scrustio: ", e$path, " is ", nrow(mm), " x ", ncol(mm), "; the manifest's shape is ",
             paste(e$shape, collapse = " x "), call. = FALSE)
      }
      if (!is.null(e$nnz) && length(mm@x) != e$nnz) {
        stop("scrustio: ", e$path, " stores ", length(mm@x), " entries; the manifest lists ", e$nnz, call. = FALSE)
      }
      dimnames(mm) <- list(.identifiers(m, e$rows), .identifiers(m, e$columns))
      mm
    },
    tsv = {
      if (isTRUE(e$header)) {
        utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                          check.names = FALSE, stringsAsFactors = FALSE)
      } else {
        utils::read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                          colClasses = "character", stringsAsFactors = FALSE)
      }
    },
    npy = {
      a <- .read_npy(path)
      if (!is.null(e$shape) && !identical(as.numeric(if (is.null(dim(a))) length(a) else dim(a)), e$shape)) {
        stop("scrustio: ", e$path, ": the array's shape differs from the manifest's", call. = FALSE)
      }
      if (!is.null(dim(a)) && length(dim(a)) == 2) dimnames(a) <- list(.identifiers(m, e$rows), .identifiers(m, e$columns))
      a
    },
    json = jsonlite::fromJSON(path, simplifyVector = FALSE),
    text = .read_lines(path),
    stop("scrustio: ", e$path, ": format ", e$format, " is not read by this loader", call. = FALSE))
  out
}
