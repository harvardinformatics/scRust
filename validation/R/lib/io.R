# Plain-file I/O for scRust's reference data (docs/PLAN.md A3 and A4), in base R plus Matrix.
# Sourced by the export, reference and instrumented scripts in validation/R/.
#
# Exactness rules:
# - Doubles written as text use sprintf("%.17g"): glibc's printf rounds correctly, and 17
#   significant digits identify every double, so a correctly rounded reader (Rust's
#   str::parse::<f64>) recovers the exact value. Integer-valued doubles print as integers.
#   R's own decimal parser (R_strtod) is NOT correctly rounded: it accumulates digits in a long
#   double and scales by powers of ten, so reading "%.17g" text back into R can be off by one
#   ULP for non-integers. Exactness checks here therefore compare the written text with
#   sprintf("%.17g") of the original (an injective map), never R-parsed values, except for
#   integer-valued numbers, which R parses exactly.
# - Dense arrays and RNG states go to .npy (format 1.0, little-endian, Fortran order for
#   matrices), written byte for byte by npy_write().

# Numbers to text ---------------------------------------------------------------------------

fmt17 <- function(x) {
  x <- as.double(x)
  out <- sprintf("%.17g", x)
  out[which(is.na(x) & !is.nan(x))] <- "NA"
  out[which(is.nan(x))] <- "NaN"
  out[which(x == Inf)] <- "Inf"
  out[which(x == -Inf)] <- "-Inf"
  out
}

# .npy (NumPy format 1.0) ---------------------------------------------------------------------

# Writes a numeric/integer/logical vector or matrix. dtype "<f8" (double) or "<i4" (integer;
# logicals become 0/1). Matrices are written in Fortran (column-major) order, R's own layout,
# with shape (nrow, ncol); vectors have shape (n,). The header is padded with spaces to a
# multiple of 64 bytes, as NumPy does.
npy_write <- function(path, x, dtype = NULL) {
  if (is.logical(x)) {
    if (anyNA(x)) stop("npy_write: NA in logical array ", path)
    storage.mode(x) <- "integer"
  }
  if (is.null(dtype)) dtype <- if (is.integer(x)) "<i4" else "<f8"
  if (dtype == "<i4") {
    if (!is.integer(x)) {
      if (any(!is.finite(x)) || any(x != round(x)) || any(abs(x) > .Machine$integer.max)) {
        stop("npy_write: values not representable as int32 in ", path)
      }
      storage.mode(x) <- "integer"
    }
    if (anyNA(x)) stop("npy_write: NA in int32 array ", path)
  } else if (dtype == "<f8") {
    storage.mode(x) <- "double"
  } else {
    stop("npy_write: unsupported dtype ", dtype)
  }
  d <- dim(x)
  shape <- if (is.null(d)) {
    sprintf("(%d,)", length(x))
  } else {
    paste0("(", paste(d, collapse = ", "), ")")
  }
  fortran <- if (!is.null(d) && length(d) > 1) "True" else "False"
  header <- sprintf("{'descr': '%s', 'fortran_order': %s, 'shape': %s, }", dtype, fortran, shape)
  # magic (6) + version (2) + header length (2) + header + '\n', padded to a multiple of 64
  total <- 10 + nchar(header) + 1
  pad <- (64 - total %% 64) %% 64
  header <- paste0(header, strrep(" ", pad), "\n")
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59, 0x01, 0x00)), con)
  writeBin(as.integer(nchar(header)), con, size = 2, endian = "little")
  writeBin(charToRaw(header), con)
  if (dtype == "<f8") {
    writeBin(as.vector(x), con, size = 8, endian = "little")
  } else {
    writeBin(as.vector(x), con, size = 4, endian = "little")
  }
  invisible(path)
}

# Reads what npy_write writes (used to verify the files and to compare runs).
npy_read <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 8)
  if (!identical(magic[1:6], as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59)))) {
    stop("not a .npy file: ", path)
  }
  hlen <- readBin(con, "integer", size = 2, signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", hlen))
  dtype <- sub(".*'descr': *'([^']+)'.*", "\\1", header)
  fortran <- grepl("'fortran_order': *True", header)
  shape_txt <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
  shape <- as.integer(strsplit(gsub(" ", "", shape_txt), ",")[[1]])
  shape <- shape[!is.na(shape)]
  n <- if (length(shape)) prod(shape) else 1
  x <- switch(dtype,
    "<f8" = readBin(con, "double", n, size = 8, endian = "little"),
    "<i4" = readBin(con, "integer", n, size = 4, endian = "little"),
    stop("npy_read: unsupported dtype ", dtype))
  if (length(shape) > 1) {
    if (fortran) {
      dim(x) <- shape
    } else {
      x <- aperm(array(x, rev(shape)))
    }
  }
  x
}

# TSV and JSON --------------------------------------------------------------------------------

# Data frame to TSV with a header line; double columns via fmt17(), integer and character
# columns as is. Fails on tabs or newlines inside fields.
tsv_write <- function(df, path) {
  cols <- lapply(df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.double(v)) return(fmt17(v))
    if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
    v <- as.character(v)
    if (any(grepl("[\t\n\r]", v))) stop("tsv_write: tab or newline in a field of ", path)
    v[is.na(v)] <- "NA"
    v
  })
  lines <- c(paste(names(df), collapse = "\t"), do.call(paste, c(cols, sep = "\t")))
  if (nrow(df) == 0) lines <- lines[1]
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

json_escape <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  paste0("\"", s, "\"")
}

# Minimal JSON writer: named lists become objects, unnamed lists and vectors of length != 1
# become arrays (wrap a length-1 vector in I() to force an array), doubles via fmt17 (NA, NaN
# and infinities become null), integers as integers, logicals as true/false, NULL as null.
to_json <- function(x, indent = 0) {
  pad <- strrep("  ", indent)
  pad1 <- strrep("  ", indent + 1)
  scalar <- function(v) {
    if (is.null(v)) return("null")
    if (is.factor(v)) v <- as.character(v)
    if (is.logical(v)) return(ifelse(is.na(v), "null", ifelse(v, "true", "false")))
    if (is.integer(v)) return(ifelse(is.na(v), "null", as.character(v)))
    if (is.double(v)) return(ifelse(is.finite(v), fmt17(v), "null"))
    ifelse(is.na(v), "null", json_escape(as.character(v)))
  }
  if (is.null(x)) return("null")
  if (is.list(x)) {
    if (length(x) == 0) return(if (is.null(names(x))) "[]" else "{}")
    nm <- names(x)
    items <- vapply(seq_along(x), function(i) to_json(x[[i]], indent + 1), "")
    if (!is.null(nm) && all(nzchar(nm))) {
      body <- paste0(pad1, json_escape(nm), ": ", items, collapse = ",\n")
      return(paste0("{\n", body, "\n", pad, "}"))
    }
    return(paste0("[\n", paste0(pad1, items, collapse = ",\n"), "\n", pad, "]"))
  }
  as_array <- inherits(x, "AsIs") || length(x) != 1
  vals <- scalar(unclass(x))
  if (as_array) return(paste0("[", paste(vals, collapse = ", "), "]"))
  vals
}

json_write <- function(x, path) {
  writeLines(to_json(x), path, useBytes = TRUE)
  invisible(path)
}

# 10x-style Matrix Market --------------------------------------------------------------------

# Writes a dgCMatrix as a Cell Ranger v3 style directory: matrix.mtx.gz (coordinate, 1-based,
# entries in column-major order exactly as stored, explicit zeros included), barcodes.tsv.gz
# and features.tsv.gz (id, name, type). The header says "integer" when every value is
# integer-valued, "real" otherwise. Values use fmt17, so integer-valued doubles print as
# integers. Written in chunks of columns to bound memory on large objects. `software` goes in
# the metadata comment line.
mtx_write_10x <- function(m, dir, feature_ids = rownames(m), chunk_nnz = 5e6,
                          software = "scRust validation/R/export_doubletfinder.R") {
  stopifnot(inherits(m, "dgCMatrix"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  x <- m@x
  if (any(!is.finite(x))) stop("mtx_write_10x: non-finite values")
  if (any(x == 0 & 1 / x < 0)) stop("mtx_write_10x: negative zeros")
  field <- if (all(x == round(x))) "integer" else "real"
  con <- gzfile(file.path(dir, "matrix.mtx.gz"), "w")
  writeLines(c(sprintf("%%%%MatrixMarket matrix coordinate %s general", field),
               sprintf("%%metadata_json: {\"software_version\": \"%s\", \"format_version\": 2}", software),
               sprintf("%d %d %d", nrow(m), ncol(m), length(x))), con)
  p <- m@p
  ncol <- ncol(m)
  j <- 1L
  while (j <= ncol) {
    # columns j..k with at most chunk_nnz entries (at least one column)
    k <- j
    while (k < ncol && p[k + 2] - p[j] <= chunk_nnz) k <- k + 1L
    lo <- p[j] + 1
    hi <- p[k + 1]
    if (hi >= lo) {
      cols <- rep.int(j:k, diff(p[j:(k + 1)]))
      writeLines(paste(m@i[lo:hi] + 1L, cols, fmt17(x[lo:hi])), con)
    }
    j <- k + 1L
  }
  close(con)
  con <- gzfile(file.path(dir, "barcodes.tsv.gz"), "w")
  writeLines(colnames(m), con)
  close(con)
  con <- gzfile(file.path(dir, "features.tsv.gz"), "w")
  writeLines(paste(feature_ids, rownames(m), "Gene Expression", sep = "\t"), con)
  close(con)
  invisible(field)
}

# Reads a directory written by mtx_write_10x back into the identical dgCMatrix (slots i, p, x,
# Dim, Dimnames), keeping explicit zeros and the stored entry order. Requires the entries to be
# in column-major order with increasing rows inside each column, as mtx_write_10x writes them.
# Values are parsed by R's scan(), exact for integer-valued counts (all workflow inputs are);
# export_doubletfinder.R checks identical() against the original object.
mtx_read_10x <- function(dir) {
  con <- gzfile(file.path(dir, "matrix.mtx.gz"), "r")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  if (!grepl("^%%MatrixMarket matrix coordinate (integer|real) general", header)) {
    stop("mtx_read_10x: unsupported header: ", header)
  }
  line <- readLines(con, n = 1)
  while (startsWith(line, "%")) line <- readLines(con, n = 1)
  dims <- as.numeric(strsplit(trimws(line), " +")[[1]])
  nr <- dims[1]; nc <- dims[2]; nnz <- dims[3]
  tr <- scan(con, what = list(0L, 0L, 0), nmax = nnz, quiet = TRUE)
  if (length(tr[[1]]) != nnz) stop("mtx_read_10x: expected ", nnz, " entries, read ", length(tr[[1]]))
  i <- tr[[1]]; j <- tr[[2]]; x <- tr[[3]]
  rm(tr)
  if (is.unsorted(j)) stop("mtx_read_10x: entries are not in column-major order")
  p <- c(0L, cumsum(tabulate(j, nbins = nc)))
  # rows strictly increasing inside each column
  if (nnz > 1) {
    same_col <- j[-1] == j[-nnz]
    if (any(i[-1][same_col] <= i[-nnz][same_col])) stop("mtx_read_10x: rows not increasing within a column")
  }
  bcon <- gzfile(file.path(dir, "barcodes.tsv.gz"), "r")
  barcodes <- readLines(bcon)
  close(bcon)
  features <- read.delim(file.path(dir, "features.tsv.gz"), header = FALSE,
                         colClasses = "character", quote = "", comment.char = "")
  if (length(barcodes) != nc || nrow(features) != nr) stop("mtx_read_10x: dimnames do not match the matrix")
  m <- new("dgCMatrix", i = i - 1L, p = as.integer(p), x = x, Dim = as.integer(c(nr, nc)),
           Dimnames = list(features[[2]], barcodes))
  invisible(validObject(m))
  m
}

# Environment, memory, RNG --------------------------------------------------------------------

# Peak resident set size of this R process so far (kernel's VmHWM), in KiB.
peak_rss_kib <- function() {
  st <- readLines("/proc/self/status")
  as.numeric(sub("^VmHWM:\\s+([0-9]+) kB$", "\\1", grep("^VmHWM:", st, value = TRUE)))
}

rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

# The environment fingerprint (plan A5b): `scrust fingerprint` (CPU, OpenBLAS core and threads,
# glibc, node, SLURM job), run as a child of this R process so it sees the same OPENBLAS_*
# variables, plus R's own view: versions, BLAS/LAPACK libraries, locale, RNG kind.
fingerprint <- function(packages = c("Seurat", "SeuratObject", "sctransform", "glmGamPoi",
                                     "irlba", "fields", "KernSmooth", "matrixStats", "Matrix",
                                     "DoubletFinder", "future", "future.apply", "Rcpp")) {
  root <- Sys.getenv("SCRUST_ROOT")
  scrust <- file.path(root, "target", "release", "scrust")
  fp <- if (nzchar(root) && file.exists(scrust)) {
    jsonlite::fromJSON(paste(system2(scrust, "fingerprint", stdout = TRUE), collapse = "\n"))
  } else {
    list(error = "scrust binary not found; build it with `cargo build --release -p scrust`")
  }
  pkg <- lapply(packages, function(p) {
    d <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
    if (is.null(d) || !is.list(d)) return(NULL)
    v <- list(version = d$Version)
    if (!is.null(d$RemoteSha)) v$remote_sha <- d$RemoteSha
    v
  })
  names(pkg) <- packages
  si <- utils::sessionInfo()
  blas <- si$BLAS
  c(fp, list(
    R = R.version.string,
    platform = R.version$platform,
    packages = pkg,
    blas_library = blas,
    lapack_library = si$LAPACK,
    lapack_version = La_version(),
    blas_md5 = if (file.exists(blas)) unname(tools::md5sum(blas)) else NULL,
    rng_kind = RNGkind(),
    locale = Sys.getlocale(),
    sizeof_longdouble = .Machine$sizeof.longdouble
  ))
}
