# Shared pieces of milestone M14's loader checks (docs/PLAN.md C3g, C4 "Loader tests", C5's K11
# and K13): the loader's library, the R-vs-R floor of validation/tolerances.toml, comparison
# statistics and grading. Requires lib/io.R.

# scrustio from the repository's own library (.rlib/, never the R env), first on the path.
ld_setup <- function(root = Sys.getenv("SCRUST_ROOT")) {
  if (!nzchar(root)) stop("SCRUST_ROOT is not set (source env.sh)")
  lib <- file.path(root, ".rlib")
  if (!dir.exists(file.path(lib, "scrustio"))) stop("scrustio is not installed in ", lib, " (R CMD INSTALL --library=.rlib loaders/R/scrustio)")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages({
    library("Matrix")
    library("scrustio")
  })
  invisible(lib)
}

# [clustering.floor.kernel] of validation/tolerances.toml: one inline table per array,
# `name = { key = value, ... }`, read into a named list of named lists (numbers as doubles,
# strings without quotes, booleans as logicals).
ld_floor <- function(path, section = "clustering.floor.kernel") {
  lines <- readLines(path, warn = FALSE)
  start <- which(trimws(lines) == paste0("[", section, "]"))
  if (length(start) != 1) stop("no [", section, "] in ", path)
  out <- list()
  for (line in lines[(start + 1):length(lines)]) {
    t <- trimws(line)
    if (startsWith(t, "[")) break
    if (!nzchar(t) || startsWith(t, "#")) next
    mm <- regmatches(t, regexec("^([A-Za-z0-9_.]+)[[:space:]]*=[[:space:]]*\\{(.*)\\}[[:space:]]*$", t))[[1]]
    if (!length(mm)) next
    body <- mm[3]
    # split on commas outside double quotes
    chars <- strsplit(body, "")[[1]]
    inq <- FALSE; cur <- ""; parts <- character()
    for (ch in chars) {
      if (ch == "\"") inq <- !inq
      if (ch == "," && !inq) { parts <- c(parts, cur); cur <- "" } else cur <- paste0(cur, ch)
    }
    parts <- c(parts, cur)
    entry <- list()
    for (p in parts) {
      kv <- regmatches(p, regexec("^[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*=[[:space:]]*(.*?)[[:space:]]*$", p))[[1]]
      if (!length(kv)) next
      v <- kv[3]
      entry[[kv[2]]] <- if (grepl("^\".*\"$", v)) substr(v, 2, nchar(v) - 1) else if (v %in% c("true", "false")) v == "true" else as.numeric(v)
    }
    out[[mm[2]]] <- entry
  }
  out
}

# Statistics of a numeric comparison of a (the loader's) with b (the reference): values, values
# differing (NA patterns included), the largest absolute difference, that relative to b's
# largest absolute value, and the fraction differing. With align = TRUE (PCA), each column of a
# is first given the sign that brings it closer to b's, as the floor was measured.
ld_num_stats <- function(a, b, align = FALSE) {
  a <- unclass(a); b <- unclass(b)
  flips <- 0L
  if (align && length(dim(a)) == 2 && identical(dim(a), dim(b))) {
    for (k in seq_len(ncol(a))) {
      if (sum((a[, k] + b[, k])^2) < sum((a[, k] - b[, k])^2)) {
        a[, k] <- -a[, k]
        flips <- flips + 1L
      }
    }
  }
  av <- as.double(a); bv <- as.double(b)
  if (length(av) != length(bv)) return(list(n = length(bv), comparable = FALSE))
  na_diff <- is.na(av) != is.na(bv)
  both <- !is.na(av) & !is.na(bv)
  d <- abs(av[both] - bv[both])
  differing <- sum(na_diff) + sum(av[both] != bv[both])
  scale <- if (any(both)) max(abs(bv[both])) else 0
  max_abs <- if (length(d)) max(d) else 0
  list(n = length(bv), comparable = TRUE, differing = differing, max_abs = max_abs,
       max_scaled = if (scale > 0) max_abs / scale else if (max_abs > 0) Inf else 0,
       fraction_differing = if (length(bv)) differing / length(bv) else 0, flipped = flips)
}

# A sparse matrix comparison: dimensions, dimnames, stored structure (i, p) and values.
ld_sparse_stats <- function(a, b) {
  a <- methods::as(a, "dgCMatrix"); b <- methods::as(b, "dgCMatrix")
  same_dim <- identical(dim(a), dim(b))
  same_names <- identical(dimnames(a), dimnames(b))
  same_struct <- same_dim && identical(a@i, b@i) && identical(a@p, b@p)
  out <- list(n = length(b@x), same_dim = same_dim, same_names = same_names, same_structure = same_struct)
  if (same_struct) out <- c(out, ld_num_stats(a@x, b@x)[c("differing", "max_abs", "max_scaled", "fraction_differing")])
  out
}

# Grades one comparison. mode "exact": identical or failing. mode "floor": identical, or, for a
# BLAS-dependent component with an entry in the floor, within every bound of that entry (keys
# starting with max_, values_differing and cells_differing); otherwise failing.
ld_grade <- function(identical, stats, floor_entry, mode) {
  if (isTRUE(identical)) return("identical")
  if (mode != "floor" || is.null(floor_entry) || is.null(stats) || !isTRUE(stats$comparable %||% TRUE)) return("differing")
  ok <- TRUE
  for (k in names(floor_entry)) {
    v <- floor_entry[[k]]
    if (!is.numeric(v)) next
    got <- switch(k,
                  max_abs = stats$max_abs, max_scaled = stats$max_scaled,
                  max_fraction_differing = stats$fraction_differing,
                  values_differing = stats$differing, cells_differing = stats$differing,
                  NULL)
    if (is.null(got)) next
    if (!is.finite(got) || got > v) ok <- FALSE
  }
  if (ok) "within floor" else "beyond floor"
}

`%||%` <- function(a, b) if (is.null(a)) b else a

ld_vmhwm_kib <- function() {
  st <- readLines("/proc/self/status")
  as.numeric(sub("^VmHWM:\\s+([0-9]+) kB$", "\\1", grep("^VmHWM:", st, value = TRUE)))
}
