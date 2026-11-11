# Comparisons for DoubletFinder's reference data (docs/PLAN.md A4, A5b and B6, milestone M2).
# Reads only plain files (no .rds), so small comparisons may run anywhere; large stage files
# are better compared in a batch job.
#
#   Rscript validation/R/compare_doubletfinder.R crosscheck <run dir> <export dir> <out.json>
#       A reference run against the saved workflow output (export_doubletfinder.R's saved/):
#       the doublet set, the singlet pANN (bitwise, compared as 17-digit text) and the column
#       names, hence pK and nExp.
#   Rscript validation/R/compare_doubletfinder.R runs <run dir A> <run dir B> <out prefix>
#       Two runs of the same object: whether the final outputs are byte-identical (cells.tsv,
#       sweep_stats.tsv, bcmvn.tsv) and the summary.json "result" sections equal; and, when both
#       have stages/, every stage array, block by block (context), bitwise: number of differing
#       elements, max absolute and relative difference. Writes <prefix>.json and
#       <prefix>_stages.tsv. Used for instrumented vs direct (same profile and CPU: must be
#       identical), the same profile on different CPU types (must be identical), and portable
#       vs production (the R-vs-R floor per stage).
#   Rscript validation/R/compare_doubletfinder.R sha256 <out.tsv> <dir> [<dir> ...]
#       sha256 and size of every file under the given directories.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))

read_tsv_chr <- function(path) read.delim(path, colClasses = "character", quote = "", comment.char = "")

crosscheck <- function(run_dir, export_dir, out) {
  run <- read_tsv_chr(file.path(run_dir, "cells.tsv"))
  res <- jsonlite::fromJSON(file.path(run_dir, "summary.json"))$result
  saved <- read_tsv_chr(file.path(export_dir, "saved", "cells.tsv"))
  saved_info <- jsonlite::fromJSON(file.path(export_dir, "saved", "saved.json"))
  saved_doublets <- readLines(file.path(export_dir, "saved", "doublets.txt"))
  run_doublets <- run$barcode[run$classification == "Doublet"]
  m <- match(saved$barcode, run$barcode)
  if (anyNA(m)) stop("saved cells missing from the run")
  pann_equal <- run$pANN[m] == saved$pANN
  d <- abs(as.numeric(run$pANN[m]) - as.numeric(saved$pANN))
  r <- list(
    run = normalizePath(run_dir),
    export = normalizePath(export_dir),
    pANN_column_run = res$pANN_column,
    pANN_column_saved = saved_info$pANN_column,
    classification_column_run = res$classification_column,
    classification_column_saved = saved_info$classification_column,
    column_names_identical = identical(res$pANN_column, saved_info$pANN_column) &&
      identical(res$classification_column, saved_info$classification_column),
    pK_run = res$pK_label,
    pK_saved = saved_info$pK_label,
    nExp_run = res$nExp,
    nExp_saved = saved_info$nExp,
    n_cells = nrow(run),
    n_doublets_run = length(run_doublets),
    n_doublets_saved = length(saved_doublets),
    doublet_set_identical = identical(sort(run_doublets), sort(saved_doublets)),
    doublets_only_in_run = length(setdiff(run_doublets, saved_doublets)),
    doublets_only_in_saved = length(setdiff(saved_doublets, run_doublets)),
    n_singlets = nrow(saved),
    n_singlet_pANN_identical = sum(pann_equal),
    singlet_pANN_max_abs_diff = if (length(d)) max(d) else 0
  )
  r$pass <- r$column_names_identical && r$doublet_set_identical && r$n_singlet_pANN_identical == r$n_singlets
  json_write(r, out)
  cat(sprintf("%s: columns %s, doublet set %s (%d), singlet pANN %d/%d identical -> %s\n",
              basename(dirname(dirname(run_dir))), r$column_names_identical, r$doublet_set_identical,
              r$n_doublets_run, r$n_singlet_pANN_identical, r$n_singlets, if (r$pass) "PASS" else "FAIL"))
  invisible(r)
}

array_block <- function(x, start, nrow) {
  if (is.null(dim(x))) return(x[start + seq_len(nrow)])
  x[start + seq_len(nrow), , drop = FALSE]
}

diff_stats <- function(a, b) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b)) {
    return(list(n = NA, n_diff = NA, max_abs = NA, max_rel = NA, max_scaled = NA, note = "shape differs"))
  }
  if (is.integer(a) && is.integer(b)) {
    nd <- sum(a != b)
    return(list(n = length(a), n_diff = nd, max_abs = if (nd) max(abs(as.double(a) - as.double(b))) else 0,
                max_rel = NA, max_scaled = NA, note = ""))
  }
  a <- as.double(a); b <- as.double(b)
  same <- (a == b & !is.na(a) & !is.na(b)) | (is.na(a) & is.na(b) & is.nan(a) == is.nan(b))
  # bitwise: also distinguish -0 from 0
  same <- same & !(a == 0 & b == 0 & (1 / a) != (1 / b))
  fin <- is.finite(a) & is.finite(b)
  ad <- abs(a - b)[fin]
  rel <- (abs(a - b) / pmax(abs(a), abs(b)))[fin & !(a == 0 & b == 0)]
  scale <- if (any(fin)) max(abs(a[fin]), abs(b[fin])) else 0
  list(n = length(a), n_diff = sum(!same), max_abs = if (length(ad)) max(ad) else 0,
       max_rel = if (length(rel)) max(rel) else 0,
       max_scaled = if (length(ad) && scale > 0) max(ad) / scale else 0,
       note = if (sum(!same & !fin)) sprintf("%d non-finite mismatches", sum(!same & !fin)) else "")
}

compare_runs <- function(a, b, prefix) {
  files <- c("cells.tsv", "sweep_stats.tsv", "bcmvn.tsv")
  same_files <- vapply(files, function(f) identical(sha256_file(file.path(a, f)), sha256_file(file.path(b, f))), logical(1))
  ra <- jsonlite::fromJSON(file.path(a, "summary.json"))$result
  rb <- jsonlite::fromJSON(file.path(b, "summary.json"))$result
  out <- list(a = normalizePath(a), b = normalizePath(b), final_files_identical = as.list(same_files),
              result_identical = identical(ra, rb))
  # Calls and pANN, element-wise, when the files differ
  ca <- read_tsv_chr(file.path(a, "cells.tsv")); cb <- read_tsv_chr(file.path(b, "cells.tsv"))
  out$n_cells <- nrow(ca)
  out$n_pANN_differ <- sum(ca$pANN != cb$pANN)
  out$max_abs_pANN_diff <- max(abs(as.numeric(ca$pANN) - as.numeric(cb$pANN)))
  out$n_calls_differ <- sum(ca$classification != cb$classification)
  out$pK <- c(ra$pK_label, rb$pK_label)
  out$nExp <- c(ra$nExp, rb$nExp)
  sa <- file.path(a, "stages"); sb <- file.path(b, "stages")
  if (dir.exists(sa) && dir.exists(sb)) {
    ia <- read_tsv_chr(file.path(sa, "index.tsv"))
    ib <- read_tsv_chr(file.path(sb, "index.tsv"))
    # blocks are matched by file, context and occurrence
    occ <- function(k) ave(seq_along(k), k, FUN = seq_along)
    key_a <- paste(ia$file, ia$ctx); key_a <- paste(key_a, occ(key_a))
    key_b <- paste(ib$file, ib$ctx); key_b <- paste(key_b, occ(key_b))
    rows <- list()
    cache <- list(file = NULL)
    for (f in unique(ia$file)) {
      if (!f %in% ib$file) {
        rows[[length(rows) + 1]] <- data.frame(file = f, ctx = "*", stage = ia$stage[ia$file == f][1],
          n = NA, n_diff = NA, max_abs = NA, max_rel = NA, max_scaled = NA, note = "missing in B")
        next
      }
      xa <- npy_read(file.path(sa, f)); xb <- npy_read(file.path(sb, f))
      for (i in which(ia$file == f)) {
        j <- match(key_a[i], key_b)
        if (is.na(j)) {
          rows[[length(rows) + 1]] <- data.frame(file = f, ctx = ia$ctx[i], stage = ia$stage[i],
            n = NA, n_diff = NA, max_abs = NA, max_rel = NA, max_scaled = NA, note = "context missing in B")
          next
        }
        ba <- array_block(xa, as.integer(ia$row_start[i]), as.integer(ia$nrow[i]))
        bb <- array_block(xb, as.integer(ib$row_start[j]), as.integer(ib$nrow[j]))
        s <- diff_stats(ba, bb)
        if (f %in% c("S3_embeddings.npy", "S3_loadings.npy") && identical(dim(ba), dim(bb))) {
          # PCA components are defined up to sign (plan B3b, B5 S3): align each column of B
          flip <- colSums(ba * bb) < 0
          if (any(flip)) {
            s <- diff_stats(ba, sweep(bb, 2, ifelse(flip, -1, 1), `*`))
            s$note <- paste0(sum(flip), " component(s) sign-flipped, aligned")
          }
        }
        rows[[length(rows) + 1]] <- data.frame(file = f, ctx = ia$ctx[i], stage = ia$stage[i],
          n = s$n, n_diff = s$n_diff, max_abs = s$max_abs, max_rel = s$max_rel, max_scaled = s$max_scaled, note = s$note)
      }
      rm(xa, xb)
    }
    tab <- do.call(rbind, rows)
    # RNG checkpoints, by label
    ra_ <- npy_read(file.path(sa, "rng.npy")); rb_ <- npy_read(file.path(sb, "rng.npy"))
    la <- read_tsv_chr(file.path(sa, "rng.tsv")); lb <- read_tsv_chr(file.path(sb, "rng.tsv"))
    out$rng_checkpoints_identical <- identical(ra_, rb_) && identical(la, lb)
    out$rng_checkpoints <- nrow(la)
    out$stage_arrays <- nrow(tab)
    out$stage_arrays_identical <- sum(tab$n_diff == 0, na.rm = TRUE)
    out$stage_arrays_differing <- sum(is.na(tab$n_diff) | tab$n_diff > 0)
    tsv_write(tab, paste0(prefix, "_stages.tsv"))
  }
  json_write(out, paste0(prefix, ".json"))
  cat(sprintf("%s vs %s: final files identical %s, result identical %s, calls differ %d, pANN differ %d%s\n",
              basename(a), basename(b), paste(same_files, collapse = "/"), out$result_identical,
              out$n_calls_differ, out$n_pANN_differ,
              if (!is.null(out$stage_arrays)) sprintf(", stage arrays identical %d/%d, RNG checkpoints identical %s",
                                                      out$stage_arrays_identical, out$stage_arrays, out$rng_checkpoints_identical) else ""))
  invisible(out)
}

sha256_dirs <- function(out, dirs) {
  rows <- lapply(dirs, function(d) {
    f <- list.files(d, recursive = TRUE, full.names = TRUE)
    data.frame(dir = d, file = substring(f, nchar(d) + 2), bytes = file.size(f),
               sha256 = vapply(f, sha256_file, ""), stringsAsFactors = FALSE)
  })
  tsv_write(do.call(rbind, rows), out)
}

a <- commandArgs(trailingOnly = TRUE)
switch(a[1],
  crosscheck = crosscheck(a[2], a[3], a[4]),
  runs = compare_runs(a[2], a[3], a[4]),
  sha256 = sha256_dirs(a[2], a[-(1:2)]),
  stop("usage: compare_doubletfinder.R crosscheck|runs|sha256 ..."))
