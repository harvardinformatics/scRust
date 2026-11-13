# Every portable-profile run against its saved object, output by output (docs/PLAN.md A5b and C4;
# M8's R-vs-R floor, docs/validation/clustering_m8.md). The saved objects were made with OpenBLAS's
# SapphireRapids kernel, so only production-profile runs can reproduce them bit for bit (A5b,
# K8); the portable runs' differences from them are measured here. Reads only plain files.
#
#   Rscript validation/R/floor_clustering.R <out.tsv>
#
# For every object of validation/clustering_datasets.tsv and every portable run in
# data/clustering/<object>/runs/*-portable-*, each output against the saved object (saved/):
# identical or not and, if not, how many values differ, the largest absolute difference and the
# largest difference relative to the saved value. The outputs: the labels (seurat_clusters and
# SCT_snn_res.0.8), the three levels files, silhouette_width, neighborhood_purity, nCount_SCT,
# nFeature_SCT, the SCT corrected counts and their log1p, scale.data, the SCTModel's feature and
# cell attributes and model.json, the variable features, the PCA's embeddings and loadings (each
# component's sign aligned first; the components that changed sign are counted) and stdev, and
# both graphs. Text values are compared as written (17 significant digits, so bitwise); the sizes
# of the differences use R's parser, which is accurate enough for that. Where scale.data was kept
# only as a digest (level light), it comes from floor_scale_data.R's measurement
# (data/clustering/<object>/floor/scale_data_portable.json), whose digests must equal the run's
# and the saved object's. Every row also carries the PCA's own largest difference (signs aligned),
# over the 30 dimensions the neighbour search, silhouette and purity read and over all 50.
#
# Writes one row per run and output to <out.tsv>, checks that the portable runs of each object
# differ from the saved object in exactly the same way (they are byte-identical to each other), and
# prints the report's tables: first the runs whose labels, levels, graphs or purities differ, then
# one row per object.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: floor_clustering.R <out.tsv>")

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
rj <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
sha_cache <- new.env()
content_sha <- function(path) {  # a Matrix Market file's content without comment lines, memoised
  key <- normalizePath(path)
  if (is.null(sha_cache[[key]])) {
    assign(key, sub(" .*", "", system(sprintf("zcat %s | grep -v '^%%' | sha256sum", shQuote(path)), intern = TRUE)), envir = sha_cache)
  }
  sha_cache[[key]]
}
same_file <- function(a, b) file.exists(a) && file.exists(b) && identical(sha256_file(a), sha256_file(b))

res_row <- function(identical, values, differing = 0, max_abs = 0, max_rel = 0, note = "") {
  list(identical = identical, values = values, differing = differing, max_abs = max_abs, max_rel = max_rel, note = note)
}
cmp_text <- function(a, b) {  # saved, run: values as written
  if (length(a) != length(b)) stop("lengths differ")
  d <- a != b
  if (!any(d)) return(res_row(TRUE, length(a)))
  x <- as.numeric(a[d])
  y <- as.numeric(b[d])
  diff <- abs(x - y)
  nz <- x != 0
  res_row(FALSE, length(a), sum(d), max(diff), if (any(nz)) max(diff[nz] / abs(x[nz])) else NA)
}
cmp_num <- function(a, b, note = "") {  # saved, run: doubles
  if (length(a) != length(b)) stop("lengths differ")
  d <- a != b
  if (!any(d)) return(res_row(TRUE, length(a), note = note))
  diff <- abs(a[d] - b[d])
  nz <- a[d] != 0
  res_row(FALSE, length(a), sum(d), max(diff), if (any(nz)) max(diff[nz] / abs(a[d][nz])) else NA, note)
}
align <- function(a, b) {  # flips b's columns to a's sign; returns b and the number flipped
  flips <- 0L
  for (k in seq_len(ncol(a))) {
    if (sum((a[, k] + b[, k])^2) < sum((a[, k] - b[, k])^2)) {
      b[, k] <- -b[, k]
      flips <- flips + 1L
    }
  }
  list(b = b, flips = flips)
}
cmp_tsv_numeric <- function(pa, pb) {  # attribute tables: every numeric column, as written
  a <- read.delim(pa, colClasses = "character", quote = "", check.names = FALSE, na.strings = character())
  b <- read.delim(pb, colClasses = "character", quote = "", check.names = FALSE, na.strings = character())
  if (!identical(names(a), names(b)) || !identical(a[[1]], b[[1]])) stop("tables differ in shape: ", pa)
  cols <- names(a)[-1]
  vals <- 0
  diffs <- list()
  for (col in cols) {
    if (all(a[[col]] %in% c("TRUE", "FALSE"))) {
      vals <- vals + nrow(a)
      if (any(a[[col]] != b[[col]])) diffs[[col]] <- res_row(FALSE, nrow(a), sum(a[[col]] != b[[col]]), NA, NA)
      next
    }
    r <- cmp_text(a[[col]], b[[col]])
    vals <- vals + r$values
    if (!r$identical) diffs[[col]] <- r
  }
  if (!length(diffs)) return(res_row(TRUE, vals))
  res_row(FALSE, vals, sum(vapply(diffs, `[[`, 0, "differing")),
          max(vapply(diffs, function(x) if (is.na(x$max_abs)) 0 else x$max_abs, 0)),
          max(vapply(diffs, function(x) if (is.na(x$max_rel)) 0 else x$max_rel, 0)),
          paste("columns", paste(names(diffs), collapse = ", ")))
}
cmp_mtx <- function(pa, pb) {
  if (identical(content_sha(pa), content_sha(pb))) {
    n <- as.numeric(strsplit(system(sprintf("zcat %s | grep -v '^%%' | head -n 1", shQuote(pa)), intern = TRUE), " +")[[1]][3])
    return(res_row(TRUE, n))
  }
  read_mtx <- function(path) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    line <- readLines(con, n = 1)
    while (startsWith(line, "%")) line <- readLines(con, n = 1)
    nnz <- as.numeric(strsplit(trimws(line), " +")[[1]][3])
    scan(con, what = list(0L, 0L, ""), nmax = nnz, quiet = TRUE)
  }
  a <- read_mtx(pa)
  b <- read_mtx(pb)
  if (!identical(a[[1]], b[[1]]) || !identical(a[[2]], b[[2]])) {
    return(res_row(FALSE, length(a[[3]]), NA, NA, NA, "the stored entries differ"))
  }
  cmp_text(a[[3]], b[[3]])
}

rows <- list()
for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  base <- file.path("data/clustering", obj)
  saved <- file.path(base, "saved")
  if (!dir.exists(saved)) next
  sf <- rj(file.path(saved, "object.json"))
  ms <- read.delim(file.path(saved, "metadata.tsv"), colClasses = "character", quote = "", check.names = FALSE, na.strings = character())
  es <- npy_read(file.path(saved, "pca", "embeddings.npy"))
  ls_ <- npy_read(file.path(saved, "pca", "loadings.npy"))
  runs <- list.dirs(file.path(base, "runs"), recursive = FALSE, full.names = FALSE)
  runs <- sort(runs[grepl("-portable-", runs) & !grepl("\\.partial\\.", runs)])
  for (r in runs) {
    rd <- file.path(base, "runs", r)
    rf <- rj(file.path(rd, "summary.json"))$result
    mr <- read.delim(file.path(rd, "metadata.tsv"), colClasses = "character", quote = "", check.names = FALSE, na.strings = character())
    if (!identical(ms$barcode, mr$barcode)) stop(obj, " ", r, ": cells differ")
    out <- list()
    out$labels <- cmp_text(c(ms$seurat_clusters, ms$`SCT_snn_res.0.8`), c(mr$seurat_clusters, mr$`SCT_snn_res.0.8`))
    lv <- c("seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt", "idents_levels.txt")
    out$levels <- res_row(all(vapply(lv, function(f) same_file(file.path(saved, f), file.path(rd, f)), logical(1))),
                          length(readLines(file.path(saved, lv[1]))))
    for (col in c("silhouette_width", "neighborhood_purity", "nCount_SCT", "nFeature_SCT")) out[[col]] <- cmp_text(ms[[col]], mr[[col]])
    out$sct_corrected_counts <- cmp_mtx(file.path(saved, "sct", "corrected_counts.mtx.gz"), file.path(rd, "sct", "corrected_counts.mtx.gz"))
    l1 <- file.path(c(saved, rd), "sct", "log1p_corrected_counts.mtx.gz")
    out$sct_log1p <- if (all(file.exists(l1))) cmp_mtx(l1[1], l1[2]) else
      res_row(identical(sf$sct_data_x_sha256, rf$sct_data_x_sha256), sf$sct_nnz, note = "compared through the digests of the stored doubles")
    sd1 <- file.path(c(saved, rd), "sct", "scale_data.npy")
    if (all(file.exists(sd1))) {
      out$sct_scale_data <- if (same_file(sd1[1], sd1[2])) res_row(TRUE, prod(unlist(sf$scale_data_dim))) else cmp_num(npy_read(sd1[1]), npy_read(sd1[2]))
    } else if (identical(sf$scale_data_sha256, rf$scale_data_sha256)) {
      out$sct_scale_data <- res_row(TRUE, prod(unlist(sf$scale_data_dim)), note = "compared through the digests of the doubles")
    } else {
      fj <- file.path(base, "floor", "scale_data_portable.json")
      if (!file.exists(fj)) stop(obj, " ", r, ": scale.data differs and ", fj, " is missing")
      m <- rj(fj)$scale_data
      if (!identical(m$run_sha256, rf$scale_data_sha256) || !identical(m$saved_sha256, sf$scale_data_sha256)) {
        stop(obj, " ", r, ": the measurement's scale.data is not the run's or the saved object's")
      }
      out$sct_scale_data <- res_row(FALSE, m$values, m$differing, m$max_abs, m$max_rel,
                                    "measured by floor_scale_data.R (the same SCTransform call; digests match the run's)")
    }
    out$sct_feature_attributes <- cmp_tsv_numeric(file.path(saved, "sct", "model", "feature_attributes.tsv"), file.path(rd, "sct", "model", "feature_attributes.tsv"))
    out$sct_cell_attributes <- cmp_tsv_numeric(file.path(saved, "sct", "model", "cell_attributes.tsv"), file.path(rd, "sct", "model", "cell_attributes.tsv"))
    out$sct_model_json <- res_row(same_file(file.path(saved, "sct", "model", "model.json"), file.path(rd, "sct", "model", "model.json")), 1)
    out$variable_features <- res_row(same_file(file.path(saved, "sct", "variable_features.tsv"), file.path(rd, "sct", "variable_features.tsv")),
                                     sf$n_variable_features)
    er <- npy_read(file.path(rd, "pca", "embeddings.npy"))
    al <- align(es, er)
    out$pca_embeddings <- cmp_num(es, al$b, sprintf("signs aligned; %d of %d components had changed sign", al$flips, ncol(es)))
    lr <- npy_read(file.path(rd, "pca", "loadings.npy"))
    all_l <- align(ls_, lr)
    out$pca_loadings <- cmp_num(ls_, all_l$b, sprintf("signs aligned; %d components had changed sign", all_l$flips))
    out$pca_stdev <- cmp_num(npy_read(file.path(saved, "pca", "stdev.npy")), npy_read(file.path(rd, "pca", "stdev.npy")))
    for (g in c("nn", "snn")) out[[paste0("graph_", g)]] <- cmp_mtx(file.path(saved, "graphs", paste0(g, ".mtx.gz")), file.path(rd, "graphs", paste0(g, ".mtx.gz")))
    pca_30 <- max(abs(es[, 1:30] - al$b[, 1:30]))
    pca_all <- max(abs(es - al$b))
    for (o in names(out)) {
      x <- out[[o]]
      rows[[length(rows) + 1]] <- data.frame(object = obj, tier = ds$tier[i], run = r, output = o,
        identical = x$identical, values = x$values, differing = x$differing, max_abs_diff = x$max_abs,
        max_rel_diff = x$max_rel, pca_max_abs_diff_dims_1_30 = pca_30, pca_max_abs_diff_all = pca_all,
        note = x$note, stringsAsFactors = FALSE)
    }
    message(obj, " ", r, ": ", sum(!vapply(out, `[[`, TRUE, "identical")), " outputs differ")
  }
}
tab <- do.call(rbind, rows)
tsv_write(tab, args[1])

# The portable runs of an object must differ from the saved object in the same way
key <- c("identical", "values", "differing", "max_abs_diff", "max_rel_diff", "pca_max_abs_diff_dims_1_30", "pca_max_abs_diff_all")
for (o in unique(tab$object)) {
  t <- tab[tab$object == o, ]
  sig <- vapply(split(t, t$run), function(d) paste(unlist(d[order(d$output), key]), collapse = "|"), "")
  if (length(unique(sig)) != 1) stop(o, ": the portable runs differ from the saved object in different ways")
}

# Tables for the report
g <- function(x, d = 2) formatC(as.numeric(x), digits = d, format = "g")
cell <- function(d) {
  if (isTRUE(d$identical)) return("=")
  if (is.na(d$differing)) return("differs")
  sprintf("%s of %s; %s (%s)", formatC(d$differing, format = "d", big.mark = ","), formatC(d$values, format = "d", big.mark = ","),
          g(d$max_abs_diff), if (is.na(d$max_rel_diff)) "-" else g(d$max_rel_diff))
}
run_desc <- function(runs) {  # "4: direct, instrumented; cascadelake, sapphirerapids" when they cross
  modes <- unique(sub("-.*", "", runs))
  archs <- unique(sub(".*-portable-", "", runs))
  if (length(runs) == length(modes) * length(archs)) {
    return(sprintf("%d: %s; %s", length(runs), paste(modes, collapse = ", "), paste(archs, collapse = ", ")))
  }
  sprintf("%d: %s", length(runs), paste(runs, collapse = ", "))
}
md_table <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
                          apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
n_runs <- length(unique(paste(tab$object, tab$run)))
first <- tab[tab$output %in% c("labels", "levels", "graph_nn", "graph_snn", "neighborhood_purity") & !tab$identical, ]
lines <- c(sprintf("Portable runs whose labels, levels, graphs or purities differ from the saved object: %d of %d (one row per object; its runs differ in the same way):",
                   length(unique(paste(first$object, first$run))), n_runs), "")
if (nrow(first)) {
  keys <- unique(first[, c("object", "output")])
  ft <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    d <- first[first$object == keys$object[k] & first$output == keys$output[k], ]
    data.frame(object = d$object[1], `portable runs` = run_desc(unique(d$run)), output = d$output[1],
               `values differing` = sprintf("%d of %s", d$differing[1], formatC(d$values[1], format = "d", big.mark = ",")),
               `max abs` = g(d$max_abs_diff[1]), `max rel` = g(d$max_rel_diff[1]),
               `PCA max abs, dims 1:30` = g(d$pca_max_abs_diff_dims_1_30[1]), check.names = FALSE)
  }))
  lines <- c(lines, md_table(ft))
}
same_everywhere <- names(which(tapply(tab$identical, tab$output, all)))
lines <- c(lines, "", sprintf("In all %d portable runs these outputs are identical to the saved object: %s.",
                              n_runs, paste(same_everywhere, collapse = ", ")))
outs <- setdiff(unique(tab$output), same_everywhere)
per <- lapply(unique(tab$object), function(o) {
  t <- tab[tab$object == o, ]
  runs <- unique(t$run)
  r1 <- t[t$run == runs[1], ]
  c(list(object = o, `portable runs` = run_desc(runs)),
    setNames(lapply(outs, function(x) cell(r1[r1$output == x, ])), outs),
    list(`PCA max abs, dims 1:30` = g(r1$pca_max_abs_diff_dims_1_30[1])))
})
pt <- do.call(rbind, lapply(per, as.data.frame, check.names = FALSE, stringsAsFactors = FALSE))
lines <- c(lines, "", "Every portable run, by object (its portable runs differ from the saved object in exactly the same way; each cell: values differing of all; largest absolute difference (largest difference relative to the saved value); = when identical):", "",
           md_table(pt))
writeLines(lines)
