# M8's checks on plain files (docs/PLAN.md C4; milestone M8): a reference run against the saved
# object it must reproduce, and two runs against each other. Reads no .rds; large objects run it
# through validation/slurm/compare_clustering.sbatch.
#
#   Rscript validation/R/compare_clustering.R <object> <a> <b> <out.json>
#
# <a> and <b> are result directories in lib/clustering.R's layout (cl_write_object): an export's
# saved/ (data/clustering/<object>/saved) or a run directory (data/clustering/<object>/runs/<run>).
# Files are compared byte for byte: text files as written; .npy files as bytes (if they differ:
# elements differing, largest absolute difference, largest difference relative to the array's
# largest absolute value; the PCA also after aligning each component's sign); Matrix Market files
# on their content without comment lines (if they differ: the entries, the structure and the
# values, values compared as written). metadata.tsv is also compared column by column. For level
# light, the log1p data and scale.data are compared through the digests in object.json (saved/) or
# summary.json's result (runs). When both are instrumented runs, every stage array is compared
# block by block (stages/index.tsv), with every RNG checkpoint, every scalar of stages.json and
# the Annoy index file; the floor statistics (C4) say how many float items, neighbour lists,
# graph entries and labels changed.
#
# Output: <out.json>, with a "summary" section (one value per check) that report_clustering_m8.R
# tabulates.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: compare_clustering.R <object> <a> <b> <out.json>")
object <- args[1]
A <- args[2]
B <- args[3]
out_path <- args[4]
suppressPackageStartupMessages(library("Matrix"))

facts_of <- function(dir) {
  if (file.exists(file.path(dir, "object.json"))) return(jsonlite::fromJSON(file.path(dir, "object.json"), simplifyVector = FALSE))
  jsonlite::fromJSON(file.path(dir, "summary.json"), simplifyVector = FALSE)$result
}
fa <- facts_of(A)
fb <- facts_of(B)

same_bytes <- function(f) {
  pa <- file.path(A, f)
  pb <- file.path(B, f)
  if (!file.exists(pa) || !file.exists(pb)) return(NA)
  unname(tools::md5sum(pa) == tools::md5sum(pb)) && identical(sha256_file(pa), sha256_file(pb))
}

mtx_content_sha <- function(path) {
  out <- system(sprintf("zcat %s | grep -v '^%%' | sha256sum", shQuote(path)), intern = TRUE)
  sub(" .*", "", out)
}

mtx_read_text <- function(path) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  line <- readLines(con, n = 1)
  while (startsWith(line, "%")) line <- readLines(con, n = 1)
  dims <- as.numeric(strsplit(trimws(line), " +")[[1]])
  tr <- scan(con, what = list(0L, 0L, ""), nmax = dims[3], quiet = TRUE)
  list(dims = dims, i = tr[[1]], j = tr[[2]], x = tr[[3]])
}

cmp_mtx <- function(f) {
  pa <- file.path(A, f)
  pb <- file.path(B, f)
  if (!file.exists(pa) || !file.exists(pb)) return(list(identical = NA))
  if (identical(mtx_content_sha(pa), mtx_content_sha(pb))) return(list(identical = TRUE))
  ma <- mtx_read_text(pa)
  mb <- mtx_read_text(pb)
  same_structure <- identical(ma$dims, mb$dims) && identical(ma$i, mb$i) && identical(ma$j, mb$j)
  res <- list(identical = FALSE, dims_a = ma$dims, dims_b = mb$dims, same_structure = same_structure)
  if (same_structure) {
    d <- ma$x != mb$x
    res$entries_differing <- sum(d)
    if (any(d)) {
      xa <- as.numeric(ma$x[d])
      xb <- as.numeric(mb$x[d])
      res$max_abs_diff <- max(abs(xa - xb))
      res$example <- head(paste(ma$i[d], ma$j[d], ma$x[d], mb$x[d]), 5)
    }
  } else {
    ka <- paste(ma$i, ma$j)
    kb <- paste(mb$i, mb$j)
    res$entries_only_a <- sum(!ka %in% kb)
    res$entries_only_b <- sum(!kb %in% ka)
    common <- intersect(ka, kb)
    res$entries_common_value_differing <- sum(ma$x[match(common, ka)] != mb$x[match(common, kb)])
  }
  res
}

cmp_npy <- function(f, signs = FALSE) {
  pa <- file.path(A, f)
  pb <- file.path(B, f)
  if (!file.exists(pa) || !file.exists(pb)) return(list(identical = NA))
  if (isTRUE(same_bytes(f))) return(list(identical = TRUE))
  a <- npy_read(pa)
  b <- npy_read(pb)
  res <- list(identical = FALSE, same_shape = identical(dim(a), dim(b)) && length(a) == length(b))
  if (!res$same_shape) return(res)
  scale <- max(abs(a))
  d <- a != b
  res$elements_differing <- sum(d)
  res$elements <- length(a)
  res$max_abs_diff <- max(abs(a - b))
  res$max_scaled_diff <- if (scale > 0) res$max_abs_diff / scale else NA
  if (signs && !is.null(dim(a))) {
    flips <- 0L
    for (k in seq_len(ncol(a))) {
      if (sum((a[, k] + b[, k])^2) < sum((a[, k] - b[, k])^2)) {
        b[, k] <- -b[, k]
        flips <- flips + 1L
      }
    }
    res$components_flipped <- flips
    res$aligned_max_abs_diff <- max(abs(a - b))
    res$aligned_max_scaled_diff <- if (scale > 0) res$aligned_max_abs_diff / scale else NA
  }
  res
}

cmp_text <- function(f) list(identical = same_bytes(f))

# metadata.tsv, column by column (doubles compared as written: 17 significant digits, exact)
meta_cmp <- function() {
  ma <- read.delim(file.path(A, "metadata.tsv"), colClasses = "character", quote = "", comment.char = "")
  mb <- read.delim(file.path(B, "metadata.tsv"), colClasses = "character", quote = "", comment.char = "")
  if (!identical(ma$barcode, mb$barcode)) return(list(same_cells = FALSE))
  cols <- setdiff(intersect(names(ma), names(mb)), "barcode")
  out <- list(same_cells = TRUE, n_cells = nrow(ma))
  for (col in cols) {
    d <- ma[[col]] != mb[[col]]
    r <- list(differing = sum(d))
    if (any(d) && col %in% c("silhouette_width", "neighborhood_purity", "nCount_SCT")) {
      xa <- as.numeric(ma[[col]][d])
      xb <- as.numeric(mb[[col]][d])
      r$max_abs_diff <- max(abs(xa - xb))
    }
    out[[col]] <- r
  }
  # labels: agreement beyond identity (adjusted Rand index), for the floor
  la <- ma$seurat_clusters
  lb <- mb$seurat_clusters
  out$seurat_clusters$ari <- if (identical(la, lb)) 1 else {
    tab <- table(la, lb)
    comb2 <- function(x) x * (x - 1) / 2
    s <- sum(comb2(tab))
    ra <- sum(comb2(rowSums(tab)))
    rb <- sum(comb2(colSums(tab)))
    e <- ra * rb / comb2(nrow(ma))
    (s - e) / ((ra + rb) / 2 - e)
  }
  out
}

level_of <- function(f) if (is.null(f$level)) NA else f$level
res <- list(object = object, a = A, b = B, level_a = level_of(fa), level_b = level_of(fb))
res$metadata <- meta_cmp()
for (f in c("seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt", "idents_levels.txt")) res$levels[[f]] <- cmp_text(f)
light <- any(c(level_of(fa), level_of(fb)) %in% c("labels"))
if (!light) {
  res$sct <- list(
    corrected_counts = cmp_mtx("sct/corrected_counts.mtx.gz"),
    log1p = cmp_mtx("sct/log1p_corrected_counts.mtx.gz"),
    features = list(identical = identical(mtx_content_sha(file.path(A, "sct/sct_features.tsv.gz")),
                                          mtx_content_sha(file.path(B, "sct/sct_features.tsv.gz")))),
    scale_data = cmp_npy("sct/scale_data.npy"),
    scale_data_features = cmp_text("sct/scale_data_features.tsv"),
    variable_features = cmp_text("sct/variable_features.tsv"),
    feature_attributes = cmp_text("sct/model/feature_attributes.tsv"),
    cell_attributes = cmp_text("sct/model/cell_attributes.tsv"),
    model = cmp_text("sct/model/model.json"),
    digests = list(
      counts_x = identical(fa$sct_counts_x_sha256, fb$sct_counts_x_sha256),
      data_x = identical(fa$sct_data_x_sha256, fb$sct_data_x_sha256),
      scale_data = identical(fa$scale_data_sha256, fb$scale_data_sha256)),
    data_is_log1p_counts = c(isTRUE(fa$sct_data_is_log1p_counts), isTRUE(fb$sct_data_is_log1p_counts)))
  res$pca <- list(embeddings = cmp_npy("pca/embeddings.npy", signs = TRUE),
                  loadings = cmp_npy("pca/loadings.npy", signs = TRUE),
                  stdev = cmp_npy("pca/stdev.npy"),
                  loadings_features = cmp_text("pca/loadings_features.tsv"),
                  json = cmp_text("pca/pca.json"))
  res$graphs <- list(nn = cmp_mtx("graphs/nn.mtx.gz"), snn = cmp_mtx("graphs/snn.mtx.gz"),
                     json = cmp_text("graphs/graphs.json"))
}

# Stage arrays, when both are instrumented runs ---------------------------------------------------
sa <- file.path(A, "stages")
sb <- file.path(B, "stages")
if (dir.exists(sa) && dir.exists(sb)) {
  ia <- read.delim(file.path(sa, "index.tsv"), colClasses = "character")
  ib <- read.delim(file.path(sb, "index.tsv"), colClasses = "character")
  st <- list(same_index = identical(ia[, c("file", "ctx", "row_start", "nrow", "ncol")],
                                    ib[, c("file", "ctx", "row_start", "nrow", "ncol")]))
  arrays <- list()
  for (f in unique(intersect(ia$file, ib$file))) {
    same <- identical(sha256_file(file.path(sa, f)), sha256_file(file.path(sb, f)))
    r <- list(identical = same, stage = ia$stage[match(f, ia$file)])
    if (!same) {
      a <- npy_read(file.path(sa, f))
      b <- npy_read(file.path(sb, f))
      if (length(a) == length(b)) {
        r$elements <- length(a)
        r$elements_differing <- sum(a != b)
        if (is.double(a)) {
          r$max_abs_diff <- max(abs(a - b))
          s <- max(abs(a))
          r$max_scaled_diff <- if (s > 0) r$max_abs_diff / s else NA
        }
        if (!is.null(dim(a)) && length(dim(a)) == 2) r$rows_differing <- sum(rowSums(a != b) > 0)
        if (f %in% c("S3_embeddings.npy", "S3_loadings.npy")) {
          for (k in seq_len(ncol(a))) if (sum((a[, k] + b[, k])^2) < sum((a[, k] - b[, k])^2)) b[, k] <- -b[, k]
          r$aligned_max_abs_diff <- max(abs(a - b))
        }
      } else {
        r$lengths <- c(length(a), length(b))
      }
    }
    arrays[[f]] <- r
  }
  st$arrays <- arrays
  rng_a <- npy_read(file.path(sa, "rng.npy"))
  rng_b <- npy_read(file.path(sb, "rng.npy"))
  la <- read.delim(file.path(sa, "rng.tsv"), colClasses = "character")
  lb <- read.delim(file.path(sb, "rng.tsv"), colClasses = "character")
  st$rng_labels_identical <- identical(la$label, lb$label)
  common <- intersect(la$label, lb$label)
  st$rng_checkpoints_identical <- vapply(common, function(l) identical(rng_a[match(l, la$label), ], rng_b[match(l, lb$label), ]), logical(1))
  ja <- jsonlite::fromJSON(file.path(sa, "stages.json"), simplifyVector = FALSE)
  jb <- jsonlite::fromJSON(file.path(sb, "stages.json"), simplifyVector = FALSE)
  keys <- union(names(ja$scalars$chain), names(jb$scalars$chain))
  st$scalars_differing <- I(keys[!vapply(keys, function(k) identical(ja$scalars$chain[[k]], jb$scalars$chain[[k]]), logical(1))])
  st$annoy_index_identical <- identical(sha256_file(file.path(sa, "annoy.ann")), sha256_file(file.path(sb, "annoy.ann")))
  # floor statistics: float items, neighbour lists, labels before GroupSingletons
  if (file.exists(file.path(sa, "K3_items.npy"))) {
    ka <- npy_read(file.path(sa, "K3_items.npy"))
    kb <- npy_read(file.path(sb, "K3_items.npy"))
    st$floor <- list(items_differing = sum(rowSums(ka != kb) > 0), items = nrow(ka),
                     item_values_differing = sum(ka != kb))
    # A component can change sign between kernels; that changes no distance. Counted again after
    # aligning each column's sign (float(-x) is -float(x), so a pure flip leaves nothing).
    flips <- 0L
    for (k in seq_len(ncol(ka))) {
      if (sum((ka[, k] + kb[, k])^2) < sum((ka[, k] - kb[, k])^2)) {
        kb[, k] <- -kb[, k]
        flips <- flips + 1L
      }
    }
    st$floor$item_components_flipped <- flips
    st$floor$items_differing_aligned <- sum(rowSums(ka != kb) > 0)
    st$floor$item_values_differing_aligned <- sum(ka != kb)
    st$floor$item_max_abs_diff_aligned <- max(abs(ka - kb))
    ca <- npy_read(file.path(sa, "K9_neighbor_count.npy"))
    cb <- npy_read(file.path(sb, "K9_neighbor_count.npy"))
    st$floor$purity_neighbor_counts_differing <- sum(ca != cb)
    # the PCA dims the neighbour search, silhouette and purity read (1:30), after aligning signs
    ea <- npy_read(file.path(sa, "S3_embeddings.npy"))[, 1:30, drop = FALSE]
    eb <- npy_read(file.path(sb, "S3_embeddings.npy"))[, 1:30, drop = FALSE]
    for (k in seq_len(ncol(ea))) if (sum((ea[, k] + eb[, k])^2) < sum((ea[, k] - eb[, k])^2)) eb[, k] <- -eb[, k]
    st$floor$embeddings_1_30_max_abs_diff_aligned <- max(abs(ea - eb))
    st$floor$embeddings_1_30_max_scaled_diff_aligned <- max(abs(ea - eb)) / max(abs(ea))
    na <- npy_read(file.path(sa, "K5_nn_index.npy"))
    nb <- npy_read(file.path(sb, "K5_nn_index.npy"))
    da <- npy_read(file.path(sa, "K5_nn_dist.npy"))
    db <- npy_read(file.path(sb, "K5_nn_dist.npy"))
    st$floor$neighbour_lists_differing <- sum(rowSums(na != nb) > 0)
    st$floor$neighbour_sets_differing <- sum(vapply(seq_len(nrow(na)), function(i) !setequal(na[i, ], nb[i, ]), logical(1)))
    st$floor$neighbour_distances_differing <- sum(rowSums(da != db) > 0)
    oa <- npy_read(file.path(sa, "K7_optimizer_ids.npy"))
    ob <- npy_read(file.path(sb, "K7_optimizer_ids.npy"))
    st$floor$optimizer_ids_differing <- sum(oa != ob)
  }
  res$stages <- st
}

# Optimizer numbers and summary results ------------------------------------------------------------
if (!is.null(fa$optimizer) && !is.null(fb$optimizer)) res$optimizer_identical <- identical(fa$optimizer, fb$optimizer)
if (!is.null(fa$rng_after) && !is.null(fb$rng_after)) {
  common <- intersect(names(fa$rng_after), names(fb$rng_after))
  res$rng_after_identical <- vapply(common, function(k) identical(fa$rng_after[[k]], fb$rng_after[[k]]), logical(1))
}
if (!file.exists(file.path(A, "object.json")) && !file.exists(file.path(B, "object.json"))) {
  res$result_identical <- identical(fa, fb)
}

# Summary ------------------------------------------------------------------------------------------
m <- res$metadata
s <- list(
  cells = m$n_cells,
  labels_identical = isTRUE(m$same_cells) && m$seurat_clusters$differing == 0 && m$`SCT_snn_res.0.8`$differing == 0,
  levels_identical = all(vapply(res$levels, function(x) isTRUE(x$identical), logical(1))),
  ari = m$seurat_clusters$ari,
  silhouette_differing = m$silhouette_width$differing,
  purity_differing = m$neighborhood_purity$differing,
  nCount_SCT_differing = m$nCount_SCT$differing,
  nFeature_SCT_differing = m$nFeature_SCT$differing)
if (!is.null(res$sct)) {
  s$sct_counts_identical <- isTRUE(res$sct$corrected_counts$identical)
  s$sct_data_identical <- if (!is.na(res$sct$log1p$identical)) isTRUE(res$sct$log1p$identical) else
    isTRUE(res$sct$digests$data_x) && all(res$sct$data_is_log1p_counts)
  s$sct_scale_data_identical <- if (!is.na(res$sct$scale_data$identical)) isTRUE(res$sct$scale_data$identical) else
    isTRUE(res$sct$digests$scale_data)
  s$sct_model_identical <- all(vapply(res$sct[c("variable_features", "feature_attributes", "cell_attributes", "model")],
                                      function(x) isTRUE(x$identical), logical(1)))
  s$pca_identical <- all(vapply(res$pca[c("embeddings", "loadings", "stdev")], function(x) isTRUE(x$identical), logical(1)))
  s$pca_embeddings_max_abs_diff <- if (isTRUE(res$pca$embeddings$identical)) 0 else res$pca$embeddings$aligned_max_abs_diff
  s$nn_graph_identical <- isTRUE(res$graphs$nn$identical)
  s$snn_graph_identical <- isTRUE(res$graphs$snn$identical)
}
if (!is.null(res$stages)) {
  s$stage_arrays <- length(res$stages$arrays)
  s$stage_arrays_identical <- sum(vapply(res$stages$arrays, function(x) isTRUE(x$identical), logical(1)))
  s$stage_arrays_differing <- I(names(res$stages$arrays)[!vapply(res$stages$arrays, function(x) isTRUE(x$identical), logical(1))])
  s$rng_checkpoints_identical <- isTRUE(res$stages$rng_labels_identical) && all(res$stages$rng_checkpoints_identical)
  s$annoy_index_identical <- res$stages$annoy_index_identical
  if (!is.null(res$stages$floor)) s <- c(s, res$stages$floor)
}
if (!is.null(res$result_identical)) s$result_identical <- res$result_identical
if (!is.null(res$optimizer_identical)) s$optimizer_identical <- res$optimizer_identical
res$summary <- s
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
json_write(res, out_path)
cat(object, basename(A), "vs", basename(B), ": labels", s$labels_identical, "silhouette differing", s$silhouette_differing,
    "purity differing", s$purity_differing, if (!is.null(s$sct_counts_identical)) paste("sct counts", s$sct_counts_identical,
    "pca", s$pca_identical, "snn", s$snn_graph_identical) else "", "\n")
