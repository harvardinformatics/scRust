# Reference runs of Seurat 5.5.1's FindNeighbors for milestone M9 (docs/PLAN.md C4, C5 K3-K6 and
# C7): the Annoy metrics M8's instrumented runs did not cover, annoy.metric = "cosine" (RcppAnnoy's
# AnnoyAngular) and "manhattan" (AnnoyManhattan), from the same embeddings; and constructed edge
# cases for all three metrics. No .rds is read. Run through validation/slurm/neighbors_references.sbatch.
#
#   Rscript validation/R/neighbors_references.R --mode real --runs <run dir>[,<run dir>...] \
#     --metrics cosine,manhattan
#   Rscript validation/R/neighbors_references.R --mode edge --out <dir>
#   Rscript validation/R/neighbors_references.R --mode bench --runs <run dir> --out <file.json>
#
# --mode real: for each instrumented run data/clustering/<object>/runs/<run>/ and each metric, the
#   run's pca/embeddings.npy, columns 1:30 (the workflow's dims), barcodes from its metadata.tsv,
#   FindNeighbors' defaults otherwise (k.param 20, n.trees 50, prune.SNN 1/15), written to
#   data/clustering/<object>/neighbors/<run>/<metric>/.
# --mode edge: the constructed inputs of nr_cases() below, each written to <out>/<case>/input/
#   (embeddings.npy, barcodes.tsv), with one reference per metric in <out>/<case>/<metric>/;
#   <out>/cases.tsv lists each case's size and parameters, and <out>/expected.tsv the bytes and
#   sha256 of every reference file, which cargo test compares with Rust's
#   (crates/neighbors/tests/edge_cases.rs).
# --mode bench: the memory and time baseline of C5's performance row for FindNeighbors alone: the
#   workflow's call, FindNeighbors(data, dims 1:30, Euclidean, defaults), on the run's embeddings,
#   timed, with VmHWM before and after; written to --out as JSON (run under /usr/bin/time -v).
#
# Each reference directory has the layout of an instrumented clustering run (lib/clustering_hooks.R),
# so `scrust validate neighbors` reads either:
#   stages/annoy.ann         the index of Seurat's own AnnoyBuildIndex(data, metric, n.trees), saved
#                            with save() (which then reloads it)
#   stages/K3_items.npy      getItemsVector(i) for every item (cells x dims, the floats as doubles)
#   stages/K5_nn_index.npy   FindNeighbors(..., return.neighbor = TRUE): Indices (1-based, int32)
#   stages/K5_nn_dist.npy    and Distances (Annoy's floats as doubles)
#   stages/stages.json       scalars (context "chain"): parameters, sizes and the checks below
#   graphs/nn.mtx.gz, graphs/snn.mtx.gz   FindNeighbors' graphs as stored (lib/clustering.R's
#                            mtx_write_plain), graphs/graphs.json their facts
#   graphs/{nn,snn}_{p,i,x}.npy           the same graphs' dgCMatrix slots, for the digests
#   summary.json             run record: parameters, timings, fingerprint
# Checks recorded in stages.json: the items equal R's own float rounding of the data; FindNeighbors'
# neighbours equal Seurat's AnnoySearch on a second, unsaved build and on the saved (reloaded) index;
# how many cells are first in, and missing from, their own lists. Where FindNeighbors stops (a single
# cell), stages.json records its error and only the index and items are written.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))     # mtx_write_plain

args <- df_args(list(mode = NULL, runs = "", metrics = "cosine,manhattan", out = ""))

suppressPackageStartupMessages(library("Seurat"))

nr_node_bytes <- function(metric, f) if (metric == "cosine") 12 + 4 * f else 16 + 4 * f

# FindNeighbors with these arguments on `data` (cells x dims, with rownames), every reference file.
nr_reference <- function(data, metric, out, k.param = 20, n.trees = 50, prune = 1 / 15, what = "") {
  stopifnot(!is.null(rownames(data)), metric %in% c("euclidean", "cosine", "manhattan"))
  dir.create(file.path(out, "stages"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "graphs"), recursive = TRUE, showWarnings = FALSE)
  n <- nrow(data)
  f <- ncol(data)
  timing <- list()
  tic <- proc.time()[["elapsed"]]
  # the index FindNeighbors builds (Seurat's own function), its items, and the file save() writes
  a <- Seurat:::AnnoyBuildIndex(data = data, metric = metric, n.trees = n.trees)
  items <- matrix(0, n, f)
  for (i in seq_len(n)) items[i, ] <- a$getItemsVector(i - 1L)
  as_float <- readBin(writeBin(as.vector(data), raw(), size = 4), "double", n = n * f, size = 4)
  dim(as_float) <- dim(data)
  npy_write(file.path(out, "stages", "K3_items.npy"), items)
  ann <- file.path(out, "stages", "annoy.ann")
  a$save(ann)
  node_bytes <- nr_node_bytes(metric, f)
  timing$index_s <- proc.time()[["elapsed"]] - tic
  scal <- list(annoy_metric = metric, annoy_f = f, annoy_n_trees_arg = n.trees, annoy_n_items = a$getNItems(),
               annoy_n_trees_reloaded = a$getNTrees(),
               annoy_items_are_float_of_embeddings = identical(unname(items), unname(as_float)),
               annoy_file_bytes = file.size(ann), annoy_node_bytes = node_bytes,
               annoy_n_nodes = file.size(ann) / node_bytes, n_cells = n, k_param = k.param,
               prune_snn = prune, what = what)
  problems <- character()
  if (!isTRUE(scal$annoy_items_are_float_of_embeddings)) problems <- c(problems, "items differ from the float rounding")
  if (scal$annoy_file_bytes %% node_bytes != 0) problems <- c(problems, "file size is not a multiple of the node size")

  warned <- character()
  fn <- function(...) {
    withCallingHandlers(
      FindNeighbors(data, k.param = k.param, prune.SNN = prune, n.trees = n.trees, annoy.metric = metric,
                    verbose = FALSE, ...),
      warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }
  tic <- proc.time()[["elapsed"]]
  nb <- tryCatch(fn(return.neighbor = TRUE), error = function(e) e)
  timing$findneighbors_neighbours_s <- proc.time()[["elapsed"]] - tic
  if (inherits(nb, "error")) {
    scal$findneighbors_error <- conditionMessage(nb)
  } else {
    idx <- unname(Indices(nb))
    dst <- unname(Distances(nb))
    if (any(idx != round(idx))) stop("non-integer neighbour ids")
    storage.mode(idx) <- "integer"
    k <- ncol(idx)
    npy_write(file.path(out, "stages", "K5_nn_index.npy"), idx)
    npy_write(file.path(out, "stages", "K5_nn_dist.npy"), dst)
    # the same search on a second, unsaved build, and on the saved (reloaded) index
    b <- Seurat:::AnnoyBuildIndex(data = data, metric = metric, n.trees = n.trees)
    s_unsaved <- Seurat:::AnnoySearch(index = b, query = data, k = k)
    s_reloaded <- Seurat:::AnnoySearch(index = a, query = data, k = k)
    scal$knn_k <- k
    scal$knn_search_k <- -1
    scal$search_equals_unsaved_build <- identical(s_unsaved$nn.idx, unname(Indices(nb)) + 0) &&
      identical(s_unsaved$nn.dists, dst)
    scal$search_equals_reloaded_index <- identical(s_reloaded$nn.idx, unname(Indices(nb)) + 0) &&
      identical(s_reloaded$nn.dists, dst)
    scal$knn_self_first <- if (k > 0) sum(idx[, 1] == seq_len(n)) else 0
    scal$knn_self_missing <- sum(!vapply(seq_len(n), function(i) any(idx[i, ] == i), logical(1)))
    if (!isTRUE(scal$search_equals_unsaved_build)) problems <- c(problems, "FindNeighbors differs from a second build's search")
    tic <- proc.time()[["elapsed"]]
    g <- tryCatch(fn(), error = function(e) e)
    timing$findneighbors_graphs_s <- proc.time()[["elapsed"]] - tic
    if (inherits(g, "error")) {
      scal$findneighbors_graphs_error <- conditionMessage(g)
      g <- list()
    }
    facts <- list()
    for (nm in intersect(c("nn", "snn"), names(g))) {
      G <- g[[nm]]
      mtx_write_plain(G, file.path(out, "graphs", paste0(nm, ".mtx.gz")),
                      comments = c(sprintf("metadata_json: {\"software_version\": \"scRust validation/R/neighbors_references.R\", \"format_version\": 2}"),
                                   sprintf("Seurat FindNeighbors graph %s as stored (annoy.metric = %s); rows and columns are the cells.", nm, metric)))
      npy_write(file.path(out, "graphs", paste0(nm, "_p.npy")), G@p)
      npy_write(file.path(out, "graphs", paste0(nm, "_i.npy")), G@i)
      npy_write(file.path(out, "graphs", paste0(nm, "_x.npy")), G@x)
      sorted <- all(vapply(seq_len(ncol(G)), function(j) {
        r <- G@i[(G@p[j] + 1):G@p[j + 1]]
        G@p[j] == G@p[j + 1] || !is.unsorted(r, strictly = TRUE)
      }, logical(1)))
      facts[[nm]] <- list(class = class(G)[1], dim = dim(G), nnz = length(G@x),
                          symmetric = Matrix::isSymmetric(unname(as(G, "dgCMatrix"))),
                          rows_sorted_within_columns = sorted,
                          diagonal_entries = sum(Matrix::diag(G) != 0),
                          distinct_values = length(unique(G@x)))
    }
    json_write(facts, file.path(out, "graphs", "graphs.json"))
  }
  scal$warnings <- I(unique(warned))
  json_write(list(level = "neighbors", contexts = I("chain"), scalars = list(chain = scal), problems = I(problems)),
             file.path(out, "stages", "stages.json"))
  invisible(list(problems = problems, timing = timing, scal = scal))
}

nr_summary <- function(out, kind, extra, res) {
  json_write(c(list(kind = kind), extra,
               list(timing = res$timing, peak_rss_kib = peak_rss_kib(),
                    start_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                    slurm_job_id = Sys.getenv("SLURM_JOB_ID"), fingerprint = fingerprint(packages = c("Seurat", "SeuratObject", "RcppAnnoy", "Matrix", "Rcpp")))),
             file.path(out, "summary.json"))
}

# The constructed inputs: each exercises a branch of Annoy or FindNeighbors that real data reaches
# rarely or not at all. Values are chosen so that the float arithmetic is exact where it matters.
nr_cases <- function() {
  cases <- list()
  add <- function(name, data, what, k.param = 20, n.trees = 50, prune = 1 / 15) {
    data <- as.matrix(data) * 1.0
    dimnames(data) <- list(sprintf("cell%04d", seq_len(nrow(data))), paste0("PC_", seq_len(ncol(data))))
    cases[[name]] <<- list(data = data, what = what, k.param = k.param, n.trees = n.trees, prune = prune)
  }
  set.seed(9091)
  base <- matrix(round(rnorm(40 * 30), 3), 40, 30)
  add("duplicates", rbind(base[rep(1:40, times = 5), ], base[rep(7, 40), ]),
      "240 cells: 40 points 5 times each, interleaved, and 40 more copies of point 7, so 45 identical cells (more than k): ties at distance 0, cells missing from their own lists, splits of identical points")
  add("lattice", as.matrix(expand.grid(-7:7, -7:7)),
      "the integer lattice -7..7 squared (f = 2, leaves of 4): exact midpoints and margins of exactly 0")
  add("two_sites", rbind(matrix(c(1, 0, 0), 70, 3, byrow = TRUE), matrix(c(-1, 0, 0), 70, 3, byrow = TRUE), cbind(0, 1:10, 0)),
      "70 copies of (1, 0, 0), 70 of (-1, 0, 0) and 10 cells on the plane between them: centroids that stay on the sites, margins of exactly 0")
  x <- c(1.5 + 2^-23, -0.25, 2, 0.75)
  add("ulp_outlier", rbind(matrix(x, 199, 4, byrow = TRUE), x + c(2^-23, 0, 0, 0)),
      "199 copies of one point and one cell a float ulp away (in a coordinate with an odd mantissa, so the midpoint rounds onto the outlier): splits that cannot separate them, retried, kept imbalanced or randomised")
  add("collinear", outer((1:150) / 7, c(0.6, -0.8, 0, 0.5, 0.25, -1)),
      "150 cells on a line through the origin (rank one, f = 6)")
  add("zero_rows", rbind(matrix(0, 20, 4), matrix(round(rnorm(80 * 4), 2), 80, 4)),
      "20 zero vectors among 80 others: Angular's two_means skips them and its distance to them is 2")
  for (n in c(1, 2, 3, 5, 32, 33)) {
    add(sprintf("n%d_f30", n), matrix(round(rnorm(n * 30), 3), n, 30),
        sprintf("%d cells, f = 30 (leaves of K = 32): %s", n, if (n <= 32) "every root is a leaf" else "the first split"))
  }
  add("n4_f2", matrix(c(0, 1, 3, 7, 2, 0, -1, 5), 4, 2), "4 cells, f = 2 (K = 4): roots are leaves")
  add("n5_f2", matrix(c(0, 1, 3, 7, 4, 2, 0, -1, 5, 6), 5, 2), "5 cells, f = 2 (K = 4): the smallest split")
  add("fewer_than_k", matrix(rnorm(12 * 30), 12, 30), "12 cells with k.param = 20: FindNeighbors lowers k to 11")
  for (n in c(64, 256, 1024)) add(sprintf("pow2_%d", n), matrix(rnorm(n * 10), n, 10), sprintf("%d cells (a power of two), f = 10", n))
  d300 <- matrix(rnorm(300 * 30), 300, 30)
  add("trees10", d300, "300 cells, n.trees = 10", n.trees = 10)
  add("k5", d300, "300 cells, k.param = 5", k.param = 5)
  add("k50", d300, "300 cells, k.param = 50", k.param = 50)
  add("prune0", d300, "300 cells, prune.SNN = 0 (nothing pruned)", prune = 0)
  add("prune02", d300, "300 cells, prune.SNN = 0.2", prune = 0.2)
  cases
}

nr_files <- c("stages/annoy.ann", "stages/K3_items.npy", "stages/K5_nn_index.npy", "stages/K5_nn_dist.npy",
              "graphs/nn_p.npy", "graphs/nn_i.npy", "graphs/nn_x.npy",
              "graphs/snn_p.npy", "graphs/snn_i.npy", "graphs/snn_x.npy")

if (args$mode == "real") {
  runs <- strsplit(args$runs, ",")[[1]]
  metrics <- strsplit(args$metrics, ",")[[1]]
  if (!length(runs)) stop("--mode real needs --runs")
  for (run in runs) {
    run <- sub("/+$", "", run)
    emb_path <- file.path(run, "pca", "embeddings.npy")
    emb <- npy_read(emb_path)
    meta <- read.delim(file.path(run, "metadata.tsv"), colClasses = "character", quote = "", comment.char = "")
    if (nrow(meta) != nrow(emb)) stop(run, ": metadata.tsv and embeddings.npy differ in cells")
    data <- emb[, 1:30, drop = FALSE]
    dimnames(data) <- list(meta$barcode, paste0("PC_", 1:30))
    for (m in metrics) {
      out <- file.path(dirname(dirname(run)), "neighbors", basename(run), m)
      res <- nr_reference(data, m, out, what = paste("embeddings of", run))
      nr_summary(out, "validation/R/neighbors_references.R --mode real",
                 list(run = run, embeddings = emb_path, embeddings_sha256 = sha256_file(emb_path),
                      metric = m, dims = "1:30", k_param = 20, n_trees = 50, prune = 1 / 15), res)
      if (length(res$problems)) stop(run, " ", m, ": ", paste(res$problems, collapse = "; "))
      cat(sprintf("%s %s: %d cells, %.1f s index, %.1f s FindNeighbors\n", run, m, nrow(data), res$timing$index_s,
                  res$timing$findneighbors_neighbours_s + res$timing$findneighbors_graphs_s))
    }
  }
} else if (args$mode == "edge") {
  if (!nzchar(args$out)) stop("--mode edge needs --out")
  cases <- nr_cases()
  rows <- list()
  exp <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    idir <- file.path(args$out, nm, "input")
    dir.create(idir, recursive = TRUE, showWarnings = FALSE)
    npy_write(file.path(idir, "embeddings.npy"), unname(cs$data))
    writeLines(rownames(cs$data), file.path(idir, "barcodes.tsv"))
    rows[[nm]] <- data.frame(case = nm, n = nrow(cs$data), f = ncol(cs$data), k_param = cs$k.param,
                             n_trees = cs$n.trees, prune = fmt17(cs$prune),
                             embeddings_sha256 = sha256_file(file.path(idir, "embeddings.npy")), what = cs$what)
    for (m in c("euclidean", "cosine", "manhattan")) {
      out <- file.path(args$out, nm, m)
      res <- nr_reference(cs$data, m, out, k.param = cs$k.param, n.trees = cs$n.trees, prune = cs$prune, what = cs$what)
      nr_summary(out, "validation/R/neighbors_references.R --mode edge",
                 list(case = nm, metric = m, n = nrow(cs$data), f = ncol(cs$data), k_param = cs$k.param,
                      n_trees = cs$n.trees, prune = cs$prune), res)
      if (length(res$problems)) stop(nm, " ", m, ": ", paste(res$problems, collapse = "; "))
      for (fl in nr_files) {
        p <- file.path(out, fl)
        if (file.exists(p)) exp[[length(exp) + 1]] <- data.frame(case = nm, metric = m, file = fl, bytes = file.size(p), sha256 = sha256_file(p))
      }
      cat(nm, m, if (is.null(res$scal$findneighbors_error)) "" else paste("(FindNeighbors:", res$scal$findneighbors_error, ")"),
          if (is.null(res$scal$findneighbors_graphs_error)) "" else paste("(graphs:", res$scal$findneighbors_graphs_error, ")"), "\n")
    }
  }
  tsv_write(do.call(rbind, rows), file.path(args$out, "cases.tsv"))
  tsv_write(do.call(rbind, exp), file.path(args$out, "expected.tsv"))
} else if (args$mode == "bench") {
  run <- sub("/+$", "", strsplit(args$runs, ",")[[1]][1])
  emb <- npy_read(file.path(run, "pca", "embeddings.npy"))
  meta <- read.delim(file.path(run, "metadata.tsv"), colClasses = "character", quote = "", comment.char = "")
  data <- emb[, 1:30, drop = FALSE]
  dimnames(data) <- list(meta$barcode, paste0("PC_", 1:30))
  rm(emb)
  invisible(gc())
  before <- peak_rss_kib()
  tic <- proc.time()[["elapsed"]]
  g <- FindNeighbors(data, verbose = FALSE)
  elapsed <- proc.time()[["elapsed"]] - tic
  after <- peak_rss_kib()
  json_write(list(run = run, n_cells = nrow(data), call = "FindNeighbors(data[, 1:30]) with its defaults",
                  elapsed_s = elapsed, vmhwm_before_kib = before, vmhwm_after_kib = after,
                  snn_nnz = length(g$snn@x), nn_nnz = length(g$nn@x),
                  start_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                  fingerprint = fingerprint(packages = c("Seurat", "SeuratObject", "RcppAnnoy", "Matrix", "Rcpp"))),
             args$out)
  cat(sprintf("FindNeighbors: %d cells, %.1f s, VmHWM %.0f -> %.0f KiB\n", nrow(data), elapsed, before, after))
} else {
  stop("--mode must be real, edge or bench")
}
