# Reference runs of Seurat 5.5.1's FindClusters for milestone M10 (docs/PLAN.md C4, C5 K7 and C7):
# the settings M8's runs did not cover, from the same saved SNN graphs, and constructed edge cases;
# each with the instrumented optimizer's per-start results. No .rds is read and no BLAS is
# involved. Run through validation/slurm/cluster_references.sbatch.
#
#   Rscript validation/R/cluster_references.R --mode settings --object <name> --tier T|P --out <dir> \
#     [--settings all|<name>,<name>,...]
#   Rscript validation/R/cluster_references.R --mode edge --out <dir>
#   Rscript validation/R/cluster_references.R --mode bench --object <name> --out <file.json>
#
# --mode settings: the saved object's SCT_snn graph (data/clustering/<object>/saved/graphs/
#   snn.mtx.gz, cells in saved/metadata.tsv's order) in a minimal Seurat object, and
#   FindClusters(obj, graph.name = "SCT_snn", ...) with each setting of SETTINGS below: the
#   workflow's call (random.seed = 12345, every other default) with one argument changed, as the
#   parameter grid changes them (C4), plus multi-resolution and singleton-label variants. Each
#   setting goes to <out>/<setting>/, the instrumented optimizer's runs to
#   <out>/<setting>/optimizer/<res.r>/.
# --mode edge: the constructed graphs of edge_cases() below, each written to <out>/<case>/input/
#   (snn.mtx.gz, barcodes.tsv), with one reference per setting in <out>/<case>/<setting>/;
#   <out>/cases.tsv lists them.
# --mode bench: the memory and time baseline of C5's performance row for FindClusters alone: the
#   workflow's call on the saved graph, timed, with VmHWM before and after; written to --out as
#   JSON (run under /usr/bin/time -v).
#
# Reading the graph exactly. The saved graphs hold their values as %.17g text, and R's decimal
# parser is not correctly rounded (lib/io.R). Each distinct text is parsed, then corrected to the
# double whose %.17g text it is (one exists, and it is within an ulp or two of R's parse); the
# matrix is rebuilt in stored order and checked to write back the same text. The workflow setting
# must then reproduce the saved object's labels and levels and M8's printed numbers, which each
# setting directory records (summary.json, run.checks).
#
# Each reference directory has the layout of an instrumented clustering run
# (lib/clustering_hooks.R), so `scrust validate cluster` reads either:
#   metadata.tsv                 barcode, each SCT_snn_res.<r> column (in call order),
#                                seurat_clusters
#   <column>_levels.txt          each column's levels; seurat_clusters_levels.txt,
#                                idents_levels.txt; cluster_ids.txt and nclusters.txt
#                                (write_cluster_metadata)
#   stages/K7_optimizer_ids.npy  RunModularityClustering's result for each call (cells x calls,
#                                0-based, numbered by size, before GroupSingletons)
#   stages/K7_final_ids.npy      GroupSingletons' result (cells x calls; -1 for "singleton")
#   stages/rng.npy, rng.tsv      .Random.seed at each GroupSingletons call's entry and exit
#   stages/stages.json           the singletons of each call, in GroupSingletons' order
#   summary.json                 parameters (FindClusters' arguments), result (per call), run
#                                (what FindClusters printed, its messages, an error if R stopped,
#                                timings, checks)
#   optimizer/<res.r>/           the instrumented optimizer (data/clustering/bin/ModularityOptimizer,
#                                validation/cpp/build.sh) on the same graph with the same settings:
#                                starts.tsv (start, iterations, modularity as %.17g),
#                                partitions.npy (starts x cells, each start's own ids), ids.txt
#                                (its result, numbered by size), stdout.txt, optimizer.json
# Where FindClusters stops with an error, summary.json records it (run.error) and nothing else is
# written but the optimizer's runs.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))     # cl_setup, write_cluster_metadata

args <- df_args(list(mode = NULL, object = "", tier = "", out = NULL, settings = "all"))
cl_setup()
suppressPackageStartupMessages(library(Matrix))
OPTIMIZER <- file.path(root, "data", "clustering", "bin", "ModularityOptimizer")
GRAPH <- "SCT_snn"

# The workflow's FindClusters call (C1) ---------------------------------------------------------
BASE <- list(modularity.fxn = 1L, resolution = 0.8, algorithm = 1L, n.start = 10L, n.iter = 10L,
             random.seed = 12345, group.singletons = TRUE)

# One argument changed at a time, as the parameter grid changes them (C4), and variants that reach
# singletons on small objects. modularity.fxn = 2 leaves the resolution unscaled: at the grid's 0.8
# nearly every cell of a test object is a singleton (on some, every one: R then stops), so a second
# setting uses resolutions that give clusters, which depend on the object's size.
settings_for <- function(tier, n_cells) {
  s <- list(
    workflow = list(),
    algorithm2 = list(algorithm = 2L),
    algorithm3 = list(algorithm = 3L),
    modularity2 = list(modularity.fxn = 2L),
    modularity2_resolutions = list(modularity.fxn = 2L,
                                   resolution = if (tier == "T") c(0.05, 0.1, 0.2) else c(0.002, 0.005, 0.01)),
    resolutions = list(resolution = c(0.2, 0.4, 0.8, 1.2, 2.0)),
    nogroup = list(group.singletons = FALSE),
    nogroup_resolutions = list(group.singletons = FALSE, resolution = c(0.8, 2.0)),
    starts1 = list(n.start = 1L),
    starts20 = list(n.start = 20L),
    iter1 = list(n.iter = 1L),
    iter50 = list(n.iter = 50L),
    seed0 = list(random.seed = 0),
    seed1 = list(random.seed = 1)
  )
  if (tier == "P") {
    # modularity.fxn = 2 on large objects: at 0.8 every cell is a singleton or nearly, and R's
    # GroupSingletons costs singletons x clusters x cells (hours); the resolutions that give
    # clusters there are below 1e-3. Kept on the objects of up to 12,000 cells, at the second
    # setting's resolutions.
    s$modularity2 <- NULL
    if (n_cells > 12000) s$modularity2_resolutions <- NULL
  }
  s
}

# Reading the graph exactly -----------------------------------------------------------------------

# The double whose %.17g text each element of `s` is.
exact_doubles <- function(s) {
  x <- as.numeric(s)
  bad <- which(is.na(x) | fmt17(x) != s)
  for (k in bad) {
    if (is.na(x[k]) || !is.finite(x[k]) || x[k] == 0) stop("cannot read ", s[k], " exactly")
    e <- 2^(floor(log2(abs(x[k]))) - 52)
    for (cand in x[k] + c(-e, e, -e / 2, e / 2, -2 * e, 2 * e)) {
      if (fmt17(cand) == s[k]) {
        x[k] <- cand
        break
      }
    }
    if (fmt17(x[k]) != s[k]) stop("cannot read ", s[k], " exactly")
  }
  x
}

# A Matrix Market graph (entries in stored order: column-major, rows increasing) as the identical
# dgCMatrix, its values exact; dimnames are the cells.
graph_read <- function(path, cells) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  if (!grepl("^%%MatrixMarket matrix coordinate (integer|real) general", header)) stop("unsupported header: ", header)
  line <- readLines(con, n = 1)
  while (startsWith(line, "%")) line <- readLines(con, n = 1)
  dims <- as.numeric(strsplit(trimws(line), " +")[[1]])
  nnz <- dims[3]
  tr <- scan(con, what = list(0L, 0L, ""), nmax = nnz, quiet = TRUE)
  i <- tr[[1]]; j <- tr[[2]]; txt <- tr[[3]]
  rm(tr)
  if (length(i) != nnz) stop("expected ", nnz, " entries, read ", length(i))
  if (is.unsorted(j)) stop("entries are not in column-major order")
  if (nnz > 1) {
    same <- j[-1] == j[-nnz]
    if (any(i[-1][same] <= i[-nnz][same])) stop("rows not increasing within a column")
  }
  u <- unique(txt)
  x <- exact_doubles(u)[match(txt, u)]
  if (!identical(fmt17(x), txt)) stop("values do not write back as read")
  p <- c(0L, cumsum(tabulate(j, nbins = dims[2])))
  if (length(cells) != dims[2] || dims[1] != dims[2]) stop("the graph does not fit the cells")
  m <- new("dgCMatrix", i = i - 1L, p = as.integer(p), x = x, Dim = as.integer(dims[1:2]),
           Dimnames = list(cells, cells))
  invisible(validObject(m))
  list(m = m, distinct_values = length(u), corrected = sum(fmt17(as.numeric(u)) != u))
}

# The edge list the instrumented optimizer reads (optimizer_clustering.R's layout): one
# "col<TAB>row<TAB>weight" line per entry with row > col, 0-based, in stored order, weights as
# %.17g; when the last cells have no edge, "n-1<TAB>n-1<TAB>0" sets the node count and adds nothing.
edges_write <- function(m, path) {
  j <- rep.int(seq_len(ncol(m)) - 1L, diff(m@p))
  i <- m@i
  keep <- i > j
  n <- ncol(m)
  lines <- paste(j[keep], i[keep], fmt17(m@x[keep]), sep = "\t")
  mx <- if (any(keep)) max(c(j[keep], i[keep])) else -1L
  padded <- mx < n - 1L && any(keep)
  if (padded) lines <- c(lines, sprintf("%d\t%d\t0", n - 1L, n - 1L))
  writeLines(lines, path)
  list(edges = sum(keep), padded = padded)
}

# Hooks: per call, the optimizer's ids and GroupSingletons' result, singletons and RNG --------------
.k7 <- new.env()
k7_reset <- function() {
  .k7$opt <- list(); .k7$final <- list(); .k7$entry <- list(); .k7$exit <- list(); .k7$singletons <- list()
}
k7_rng <- function() {
  s <- rng_state()
  if (is.null(s)) rep(NA_integer_, 626L) else s
}
k7_rmc_exit <- function(rv) .k7$opt[[length(.k7$opt) + 1L]] <- as.integer(rv)
k7_gs_entry <- function(env) .k7$entry[[length(.k7$entry) + 1L]] <- k7_rng()
k7_gs_exit <- function(env, rv) {
  .k7$final[[length(.k7$final) + 1L]] <- as.character(rv)
  .k7$singletons[[length(.k7$singletons) + 1L]] <- as.character(env$singletons)
  .k7$exit[[length(.k7$exit) + 1L]] <- k7_rng()
}
suppressMessages({
  trace("RunModularityClustering", where = asNamespace("Seurat"), print = FALSE,
        exit = quote(k7_rmc_exit(returnValue())))
  trace("GroupSingletons", where = asNamespace("Seurat"), print = FALSE,
        tracer = quote(k7_gs_entry(environment())), exit = quote(k7_gs_exit(environment(), returnValue())))
})

# A minimal Seurat object holding the graph as SCT_snn (FindClusters.Seurat reads nothing else).
make_object <- function(g) {
  cells <- colnames(g)
  counts <- sparseMatrix(i = rep(1L, length(cells)), j = seq_along(cells), x = rep(1, length(cells)),
                         dims = c(1L, length(cells)), dimnames = list("feature1", cells))
  obj <- suppressWarnings(CreateSeuratObject(counts = counts))
  gr <- as.Graph(g)
  DefaultAssay(gr) <- "RNA"
  obj[[GRAPH]] <- gr
  if (!identical(obj[[GRAPH]]@x, g@x) || !identical(obj[[GRAPH]]@i, g@i) || !identical(obj[[GRAPH]]@p, g@p)) {
    stop("the Graph differs from the matrix")
  }
  obj
}

# The numbers FindClusters printed, per optimizer run.
printed_blocks <- function(printed) {
  starts <- grep("^Modularity Optimizer version", printed)
  lapply(seq_along(starts), function(k) printed[starts[k]:(if (k < length(starts)) starts[k + 1] - 1L else length(printed))])
}
block_numbers <- function(b) {
  get <- function(re) {
    h <- grep(re, b, value = TRUE)
    if (length(h)) sub(re, "\\1", h[1]) else NA_character_
  }
  list(nodes = as.integer(get("^Number of nodes: ([0-9]+)$")),
       edges = as.integer(get("^Number of edges: ([0-9]+)$")),
       max_modularity = get("^(?:Maximum modularity in [0-9]+ random starts|Modularity): (-?[0-9.]+)$"),
       communities = as.integer(get("^Number of communities: ([0-9]+)$")))
}

# The instrumented optimizer for one call: modularity function, resolution (as %.17g), algorithm,
# starts, iterations, seed (FindClusters' int; the program takes it as unsigned long long, which a
# non-negative seed fits).
run_optimizer <- function(edges, n, s, resolution, dir, r_ids) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ids_path <- file.path(dir, "ids.txt")
  # a negative seed as text: std::stoull negates it modulo 2^64, as R's int becomes a uint64_t
  seed <- as.integer(s$random.seed)
  t1 <- proc.time()[["elapsed"]]
  out <- system2(OPTIMIZER, c(shQuote(edges), shQuote(ids_path), s$modularity.fxn, fmt17(resolution),
                              s$algorithm, s$n.start, s$n.iter, seed, "1"), stdout = TRUE, stderr = TRUE)
  secs <- proc.time()[["elapsed"]] - t1
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop("the optimizer failed")
  starts <- grep("^SCRUST_START ", out, value = TRUE)
  writeLines(out[!startsWith(out, "SCRUST_START ")], file.path(dir, "stdout.txt"))
  parts <- strsplit(starts, " ", fixed = TRUE)
  writeLines(c("start\titerations\tmodularity",
               vapply(parts, function(p) paste(p[2], p[3], p[4], sep = "\t"), "")),
             file.path(dir, "starts.tsv"))
  P <- do.call(rbind, lapply(parts, function(p) as.integer(p[-(1:4)])))
  if (ncol(P) != n) stop("partitions have ", ncol(P), " nodes, not ", n)
  npy_write(file.path(dir, "partitions.npy"), P, "<i4")
  ids <- as.integer(readLines(ids_path))
  json_write(list(starts = length(parts), seconds = secs,
                  ids_equal_R = if (is.null(r_ids)) NULL else identical(ids, r_ids),
                  printed = I(out[!startsWith(out, "SCRUST_START ")][grepl("^(Number|Maximum|Modularity:)", out[!startsWith(out, "SCRUST_START ")])])),
             file.path(dir, "optimizer.json"))
  invisible(identical(ids, r_ids))
}

# One FindClusters call with setting `s` on graph `g`, written to `dir`. `check(obj, printed, msgs)`
# may return checks for summary.json (run.checks); it is called when FindClusters succeeded.
run_setting <- function(g, s, name, dir, edges, extra = list(), check = NULL) {
  if (file.exists(file.path(dir, "summary.json"))) stop(dir, " exists")
  dir.create(file.path(dir, "stages"), recursive = TRUE, showWarnings = FALSE)
  s <- modifyList(BASE, s)
  obj <- make_object(g)
  cells <- colnames(g)
  k7_reset()
  set.seed(20261114)  # a defined stream before the calls, so "no draw" is visible
  msgs <- character()
  err <- NULL
  cl <- as.call(c(list(as.name("FindClusters"), object = as.name("obj"), graph.name = GRAPH, verbose = TRUE), s))
  t0 <- proc.time()[["elapsed"]]
  printed <- utils::capture.output(
    tryCatch(withCallingHandlers(obj2 <- eval(cl),
                                 message = function(m) {
                                   msgs <<- c(msgs, conditionMessage(m))
                                   invokeRestart("muffleMessage")
                                 }),
             error = function(e) err <<- conditionMessage(e)))
  secs <- proc.time()[["elapsed"]] - t0
  res_names <- paste0("res.", s$resolution)
  columns <- paste(GRAPH, res_names, sep = "_")
  params <- list(graph_name = GRAPH, modularity_fxn = as.integer(s$modularity.fxn), resolution = I(as.double(s$resolution)),
                 algorithm = as.integer(s$algorithm), n_start = as.integer(s$n.start), n_iter = as.integer(s$n.iter),
                 random_seed = as.integer(s$random.seed), group_singletons = isTRUE(s$group.singletons))
  blocks <- printed_blocks(printed)
  per_call <- lapply(seq_along(.k7$opt), function(k) {
    c(list(name = res_names[k], column = columns[k], communities_optimizer = length(unique(.k7$opt[[k]])),
           singletons = I(if (k <= length(.k7$singletons)) .k7$singletons[[k]] else character()),
           drew = if (k <= length(.k7$exit)) !identical(.k7$entry[[k]], .k7$exit[[k]]) else NA),
      if (k <= length(blocks)) block_numbers(blocks[[k]]) else list())
  })
  checks <- list()
  result <- list(calls = per_call)
  if (is.null(err)) {
    md <- obj2@meta.data
    if (!identical(rownames(md), cells)) stop("meta.data rows differ from the cells")
    if (!all(columns %in% colnames(md))) stop("columns missing: ", paste(setdiff(columns, colnames(md)), collapse = ", "))
    meta <- data.frame(barcode = cells, stringsAsFactors = FALSE)
    for (col in columns) meta[[col]] <- as.character(md[[col]])
    meta$seurat_clusters <- as.character(md$seurat_clusters)
    tsv_write(meta, file.path(dir, "metadata.tsv"))
    for (col in c(columns, "seurat_clusters")) writeLines(levels(md[[col]]), file.path(dir, paste0(col, "_levels.txt")))
    writeLines(levels(Idents(obj2)), file.path(dir, "idents_levels.txt"))
    write_cluster_metadata(obj2, file.path(dir, "nclusters.txt"), file.path(dir, "cluster_ids.txt"))
    O <- do.call(cbind, .k7$opt)
    npy_write(file.path(dir, "stages", "K7_optimizer_ids.npy"), O, "<i4")
    Fn <- do.call(cbind, lapply(.k7$final, function(v) ifelse(v == "singleton", -1L, suppressWarnings(as.integer(v)))))
    npy_write(file.path(dir, "stages", "K7_final_ids.npy"), Fn, "<i4")
    R <- do.call(rbind, unlist(lapply(seq_along(.k7$entry), function(k) list(.k7$entry[[k]], .k7$exit[[k]])), recursive = FALSE))
    npy_write(file.path(dir, "stages", "rng.npy"), R, "<i4")
    tsv_write(data.frame(row = seq_len(nrow(R)) - 1L, ctx = "chain", stage = "K7",
                         label = rep(c("groupsingletons_entry", "groupsingletons_exit"), length(.k7$entry))),
              file.path(dir, "stages", "rng.tsv"))
    sc <- list(n_calls = length(.k7$opt), singletons_by_call = lapply(.k7$singletons, I))
    if (length(.k7$singletons) == 1) sc$singletons <- I(.k7$singletons[[1]])
    json_write(list(level = "k7", contexts = "chain", scalars = list(chain = sc)), file.path(dir, "stages", "stages.json"))
    result$columns <- I(columns)
    result$seurat_clusters_levels <- I(levels(md$seurat_clusters))
    if (!is.null(check)) checks <- check(obj2, printed, msgs)
  } else {
    unlink(file.path(dir, "stages"), recursive = TRUE)
  }
  # the instrumented optimizer, for every call R's optimizer finished
  opt_ok <- NA
  if (!is.null(edges) && length(.k7$opt)) {
    opt_ok <- all(vapply(seq_along(.k7$opt), function(k)
      run_optimizer(edges, length(cells), s, s$resolution[k], file.path(dir, "optimizer", res_names[k]), .k7$opt[[k]]), NA))
  }
  checks$optimizer_equals_R <- opt_ok
  run <- c(list(setting = name, printed = I(printed), messages = I(msgs), error = err, seconds = secs,
                peak_rss_kib = peak_rss_kib(), checks = checks), extra)
  json_write(list(parameters = params, result = result, run = run), file.path(dir, "summary.json"))
  invisible(list(error = err, optimizer_equal = opt_ok, checks = checks))
}

# --mode settings -----------------------------------------------------------------------------------
if (args$mode == "settings") {
  obj_dir <- file.path(root, "data", "clustering", args$object)
  saved_meta <- read.delim(file.path(obj_dir, "saved", "metadata.tsv"), colClasses = "character", quote = "")
  cells <- saved_meta$barcode
  rd <- graph_read(file.path(obj_dir, "saved", "graphs", "snn.mtx.gz"), cells)
  g <- rd$m
  dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
  edges <- file.path(args$out, "edges.tsv")
  ew <- edges_write(g, edges)
  sets <- settings_for(args$tier, length(cells))
  wanted <- if (args$settings == "all") names(sets) else strsplit(args$settings, ",")[[1]]
  if (!all(wanted %in% names(sets))) stop("unknown settings: ", paste(setdiff(wanted, names(sets)), collapse = ", "))
  # M8's printed numbers for this object (the production instrumented run)
  m8 <- jsonlite::fromJSON(file.path(obj_dir, "runs", "instrumented-production-sapphirerapids", "summary.json"), simplifyVector = FALSE)
  # is SNN[i.cells, j.cells] a plain vector, so that GroupSingletons takes mean()?
  sub <- g[cells[1], cells[1:min(5, length(cells))]]
  facts <- list(object = args$object, tier = args$tier, graph = file.path("data", "clustering", args$object, "saved", "graphs", "snn.mtx.gz"),
                n_cells = length(cells), graph_nnz = length(g@x), graph_distinct_values = rd$distinct_values,
                graph_values_corrected = rd$corrected, edges = ew$edges, edges_padded = ew$padded,
                subsnn_is_object = is.object(sub))
  wf_check <- function(obj2, printed, msgs) {
    md <- obj2@meta.data
    keep <- function(p) p[!grepl("^Elapsed time", p)]
    list(labels_equal_saved = identical(as.character(md$SCT_snn_res.0.8), saved_meta$SCT_snn_res.0.8) &&
           identical(as.character(md$seurat_clusters), saved_meta$seurat_clusters),
         levels_equal_saved = identical(levels(md$seurat_clusters), readLines(file.path(obj_dir, "saved", "seurat_clusters_levels.txt"))) &&
           identical(levels(md$SCT_snn_res.0.8), readLines(file.path(obj_dir, "saved", "SCT_snn_res.0.8_levels.txt"))),
         printed_equal_m8 = identical(keep(printed), keep(unlist(m8$run$printed))),
         messages_equal_m8 = identical(trimws(msgs), trimws(unlist(m8$run$messages))))
  }
  for (nm in wanted) {
    cat(sprintf("%s %s ...\n", args$object, nm))
    r <- run_setting(g, sets[[nm]], nm, file.path(args$out, nm), edges, list(facts = facts),
                     check = if (nm == "workflow") wf_check else NULL)
    if (nm == "workflow") {
      bad <- names(r$checks)[!vapply(r$checks, isTRUE, NA)]
      if (!is.null(r$error) || length(bad)) stop("the workflow setting does not reproduce the saved object: ", paste(c(r$error, bad), collapse = ", "))
    }
  }
  unlink(edges)
  cat("done\n")
}

# --mode edge -----------------------------------------------------------------------------------
# A symmetric graph with unit diagonal (as SNN graphs have) from (i, j, w) pairs, 1-based, i != j.
sym_graph <- function(n, i, j, w, diag = TRUE) {
  ii <- c(i, j, if (diag) seq_len(n))
  jj <- c(j, i, if (diag) seq_len(n))
  xx <- c(w, w, if (diag) rep(1, n))
  m <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  m <- as(m, "generalMatrix")
  m <- as(m, "CsparseMatrix")
  m <- as(m, "dMatrix")
  dimnames(m) <- list(sprintf("cell%04d", seq_len(n)), sprintf("cell%04d", seq_len(n)))
  m
}
# A clique of the given cells, weights from a fixed cycle (so that sums are not all equal).
clique <- function(cells, ws = c(0.9, 0.6, 0.75, 0.45, 0.8, 0.55)) {
  pr <- t(utils::combn(cells, 2))
  list(i = pr[, 1], j = pr[, 2], w = rep_len(ws, nrow(pr)))
}
bind_edges <- function(...) {
  l <- list(...)
  list(i = unlist(lapply(l, `[[`, "i")), j = unlist(lapply(l, `[[`, "j")), w = unlist(lapply(l, `[[`, "w")))
}
# An SNN graph as FindNeighbors builds it (k = 20, Annoy, 1/15) from a Gaussian mixture.
snn_of <- function(n, dims = 10, centers = 4, seed = 1, k = 20, sd = 1) {
  set.seed(seed)
  mu <- matrix(rnorm(centers * dims, sd = 4), centers, dims)
  lab <- rep_len(seq_len(centers), n)
  x <- mu[lab, , drop = FALSE] + matrix(rnorm(n * dims, sd = sd), n, dims)
  rownames(x) <- sprintf("cell%04d", seq_len(n))
  gr <- suppressWarnings(FindNeighbors(x, k.param = min(k, n - 1), verbose = FALSE))
  m <- as(gr$snn, "dgCMatrix")
  m
}

edge_cases <- function() {
  cases <- list()
  add <- function(name, g, settings, what) cases[[name]] <<- list(g = g, settings = settings, what = what)
  wf <- list(workflow = list())
  # isolated cells: no SNN edge, so the optimizer leaves them alone and GroupSingletons finds them
  # tied at 0 with every cluster; R's draw places them
  e <- bind_edges(clique(c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)), clique(c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)))
  add("isolated_two_clusters", sym_graph(30, e$i, e$j, e$w),
      c(wf, list(nogroup = list(group.singletons = FALSE), seed1 = list(random.seed = 1))),
      "two weighted cliques of 12 interleaved, 6 isolated cells at the end: ties of 2")
  e <- do.call(bind_edges, lapply(0:9, function(c) clique(c * 6 + 1:6)))
  add("isolated_ten_clusters", sym_graph(65, e$i, e$j, e$w),
      c(wf, list(nogroup = list(group.singletons = FALSE))),
      "ten cliques of 6 and 5 isolated cells: ties of 10, which sample() resolves to the 9th")
  e <- do.call(bind_edges, lapply(0:11, function(c) clique(c * 5 + 1:5)))
  add("isolated_twelve_clusters", sym_graph(64, c(e$i), c(e$j), c(e$w)),
      wf, "twelve cliques of 5, 4 isolated cells; 64 nodes, a power of two")
  # singletons whose order matters: with modularity function 2 at 0.8 a node joins a cluster only
  # for an edge weight above 0.8 per member, so cells 5 and 6 stay alone; cell 5 joins A, and then
  # cell 6, tied to A only through cell 5, joins A too (without cell 5 there, it would join B)
  add("singleton_order", sym_graph(6, c(1, 3, 5, 5, 6), c(2, 4, 1, 6, 3), c(1, 1, 0.5, 0.75, 0.4)),
      list(modularity2 = list(modularity.fxn = 2L), modularity2_nogroup = list(modularity.fxn = 2L, group.singletons = FALSE)),
      "pairs A = {1, 2}, B = {3, 4}; singletons 5 (to A 0.5, to 6 0.75) and 6 (to B 0.4)")
  add("singleton_order_reversed", sym_graph(6, c(1, 3, 6, 6, 5), c(2, 4, 1, 5, 3), c(1, 1, 0.5, 0.75, 0.4)),
      list(modularity2 = list(modularity.fxn = 2L)),
      "as singleton_order with cells 5 and 6 swapped: the singleton tied to B through the other comes first")
  # networks of one and two cells, and disconnected components
  add("one_cell", sym_graph(1, integer(), integer(), numeric()), wf, "a single cell: no network (R stops)")
  add("two_cells", sym_graph(2, 1, 2, 0.5), wf, "two connected cells")
  add("two_cells_apart", sym_graph(2, integer(), integer(), numeric()), wf, "two cells without an edge: no network (R stops)")
  e <- bind_edges(clique(1:7), clique(8:12), clique(13:17), list(i = 19, j = 20, w = 0.3))
  add("components", sym_graph(20, e$i, e$j, e$w), c(wf, list(algorithm2 = list(algorithm = 2L), algorithm3 = list(algorithm = 3L))),
      "cliques of 7, 5 and 5 (equal sizes keep their order), an isolated cell, a pair")
  # equal gains and equal-sized clusters
  add("ring8", sym_graph(8, 1:8, c(2:8, 1), rep(0.5, 8)), c(wf, list(algorithm3 = list(algorithm = 3L), resolution2 = list(resolution = 2))),
      "a ring of 8 equal edges: every gain ties")
  lat <- expand.grid(r = 1:6, c = 1:6)
  idx <- function(r, c) (c - 1) * 6 + r
  hi <- c(idx(lat$r[lat$r < 6], lat$c[lat$r < 6]), idx(lat$r[lat$c < 6], lat$c[lat$c < 6]))
  hj <- c(idx(lat$r[lat$r < 6] + 1, lat$c[lat$r < 6]), idx(lat$r[lat$c < 6], lat$c[lat$c < 6] + 1))
  add("lattice6", sym_graph(36, hi, hj, rep(0.25, length(hi))), c(wf, list(starts1 = list(n.start = 1L))),
      "a 6 x 6 lattice of equal edges")
  # everything a singleton: R's GroupSingletons stops (sample() of nothing); without grouping it does not
  add("all_singletons", sym_graph(10, 1:9, 2:10, rep(0.3, 9)),
      list(modularity2 = list(modularity.fxn = 2L), modularity2_nogroup = list(modularity.fxn = 2L, group.singletons = FALSE)),
      "a path of weak edges: with modularity function 2 at 0.8 no node joins another")
  # SNN graphs as FindNeighbors builds them: powers of two (nextInt's other branch) and their neighbours
  for (n in c(63, 64, 65, 128, 256, 1024)) {
    add(sprintf("snn%d", n), snn_of(n, seed = n),
        c(wf, list(resolutions = list(resolution = c(0.3, 1, 3)), algorithm2 = list(algorithm = 2L), algorithm3 = list(algorithm = 3L))),
        sprintf("an SNN graph of %d cells from a 4-component Gaussian mixture (FindNeighbors, k = 20)", n))
  }
  # many singletons with edges, merged in turn (high resolution), labels sorted as text without grouping
  g300 <- snn_of(300, seed = 300, centers = 6)
  add("snn300_high_resolution", g300,
      list(resolution5 = list(resolution = 5), nogroup2 = list(group.singletons = FALSE, resolution = 2),
           order = list(resolution = c(1.2, 0.4, 2.0)), modularity2 = list(modularity.fxn = 2L, resolution = c(0.05, 0.2, 0.5))),
      "300 cells, 6 components: many singletons with edges")
  add("seeds", snn_of(200, seed = 7),
      list(seed_minus1 = list(random.seed = -1), seed_max = list(random.seed = 2147483647), seed_minus12345 = list(random.seed = -12345)),
      "negative and extreme seeds (the int converts to uint64_t by sign extension)")
  cases
}

if (args$mode == "edge") {
  if (file.exists(args$out)) stop(args$out, " exists")
  dir.create(args$out, recursive = TRUE)
  cases <- edge_cases()
  rows <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    g <- cs$g
    idir <- file.path(args$out, nm, "input")
    dir.create(idir, recursive = TRUE)
    # the graph as scrust reads it: snn.mtx.gz in stored order, values %.17g, and the cells
    con <- gzfile(file.path(idir, "snn.mtx.gz"), "w")
    jj <- rep.int(seq_len(ncol(g)), diff(g@p))
    writeLines(c("%%MatrixMarket matrix coordinate real general",
                 sprintf("%% scRust validation/R/cluster_references.R edge case %s", nm),
                 sprintf("%d %d %d", nrow(g), ncol(g), length(g@x)),
                 if (length(g@x)) paste(g@i + 1L, jj, fmt17(g@x))), con)
    close(con)
    writeLines(colnames(g), file.path(idir, "barcodes.tsv"))
    back <- graph_read(file.path(idir, "snn.mtx.gz"), colnames(g))$m
    if (!identical(back@x, g@x) || !identical(back@i, g@i) || !identical(back@p, g@p)) stop(nm, ": the written graph does not read back")
    edges <- file.path(args$out, nm, "edges.tsv")
    ew <- edges_write(g, edges)
    if (ew$edges == 0) edges <- NULL
    for (sn in names(cs$settings)) {
      cat(sprintf("%s %s ...\n", nm, sn))
      r <- run_setting(g, cs$settings[[sn]], sn, file.path(args$out, nm, sn), edges,
                       list(facts = list(case = nm, what = cs$what, n_cells = ncol(g), edges = ew$edges)))
      s <- modifyList(BASE, cs$settings[[sn]])
      rows[[length(rows) + 1]] <- data.frame(case = nm, setting = sn, n_cells = ncol(g), edges = ew$edges,
        modularity_fxn = s$modularity.fxn, resolution = paste(fmt17(s$resolution), collapse = ","),
        algorithm = s$algorithm, n_start = s$n.start, n_iter = s$n.iter, random_seed = as.integer(s$random.seed),
        group_singletons = isTRUE(s$group.singletons), error = if (is.null(r$error)) "" else r$error,
        optimizer_equals_R = if (is.na(r$optimizer_equal)) "" else as.character(r$optimizer_equal),
        what = cs$what, stringsAsFactors = FALSE)
    }
    if (!is.null(edges)) unlink(edges)
  }
  tsv_write(do.call(rbind, rows), file.path(args$out, "cases.tsv"))
  cat("done\n")
}

# --mode bench --------------------------------------------------------------------------------
if (args$mode == "bench") {
  obj_dir <- file.path(root, "data", "clustering", args$object)
  saved_meta <- read.delim(file.path(obj_dir, "saved", "metadata.tsv"), colClasses = "character", quote = "")
  t0 <- proc.time()[["elapsed"]]
  g <- graph_read(file.path(obj_dir, "saved", "graphs", "snn.mtx.gz"), saved_meta$barcode)$m
  obj <- make_object(g)
  rm(g)
  invisible(gc())
  load_s <- proc.time()[["elapsed"]] - t0
  hwm0 <- peak_rss_kib()
  rss0 <- as.numeric(sub("^VmRSS:\\s+([0-9]+) kB$", "\\1", grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)))
  t1 <- proc.time()[["elapsed"]]
  obj <- FindClusters(obj, graph.name = GRAPH, random.seed = 12345)
  fc_s <- proc.time()[["elapsed"]] - t1
  hwm1 <- peak_rss_kib()
  ok <- identical(as.character(obj$seurat_clusters), saved_meta$seurat_clusters)
  json_write(list(object = args$object, n_cells = ncol(obj), load_s = load_s, findclusters_s = fc_s,
                  vmhwm_before_kib = hwm0, vmrss_before_kib = rss0, vmhwm_after_kib = hwm1,
                  labels_equal_saved = ok, slurm_job_id = Sys.getenv("SLURM_JOB_ID"), host = Sys.info()[["nodename"]]),
             args$out)
  cat("done\n")
}
