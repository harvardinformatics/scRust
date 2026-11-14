# Writes the FindClusters part of validation/manifest.toml (docs/PLAN.md A4, C4 and C7; milestone
# M10) from the files on disk: the R references M10 made with validation/R/cluster_references.R,
# that is FindClusters with other settings on each object's saved SCT_snn graph
# (data/clustering/<object>/cluster_refs/<setting>/), the constructed edge cases
# (data/clustering/edge/), R's benchmark runs, and the committed edge fixture
# (tests/fixtures/clustering/edge/). M8's instrumented runs, the workflow-setting references, are in
# the clustering section already. Reads only plain files, and writes one SHA256SUMS file into each
# reference directory (data/, gitignored), which the manifest then lists:
#
#   Rscript validation/R/manifest_cluster.R
#
# The section runs from the line "# FindClusters (plan C4 and C7) " to the line before
# "# Neighbours (plan C4 and C7) ": manifest_neighbors.R regenerates from that line to the
# clustering section, manifest_clustering.R from there to the DoubletFinder section, and
# manifest_doubletfinder.R from there on, so this section sits before all three. Anything else in
# the manifest is kept as is.
#
# The format is compact: the 518 reference directories hold about 10,000 files, so each directory's
# files are listed in one SHA256SUMS (the format of `sha256sum`, so `sha256sum -c SHA256SUMS` in the
# directory checks every file), itself listed as a [[file]] with its own sha256. Each directory under
# cluster_refs/ also gets a [[run]] table with its file count and bytes; the edge cases share one
# [[run]], and each setting's SHA256SUMS [[file]] gives its directory's file count and bytes as
# dir_files and dir_bytes (its own bytes are the SHA256SUMS file's). The edge cases' inputs and the
# committed fixture are listed file by file.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
manifest <- "validation/manifest.toml"

toml_str <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")
toml_val <- function(v) {
  if (is.logical(v)) return(tolower(as.character(v)))
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  toml_str(v)
}
table_lines <- function(kind, fields) {
  fields <- fields[!vapply(fields, function(v) is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v)), logical(1))]
  c(sprintf("[[%s]]", kind),
    vapply(names(fields), function(k) {
      v <- fields[[k]]
      if (length(v) > 1) return(paste0(k, " = [", paste(vapply(v, toml_val, ""), collapse = ", "), "]"))
      paste0(k, " = ", toml_val(v))
    }, ""),
    "")
}
time_v <- function(path) {
  if (length(path) != 1 || !file.exists(path)) return(list())
  t <- readLines(path)
  list(job_max_rss_kib = as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE))),
       job_elapsed = sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE)))
}
call_of <- function(p) {
  sprintf("FindClusters(obj, graph.name = \"SCT_snn\", modularity.fxn = %d, resolution = c(%s), algorithm = %d, n.start = %d, n.iter = %d, random.seed = %d, group.singletons = %s)",
          p$modularity_fxn, paste(vapply(p$resolution, fmt17, ""), collapse = ", "), p$algorithm, p$n_start, p$n_iter,
          p$random_seed, if (isTRUE(p$group_singletons)) "TRUE" else "FALSE")
}

files_out <- character()
runs_out <- character()
n_files <- 0
n_listed <- 0
add_file <- function(fields) {
  files_out <<- c(files_out, table_lines("file", fields))
  n_files <<- n_files + 1
}
# SHA256SUMS of a directory's files (recursively, sorted), written into it; returns its fields
sums_of <- function(dir) {
  fs <- sort(setdiff(list.files(dir, recursive = TRUE), "SHA256SUMS"))
  sums <- vapply(file.path(dir, fs), sha256_file, "")
  writeLines(paste0(sums, "  ", fs), file.path(dir, "SHA256SUMS"))
  n_listed <<- n_listed + length(fs)
  list(files = length(fs), bytes = sum(file.size(file.path(dir, fs))))
}

# FindClusters with other settings, on each object's saved graph ---------------------------------
for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  cr <- file.path("data", "clustering", obj, "cluster_refs")
  if (!dir.exists(cr)) next
  graph <- file.path("data", "clustering", obj, "saved", "graphs", "snn.mtx.gz")
  graph_sha <- sha256_file(graph)
  tv <- list.files(cr, pattern = "^time_v\\.[0-9]+\\.txt$", full.names = TRUE)
  job <- sub("^time_v\\.([0-9]+)\\.txt$", "\\1", basename(tv))
  for (s in sort(list.dirs(cr, recursive = FALSE, full.names = FALSE))) {
    d <- file.path(cr, s)
    sm <- jsonlite::fromJSON(file.path(d, "summary.json"), simplifyVector = FALSE)
    id <- paste0("clustering/", obj, "/cluster_refs/", s)
    sz <- sums_of(d)
    runs_out <- c(runs_out, table_lines("run", c(list(
      id = id, tool = "clustering", kind = "findclusters-reference", stage = "K7", dataset = obj, tier = ds$tier[i],
      path = d, setting = s, call = call_of(sm$parameters), graph = graph, graph_sha256 = graph_sha,
      r_error = sm$run$error, optimizer_equals_R = sm$run$checks$optimizer_equals_R,
      reproduces_saved_object = if (s == "workflow") all(unlist(sm$run$checks[c("labels_equal_saved", "levels_equal_saved", "printed_equal_m8", "messages_equal_m8")])) else NULL,
      recipe = sprintf("validation/slurm/cluster_references.sbatch settings %s (validation/R/cluster_references.R --mode settings)", obj),
      slurm_job_id = if (length(job) == 1) job else NULL, env = "envs/r-reference", r_seconds = sm$run$seconds),
      sz, time_v(tv))))
    add_file(list(path = file.path(d, "SHA256SUMS"), sha256 = sha256_file(file.path(d, "SHA256SUMS")),
                  bytes = file.size(file.path(d, "SHA256SUMS")), role = "reference", stage = "K7", run = id))
  }
  b <- file.path(cr, "bench.json")
  if (file.exists(b)) {
    bj <- jsonlite::fromJSON(b, simplifyVector = FALSE)
    id <- paste0("clustering/", obj, "/cluster_refs/bench")
    runs_out <- c(runs_out, table_lines("run", c(list(
      id = id, tool = "clustering", kind = "findclusters-benchmark", dataset = obj, tier = ds$tier[i], path = b,
      call = "FindClusters(obj, graph.name = \"SCT_snn\", random.seed = 12345)", graph = graph, graph_sha256 = graph_sha,
      findclusters_s = bj$findclusters_s, labels_equal_saved = bj$labels_equal_saved,
      recipe = sprintf("validation/slurm/cluster_references.sbatch bench %s", obj), slurm_job_id = bj$slurm_job_id,
      node = bj$host, env = "envs/r-reference"),
      time_v(file.path(cr, "bench_time_v.txt")))))
    for (f in c(b, file.path(cr, "bench_time_v.txt"))) {
      add_file(list(path = f, sha256 = sha256_file(f), bytes = file.size(f), role = "run-record", run = id))
    }
  }
}

# Constructed edge cases -----------------------------------------------------------------------------
edge <- file.path("data", "clustering", "edge")
if (dir.exists(edge)) {
  cases <- read.delim(file.path(edge, "cases.tsv"), colClasses = "character", quote = "")
  tv <- list.files("data/clustering", pattern = "^edge_time_v\\.[0-9]+\\.txt$", full.names = TRUE)
  job <- sub("^edge_time_v\\.([0-9]+)\\.txt$", "\\1", basename(tv))
  runs_out <- c(runs_out, table_lines("run", c(list(
    id = "clustering/edge", tool = "clustering", kind = "findclusters-edge-cases", stage = "K7", path = edge,
    cases = length(unique(cases$case)), settings = nrow(cases),
    recipe = "validation/slurm/cluster_references.sbatch edge (validation/R/cluster_references.R --mode edge)",
    slurm_job_id = if (length(job) == 1) job else NULL, env = "envs/r-reference",
    note = "cases.tsv lists each case and setting; each <case>/<setting>/ holds a SHA256SUMS of its files"),
    time_v(tv))))
  add_file(list(path = file.path(edge, "cases.tsv"), sha256 = sha256_file(file.path(edge, "cases.tsv")),
                bytes = file.size(file.path(edge, "cases.tsv")), role = "reference", stage = "K7", run = "clustering/edge"))
  for (cs in unique(cases$case)) {
    for (f in sort(list.files(file.path(edge, cs, "input")))) {
      p <- file.path(edge, cs, "input", f)
      add_file(list(path = p, sha256 = sha256_file(p), bytes = file.size(p), role = "input", run = "clustering/edge"))
    }
  }
  for (k in seq_len(nrow(cases))) {
    d <- file.path(edge, cases$case[k], cases$setting[k])
    sz <- sums_of(d)
    add_file(list(path = file.path(d, "SHA256SUMS"), sha256 = sha256_file(file.path(d, "SHA256SUMS")),
                  bytes = file.size(file.path(d, "SHA256SUMS")), role = "reference", stage = "K7", run = "clustering/edge",
                  r_error = if (nzchar(cases$error[k])) cases$error[k] else NULL,
                  dir_files = sz$files, dir_bytes = sz$bytes))
  }
}

# Committed fixture ----------------------------------------------------------------------------------
fx <- file.path("tests", "fixtures", "clustering", "edge")
if (dir.exists(fx)) {
  for (f in sort(list.files(fx, recursive = TRUE))) {
    p <- file.path(fx, f)
    src <- if (f == "cases.tsv") file.path(edge, f) else if (grepl("/", f)) file.path(edge, dirname(f), "input", basename(f)) else NA
    if (!is.na(src) && (!file.exists(src) || !identical(sha256_file(p), sha256_file(src)))) stop("fixture ", p, " is not a copy of ", src)
    add_file(list(path = p, sha256 = sha256_file(p), bytes = file.size(p), role = "fixture",
                  source = if (is.na(src)) "validation/slurm/fixtures_cluster.sh (digests and printed lines of data/clustering/edge/)" else src))
  }
}

header <- c(
  "# FindClusters (plan C4 and C7) --------------------------------------------------------------",
  "#",
  "# Milestone M10: R references for the cluster crate (FindClusters: the modularity optimizer and",
  "# GroupSingletons) beyond M8's instrumented runs, which are the workflow-setting references",
  "# (clustering section, stage K7). Generated by validation/R/manifest_cluster.R from the files on",
  "# disk (do not edit by hand). Results: docs/validation/clustering_m10.md.",
  "#",
  "#   data/clustering/<object>/cluster_refs/<setting>/   FindClusters on the saved SCT_snn graph with",
  "#               the workflow's call and one argument changed (the parameter grid's clustering",
  "#               settings), in the layout of an instrumented run: metadata.tsv, the levels files,",
  "#               cluster_ids.txt, nclusters.txt, stages/ (K7_optimizer_ids, K7_final_ids, rng),",
  "#               summary.json (FindClusters' arguments, what it printed, R's error if it stopped), and",
  "#               optimizer/<res.r>/ (the instrumented optimizer: starts.tsv, partitions.npy, ids.txt)",
  "#   data/clustering/<object>/cluster_refs/bench.json   FindClusters alone, timed (the benchmark)",
  "#   data/clustering/edge/   constructed edge cases (cases.tsv), each with its graph (<case>/input/)",
  "#               and one reference per setting in the same layout",
  "#   tests/fixtures/clustering/edge/   the committed part: cases.tsv, each case's graph, expected.tsv",
  "#               (the bytes and sha256 of every file R wrote) and printed.tsv, for",
  "#               crates/cluster/tests/edge_cases.rs",
  "#",
  "# Compact form: one [[run]] table per reference directory under cluster_refs/, with its file count",
  "# and bytes, and one [[file]] for the SHA256SUMS written into it (`sha256sum -c SHA256SUMS` in the",
  "# directory checks every file). The edge cases share one [[run]], and each setting's SHA256SUMS",
  "# [[file]] gives its directory's file count and bytes as dir_files and dir_bytes. No BLAS is",
  "# involved and no .rds is read: the references depend on the graph only.",
  "# Rust's validation runs (data/clustering/<object>/rust/m10/) are not listed; their results are in",
  "# the report.",
  "")

old <- readLines(manifest)
nb_mark <- grep("^# Neighbours \\(plan C4 and C7\\) ", old)
if (length(nb_mark) != 1) stop("neighbours marker line not found in ", manifest)
fc_mark <- grep("^# FindClusters \\(plan C4 and C7\\) ", old)
before <- if (length(fc_mark) == 1) old[seq_len(fc_mark - 1)] else old[seq_len(nb_mark - 1)]
after <- old[nb_mark:length(old)]
writeLines(c(before, header, runs_out, files_out, after), manifest)
cat("wrote", n_files, "file entries (", n_listed, "files in SHA256SUMS ) and", length(grep("^\\[\\[run\\]\\]", runs_out)), "run entries to", manifest, "\n")
