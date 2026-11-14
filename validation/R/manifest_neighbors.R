# Writes the neighbour part of validation/manifest.toml (docs/PLAN.md A4, C4 and C7; milestone M9)
# from the files on disk: the R references M9 made with validation/R/neighbors_references.R, that is
# FindNeighbors with annoy.metric = "cosine" and "manhattan" on instrumented runs' own embeddings
# (data/clustering/<object>/neighbors/<run>/<metric>/) and the constructed edge cases for all three
# metrics (data/neighbors/edge/), and the committed edge fixture (tests/fixtures/neighbors/edge/).
# The Euclidean references are M8's instrumented runs, already in the clustering section. Reads only
# plain files.
#
#   Rscript validation/R/manifest_neighbors.R
#
# The section runs from the line "# Neighbours (plan C4 and C7) " to the line before
# "# Clustering (plan C4) ": manifest_clustering.R regenerates from that line to the DoubletFinder
# section, and manifest_doubletfinder.R from there on, so this section sits before both. Anything
# else in the manifest is kept as is. The format is the clustering section's compact one: one
# [[run]] table per reference directory, one [[file]] table per file.

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
      if (length(v) > 1) return(paste0(k, " = [", paste(vapply(v, toml_str, ""), collapse = ", "), "]"))
      paste0(k, " = ", toml_val(v))
    }, ""),
    "")
}
stage_of <- function(rel) {
  f <- basename(rel)
  if (f == "annoy.ann") return("K4")
  if (f == "K3_items.npy") return("K3")
  if (startsWith(f, "K5_")) return("K5")
  if (startsWith(rel, "graphs/")) return("K6")
  if (f == "stages.json") return("K3-K6")
  NA
}
role_of <- function(rel) if (basename(rel) == "summary.json") "run-record" else "reference"
fp_fields <- function(fp) {
  if (is.null(fp) || !is.null(fp$error)) return(list())
  list(node = fp$host, cpu_model = fp$cpu_model, glibc = fp$glibc, R = fp$R,
       RcppAnnoy = fp$packages$RcppAnnoy$version, Seurat = fp$packages$Seurat$version)
}
time_v <- function(path) {
  if (!file.exists(path)) return(list())
  t <- readLines(path)
  list(job_max_rss_kib = as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE))),
       job_elapsed = sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE)))
}

files_out <- character()
runs_out <- character()
n_files <- 0
add_dir <- function(dir, run_id, stage_files = TRUE) {
  for (f in sort(list.files(dir, recursive = TRUE))) {
    p <- file.path(dir, f)
    files_out <<- c(files_out, table_lines("file", list(path = p, sha256 = sha256_file(p), bytes = file.size(p),
                                                        role = role_of(f), stage = if (stage_files) stage_of(f) else NA,
                                                        run = run_id)))
    n_files <<- n_files + 1
  }
}

# Angular and Manhattan references on the instrumented runs' embeddings ------------------------------
for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  nb <- file.path("data", "clustering", obj, "neighbors")
  if (!dir.exists(nb)) next
  for (run in sort(list.dirs(nb, recursive = FALSE, full.names = FALSE))) {
    for (metric in sort(list.dirs(file.path(nb, run), recursive = FALSE, full.names = FALSE))) {
      d <- file.path(nb, run, metric)
      sm <- jsonlite::fromJSON(file.path(d, "summary.json"), simplifyVector = FALSE)
      st <- jsonlite::fromJSON(file.path(d, "stages", "stages.json"), simplifyVector = FALSE)$scalars$chain
      id <- paste0("clustering/", obj, "/neighbors/", run, "/", metric)
      runs_out <- c(runs_out, table_lines("run", c(list(
        id = id, tool = "clustering", kind = "neighbors-reference", dataset = obj, tier = ds$tier[i],
        path = d, metric = metric, annoy_index = if (metric == "cosine") "AnnoyAngular" else if (metric == "manhattan") "AnnoyManhattan" else "AnnoyEuclidean",
        embeddings = sm$embeddings, embeddings_sha256 = sm$embeddings_sha256,
        profile = sub("^instrumented-([a-z]+)-.*$", "\\1", run),
        call = "FindNeighbors(embeddings[, 1:30], k.param = 20, n.trees = 50, prune.SNN = 1/15, annoy.metric = metric)",
        recipe = sprintf("validation/slurm/neighbors_references.sbatch real %s/%s (validation/R/neighbors_references.R --mode real)", obj, run),
        slurm_job_id = sm$slurm_job_id, env = "envs/r-reference",
        findneighbors_equals_unsaved_build = isTRUE(st$search_equals_unsaved_build),
        findneighbors_equals_reloaded_index = isTRUE(st$search_equals_reloaded_index)),
        fp_fields(sm$fingerprint),
        time_v(sprintf("data/neighbors/real_time_v.%s.txt", sm$slurm_job_id)))))
      add_dir(d, id)
    }
  }
}

# Constructed edge cases -----------------------------------------------------------------------------
edge <- file.path("data", "neighbors", "edge")
if (dir.exists(edge)) {
  cases <- read.delim(file.path(edge, "cases.tsv"), colClasses = "character", quote = "")
  first <- jsonlite::fromJSON(file.path(edge, cases$case[1], "euclidean", "summary.json"), simplifyVector = FALSE)
  runs_out <- c(runs_out, table_lines("run", c(list(
    id = "neighbors/edge", tool = "clustering", kind = "neighbors-edge-cases", path = edge,
    cases = nrow(cases), metrics = c("euclidean", "cosine", "manhattan"),
    recipe = "validation/slurm/neighbors_references.sbatch edge (validation/R/neighbors_references.R --mode edge)",
    slurm_job_id = first$slurm_job_id, env = "envs/r-reference",
    note = "cases.tsv describes each case; expected.tsv lists the bytes and sha256 of every reference file"),
    fp_fields(first$fingerprint),
    time_v(sprintf("data/neighbors/edge_time_v.%s.txt", first$slurm_job_id)))))
  for (f in c("cases.tsv", "expected.tsv")) {
    p <- file.path(edge, f)
    files_out <- c(files_out, table_lines("file", list(path = p, sha256 = sha256_file(p), bytes = file.size(p),
                                                       role = "reference", stage = "K3-K6", run = "neighbors/edge")))
    n_files <- n_files + 1
  }
  for (cs in cases$case) {
    for (f in sort(list.files(file.path(edge, cs, "input")))) {
      p <- file.path(edge, cs, "input", f)
      files_out <- c(files_out, table_lines("file", list(path = p, sha256 = sha256_file(p), bytes = file.size(p),
                                                         role = "input", run = "neighbors/edge")))
      n_files <- n_files + 1
    }
    for (m in c("euclidean", "cosine", "manhattan")) add_dir(file.path(edge, cs, m), "neighbors/edge")
  }
}

# Committed fixture ----------------------------------------------------------------------------------
fx <- file.path("tests", "fixtures", "neighbors", "edge")
if (dir.exists(fx)) {
  cases <- read.delim(file.path(edge, "cases.tsv"), colClasses = "character", quote = "")
  for (f in sort(list.files(fx, recursive = TRUE))) {
    p <- file.path(fx, f)
    src <- if (f %in% c("cases.tsv", "expected.tsv")) {
      file.path(edge, f)
    } else {
      sha <- sub("^inputs/([0-9a-f]+)\\.npy$", "\\1", f)
      file.path(edge, cases$case[match(sha, cases$embeddings_sha256)], "input", "embeddings.npy")
    }
    if (!file.exists(src) || !identical(sha256_file(p), sha256_file(src))) stop("fixture ", p, " is not a copy of ", src)
    files_out <- c(files_out, table_lines("file", list(path = p, sha256 = sha256_file(p), bytes = file.size(p),
                                                       role = "fixture", source = src)))
    n_files <- n_files + 1
  }
}

header <- c(
  "# Neighbours (plan C4 and C7) ----------------------------------------------------------------",
  "#",
  "# Milestone M9: R references for the neighbors crate (FindNeighbors: Annoy and the SNN graph) beyond",
  "# M8's instrumented runs, which are the Euclidean references (clustering section, stages K3-K6).",
  "# Generated by validation/R/manifest_neighbors.R from the files on disk (do not edit by hand).",
  "# Results: docs/validation/clustering_m9.md.",
  "#",
  "#   data/clustering/<object>/neighbors/<run>/<metric>/   FindNeighbors with annoy.metric = cosine",
  "#               (AnnoyAngular) and manhattan on runs/<run>/pca/embeddings.npy, columns 1:30, in the",
  "#               layout of an instrumented run: stages/ (annoy.ann, K3_items, K5_nn_index, K5_nn_dist,",
  "#               stages.json), graphs/ (nn and snn as Matrix Market and as dgCMatrix slots), summary.json",
  "#   data/neighbors/edge/   constructed edge cases (cases.tsv), each with its input (<case>/input/) and",
  "#               one reference per metric (<case>/<metric>/, the same layout); expected.tsv lists the",
  "#               bytes and sha256 of every reference file",
  "#   tests/fixtures/neighbors/edge/   the committed part: cases.tsv, expected.tsv and every distinct",
  "#               input once (inputs/<sha256>.npy), for crates/neighbors/tests/edge_cases.rs",
  "#",
  "# No BLAS is involved and no .rds is read: the references depend on the embeddings only, not on the",
  "# profile or the CPU type. Rust's validation runs (data/clustering/<object>/rust/m9/,",
  "# data/neighbors/edge-rust/) are not listed; their results are in the report.",
  "")

old <- readLines(manifest)
cl_mark <- grep("^# Clustering \\(plan C4\\) ", old)
if (length(cl_mark) != 1) stop("clustering marker line not found in ", manifest)
nb_mark <- grep("^# Neighbours \\(plan C4 and C7\\) ", old)
before <- if (length(nb_mark) == 1) old[seq_len(nb_mark - 1)] else old[seq_len(cl_mark - 1)]
after <- old[cl_mark:length(old)]
writeLines(c(before, header, runs_out, files_out, after), manifest)
cat("wrote", n_files, "file entries and", length(grep("^\\[\\[run\\]\\]", runs_out)), "run entries to", manifest, "\n")
