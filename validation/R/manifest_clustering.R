# Writes the clustering part of validation/manifest.toml (docs/PLAN.md A4, A5b and C4; milestone M8)
# from the files on disk: every export, reference run, marker check, instrumented-optimizer run and
# comparison under data/clustering/<object>/ for the objects of validation/clustering_datasets.tsv,
# the K0 vectors (data/clustering/k0/), the C++ build record (data/clustering/bin/build.json) and the
# committed fixtures under tests/fixtures/clustering/. Reads only plain files.
#
#   Rscript validation/R/manifest_clustering.R
#
# The section runs from the line "# Clustering (plan C4) " to the line before
# "# DoubletFinder (plan B4) " (manifest_doubletfinder.R regenerates everything from that line on,
# so this section sits before it). Anything else in the manifest is kept as is.
#
# Format (a compact form of the manifest's [[file]] schema, documented in the section's header):
# one [[run]] table per run, export, marker check, optimizer run or comparison set, holding what
# its files share (recipe, profile, node, CPU, OpenBLAS core, job, peak RSS, wall time), and one
# [[file]] table per distinct file: path, sha256, bytes, role, stage, the run it belongs to, and
# identical_on, the other runs of the same object whose file at the same relative path has the
# same sha256 (those copies are not listed again).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
manifest <- "validation/manifest.toml"
PIPELINE <- "/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing"

toml_str <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")
toml_val <- function(v) {
  if (is.logical(v)) return(tolower(as.character(v)))
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  toml_str(v)
}
table_lines <- function(kind, fields, arrays = character()) {
  fields <- fields[!vapply(fields, function(v) is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v)), logical(1))]
  c(sprintf("[[%s]]", kind),
    vapply(names(fields), function(k) {
      v <- fields[[k]]
      if (k %in% arrays || length(v) > 1) return(paste0(k, " = [", paste(vapply(v, toml_str, ""), collapse = ", "), "]"))
      paste0(k, " = ", toml_val(v))
    }, ""),
    "")
}
hash_cache <- new.env()
hash_of <- function(p) {
  key <- normalizePath(p)
  if (!is.null(hash_cache[[key]])) return(hash_cache[[key]])
  h <- sha256_file(p)
  assign(key, h, envir = hash_cache)
  h
}
parse_time_v <- function(path) {
  if (!file.exists(path)) return(list())
  t <- readLines(path)
  list(max_rss_kib = as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE))),
       elapsed = sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE)))
}
stage_of <- function(rel) {
  f <- basename(rel)
  if (grepl("^stages/vst/", rel)) return("K1 (S2a-S2d: vst's result)")
  if (grepl("^stages/", rel)) {
    if (f == "annoy.ann") return("K4")
    m <- regmatches(f, regexpr("^(S2[a-d]|S3|K[0-9]+)", f))
    if (length(m)) return(m)
    return("")
  }
  if (grepl("^sct/", rel)) return(if (grepl("corrected_counts", rel)) "K10" else "K1")
  if (grepl("^pca/", rel)) return("K2")
  if (grepl("^graphs/", rel)) return("K6")
  if (rel %in% c("metadata.tsv", "seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt", "idents_levels.txt",
                 "cluster_ids.txt", "nclusters.txt")) return("K7, K8, K9")
  ""
}
role_of <- function(kind, rel) {
  if (rel %in% c("summary.json", "time_v.txt", "fingerprint_wrapper.json", "export_run.json") ||
      grepl("time_v\\.txt$|stdout\\.txt$", rel)) return("run-record")
  switch(kind, export = if (startsWith(rel, "input/")) "input" else "saved-output",
         markers = "reference", optimizer = "reference", compare = "check", "reference")
}
fp_fields <- function(fp) {
  if (is.null(fp) || !is.null(fp$error)) return(list())
  list(node = fp$host, cpu_model = fp$cpu_model, openblas_core = fp$openblas_corename,
       openblas_threads = fp$openblas_num_threads, glibc = fp$glibc, R = fp$R)
}

RUN_ORDER <- c("instrumented-portable-sapphirerapids", "direct-portable-sapphirerapids",
               "instrumented-production-sapphirerapids", "direct-production-sapphirerapids",
               "direct-production-sapphirerapids-noumap", "instrumented-production-sapphirerapids-v5",
               "instrumented-portable-cascadelake", "direct-portable-cascadelake")

files_out <- character()
runs_out <- character()
n_files <- 0
n_copies <- 0

for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  tier <- ds$tier[i]
  base <- file.path("data", "clustering", obj)
  if (!dir.exists(file.path(base, "saved"))) next
  # Every unit of this object: export, runs, markers, optimizer, comparisons
  units <- list()
  ex <- jsonlite::fromJSON(file.path(base, "input", "export.json"), simplifyVector = FALSE)
  er <- jsonlite::fromJSON(file.path(base, "export_run.json"), simplifyVector = FALSE)
  units[[length(units) + 1]] <- list(id = paste0("clustering/", obj, "/export"), kind = "export", dir = base,
    files = c("export_run.json", file.path("input", c("assay.txt", "export.json")),
              if (is.null(ex$counts_from)) file.path("input", c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")),
              file.path("saved", list.files(file.path(base, "saved"), recursive = TRUE))),
    fields = c(list(tool = "clustering", kind = "export", dataset = obj, tier = tier,
                    recipe = sprintf("validation/slurm/export_clustering.sbatch %s (validation/R/export_clustering.R --level %s)", obj, ds$level[i]),
                    source = ex$source, source_sha256 = ex$source_sha256,
                    counts_from = ex$counts_from, slurm_job_id = er$slurm_job_id,
                    peak_rss_kib = er$peak_rss_kib, env = "envs/r-reference",
                    note = if (!is.null(ex$counts_from)) "input/matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz are relative links to the counts_from export, listed in the DoubletFinder section"),
               fp_fields(er$fingerprint)))
  run_dirs <- list.dirs(file.path(base, "runs"), recursive = FALSE, full.names = FALSE)
  run_dirs <- run_dirs[!grepl("\\.partial\\.", run_dirs)]
  run_dirs <- c(intersect(RUN_ORDER, run_dirs), setdiff(sort(run_dirs), RUN_ORDER))
  for (r in run_dirs) {
    rd <- file.path(base, "runs", r)
    s <- jsonlite::fromJSON(file.path(rd, "summary.json"), simplifyVector = FALSE)
    tv <- parse_time_v(file.path(rd, "time_v.txt"))
    mode <- sub("-.*", "", r)
    units[[length(units) + 1]] <- list(id = paste0("clustering/", obj, "/", r), kind = mode, dir = rd,
      files = list.files(rd, recursive = TRUE),
      fields = c(list(tool = "clustering", kind = mode, dataset = obj, tier = tier, path = rd,
                      profile = s$run$profile, umap = s$run$umap, assay = s$run$assay, level = s$run$level,
                      recipe = sprintf("%svalidation/slurm/reference_clustering.sbatch %s %s %s",
                                       paste0(if (!isTRUE(s$run$umap)) "CL_UMAP=no " else "",
                                              if (grepl("-v[35]$", r)) paste0("CL_ASSAY=", sub(".*-", "", r), " ") else ""),
                                       s$run$profile, mode, obj),
                      slurm_job_id = s$run$slurm_job_id, env = "envs/r-reference"),
                 fp_fields(s$run$fingerprint),
                 list(max_rss_kib = tv$max_rss_kib, elapsed = tv$elapsed)))
  }
  if (file.exists(file.path(base, "markers", "markers.json"))) {
    mj <- jsonlite::fromJSON(file.path(base, "markers", "markers.json"), simplifyVector = FALSE)
    tv <- parse_time_v(file.path(base, "markers", "time_v.txt"))
    units[[length(units) + 1]] <- list(id = paste0("clustering/", obj, "/markers"), kind = "markers", dir = file.path(base, "markers"),
      files = list.files(file.path(base, "markers"), recursive = TRUE),
      fields = c(list(tool = "clustering", kind = "markers", dataset = obj, tier = tier,
                      recipe = sprintf("validation/slurm/markers_clustering.sbatch %s (validation/R/markers_clustering.R)", obj),
                      source = mj$rds, source_sha256 = mj$rds_sha256, saved_markers = mj$saved_markers,
                      saved_markers_sha256 = mj$saved_markers_sha256, reference = mj$reference,
                      slurm_job_id = mj$slurm_job_id, env = "envs/r-reference"),
                 fp_fields(mj$fingerprint), list(max_rss_kib = tv$max_rss_kib, elapsed = tv$elapsed)))
  }
  if (file.exists(file.path(base, "optimizer", "optimizer.json"))) {
    units[[length(units) + 1]] <- list(id = paste0("clustering/", obj, "/optimizer"), kind = "optimizer", dir = file.path(base, "optimizer"),
      files = list.files(file.path(base, "optimizer"), recursive = TRUE),
      fields = list(tool = "clustering", kind = "optimizer", dataset = obj, tier = tier,
                    recipe = sprintf("validation/slurm/optimizer_clustering.sbatch %s (validation/R/optimizer_clustering.R; data/clustering/bin/ModularityOptimizer from validation/cpp/build.sh)", obj),
                    source = file.path(base, "saved", "graphs", "snn.mtx.gz"), env = "envs/r-reference (compiler and Makeconf flags)"))
  }
  cmpdir <- file.path("data", "clustering", "compare", obj)
  if (dir.exists(cmpdir)) {
    units[[length(units) + 1]] <- list(id = paste0("clustering/", obj, "/compare"), kind = "compare", dir = cmpdir,
      files = list.files(cmpdir),
      fields = list(tool = "clustering", kind = "compare", dataset = obj, tier = tier,
                    recipe = sprintf("validation/slurm/compare_clustering.sbatch %s (validation/R/compare_clustering.R)", obj),
                    env = "envs/r-reference (plain files)"))
  }

  # Files: each listed once; later units' byte-identical copies at the same relative path go to
  # identical_on of the first unit's entry. Export files under saved/ map to a run's top level.
  seen <- list()   # key: relative path within a run (saved/ prefix dropped) -> entry
  entries <- list()
  for (u in units) {
    runs_out <- c(runs_out, table_lines("run", c(list(id = u$id), u$fields)))
    for (f in sort(u$files)) {
      p <- file.path(u$dir, f)
      if (!file.exists(p)) next
      rel <- if (u$kind == "export" && startsWith(f, "saved/")) sub("^saved/", "", f) else f
      h <- hash_of(p)
      key <- paste(if (u$kind %in% c("export", "direct", "instrumented")) "result" else u$kind, rel, h)
      role <- role_of(u$kind, f)
      if (!is.null(seen[[key]]) && role != "run-record") {
        k <- seen[[key]]
        entries[[k]]$identical_on <- c(entries[[k]]$identical_on, u$id)
        n_copies <- n_copies + 1
        next
      }
      entries[[length(entries) + 1]] <- list(path = p, sha256 = h, bytes = file.size(p), role = role,
                                             stage = if (u$kind %in% c("export", "direct", "instrumented")) stage_of(rel) else "",
                                             run = u$id, identical_on = character())
      seen[[key]] <- length(entries)
    }
  }
  for (e in entries) {
    if (!nzchar(e$stage)) e$stage <- NULL
    if (!length(e$identical_on)) e$identical_on <- NULL
    files_out <- c(files_out, table_lines("file", e, arrays = "identical_on"))
    n_files <- n_files + 1
  }
}

# K0 vectors and the C++ build ------------------------------------------------------------------------
k0 <- file.path("data", "clustering", "k0")
bj <- jsonlite::fromJSON(file.path("data", "clustering", "bin", "build.json"), simplifyVector = FALSE)
runs_out <- c(runs_out, table_lines("run", list(
  id = "clustering/k0", tool = "clustering", kind = "k0", stage = "K0",
  recipe = "bash validation/cpp/build.sh; data/clustering/bin/javarandom_vectors and kiss64_vectors <dir>; validation/R/k0_vectors.R <dir>",
  compiler = bj$compiler_version, cxxflags = bj$cxxflags, seurat_tarball_sha256 = bj$seurat_tarball_sha256,
  kissrandom_h_sha256 = bj$`kissrandom.h_sha256`, env = "envs/r-reference (compiler, Makeconf flags, R 4.5.3)")))
for (p in c(file.path(k0, sort(list.files(k0))), file.path("data", "clustering", "bin", "build.json"))) {
  files_out <- c(files_out, table_lines("file", list(path = p, sha256 = hash_of(p), bytes = file.size(p),
                                                     role = if (basename(p) == "build.json") "run-record" else "reference",
                                                     stage = "K0", run = "clustering/k0")))
  n_files <- n_files + 1
}

# Committed fixtures ------------------------------------------------------------------------------------
fx_root <- file.path("tests", "fixtures", "clustering")
source_of <- function(rel) {
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  if (parts[1] == "k0") return(file.path(k0, parts[2]))
  obj <- parts[1]
  sub_ <- paste(parts[-1], collapse = "/")
  map <- c("input/" = file.path("data/clustering", obj, "input/"),
           "saved/" = file.path("data/clustering", obj, "saved/"),
           "reference/" = file.path("data/clustering", obj, "runs/instrumented-portable-sapphirerapids/"),
           "optimizer/" = file.path("data/clustering", obj, "optimizer/"),
           "markers/" = file.path("data/clustering", obj, "markers/"))
  for (m in names(map)) {
    if (!startsWith(sub_, m)) next
    cand <- paste0(map[[m]], substring(sub_, nchar(m) + 1))
    if (file.exists(cand)) return(cand)
    return(sub("\\.gz$", "", cand))  # gzipped in the fixture only
  }
  NA
}
if (dir.exists(fx_root)) {
  for (f in sort(list.files(fx_root, recursive = TRUE))) {
    p <- file.path(fx_root, f)
    if (f == "README.md") next
    src <- source_of(f)
    gz <- grepl("\\.gz$", f) && !is.na(src) && !grepl("\\.gz$", src) && file.exists(src)
    ok <- !is.na(src) && file.exists(src) &&
      if (gz) identical(digest::digest(memDecompress(readBin(p, "raw", file.size(p)), "gzip"), algo = "sha256", serialize = FALSE), hash_of(src))
      else identical(hash_of(p), hash_of(src))
    if (!ok) stop("fixture ", p, " is not a copy of ", src)
    files_out <- c(files_out, table_lines("file", list(path = p, sha256 = hash_of(p), bytes = file.size(p), role = "fixture",
                                                       source = src, note = if (gz) "gzip of source" else NULL)))
    n_files <- n_files + 1
  }
}

header <- c(
  "# Clustering (plan C4) ------------------------------------------------------------------------",
  "#",
  "# Milestone M8: the R-side reference data of Step 2, the clustering crate. Generated by",
  "# validation/R/manifest_clustering.R from the files on disk (do not edit by hand). Objects and their",
  "# pipeline paths: validation/clustering_datasets.tsv (tiers T and P); the full panel's 840 objects:",
  "# validation/clustering_panel.tsv. Results and checks: docs/validation/clustering_m8.md.",
  "#",
  "# Layout under data/clustering/<object>/ (gitignored; regenerate from the recipes):",
  "#   input/      the chain's input, the object's RNA counts as saved (10x Matrix Market; for the",
  "#               DoubletFinder inputs and outputs, relative links to their M2 and M6b exports);",
  "#               assay.txt (v3 or v5); export.json",
  "#   saved/      what the chain left in the saved object (export_clustering.R): metadata.tsv, the",
  "#               levels, sct/ (corrected counts, model, scale.data at levels full and standard,",
  "#               log1p at full), pca/, graphs/, object.json, markergenes.csv (the workflow's table)",
  "#   runs/<mode>-<profile>-<cpu arch>[-noumap][-v5]/   reference runs, in saved/'s layout, with",
  "#               cluster_ids.txt, nclusters.txt, summary.json, time_v.txt; instrumented runs add",
  "#               stages/ (index.tsv explains every array; conventions in",
  "#               validation/R/lib/clustering_hooks.R)",
  "#   markers/    the marker steps on the saved .rds and reference_markergenes.csv",
  "#   optimizer/  the instrumented modularity optimizer on saved/graphs/snn.mtx.gz",
  "# data/clustering/compare/<object>/: compare_clustering.R's checks. data/clustering/k0/: the K0",
  "# vectors. tests/fixtures/clustering/: the committed subset cargo test uses (its README).",
  "#",
  "# Compact entries: each [[run]] table (ids start with 'clustering/') holds what its files share (recipe, kind, profile, node,",
  "# CPU model, OpenBLAS core and threads, glibc, R, SLURM job, peak RSS and wall time from",
  "# /usr/bin/time -v); each [[file]] table gives path, sha256, bytes, role, stage (the plan's C5",
  "# check it serves: K1 SCTransform, K2 PCA, K3 items, K4 index, K5 neighbours, K6 graphs,",
  "# K7 clustering, K8 silhouette, K9 purity, K10 corrected counts), its run, and identical_on: the",
  "# other runs and exports of the same object whose file at the same relative path (saved/ counted",
  "# as a run's top level) has the same sha256. Those copies are not listed again. Run records",
  "# (summary.json, time_v.txt, fingerprint_wrapper.json) differ between runs by construction.",
  "")

old <- readLines(manifest)
df_mark <- grep("^# DoubletFinder \\(plan B4\\) ", old)
if (length(df_mark) != 1) stop("DoubletFinder marker line not found in ", manifest)
cl_mark <- grep("^# Clustering \\(plan C4\\) ", old)
before <- if (length(cl_mark) == 1) old[seq_len(cl_mark - 1)] else old[seq_len(df_mark - 1)]
after <- old[df_mark:length(old)]
writeLines(c(before, header, files_out, "# Runs (clustering, M8) --------------------------------------------------------------------", "",
             runs_out, after), manifest)
cat("wrote", n_files, "file entries (", n_copies, "byte-identical copies recorded in identical_on ) and",
    length(grep("^\\[\\[run\\]\\]", runs_out)), "run entries to", manifest, "\n")
