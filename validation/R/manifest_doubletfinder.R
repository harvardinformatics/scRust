# Writes the DoubletFinder part of validation/manifest.toml (docs/PLAN.md A4 and A5b, milestone
# M2) from the files on disk: every export, reference run and committed fixture under
# data/doubletfinder/ and tests/fixtures/doubletfinder/. Reads only plain files.
#
#   Rscript validation/R/manifest_doubletfinder.R
#
# The manifest's text before the line "# DoubletFinder (plan B4) " is kept as is; everything
# from that line on is regenerated. Replica runs (the portable profile on other CPU types,
# made to show the portable references are CPU-independent) get a [[run]] table each; their
# files are not listed again when they are byte-identical to the primary run's, which instead
# records the replicas in identical_on.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)

ds <- read.delim("validation/doubletfinder_datasets.tsv", colClasses = "character", quote = "")
pipeline <- "scRNAseq-preprocessing"
manifest <- "validation/manifest.toml"
primary_arch <- "sapphirerapids"

toml_str <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")
toml_val <- function(v) {
  if (is.logical(v)) return(tolower(as.character(v)))
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  if (length(v) != 1) return(paste0("[", paste(vapply(v, toml_str, ""), collapse = ", "), "]"))
  toml_str(v)
}
table_lines <- function(kind, fields) {
  fields <- fields[!vapply(fields, function(v) is.null(v) || (length(v) == 1 && is.na(v)), logical(1))]
  c(sprintf("[[%s]]", kind),
    vapply(names(fields), function(k) {
      v <- fields[[k]]
      if (is.character(v) && length(v) > 1 || (is.character(v) && k == "identical_on")) {
        return(paste0(k, " = [", paste(vapply(v, toml_str, ""), collapse = ", "), "]"))
      }
      paste0(k, " = ", toml_val(v))
    }, ""),
    "")
}

fp_string <- function(fp, job) {
  if (is.null(fp) || !is.null(fp$error)) return("unknown")
  sprintf("%s; OpenBLAS %s, %s thread(s); glibc %s; R %s; %s, job %s",
          fp$cpu_model, fp$openblas_corename, fp$openblas_num_threads, fp$glibc,
          sub("^R version ([0-9.]+).*", "\\1", fp$R), sub("\\..*", "", fp$host), job)
}

parse_time_v <- function(path) {
  if (!file.exists(path)) return(list(max_rss_kib = NA, elapsed = NA))
  l <- readLines(path)
  rss <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", l, value = TRUE)))
  el <- sub(".*: ", "", grep("Elapsed \\(wall clock\\)", l, value = TRUE))
  list(max_rss_kib = rss, elapsed = el)
}

stage_of <- function(rel) {
  b <- basename(rel)
  if (grepl("^S[0-9][a-d]?_", b)) return(sub("^(S[0-9][a-d]?)_.*", "\\1", b))
  switch(b,
    "cells.tsv" = "S7", "summary.json" = "S7", "sweep_stats.tsv" = "S6", "bcmvn.tsv" = "S6",
    "find_pK.pdf" = "S6", "index.tsv" = "S1-S7", "rng.npy" = "S1-S7", "rng.tsv" = "S1-S7",
    "stages.json" = "S1-S7", "")
}

file_hash <- new.env()
hash_of <- function(path) {
  if (is.null(file_hash[[path]])) file_hash[[path]] <- sha256_file(path)
  file_hash[[path]]
}

# Files whose bytes depend on the run (time, node, job, plot timestamps) and are therefore never
# compared by checksum across runs.
run_specific <- c("summary.json", "find_pK.pdf", "time_v.txt", "fingerprint_wrapper.json",
                  "export_run.json", "export_time_v.txt", "export_singlets_run.json",
                  "export_singlets_time_v.txt")

out <- character()
runs_out <- character()
n_files <- 0

for (k in seq_len(nrow(ds))) {
  obj <- ds$object[k]
  base <- file.path("data", "doubletfinder", obj)
  if (!dir.exists(base)) next
  tier <- ds$tier[k]
  src_input <- file.path(pipeline, ds$input[k])
  src_saved <- file.path(pipeline, ds$saved[k])
  exp_info <- if (file.exists(file.path(base, "input", "export.json"))) jsonlite::fromJSON(file.path(base, "input", "export.json")) else list()
  sav_info <- if (file.exists(file.path(base, "saved", "saved.json"))) jsonlite::fromJSON(file.path(base, "saved", "saved.json")) else list()
  exp_run <- if (file.exists(file.path(base, "export_run.json"))) jsonlite::fromJSON(file.path(base, "export_run.json")) else NULL
  exp_fp <- if (is.null(exp_run)) "unknown" else fp_string(exp_run$fingerprint, exp_run$slurm_job_id)
  exp_recipe <- sprintf("validation/slurm/export_doubletfinder.sbatch %s (validation/R/export_doubletfinder.R)", obj)
  # M6b: the saved output's RNA counts (the singlets'), exported into saved/singlets/
  sgl_run <- if (file.exists(file.path(base, "export_singlets_run.json"))) jsonlite::fromJSON(file.path(base, "export_singlets_run.json")) else NULL
  sgl_fp <- if (is.null(sgl_run)) "unknown" else fp_string(sgl_run$fingerprint, sgl_run$slurm_job_id)
  sgl_recipe <- sprintf("validation/slurm/export_doubletfinder_singlets.sbatch %s (validation/R/export_doubletfinder_singlets.R)", obj)

  out <- c(out, sprintf("# %s (tier %s, %s) %s", obj, tier, ds$role[k], strrep("-", max(3, 80 - nchar(obj) - nchar(tier) - nchar(ds$role[k]) - 14))), "")
  for (sub in c("input", "saved")) {
    d <- file.path(base, sub)
    if (!dir.exists(d)) next
    src <- if (sub == "input") sprintf("%s (sha256 %s)", src_input, exp_info$source_sha256) else sprintf("%s (sha256 %s)", src_saved, sav_info$source_sha256)
    for (f in sort(list.files(d, recursive = TRUE))) {
      p <- file.path(d, f)
      singlets <- startsWith(f, "singlets/")
      out <- c(out, table_lines("file", list(
        path = p, sha256 = hash_of(p), bytes = file.size(p), tool = "doubletfinder", dataset = obj,
        tier = tier, role = if (sub == "input") "input" else "saved-output", stage = "", profile = "none (no BLAS)",
        source = src, recipe = if (singlets) sgl_recipe else exp_recipe, env = "envs/r-reference",
        fingerprint = if (singlets) sgl_fp else exp_fp)))
      n_files <- n_files + 1
    }
  }
  for (f in c("export_run.json", "export_time_v.txt", "export_singlets_run.json", "export_singlets_time_v.txt")) {
    p <- file.path(base, f)
    if (!file.exists(p)) next
    singlets <- grepl("singlets", f, fixed = TRUE)
    out <- c(out, table_lines("file", list(
      path = p, sha256 = hash_of(p), bytes = file.size(p), tool = "doubletfinder", dataset = obj,
      tier = tier, role = "run-record", stage = "", profile = "none (no BLAS)",
      source = if (singlets) src_saved else src_input,
      recipe = if (singlets) sgl_recipe else exp_recipe, env = "envs/r-reference",
      fingerprint = if (singlets) sgl_fp else exp_fp,
      note = "run-specific (timings, node, job); not compared by checksum")))
    n_files <- n_files + 1
  }

  run_dirs <- sort(list.dirs(file.path(base, "runs"), recursive = FALSE))
  run_dirs <- run_dirs[!grepl("\\.partial\\.", run_dirs)]
  for (rd in run_dirs) {
    run <- basename(rd)
    parts <- strsplit(run, "-", fixed = TRUE)[[1]]
    mode <- parts[1]; profile <- parts[2]; arch <- paste(parts[-(1:2)], collapse = "-")
    summ <- jsonlite::fromJSON(file.path(rd, "summary.json"))
    fp <- summ$run$fingerprint
    job <- summ$run$slurm_job_id
    tv <- parse_time_v(file.path(rd, "time_v.txt"))
    level <- if (mode == "instrumented") jsonlite::fromJSON(file.path(rd, "stages", "stages.json"))$level else NA
    recipe <- sprintf("validation/slurm/reference_doubletfinder.sbatch %s %s %s (validation/R/%s_doubletfinder.R%s)",
                      profile, mode, obj, if (mode == "direct") "reference" else "instrumented",
                      if (mode == "instrumented") paste0(" --level ", level) else "")
    files <- sort(list.files(rd, recursive = TRUE))
    primary <- arch == primary_arch
    # replicas: same mode and profile on other CPU types
    replicas <- run_dirs[basename(run_dirs) != run & startsWith(basename(run_dirs), paste0(mode, "-", profile, "-"))]
    rep_identical <- character()
    if (primary) {
      for (r2 in replicas) {
        f2 <- sort(list.files(r2, recursive = TRUE))
        cmp <- setdiff(files, run_specific)
        cmp <- cmp[!basename(cmp) %in% run_specific]
        same <- identical(setdiff(f2[!basename(f2) %in% run_specific], character()), cmp) &&
          all(vapply(cmp, function(f) identical(hash_of(file.path(rd, f)), hash_of(file.path(r2, f))), logical(1)))
        s2 <- jsonlite::fromJSON(file.path(r2, "summary.json"))
        same <- same && identical(s2$result, summ$result)
        lab <- sprintf("%s (%s, %s, job %s)", sub(".*-", "", basename(r2)), s2$run$fingerprint$cpu_model,
                       sub("\\..*", "", s2$run$fingerprint$host), s2$run$slurm_job_id)
        if (same) rep_identical <- c(rep_identical, lab) else stop("replica ", r2, " differs from ", rd)
      }
    }
    runs_out <- c(runs_out, table_lines("run", list(
      id = paste(obj, run, sep = "/"), path = rd, tool = "doubletfinder", dataset = obj, tier = tier,
      mode = mode, profile = profile, cpu_arch = arch, level = level, slurm_job_id = job,
      node = fp$host, cpu_model = fp$cpu_model, openblas_core = fp$openblas_corename,
      openblas_threads = fp$openblas_num_threads, glibc = fp$glibc, R = fp$R,
      max_rss_kib = tv$max_rss_kib, elapsed = tv$elapsed,
      replica_of = if (!primary && file.exists(file.path(base, "runs", paste(mode, profile, primary_arch, sep = "-")))) paste(obj, paste(mode, profile, primary_arch, sep = "-"), sep = "/") else NA,
      pK = summ$result$pK_label, nExp = summ$result$nExp)))
    if (!primary) next   # replica files are identical to the primary run's (checked above)
    for (f in files) {
      p <- file.path(rd, f)
      rs <- basename(f) %in% run_specific
      out <- c(out, table_lines("file", list(
        path = p, sha256 = hash_of(p), bytes = file.size(p), tool = "doubletfinder", dataset = obj,
        tier = tier, role = if (rs && !basename(f) %in% c("summary.json")) "run-record" else "reference",
        stage = stage_of(f), profile = profile,
        source = file.path(base, "input"), recipe = recipe, env = "envs/r-reference",
        fingerprint = fp_string(fp, job),
        identical_on = if (!rs && length(rep_identical)) rep_identical else NULL,
        note = if (rs) "run-specific bytes (timings, node, job or plot timestamp); summary.json's \"result\" section is deterministic" else NULL)))
      n_files <- n_files + 1
    }
  }

  # M2 checks (compare_doubletfinder.sbatch): cross-checks and run-against-run comparisons
  cmpdir <- file.path("data", "doubletfinder", "compare", obj)
  if (dir.exists(cmpdir)) {
    for (f in sort(list.files(cmpdir))) {
      p <- file.path(cmpdir, f)
      out <- c(out, table_lines("file", list(
        path = p, sha256 = hash_of(p), bytes = file.size(p), tool = "doubletfinder", dataset = obj,
        tier = tier, role = "check", stage = "", profile = "none (plain files)",
        source = file.path(base, "runs"),
        recipe = sprintf("validation/slurm/compare_doubletfinder.sbatch %s (validation/R/compare_doubletfinder.R)", obj),
        env = "envs/r-reference", fingerprint = "plain-file comparison, no BLAS")))
      n_files <- n_files + 1
    }
  }

  # committed fixtures for cargo test
  fx <- file.path("tests", "fixtures", "doubletfinder", obj)
  if (dir.exists(fx)) {
    for (f in sort(list.files(fx, recursive = TRUE))) {
      p <- file.path(fx, f)
      top <- strsplit(f, "/", fixed = TRUE)[[1]][1]
      orig <- switch(top,
        input = file.path(base, f), saved = file.path(base, f),
        reference = file.path(base, "runs", paste0("direct-portable-", primary_arch), sub("^reference/", "", f)),
        stages = file.path(base, "runs", paste0("instrumented-portable-", primary_arch), f), NA)
      if (is.na(orig) || !file.exists(orig) || !identical(hash_of(p), hash_of(orig))) stop("fixture ", p, " is not a copy of ", orig)
      out <- c(out, table_lines("file", list(
        path = p, sha256 = hash_of(p), bytes = file.size(p), tool = "doubletfinder", dataset = obj,
        tier = tier, role = switch(top, input = "input", saved = "saved-output", "reference"),
        stage = if (top %in% c("input", "saved")) "" else stage_of(f),
        profile = if (top %in% c("input", "saved")) "none (no BLAS)" else "portable",
        source = orig, recipe = "copy of the data/ file named in source", env = "envs/r-reference",
        fingerprint = "as the source file's entry")))
      n_files <- n_files + 1
    }
  }
}

header <- c(
  "# DoubletFinder (plan B4) ---------------------------------------------------------------------",
  "#",
  "# Milestone M2: exports, direct and instrumented reference runs, generated by",
  "# validation/R/manifest_doubletfinder.R from the files on disk (do not edit by hand).",
  "# Datasets and their pipeline paths: validation/doubletfinder_datasets.tsv. Results and checks:",
  "# docs/validation/doubletfinder_m2.md.",
  "#",
  "# Layout under data/doubletfinder/<object>/ (gitignored; regenerate from the recipes):",
  "#   input/   the counts DoubletFinder reads (RNA \"counts\" layer) as 10x Matrix Market, and",
  "#            seurat_clusters (export_doubletfinder.R); identical to the object's matrix, bit for bit",
  "#   saved/   the saved workflow output's compared fields: kept barcodes and singlet pANN, the",
  "#            doublet set, column names (pK, nExp); saved/singlets/ (M6b): its RNA counts, the",
  "#            singlets' (export_doubletfinder_singlets.R)",
  "#   runs/<mode>-<profile>-<cpu arch>/   reference runs: cells.tsv, sweep_stats.tsv, bcmvn.tsv,",
  "#            summary.json, find_pK.pdf, time_v.txt; instrumented runs add stages/ (index.tsv",
  "#            explains every array; conventions in validation/R/lib/doubletfinder_hooks.R)",
  "# tests/fixtures/doubletfinder/<object>/ holds the tier-T subset cargo test uses: input/, saved/,",
  "# reference/ (portable direct run) and stages/ (portable instrumented run, without the",
  "# per-context residual-sized matrices and S4 position matrices, which stay in data/).",
  "#",
  "# Fields beyond the schema above: role \"run-record\" (logs, timings, plots) and \"check\" (the",
  "# outputs of compare_doubletfinder.R under data/doubletfinder/compare/); note; and",
  "# identical_on, the replica runs (portable profile, other CPU types) whose file had the same",
  "# sha256. Files listed as run-specific differ between runs by construction (timestamps, node,",
  "# job). The [[run]] tables record every job: node, CPU, OpenBLAS core, peak RSS and wall time",
  "# from /usr/bin/time -v.",
  "")

old <- readLines(manifest)
cut <- grep("^# DoubletFinder \\(plan B4\\) ", old)
if (length(cut) != 1) stop("marker line not found in ", manifest)
keep <- old[seq_len(cut - 1)]
writeLines(c(keep, header, out, "# Reference runs (DoubletFinder, M2) ----------------------------------------------------------", "", runs_out), manifest)
cat("wrote", n_files, "file entries and", length(grep("^\\[\\[run\\]\\]", runs_out)), "run entries to", manifest, "\n")
