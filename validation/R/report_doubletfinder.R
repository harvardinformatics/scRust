# Markdown tables for docs/validation/doubletfinder_m2.md from the M2 outputs (comparison JSON
# and TSV files, run summaries, /usr/bin/time -v and sacct). Reads only plain files and sacct.
#
#   Rscript validation/R/report_doubletfinder.R > <tables.md>

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
setwd(normalizePath(file.path(here, "..", "..")))

ds <- read.delim("validation/doubletfinder_datasets.tsv", colClasses = "character", quote = "")
js <- function(p) if (file.exists(p)) jsonlite::fromJSON(p) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
sacct_job <- function(job) {
  if (is.null(job) || !nzchar(job)) return(NULL)
  x <- system2("sacct", c("-j", job, "--format=JobID,NodeList,Elapsed,MaxRSS,State", "-P", "-n"), stdout = TRUE)
  x <- strsplit(x, "|", fixed = TRUE)
  main <- x[[1]]
  rss <- vapply(x, function(r) r[4], "")
  rss_kib <- suppressWarnings(max(vapply(rss, function(s) {
    if (!nzchar(s)) return(NA_real_)
    u <- substring(s, nchar(s)); v <- as.numeric(substring(s, 1, nchar(s) - 1))
    switch(u, K = v, M = v * 1024, G = v * 1024^2, as.numeric(s) / 1024)
  }, 0), na.rm = TRUE))
  list(node = main[2], elapsed = main[3], max_rss_kib = rss_kib, state = main[5])
}
gib <- function(kib) sprintf("%.1f", kib / 1024^2)
time_v <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p)
  list(rss = as.numeric(sub(".*: ", "", grep("Maximum resident set size", l, value = TRUE))),
       elapsed = sub(".*: ", "", grep("Elapsed \\(wall clock\\)", l, value = TRUE)))
}
hms_to_min <- function(s) {
  p <- as.numeric(strsplit(s, ":")[[1]])
  if (length(p) == 3) p[1] * 60 + p[2] + p[3] / 60 else p[1] + p[2] / 60
}

# 0. Objects ------------------------------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(ds))) {
  o <- ds$object[i]
  e <- js(file.path("data/doubletfinder", o, "input", "export.json"))
  s <- js(file.path("data/doubletfinder", o, "saved", "saved.json"))
  if (is.null(e)) next
  rows[[length(rows) + 1]] <- data.frame(
    object = o, tier = ds$tier[i], features = format(e$n_features, big.mark = ","),
    cells = format(e$n_cells, big.mark = ","), entries = format(e$nnz, big.mark = ","),
    `explicit zeros` = format(e$explicit_zeros, big.mark = ","), clusters = length(e$cluster_levels),
    saved = sprintf("`%s` (%s, %s)", s$pANN_column, s$pK_label, s$nExp_label), check.names = FALSE)
}
cat("## Objects\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")
cat("\n")

# 1. Cross-check against the saved workflow outputs -------------------------------------------
rows <- list()
for (o in ds$object) {
  for (f in Sys.glob(file.path("data/doubletfinder/compare", o, "crosscheck_*.json"))) {
    r <- js(f)
    run <- sub("^crosscheck_(.*)\\.json$", "\\1", basename(f))
    rows[[length(rows) + 1]] <- data.frame(
      object = o, run = run, pK = r$pK_run, nExp = r$nExp_run,
      `columns` = if (isTRUE(r$column_names_identical)) "identical" else "DIFFER",
      `doublet set` = sprintf("%s (%d)", if (isTRUE(r$doublet_set_identical)) "identical" else "DIFFERS", r$n_doublets_run),
      `singlet pANN identical` = sprintf("%d / %d", r$n_singlet_pANN_identical, r$n_singlets),
      result = if (isTRUE(r$pass)) "pass" else "FAIL", check.names = FALSE)
  }
}
cat("## Cross-check against the saved workflow outputs\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")

# 2. Instrumented vs direct, CPU types ----------------------------------------------------------
rows <- list()
for (o in ds$object) {
  for (f in Sys.glob(file.path("data/doubletfinder/compare", o, "*.json"))) {
    b <- basename(f)
    if (startsWith(b, "crosscheck_")) next
    r <- js(f)
    rows[[length(rows) + 1]] <- data.frame(
      object = o, comparison = sub("\\.json$", "", b),
      `final files identical` = paste(unlist(r$final_files_identical), collapse = "/"),
      `result identical` = r$result_identical,
      `stage arrays identical` = if (!is.null(r$stage_arrays)) sprintf("%d / %d", r$stage_arrays_identical, r$stage_arrays) else "",
      `RNG checkpoints identical` = if (!is.null(r$rng_checkpoints_identical)) as.character(r$rng_checkpoints_identical) else "",
      `calls differ` = r$n_calls_differ, `pANN differ` = r$n_pANN_differ, check.names = FALSE)
  }
}
cat("\n## Run against run\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")

# 3. Floor: portable vs production, per stage array --------------------------------------------
blocks <- list()
for (o in ds$object) {
  f <- file.path("data/doubletfinder/compare", o, "portable_vs_production_stages.tsv")
  if (!file.exists(f)) next
  t <- read.delim(f, quote = "", comment.char = "", colClasses = c(note = "character"))
  t$object <- o
  t$quantity <- sub("\\.npy$", "", sub("\\.(sweep_pN[0-9.]+|final)\\.npy$", "", t$file))
  blocks[[o]] <- t
}
if (length(blocks)) {
  B <- do.call(rbind, blocks)
  B$note[is.na(B$note)] <- ""
  agg <- do.call(rbind, lapply(split(B, B$quantity), function(s) {
    d <- s[!is.na(s$n_diff) & s$n_diff > 0, ]
    flips <- suppressWarnings(as.integer(sub(" component.*", "", s$note[grepl("sign-flipped", s$note)])))
    data.frame(stage = s$stage[1], quantity = s$quantity[1],
      objects = length(unique(s$object)), `objects with differences` = length(unique(d$object)),
      `elements differing` = sprintf("%d / %d", sum(s$n_diff, na.rm = TRUE), sum(s$n, na.rm = TRUE)),
      `max abs diff` = if (nrow(d)) sprintf("%.3g", max(d$max_abs)) else "0",
      `max diff / max abs value` = if (nrow(d)) sprintf("%.3g", max(d$max_scaled)) else "0",
      note = if (length(flips)) sprintf("after aligning signs; %d of %d blocks had 1 to %d components flipped",
                                        length(flips), nrow(s), max(flips)) else "",
      check.names = FALSE)
  }))
  agg <- agg[order(agg$stage, agg$quantity), ]
  cat("\n## R-vs-R floor: OPENBLAS_CORETYPE=Haswell against SapphireRapids, same CPU type\n\n")
  cat(md_table(agg), sep = "\n")
}

# 4. Memory and time ----------------------------------------------------------------------------
# DoubletFinder's own cost: the direct run in the production profile (the workflow's), against
# the workflow job (which also ran Seurat after DoubletFinder; sacct samples every 30 s). The
# instrumented run's cost is shown for the overhead of the dumps.
rows <- list()
for (i in seq_len(nrow(ds))) {
  o <- ds$object[i]
  wf <- sacct_job(ds$workflow_job[i])
  rd <- file.path("data/doubletfinder", o, "runs", "direct-production-sapphirerapids")
  if (!dir.exists(rd)) next
  s <- js(file.path(rd, "summary.json"))
  tv <- time_v(file.path(rd, "time_v.txt"))
  sj <- sacct_job(s$run$slurm_job_id)
  ph <- s$run$phases
  inst <- Filter(dir.exists, file.path("data/doubletfinder", o, "runs",
                 c("instrumented-production-sapphirerapids", "instrumented-portable-sapphirerapids")))
  itv <- if (length(inst)) time_v(file.path(inst[1], "time_v.txt")) else NULL
  rows[[length(rows) + 1]] <- data.frame(
    object = o, job = s$run$slurm_job_id,
    `peak RSS GiB` = gib(tv$rss), `sacct MaxRSS GiB` = gib(sj$max_rss_kib),
    `wall min` = sprintf("%.1f", hms_to_min(tv$elapsed)),
    `load min` = sprintf("%.1f", s$run$load_input_s / 60),
    `paramSweep min` = sprintf("%.1f", ph$paramSweep$elapsed_s / 60),
    `doubletFinder min` = sprintf("%.1f", (ph$doubletFinder$elapsed_s - ph$sweep_cleanup$elapsed_s) / 60),
    `peak after sweep GiB` = gib(ph$sweep_cleanup$peak_rss_kib),
    `workflow job` = ds$workflow_job[i], `workflow MaxRSS GiB` = gib(wf$max_rss_kib),
    `workflow min` = sprintf("%.1f", hms_to_min(wf$elapsed)),
    `instrumented peak GiB / min` = if (is.null(itv)) "" else sprintf("%s / %.1f (%s)", gib(itv$rss), hms_to_min(itv$elapsed), sub(".*instrumented-([a-z]+)-.*", "\\1", inst[1])),
    check.names = FALSE)
}
cat("\n## Memory and time\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")
