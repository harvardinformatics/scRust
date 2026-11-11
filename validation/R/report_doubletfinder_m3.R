# Markdown tables for docs/validation/doubletfinder_m3.md from the M3 outputs: the Rust runs on
# R's own embeddings (data/doubletfinder/<object>/rust/m3/<run>/, with their scrust validate
# JSON), the direct R runs and the saved workflow outputs, and the benchmark
# (data/doubletfinder/<object>/rust/m3-benchmark/<job>/). Reads only plain files.
#
#   Rscript validation/R/report_doubletfinder_m3.R [<benchmark job id>] > <tables.md>

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
setwd(normalizePath(file.path(here, "..", "..")))
args <- commandArgs(trailingOnly = TRUE)

ds <- read.delim("validation/doubletfinder_datasets.tsv", colClasses = "character", quote = "")
js <- function(p) if (file.exists(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
time_v <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p)
  g <- function(pat) sub(".*: ", "", grep(pat, l, value = TRUE, fixed = TRUE))
  el <- g("Elapsed (wall clock) time")
  pt <- as.numeric(strsplit(el, ":")[[1]])
  secs <- if (length(pt) == 3) pt[1] * 3600 + pt[2] * 60 + pt[3] else pt[1] * 60 + pt[2]
  list(rss_kib = as.numeric(g("Maximum resident set size")), wall_s = secs,
       user_s = as.numeric(g("User time (seconds)")), cpu = g("Percent of CPU this job got"))
}
fmt_s <- function(s) if (s >= 120) sprintf("%.1f min", s / 60) else sprintf("%.1f s", s)
fmt_mem <- function(kib) if (kib >= 1024^2) sprintf("%.1f GiB", kib / 1024^2) else sprintf("%.0f MiB", kib / 1024)
check_of <- function(v, name) {
  for (c in v$checks) if (identical(c$name, name)) return(c)
  NULL
}
cell <- function(c) {
  if (is.null(c)) return("")
  if (identical(c$status, "identical") && grepl("no subsample", c$detail)) return("= (none)")
  switch(c$status,
         identical = if (c$n > 0) sprintf("= (%s)", num(c$n)) else "=",
         differ = sprintf("**%s of %s differ**", num(c$n_differ), num(c$n)),
         missing = "not dumped", `n/a` = "n/a", c$status)
}

# 1. Rust on R's embeddings against the same R run -----------------------------------------------
rows <- list()
pos_rows <- list()
for (o in ds$object) {
  for (d in sort(Sys.glob(file.path("data/doubletfinder", o, "rust/m3/instrumented-*")))) {
    v <- js(file.path(d, "validate.json"))
    if (is.null(v)) next
    run <- basename(d)
    res <- Filter(function(c) identical(c$stage, "result"), v$checks)
    n_res_id <- sum(vapply(res, function(c) identical(c$status, "identical"), TRUE))
    files <- vapply(c("cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"), function(f) {
      c <- check_of(v, f); if (!is.null(c) && identical(c$status, "identical")) "=" else "**differ**"
    }, "")
    rows[[length(rows) + 1]] <- data.frame(
      object = o, run = sub("^instrumented-", "", run),
      `S1 subsample` = cell(check_of(v, "S1_subsample (from R's RNG state)")),
      `S1 doublets` = cell(check_of(v, "S1_doublets (from R's RNG states)")),
      `S4 digests` = cell(check_of(v, "S4_digest.npy")),
      `S5 sweep` = cell(check_of(v, "S5_pann_sweep.npy")),
      `S5 final` = cell(check_of(v, "S5_pann_final.npy")),
      `S6 bkde` = cell(check_of(v, "S6_bkde.npy")),
      `S6 gkde` = cell(check_of(v, "S6_gkde.npy")),
      `S6 BCreal (R's S5)` = cell(check_of(v, "BCreal (teacher-forced on R's S5)")),
      `S7 calls` = cell(check_of(v, "S7_doublet.npy")),
      `S7 (R's S5)` = cell(check_of(v, "S7_doublet (teacher-forced on R's S5)")),
      `cells / sweep_stats / bcmvn` = paste(files, collapse = " / "),
      `result fields` = sprintf("%d / %d", n_res_id, length(res) - 1L),
      pass = if (isTRUE(v$pass)) "yes" else "**no**",
      check.names = FALSE)
    pc <- Filter(function(c) startsWith(c$name, "S4_positions."), v$checks)
    if (length(pc)) {
      n <- sum(vapply(pc, function(c) as.numeric(c$n), 0))
      nd <- sum(vapply(pc, function(c) as.numeric(c$n_differ), 0))
      pos_rows[[length(pos_rows) + 1]] <- data.frame(object = o, run = sub("^instrumented-", "", run),
        contexts = length(pc), positions = num(n), differing = num(nd), check.names = FALSE)
    }
  }
}
cat("## Rust on R's embeddings, against the same R run\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")
cat("\n\n## S4 order positions (tier T, dumped in full)\n\n")
cat(md_table(do.call(rbind, pos_rows)), sep = "\n")
cat("\n")

# 2. Against every direct R run and the saved workflow output -------------------------------------
read_tsv <- function(p) read.delim(p, colClasses = "character", quote = "", comment.char = "")
rows <- list()
for (o in ds$object) {
  rust_runs <- sort(Sys.glob(file.path("data/doubletfinder", o, "rust/m3/instrumented-*")))
  if (!length(rust_runs)) next
  direct <- sort(Sys.glob(file.path("data/doubletfinder", o, "runs/direct-*")))
  direct <- direct[!grepl("\\.partial\\.", direct)]
  saved <- file.path("data/doubletfinder", o, "saved")
  sj <- js(file.path(saved, "saved.json"))
  sc <- read_tsv(file.path(saved, "cells.tsv"))
  sd <- readLines(file.path(saved, "doublets.txt"))
  for (d in rust_runs) {
    same <- vapply(direct, function(r) all(vapply(c("cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"), function(f)
      identical(unname(tools::md5sum(file.path(d, f))), unname(tools::md5sum(file.path(r, f)))), TRUE)), TRUE)
    cells <- read_tsv(file.path(d, "cells.tsv"))
    res <- js(file.path(d, "summary.json"))$result
    singlets <- cells[cells$classification == "Singlet", ]
    dbl <- cells$barcode[cells$classification == "Doublet"]
    rows[[length(rows) + 1]] <- data.frame(
      object = o, `Rust run (embeddings of)` = sub("^instrumented-", "", basename(d)),
      `direct R runs with identical files` = sprintf("%d / %d", sum(same), length(same)),
      `saved: doublet set` = if (setequal(dbl, sd) && length(dbl) == length(sd)) sprintf("= (%s)", num(length(sd))) else "**differs**",
      `saved: singlet pANN` = if (identical(singlets$barcode, sc$barcode) && identical(singlets$pANN, sc$pANN)) sprintf("= (%s)", num(nrow(sc))) else "**differs**",
      `saved: columns` = if (identical(res$pANN_column, sj$pANN_column) && identical(res$classification_column, sj$classification_column)) sprintf("`%s`", res$pANN_column) else "**differ**",
      check.names = FALSE)
  }
}
cat("\n## Against the direct R runs and the saved workflow output\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")
cat("\n")

# 3. Time and memory of the validation runs (16 threads, with stage dumps) -------------------------
rows <- list()
for (o in ds$object) {
  for (d in sort(Sys.glob(file.path("data/doubletfinder", o, "rust/m3/instrumented-*")))) {
    t <- time_v(file.path(d, "time_v.txt"))
    s <- js(file.path(d, "summary.json"))
    if (is.null(t) || is.null(s)) next
    rows[[length(rows) + 1]] <- data.frame(object = o, run = sub("^instrumented-", "", basename(d)),
      `real cells` = num(s$result$n_cells), `with doublets` = num(s$result$n_cells_with_doublets),
      k = num(s$result$k), threads = s$run$threads, dumps = s$run$dump_stages,
      wall = fmt_s(t$wall_s), `peak RSS` = fmt_mem(t$rss_kib), check.names = FALSE)
  }
}
cat("\n## Time and memory of the validation runs\n\n")
cat(md_table(do.call(rbind, rows)), sep = "\n")
cat("\n")

# 4. Benchmark ------------------------------------------------------------------------------------
bench <- if (length(args)) file.path("data/doubletfinder/soupx_emptydrops_sscovelli/rust/m3-benchmark", args[1]) else
  tail(sort(Sys.glob("data/doubletfinder/soupx_emptydrops_sscovelli/rust/m3-benchmark/*")), 1)
if (length(bench) && dir.exists(bench)) {
  rows <- list()
  for (d in sort(Sys.glob(file.path(bench, "*_t*")))) {
    t <- time_v(file.path(d, "time_v.txt"))
    if (is.null(t)) next
    kind <- sub("_t[0-9]+$", "", basename(d))
    threads <- as.integer(sub(".*_t", "", basename(d)))
    s <- js(file.path(d, "summary.json"))
    phase <- function(p) if (!is.null(s$run$phases[[p]])) s$run$phases[[p]]$elapsed_s else NA
    rows[[length(rows) + 1]] <- data.frame(
      what = if (kind == "knn") "kNN: final run (79,409 queries, 105,879 cells, k + 1 = 20,118)" else "whole core (sweep, summary, pK, nExp, final run, outputs)",
      threads = threads, wall = fmt_s(t$wall_s), `CPU time` = fmt_s(t$user_s), `CPU %` = t$cpu,
      `peak RSS` = fmt_mem(t$rss_kib),
      `sweep (6 pN)` = if (kind == "core") fmt_s(phase("paramSweep") - phase("load")) else "",
      `final run` = if (kind == "core") fmt_s(phase("doubletFinder") - phase("find.pK")) else "",
      check.names = FALSE)
  }
  rows <- do.call(rbind, rows)
  rows <- rows[order(rows$what != rows$what[grepl("^kNN", rows$what)][1], rows$threads), ]
  cat(sprintf("\n## Benchmark (%s)\n\n", bench))
  cat(md_table(rows), sep = "\n")
  r <- "data/doubletfinder/soupx_emptydrops_sscovelli/runs/direct-production-sapphirerapids"
  t <- time_v(file.path(r, "time_v.txt"))
  s <- js(file.path(r, "summary.json"))
  if (!is.null(t) && !is.null(s)) {
    ph <- s$run$phases
    cat(sprintf("\n\nR, direct production run (one thread): wall %s, peak RSS %s; paramSweep %s, doubletFinder %s (both including SCTransform and PCA); peak after the sweep %s.\n",
                fmt_s(t$wall_s), fmt_mem(t$rss_kib), fmt_s(ph$paramSweep$elapsed_s),
                fmt_s(ph$doubletFinder$elapsed_s - ph$sweep_cleanup$elapsed_s),
                fmt_mem(ph$sweep_cleanup$peak_rss_kib)))
  }
}
