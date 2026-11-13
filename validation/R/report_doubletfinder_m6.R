# Markdown tables for docs/validation/doubletfinder.md (milestone M6) from the M6 outputs:
# - the validation (validation/slurm/m6_doubletfinder.sbatch), data/doubletfinder/<object>/rust/
#   <label>/<profile>/: vs_<R run>.json (scrust validate doubletfinder against each reference
#   run), vs_saved.json (the saved workflow output; from M6b on, singlets/ against the saved
#   object's RNA counts too), run/summary.json, time_v.txt;
# - the benchmark (validation/slurm/m6_benchmark.sbatch), data/doubletfinder/<object>/rust/
#   <bench label>/t<threads>/, and the direct R runs it compares with (runs/direct-production-
#   sapphirerapids: time_v.txt, summary.json).
# Reads only plain files.
#
#   Rscript validation/R/report_doubletfinder_m6.R <label> <bench label> > <tables.md>

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
setwd(normalizePath(file.path(here, "..", "..")))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: report_doubletfinder_m6.R <label> <bench label>")
label <- args[[1]]
bench <- args[[2]]

objects <- c("soupx_tenx_test", "soupx_emptydrops_test", "decontx_raw_tenx_test", "decontx_raw_emptydrops_test",
             "cellbender_fromraw_test", "soupx_tenx_SRR18256555", "soupx_tenx_10xmousebrain",
             "soupx_emptydrops_4500871lane", "soupx_emptydrops_sscovelli", "cellbender_fromraw_ghirsutum2")
bench_objects <- c("soupx_emptydrops_sscovelli", "cellbender_fromraw_ghirsutum2")
js <- function(p) if (file.exists(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
out <- function(lines) cat(lines, sep = "\n")
num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
time_v <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p)
  g <- function(pat) sub(".*: ", "", grep(pat, l, value = TRUE, fixed = TRUE))
  pt <- as.numeric(strsplit(g("Elapsed (wall clock) time"), ":")[[1]])
  secs <- if (length(pt) == 3) pt[1] * 3600 + pt[2] * 60 + pt[3] else pt[1] * 60 + pt[2]
  list(rss_kib = as.numeric(g("Maximum resident set size")), wall_s = secs,
       cpu_s = as.numeric(g("User time (seconds)")) + as.numeric(g("System time (seconds)")))
}
fmt_s <- function(s) {
  if (is.null(s) || is.na(s)) return("")
  if (s >= 3600) sprintf("%d:%02d:%04.1f", as.integer(s %/% 3600), as.integer((s %% 3600) %/% 60), s %% 60)
  else if (s >= 120) sprintf("%d:%04.1f", as.integer(s %/% 60), s %% 60) else sprintf("%.1f s", s)
}
fmt_gib <- function(kib) sprintf("%.2f", kib / 1024^2)
cnt <- function(ok, tot) if (tot == 0) "-" else if (ok == tot) sprintf("= (%s)", num(tot)) else sprintf("**%s of %s**", num(ok), num(tot))

# Check categories ---------------------------------------------------------------------------
kind <- function(c) {
  nm <- c$name
  if (c$stage == "result") return("result")
  if (grepl("\\.tsv$", nm)) return("file")
  if (grepl(": RNG ", nm, fixed = TRUE) || nm == "checkpoint sequence") return("rng")
  if (grepl("^[^:]+: S[0-9]", nm) || grepl("^S[0-9]", nm) || grepl("teacher-forced|from R's RNG", nm)) return("array")
  if (c$stage %in% c("stages", "saved")) return(c$stage)
  "scalar"
}
summarise <- function(v) {
  ch <- v$checks
  k <- vapply(ch, kind, character(1))
  st <- vapply(ch, function(c) c$status, character(1))
  stage <- vapply(ch, function(c) c$stage, character(1))
  n <- vapply(ch, function(c) as.numeric(c$n), numeric(1))
  nd <- vapply(ch, function(c) as.numeric(c$n_differ), numeric(1))
  cmp <- st %in% c("identical", "differ")
  list(k = k, st = st, stage = stage, n = n, nd = nd, cmp = cmp,
       checks_ok = sum(st == "identical"), checks = sum(cmp), values = sum(n[cmp]), differ = sum(nd[cmp]),
       rng_ok = sum(k == "rng" & st == "identical" & grepl(": RNG ", vapply(ch, function(c) c$name, ""), fixed = TRUE)),
       rng = sum(k == "rng" & cmp & grepl(": RNG ", vapply(ch, function(c) c$name, ""), fixed = TRUE)),
       na = sum(st == "n/a"), missing = sum(st == "missing"), pass = isTRUE(v$pass))
}
stage_cell <- function(s, stages) {
  sel <- s$cmp & s$k %in% c("array", "scalar") & s$stage %in% stages
  cnt(sum(s$st[sel] == "identical"), sum(sel))
}

profiles_of <- function(obj) {
  d <- file.path("data/doubletfinder", obj, "rust", label)
  p <- list.files(d)
  p[p %in% c("portable", "production")]
}
runs_of <- function(obj, profile, mode) {
  r <- list.files(file.path("data/doubletfinder", obj, "runs"))
  r <- r[!grepl("partial", r)]
  r[startsWith(r, paste0(mode, "-", profile, "-"))]
}

# 1. Instrumented runs -------------------------------------------------------------------------
tot <- list()
add_tot <- function(what, s) {
  t <- tot[[what]]
  if (is.null(t)) t <- c(comparisons = 0, passed = 0, checks_ok = 0, checks = 0, values = 0, differ = 0, rng_ok = 0, rng = 0, na = 0, missing = 0)
  t <- t + c(1, s$pass, s$checks_ok, s$checks, s$values, s$differ, s$rng_ok, s$rng, s$na, s$missing)
  tot[[what]] <<- t
}
rows <- list()
for (obj in objects) for (p in profiles_of(obj)) {
  d <- file.path("data/doubletfinder", obj, "rust", label, p)
  lvl <- js(file.path(d, "run", "stages", "stages.json"))$level
  for (r in runs_of(obj, p, "instrumented")) {
    v <- js(file.path(d, paste0("vs_", r, ".json")))
    if (is.null(v)) next
    s <- summarise(v)
    add_tot("instrumented", s)
    rows[[length(rows) + 1]] <- data.frame(
      object = obj, `R run` = r, dumps = lvl,
      S1 = stage_cell(s, "S1"), S2a = stage_cell(s, "S2a"), S2b = stage_cell(s, "S2b"), S2c = stage_cell(s, "S2c"),
      S2d = stage_cell(s, "S2d"), S3 = stage_cell(s, "S3"), S4 = stage_cell(s, "S4"), S5 = stage_cell(s, "S5"),
      S6 = stage_cell(s, "S6"), S7 = stage_cell(s, "S7"),
      `RNG checkpoints` = cnt(s$rng_ok, s$rng),
      files = cnt(sum(s$k == "file" & s$st == "identical"), sum(s$k == "file" & s$cmp)),
      result = cnt(sum(s$k == "result" & s$st == "identical"), sum(s$k == "result" & s$cmp)),
      `values compared` = num(s$values), differing = num(s$differ), `n/a` = s$na, missing = s$missing,
      pass = if (s$pass && s$missing == 0) "yes" else "**no**", check.names = FALSE)
  }
}
out(c("### Instrumented runs", ""))
out(md_table(do.call(rbind, rows)))

# 2. Direct runs ---------------------------------------------------------------------------------
rows <- list()
for (obj in objects) for (p in profiles_of(obj)) {
  d <- file.path("data/doubletfinder", obj, "rust", label, p)
  for (r in runs_of(obj, p, "direct")) {
    v <- js(file.path(d, paste0("vs_", r, ".json")))
    if (is.null(v)) next
    s <- summarise(v)
    add_tot("direct", s)
    f <- function(file) {
      c <- Filter(function(c) c$name == file, v$checks)
      if (!length(c)) "-" else if (c[[1]]$status == "identical") "=" else "**differs**"
    }
    rows[[length(rows) + 1]] <- data.frame(
      object = obj, `R run` = r, `Rust run` = p, `cells.tsv` = f("cells.tsv"), `sweep_stats.tsv` = f("sweep_stats.tsv"),
      `bcmvn.tsv` = f("bcmvn.tsv"), `result fields` = cnt(sum(s$k == "result" & s$st == "identical"), sum(s$k == "result" & s$cmp)),
      pass = if (s$pass) "yes" else "**no**", check.names = FALSE)
  }
}
out(c("", "### Direct runs", ""))
out(md_table(do.call(rbind, rows)))

# 3. Saved workflow outputs ---------------------------------------------------------------------
rows <- list()
for (obj in objects) for (p in profiles_of(obj)) {
  d <- file.path("data/doubletfinder", obj, "rust", label, p)
  v <- js(file.path(d, "vs_saved.json"))
  if (is.null(v)) next
  s <- summarise(v)
  add_tot("saved", s)
  g <- function(name) Filter(function(c) c$name == name, v$checks)[[1]]
  cell <- function(name) { c <- g(name); if (c$status == "identical") sprintf("= (%s)", num(c$n)) else sprintf("**%s of %s differ**", num(c$n_differ), num(c$n)) }
  res <- js(file.path(d, "run", "summary.json"))$result
  # M6b: singlets/ against the saved object's RNA counts (barcodes, feature lines, matrix)
  sgl <- Filter(function(c) startsWith(c$name, "singlets "), v$checks)
  singlets <- if (length(sgl) == 0) "-" else if (all(vapply(sgl, function(c) c$status == "identical", logical(1)))) {
    sprintf("= (%s entries)", num(g("singlets matrix")$n))
  } else "**differ**"
  rows[[length(rows) + 1]] <- data.frame(
    object = obj, `Rust run` = p, `doublet set` = cell("doublet set"), `singlet pANN` = cell("singlet pANN"),
    pK = if (g("pK_label")$status == "identical") res$pK_label else "**differs**",
    nExp = if (g("nExp")$status == "identical") num(res$nExp) else "**differs**",
    columns = if (g("pANN_column")$status == "identical" && g("classification_column")$status == "identical") sprintf("`%s`", res$pANN_column) else "**differ**",
    `singlets' counts` = singlets,
    pass = if (s$pass) "yes" else "**no**", check.names = FALSE)
}
out(c("", "### Saved workflow outputs", ""))
out(md_table(do.call(rbind, rows)))

# 4. Totals --------------------------------------------------------------------------------------
rows <- lapply(names(tot), function(w) {
  t <- tot[[w]]
  data.frame(comparisons = w, runs = num(t[["comparisons"]]), passed = num(t[["passed"]]),
             `checks identical` = sprintf("%s of %s", num(t[["checks_ok"]]), num(t[["checks"]])),
             `values compared` = num(t[["values"]]), `values differing` = num(t[["differ"]]),
             `RNG checkpoints identical` = if (t[["rng"]] > 0) sprintf("%s of %s", num(t[["rng_ok"]]), num(t[["rng"]])) else "-",
             `not applicable` = num(t[["na"]]), missing = num(t[["missing"]]), check.names = FALSE)
})
out(c("", "### Totals", ""))
out(md_table(do.call(rbind, rows)))

# 5. The validation runs' own cost -------------------------------------------------------------------
rows <- list()
for (obj in objects) for (p in profiles_of(obj)) {
  d <- file.path("data/doubletfinder", obj, "rust", label, p)
  tv <- time_v(file.path(d, "time_v.txt"))
  sm <- js(file.path(d, "run", "summary.json"))
  if (is.null(tv) || is.null(sm)) next
  rows[[length(rows) + 1]] <- data.frame(
    object = obj, profile = p, cells = num(sm$result$n_cells), `with doublets` = num(sm$result$n_cells_with_doublets),
    k = num(sm$result$k), dumps = sm$run$dump_stages, threads = sm$run$threads, wall = fmt_s(tv$wall_s),
    `peak RSS GiB` = fmt_gib(tv$rss_kib), check.names = FALSE)
}
out(c("", "### Validation runs (scrust doubletfinder with stage dumps)", ""))
out(md_table(do.call(rbind, rows)))

# 6. Benchmark ---------------------------------------------------------------------------------------
phase <- function(ph, name) if (is.null(ph[[name]])) NA else ph[[name]]$elapsed_s
hwm <- function(ph, name) if (is.null(ph[[name]])) NA else ph[[name]]$peak_rss_kib
rows <- list()
mem <- list()
for (obj in bench_objects) {
  ref <- file.path("data/doubletfinder", obj, "runs", "direct-production-sapphirerapids")
  tv <- time_v(file.path(ref, "time_v.txt"))
  sm <- js(file.path(ref, "summary.json"))
  ph <- sm$run$phases
  load_s <- sm$run$load_input_s
  sweep_s <- phase(ph, "paramSweep")
  rows[[length(rows) + 1]] <- data.frame(
    object = obj, side = "R (direct run, M2)", threads = 1, wall = fmt_s(tv$wall_s), `CPU time` = fmt_s(tv$cpu_s),
    `peak RSS GiB` = fmt_gib(tv$rss_kib), load = fmt_s(load_s), `paramSweep` = fmt_s(sweep_s),
    `summary, pK` = fmt_s(phase(ph, "find.pK") - sweep_s),
    `doubletFinder` = fmt_s(phase(ph, "doubletFinder") - phase(ph, "sweep_cleanup")),
    `vs R` = "", check.names = FALSE)
  for (t in c(1, 32)) {
    d <- file.path("data/doubletfinder", obj, "rust", bench, paste0("t", t))
    tvr <- time_v(file.path(d, "time_v.txt"))
    smr <- js(file.path(d, "run", "summary.json"))
    v <- js(file.path(d, "vs_direct.json"))
    if (is.null(tvr) || is.null(smr)) next
    p <- smr$run$phases
    s <- summarise(v)
    rows[[length(rows) + 1]] <- data.frame(
      object = obj, side = "Rust", threads = t, wall = fmt_s(tvr$wall_s), `CPU time` = fmt_s(tvr$cpu_s),
      `peak RSS GiB` = fmt_gib(tvr$rss_kib), load = fmt_s(phase(p, "load")),
      `paramSweep` = fmt_s(phase(p, "paramSweep") - phase(p, "load")),
      `summary, pK` = fmt_s(phase(p, "find.pK") - phase(p, "paramSweep")),
      `doubletFinder` = fmt_s(phase(p, "doubletFinder") - phase(p, "find.pK")),
      `vs R` = if (s$pass) sprintf("identical: files %d/3, result %d/20", sum(s$k == "file" & s$st == "identical"), sum(s$k == "result" & s$st == "identical")) else "**differs**",
      check.names = FALSE)
    fin <- function(n) phase(p, paste0("final: ", n))
    mem[[length(mem) + 1]] <- data.frame(
      object = obj, threads = t, `after load` = fmt_gib(hwm(p, "load")), `after the sweep` = fmt_gib(hwm(p, "find.pK")),
      `final SCTransform` = sprintf("%s (%s GiB)", fmt_s(fin("SCTransform") - phase(p, "find.pK")), fmt_gib(hwm(p, "final: SCTransform"))),
      `final RunPCA` = sprintf("%s (%s GiB)", fmt_s(fin("RunPCA") - fin("SCTransform")), fmt_gib(hwm(p, "final: RunPCA"))),
      `final kNN, pANN, calls` = sprintf("%s (%s GiB)", fmt_s(phase(p, "doubletFinder") - fin("RunPCA")), fmt_gib(hwm(p, "doubletFinder"))),
      `write` = fmt_s(phase(p, "write") - phase(p, "doubletFinder")), check.names = FALSE)
  }
}
out(c("", "### Benchmark", ""))
out(md_table(do.call(rbind, rows)))
out(c("", "### Rust: time and peak RSS (high-water mark) by phase", ""))
out(md_table(do.call(rbind, mem)))
