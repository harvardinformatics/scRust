# Tables of docs/validation/clustering_m10.md (milestone M10) from the files the M10 jobs wrote. Reads
# only plain files (run it through srun, with the hermetic flags of validation/R/README.md):
#
#   Rscript validation/R/report_clustering_m10.R docs/validation/clustering_m10.tsv > tables.md
#
# Inputs: validation/clustering_datasets.tsv; data/clustering/<object>/rust/m10/<run>/ and
# rust/m10/refs/<setting>/ (validation/slurm/m10_cluster.sbatch: summary.json, validate.json,
# time_v.txt, or error_check.txt where R stopped); data/clustering/<object>/cluster_refs/<setting>/
# (cluster_references.sbatch settings: R's summary.json); data/clustering/edge/cases.tsv;
# data/clustering/<object>/rust/m10/bench/ (m10_benchmark.sbatch); M8's direct runs' summary.json
# for R's FindClusters step in the chain. Writes one row per comparison to the TSV given, and the
# markdown tables to stdout.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
setwd(normalizePath(file.path(here, "..", "..")))
out_tsv <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_tsv)) stop("usage: report_clustering_m10.R <out.tsv>")

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
rd <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
fmtn <- function(x) formatC(x, format = "d", big.mark = ",")
md_table <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) cat("| ", paste(vapply(df[i, ], as.character, ""), collapse = " | "), " |\n", sep = "")
  cat("\n")
}
time_v <- function(path) {
  if (!file.exists(path)) return(list(max_rss_kib = NA, elapsed = NA))
  t <- readLines(path)
  list(max_rss_kib = as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE))),
       elapsed = sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE)))
}
checks_like <- function(v, pattern) Filter(function(c) grepl(pattern, c$name), v$checks)
sum_n <- function(cs) sum(vapply(cs, function(c) as.numeric(c$n), 0))
sum_d <- function(cs) sum(vapply(cs, function(c) as.numeric(c$n_differ), 0))

# One row per comparison -----------------------------------------------------------------------------
rows <- list()
add_row <- function(obj, tier, kind, name, d, rref) {
  err_check <- file.path(d, "error_check.txt")
  if (file.exists(err_check)) {
    txt <- readLines(err_check)
    rows[[length(rows) + 1]] <<- data.frame(object = obj, tier = tier, kind = kind, run = name, cells = NA, resolutions = NA,
      modularity_fxn = NA, algorithm = NA, n_start = NA, n_iter = NA, seed = NA, group_singletons = NA,
      pass = startsWith(txt[1], "identical"), checks = 1, checks_differ = as.integer(!startsWith(txt[1], "identical")),
      starts = 0, starts_differ = 0, partition_values = 0, partition_differ = 0, ids = 0, ids_differ = 0, labels = 0, labels_differ = 0,
      singletons = NA, rng_compared = 0, rng_differ = 0, draws = NA, pow2_draws = NA, communities = NA, clusters = NA,
      max_modularity = NA, best_start = NA, distinct_modularities = NA, log = NA,
      rust_s = NA, r_s = if (!is.null(rref)) rref$run$seconds else NA, max_rss_kib = NA,
      r_error = sub('^identical: R and Rust stop with "(.*)" \\(exit.*$', "\\1", txt[1]), stringsAsFactors = FALSE)
    return(invisible())
  }
  v <- rd(file.path(d, "validate.json"))
  s <- rd(file.path(d, "summary.json"))
  p <- s$run$parameters
  res <- s$result$resolutions
  st <- checks_like(v, "every start \\(iterations")
  pa <- checks_like(v, "every start's partition")
  ids <- checks_like(v, "optimizer ids$")
  lab <- checks_like(v, ": labels$")
  rng <- checks_like(v, "Random.seed")
  rng_drawn <- Filter(function(c) startsWith(c$detail, "Rust drew;"), rng)
  log <- checks_like(v, "logged numbers")
  q <- vapply(res, function(r) as.character(r$max_modularity_printed), "")
  dm <- vapply(res, function(r) length(unique(vapply(r$start_modularity, fmt17, ""))), 0)
  tv <- time_v(file.path(d, "time_v.txt"))
  rows[[length(rows) + 1]] <<- data.frame(object = obj, tier = tier, kind = kind, run = name, cells = s$result$n_cells,
    resolutions = length(res), modularity_fxn = p$modularity_fxn, algorithm = p$algorithm, n_start = p$n_start,
    n_iter = p$n_iter, seed = p$random_seed, group_singletons = p$group_singletons,
    pass = isTRUE(v$pass), checks = length(v$checks), checks_differ = sum(vapply(v$checks, function(c) c$status %in% c("differ", "missing"), NA)),
    starts = sum_n(st), starts_differ = sum_d(st), partition_values = sum_n(pa), partition_differ = sum_d(pa),
    ids = sum_n(ids), ids_differ = sum_d(ids), labels = sum_n(lab), labels_differ = sum_d(lab),
    singletons = sum(vapply(res, function(r) as.numeric(r$n_singletons), 0)),
    rng_compared = length(rng_drawn), rng_differ = sum_d(rng),
    draws = sum(vapply(res, function(r) as.numeric(r$next_int_draws), 0)),
    pow2_draws = sum(vapply(res, function(r) as.numeric(r$next_int_power_of_two_draws), 0)),
    communities = paste(vapply(res, function(r) as.character(r$communities), ""), collapse = ","),
    clusters = paste(vapply(res, function(r) as.character(length(r$levels)), ""), collapse = ","),
    max_modularity = paste(q, collapse = ","), best_start = paste(vapply(res, function(r) as.character(r$best_start), ""), collapse = ","),
    distinct_modularities = paste(dm, collapse = ","), log = if (length(log)) log[[1]]$status else NA,
    rust_s = s$run$resolutions_s, r_s = if (!is.null(rref)) rref$run$seconds else NA, max_rss_kib = tv$max_rss_kib,
    r_error = "", stringsAsFactors = FALSE)
}
for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  base <- file.path("data", "clustering", obj, "rust", "m10")
  if (!dir.exists(base)) next
  for (run in sort(list.dirs(base, recursive = FALSE, full.names = FALSE))) {
    if (!startsWith(run, "instrumented-")) next
    add_row(obj, ds$tier[i], "instrumented", run, file.path(base, run), NULL)
  }
  refs <- file.path(base, "refs")
  if (dir.exists(refs)) {
    for (s in sort(list.dirs(refs, recursive = FALSE, full.names = FALSE))) {
      if (grepl("\\.t1$", s)) next
      rref <- rd(file.path("data", "clustering", obj, "cluster_refs", s, "summary.json"))
      add_row(obj, ds$tier[i], "reference", s, file.path(refs, s), rref)
    }
  }
}
tab <- do.call(rbind, rows)
tsv_write(tab, out_tsv)

# Summaries --------------------------------------------------------------------------------------
cat("## Every comparison, by tier and kind\n\n")
agg <- do.call(rbind, lapply(split(tab, list(tab$tier, tab$kind), drop = TRUE), function(t) {
  data.frame(tier = t$tier[1], kind = t$kind[1], comparisons = nrow(t), identical = sum(t$pass),
             objects = length(unique(t$object)), resolutions = fmtn(sum(t$resolutions, na.rm = TRUE)),
             starts = fmtn(sum(t$starts)), partition_values = fmtn(sum(t$partition_values)),
             optimizer_ids = fmtn(sum(t$ids)), labels = fmtn(sum(t$labels)),
             singletons = fmtn(sum(t$singletons, na.rm = TRUE)), rng_states = sum(t$rng_compared),
             R_stopped = sum(nzchar(t$r_error)), differing = sum(t$checks_differ))
}))
md_table(agg[order(agg$tier, agg$kind), ])

cat("## Tier P, instrumented runs\n\n")
p <- tab[tab$tier == "P" & tab$kind == "instrumented", ]
md_table(data.frame(object = p$object, run = sub("^instrumented-", "", p$run), cells = fmtn(p$cells),
                    Q = p$max_modularity, communities = p$communities, singletons = p$singletons, clusters = p$clusters,
                    best_start = p$best_start, distinct_modularities = p$distinct_modularities,
                    K7 = ifelse(p$pass, "=", "DIFFER"), log = ifelse(p$log == "identical", "=", p$log), rust_s = sprintf("%.1f", p$rust_s)))

cat("## References, by setting\n\n")
r <- tab[tab$kind == "reference", ]
bys <- do.call(rbind, lapply(split(r, r$run), function(t) {
  data.frame(setting = t$run[1], comparisons = nrow(t), identical = sum(t$pass), T = sum(t$tier == "T"), P = sum(t$tier == "P"),
             resolutions = sum(t$resolutions, na.rm = TRUE), singletons = fmtn(sum(t$singletons, na.rm = TRUE)),
             rng_states = sum(t$rng_compared), R_stopped = sum(nzchar(t$r_error)),
             starts = fmtn(sum(t$starts)), rust_s = sprintf("%.1f", sum(t$rust_s, na.rm = TRUE)), R_s = sprintf("%.1f", sum(t$r_s, na.rm = TRUE)))
}))
md_table(bys)

cat("## Tier P references\n\n")
rp <- r[r$tier == "P", ]
cell <- function(t) ifelse(nzchar(t$r_error), "stops", paste0(ifelse(t$pass, "=", "DIFFER"), " (", t$clusters, "; ", t$singletons, ")"))
objs <- unique(rp$object)
sets <- unique(rp$run)
m <- data.frame(setting = sets, stringsAsFactors = FALSE)
for (o in objs) m[[o]] <- vapply(sets, function(s) { t <- rp[rp$object == o & rp$run == s, ]; if (nrow(t)) cell(t) else "" }, "")
md_table(m)

cat("## nextInt's draws\n\n")
nx <- tab[!is.na(tab$draws), ]
md_table(do.call(rbind, lapply(split(nx, nx$tier), function(t) data.frame(tier = t$tier[1], comparisons = nrow(t),
  draws = fmtn(sum(t$draws)), power_of_two = fmtn(sum(t$pow2_draws)),
  runs_with_power_of_two = sum(t$pow2_draws > 0)))))

# Benchmark ------------------------------------------------------------------------------------
cat("## Benchmark\n\n")
brow <- list()
for (obj in c("soupx_emptydrops_sscovelli", "filtered_seurat_emptydrops_sscovelli")) {
  b <- file.path("data", "clustering", obj, "rust", "m10", "bench")
  if (!dir.exists(b)) next
  rj <- rd(file.path(b, "r.json"))
  rt <- time_v(file.path(b, "r_time_v.txt"))
  brow[[length(brow) + 1]] <- data.frame(object = obj, side = "R", run = "FindClusters, workflow (1 resolution)", threads = 1,
    clustering_s = sprintf("%.1f", rj$findclusters_s), process = rt$elapsed, peak_rss_mib = fmtn(round(rt$max_rss_kib / 1024)),
    note = sprintf("VmHWM %s -> %s MiB across the call", fmtn(round(rj$vmhwm_before_kib / 1024)), fmtn(round(rj$vmhwm_after_kib / 1024))))
  rres <- rd(file.path("data", "clustering", obj, "cluster_refs", "resolutions", "summary.json"))
  brow[[length(brow) + 1]] <- data.frame(object = obj, side = "R", run = "FindClusters, 5 resolutions", threads = 1,
    clustering_s = sprintf("%.1f", rres$run$seconds), process = "", peak_rss_mib = "", note = "the resolutions reference (cluster_references.sbatch, another node)")
  for (r in c("workflow_t1", "resolutions_t1", grep("^resolutions_t[0-9]+$", list.files(b), value = TRUE))) {
    if (r %in% vapply(brow, function(x) paste(x$object, x$run), "")) next
    d <- file.path(b, r)
    if (!dir.exists(d)) next
    s <- rd(file.path(d, "summary.json"))
    tv <- time_v(file.path(b, paste0(r, ".time_v.txt")))
    brow[[length(brow) + 1]] <- data.frame(object = obj, side = "Rust", run = if (startsWith(r, "workflow")) "scrust cluster, workflow (1 resolution)" else "scrust cluster, 5 resolutions",
      threads = s$run$threads, clustering_s = sprintf("%.1f", s$run$resolutions_s), process = tv$elapsed,
      peak_rss_mib = fmtn(round(tv$max_rss_kib / 1024)), note = "reading the graph and writing every output included")
  }
}
bt <- unique(do.call(rbind, brow))
md_table(bt)
m8 <- list()
for (obj in c("soupx_emptydrops_sscovelli", "filtered_seurat_emptydrops_sscovelli")) {
  for (run in c("direct-production-sapphirerapids", "direct-production-sapphirerapids-noumap")) {
    ph <- rd(file.path("data", "clustering", obj, "runs", run, "summary.json"))$run$phases
    m8[[length(m8) + 1]] <- data.frame(object = obj, run = run, findclusters_s = sprintf("%.1f", ph$FindClusters$elapsed_s - ph$FindNeighbors$elapsed_s))
  }
}
cat("R's FindClusters step in M8's chain (FindClusters' cumulative time minus FindNeighbors'):\n\n")
md_table(do.call(rbind, m8))
