# Tables of docs/validation/clustering_m9.md (milestone M9) from the files the M9 jobs wrote. Reads
# only plain files, so a login node is enough (with the hermetic flags of validation/R/README.md):
#
#   Rscript validation/R/report_clustering_m9.R docs/validation/clustering_m9.tsv > tables.md
#
# Inputs: validation/clustering_datasets.tsv; data/clustering/<object>/rust/m9/<run>/<metric>/
# (validation/slurm/m9_neighbors.sbatch: summary.json, validate.json, time_v.txt);
# data/neighbors/edge/cases.tsv and data/neighbors/edge-rust/ (m9_edge.sbatch);
# data/clustering/<object>/rust/m9-bench/ (m9_benchmark.sbatch); M8's direct runs' summary.json for
# R's FindNeighbors step in the chain. Writes one row per comparison to the TSV given, and the
# markdown tables to stdout.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
setwd(normalizePath(file.path(here, "..", "..")))
out_tsv <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_tsv)) stop("usage: report_clustering_m9.R <out.tsv>")

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
rd <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
fmtn <- function(x) formatC(x, format = "d", big.mark = ",")
md_table <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) cat("| ", paste(vapply(df[i, ], as.character, ""), collapse = " | "), " |\n", sep = "")
  cat("\n")
}
check_of <- function(v, stage, name) {
  for (c in v$checks) if (identical(c$stage, stage) && identical(c$name, name)) return(c)
  NULL
}
time_v <- function(path) {
  if (!file.exists(path)) return(list(max_rss_kib = NA, elapsed = NA))
  t <- readLines(path)
  list(max_rss_kib = as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE))),
       elapsed = sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE)))
}

# One row per comparison -----------------------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(ds))) {
  obj <- ds$object[i]
  base <- file.path("data", "clustering", obj, "rust", "m9")
  if (!dir.exists(base)) next
  for (run in sort(list.dirs(base, recursive = FALSE, full.names = FALSE))) {
    for (metric in sort(list.dirs(file.path(base, run), recursive = FALSE, full.names = FALSE))) {
      d <- file.path(base, run, metric)
      v <- rd(file.path(d, "validate.json"))
      s <- rd(file.path(d, "summary.json"))
      r <- s$result
      b <- r$annoy_build
      k3 <- check_of(v, "K3", "items"); k4 <- check_of(v, "K4", "annoy.ann")
      i5 <- check_of(v, "K5", "nn_index"); d5 <- check_of(v, "K5", "nn_dist")
      nn <- check_of(v, "K6", "SCT_nn"); snn <- check_of(v, "K6", "SCT_snn")
      pr <- check_of(v, "params", "metric, dims, n.trees, k, search.k")
      tv <- time_v(file.path(d, "time_v.txt"))
      rows[[length(rows) + 1]] <- data.frame(
        object = obj, tier = ds$tier[i], run = run, metric = metric, cells = r$n_cells, k = r$k,
        pass = isTRUE(v$pass), params = pr$status, K3_values = k3$n, K3_differ = k3$n_differ,
        K4_nodes = k4$n, K4_differ = k4$n_differ, node_bytes = r$annoy_node_bytes,
        K5_ids = i5$n, K5_ids_differ = i5$n_differ, K5_dists = d5$n, K5_dists_differ = d5$n_differ,
        K6_nn_entries = nn$n, K6_nn_differ = nn$n_differ, K6_snn_entries = snn$n, K6_snn_differ = snn$n_differ,
        split_nodes = b$split_nodes, leaves = b$leaves, attempt2 = b$accepted_at_attempt[[2]],
        attempt3 = b$accepted_at_attempt[[3]], kept_imbalanced = b$kept_imbalanced,
        random_fallbacks = b$random_fallbacks, zero_margins = b$zero_margins, draws = b$draws,
        max_depth = b$max_depth, self_first = r$knn_self_first, self_missing = r$knn_self_missing,
        snn_distinct_values = r$snn_distinct_values, threads = s$run$threads,
        build_s = s$run$build_s, search_s = s$run$search_s, graphs_s = s$run$graphs_s,
        max_rss_kib = tv$max_rss_kib, stringsAsFactors = FALSE)
    }
  }
}
tab <- do.call(rbind, rows)
tsv_write(tab, out_tsv)

cat("## Every comparison, by tier and metric\n\n")
agg <- do.call(rbind, lapply(split(tab, list(tab$tier, tab$metric), drop = TRUE), function(t) {
  data.frame(tier = t$tier[1], metric = t$metric[1], comparisons = nrow(t), identical = sum(t$pass),
             objects = length(unique(t$object)), cells = fmtn(sum(t$cells)),
             `items (K3)` = fmtn(sum(t$K3_values)), `index nodes (K4)` = fmtn(sum(t$K4_nodes)),
             `neighbour ids (K5)` = fmtn(sum(t$K5_ids)), `distances (K5)` = fmtn(sum(t$K5_dists)),
             `nn entries (K6)` = fmtn(sum(t$K6_nn_entries)), `snn entries (K6)` = fmtn(sum(t$K6_snn_entries)),
             differing = sum(t$K3_differ + t$K4_differ + t$K5_ids_differ + t$K5_dists_differ + t$K6_nn_differ + t$K6_snn_differ),
             check.names = FALSE)
}))
agg <- agg[order(agg$tier == "P", match(agg$metric, c("euclidean", "cosine", "manhattan"))), ]
md_table(agg)

cat("## Tier P, run by run\n\n")
tp <- tab[tab$tier == "P", ]
tp <- tp[order(match(tp$object, ds$object), tp$metric, tp$run), ]
md_table(data.frame(object = tp$object, run = sub("^instrumented-", "", tp$run), metric = tp$metric,
                    cells = fmtn(tp$cells), `index nodes` = fmtn(tp$K4_nodes), `node bytes` = tp$node_bytes,
                    `K3-K6` = ifelse(tp$pass, "=", "DIFFER"), `snn entries` = fmtn(tp$K6_snn_entries),
                    `self missing` = tp$self_missing, `zero margins` = fmtn(tp$zero_margins),
                    `retried` = tp$attempt2 + tp$attempt3, `kept imbalanced` = tp$kept_imbalanced,
                    `random sides` = tp$random_fallbacks, depth = tp$max_depth, check.names = FALSE))

cat("## The rare branches on real data (all comparisons)\n\n")
rb <- do.call(rbind, lapply(split(tab, list(tab$tier, tab$metric), drop = TRUE), function(t) {
  data.frame(tier = t$tier[1], metric = t$metric[1], comparisons = nrow(t),
             `split nodes` = fmtn(sum(t$split_nodes)), `zero margins` = fmtn(sum(t$zero_margins)),
             `in runs` = sum(t$zero_margins > 0), `retried splits` = fmtn(sum(t$attempt2 + t$attempt3)),
             `kept imbalanced` = fmtn(sum(t$kept_imbalanced)), `random sides` = fmtn(sum(t$random_fallbacks)),
             `cells missing from their own list` = fmtn(sum(t$self_missing)), check.names = FALSE)
}))
rb <- rb[order(rb$tier == "P", match(rb$metric, c("euclidean", "cosine", "manhattan"))), ]
md_table(rb)

# Edge cases -------------------------------------------------------------------------------------------
edge <- "data/neighbors/edge"
er <- "data/neighbors/edge-rust"
if (file.exists(file.path(er, "status.tsv"))) {
  cases <- read.delim(file.path(edge, "cases.tsv"), colClasses = "character", quote = "")
  st <- read.delim(file.path(er, "status.tsv"), colClasses = "character", quote = "")
  cat("## Constructed edge cases\n\n")
  ed <- list()
  for (j in seq_len(nrow(cases))) {
    cs <- cases$case[j]
    cell <- function(m) {
      s <- st[st$case == cs & st$metric == m, ]
      if (s$rust == "refused") return(if (s$r_graphs == "no") "refused, as R" else "REFUSED")
      sm <- rd(file.path(er, cs, m, "summary.json"))$result
      b <- sm$annoy_build
      bits <- c(if (b$zero_margins > 0) paste0("0-margins ", fmtn(b$zero_margins)),
                if (b$accepted_at_attempt[[2]] + b$accepted_at_attempt[[3]] > 0) paste0("retried ", b$accepted_at_attempt[[2]] + b$accepted_at_attempt[[3]]),
                if (b$kept_imbalanced > 0) paste0("kept ", b$kept_imbalanced),
                if (b$random_fallbacks > 0) paste0("random ", b$random_fallbacks),
                if (sm$knn_self_missing > 0) paste0("self missing ", sm$knn_self_missing))
      paste0(if (s$validate == "pass") "=" else "DIFFER", if (length(bits)) paste0(" (", paste(bits, collapse = ", "), ")"))
    }
    ed[[j]] <- data.frame(case = cs, n = cases$n[j], f = cases$f[j], k = cases$k_param[j], trees = cases$n_trees[j],
                          prune = cases$prune[j], euclidean = cell("euclidean"), cosine = cell("cosine"),
                          manhattan = cell("manhattan"), check.names = FALSE)
  }
  md_table(do.call(rbind, ed))
  cat("What each case is for:\n\n")
  for (j in seq_len(nrow(cases))) cat(sprintf("- `%s`: %s\n", cases$case[j], cases$what[j]))
  cat("\n")
}

# Benchmark --------------------------------------------------------------------------------------------
cat("## Benchmark\n\n")
bt <- list()
for (obj in c("soupx_emptydrops_sscovelli", "filtered_seurat_emptydrops_sscovelli")) {
  b <- file.path("data", "clustering", obj, "rust", "m9-bench")
  if (!dir.exists(b)) next
  rj <- rd(file.path(b, "r", "findneighbors.json"))
  tv <- time_v(file.path(b, "r", "time_v.txt"))
  chain <- vapply(c("direct-production-sapphirerapids", "direct-production-sapphirerapids-noumap"), function(run) {
    st <- rd(file.path("data", "clustering", obj, "runs", run, "summary.json"))$run$phases
    nm <- names(st)
    k <- match("FindNeighbors", nm)
    st[[k]]$elapsed_s - st[[k - 1]]$elapsed_s
  }, 0)
  bt[[length(bt) + 1]] <- data.frame(object = obj, side = "R", threads = 1, `FindNeighbors s` = sprintf("%.1f", rj$elapsed_s),
                                     `build s` = "", `search s` = "", `graphs s` = "", `whole process` = tv$elapsed,
                                     `peak RSS MiB` = fmtn(round(tv$max_rss_kib / 1024)),
                                     note = sprintf("VmHWM %s -> %s MiB across the call; in M8's chain %.1f s and %.1f s (with and without RunUMAP)",
                                                    fmtn(round(rj$vmhwm_before_kib / 1024)), fmtn(round(rj$vmhwm_after_kib / 1024)), chain[1], chain[2]),
                                     check.names = FALSE)
  for (d in sort(grep("^t[0-9]+$", list.files(b), value = TRUE))) {
    s <- rd(file.path(b, d, "summary.json"))$run
    tv <- time_v(file.path(b, d, "time_v.txt"))
    bt[[length(bt) + 1]] <- data.frame(object = obj, side = "Rust", threads = s$threads,
                                       `FindNeighbors s` = sprintf("%.1f", s$build_s + s$search_s + s$graphs_s),
                                       `build s` = sprintf("%.2f", s$build_s), `search s` = sprintf("%.2f", s$search_s),
                                       `graphs s` = sprintf("%.2f", s$graphs_s), `whole process` = tv$elapsed,
                                       `peak RSS MiB` = fmtn(round(tv$max_rss_kib / 1024)),
                                       note = "reading the embeddings and writing every output included in the process", check.names = FALSE)
  }
  bt <- c(bt)
  cat(sprintf("%s: %s\n\n", obj, paste(readLines(file.path(b, "threads.txt")), collapse = "; ")))
}
bt <- do.call(rbind, bt)
bt <- bt[order(bt$object, bt$side, bt$threads), ]
bt$threads <- as.character(bt$threads)
md_table(bt)
