# Tables of the M8 report (docs/validation/clustering_m8.md; docs/PLAN.md C4 and C6) from the
# plain files M8 wrote under data/clustering/: the exports (input/export.json, saved/object.json),
# the reference runs (summary.json, time_v.txt), the comparisons (compare/<object>/*.json), the
# marker checks (markers/markers.json), the instrumented optimizer (optimizer/optimizer.json), and
# validation/clustering_datasets.tsv and validation/clustering_panel.tsv. Reads no .rds.
#
#   Rscript validation/R/report_clustering_m8.R [<runs.tsv>]
#
# Prints the report's tables as markdown; with an argument, also writes one row per comparison
# (object, comparison, and every summary value) as TSV (docs/validation/clustering_m8.tsv).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)
args <- commandArgs(trailingOnly = TRUE)

ds <- read.delim("validation/clustering_datasets.tsv", colClasses = "character", quote = "")
panel <- read.delim("validation/clustering_panel.tsv", colClasses = "character", quote = "")
rj <- function(p) if (file.exists(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
fmtn <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
fmtg <- function(x, d = 3) if (is.null(x) || length(x) == 0 || is.na(x)) "" else formatC(as.numeric(x), digits = d, format = "g")
parse_time_v <- function(path) {
  if (!file.exists(path)) return(list(max_rss_kib = NA, elapsed = NA, elapsed_s = NA))
  t <- readLines(path)
  rss <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", t, value = TRUE)))
  el <- sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", t, value = TRUE))
  parts <- as.numeric(strsplit(el, ":")[[1]])
  secs <- sum(parts * 60^(rev(seq_along(parts)) - 1))
  list(max_rss_kib = rss, elapsed = el, elapsed_s = secs)
}
gib <- function(kib) sprintf("%.1f", as.numeric(kib) / 1024^2)
mins <- function(s) sprintf("%.1f", as.numeric(s) / 60)
out <- character()
emit <- function(...) out <<- c(out, ...)

# Objects ------------------------------------------------------------------------------------------
rows <- lapply(seq_len(nrow(ds)), function(i) {
  o <- ds$object[i]
  ex <- rj(file.path("data/clustering", o, "input/export.json"))
  sv <- rj(file.path("data/clustering", o, "saved/object.json"))
  pr <- rj(file.path("data/clustering", o, "runs/direct-production-sapphirerapids/summary.json"))
  if (is.null(ex) || is.null(sv)) return(NULL)
  data.frame(object = o, tier = ds$tier[i], role = ds$role[i], assay = ex$assay,
             features = fmtn(ex$n_features), cells = fmtn(ex$n_cells), entries = fmtn(ex$nnz),
             explicit_zeros = fmtn(ex$explicit_zeros), sct_features = fmtn(sv$sct_n_features),
             variable_features = fmtn(sv$n_variable_features), clusters = sv$n_clusters,
             singletons = if (is.null(pr)) "" else pr$result$optimizer$singletons,
             counts = if (is.null(ex$counts_from)) "exported" else paste0("reused (", sub("^data/doubletfinder/", "", ex$counts_from), ")"),
             stringsAsFactors = FALSE)
})
emit("## Objects", "", md_table(do.call(rbind, rows)), "")

# Comparisons ------------------------------------------------------------------------------------------
cmp_rows <- list()
for (o in ds$object) {
  dir <- file.path("data/clustering/compare", o)
  if (!dir.exists(dir)) next
  for (f in list.files(dir, pattern = "\\.json$")) {
    j <- rj(file.path(dir, f))
    s <- j$summary
    flat <- lapply(s, function(v) if (is.null(v)) NA else if (length(v) > 1 || is.list(v)) paste(unlist(v), collapse = ",") else v)
    cmp_rows[[length(cmp_rows) + 1]] <- c(list(object = o, tier = ds$tier[ds$object == o], comparison = sub("\\.json$", "", f)), flat)
  }
}
all_names <- unique(unlist(lapply(cmp_rows, names)))
cmp <- do.call(rbind, lapply(cmp_rows, function(r) {
  r <- r[all_names]
  names(r) <- all_names
  r[vapply(r, is.null, logical(1))] <- NA
  as.data.frame(r, stringsAsFactors = FALSE)
}))
if (length(args) >= 1) tsv_write(cmp, args[1])
is_true <- function(v) !is.na(v) & as.character(v) %in% c("TRUE", "true")

vs_saved <- cmp[grepl("_vs_saved$", cmp$comparison), ]
vs_saved$run <- sub("_vs_saved$", "", vs_saved$comparison)
all_identical <- function(r) {
  keys <- intersect(c("labels_identical", "levels_identical", "sct_counts_identical", "sct_data_identical",
                      "sct_scale_data_identical", "sct_model_identical", "pca_identical", "nn_graph_identical",
                      "snn_graph_identical"), names(r))
  all(is_true(unlist(r[keys]))) && as.numeric(r$silhouette_differing) == 0 && as.numeric(r$purity_differing) == 0 &&
    as.numeric(r$nCount_SCT_differing) == 0 && as.numeric(r$nFeature_SCT_differing) == 0
}
vs_saved$everything <- vapply(seq_len(nrow(vs_saved)), function(i) all_identical(vs_saved[i, ]), logical(1))
agg <- do.call(rbind, lapply(split(vs_saved, list(vs_saved$tier, vs_saved$run), drop = TRUE), function(d) {
  data.frame(tier = d$tier[1], run = d$run[1], objects = nrow(d),
             `every file identical` = sum(d$everything),
             `labels and levels` = sum(is_true(d$labels_identical) & is_true(d$levels_identical)),
             silhouette = sum(as.numeric(d$silhouette_differing) == 0),
             purity = sum(as.numeric(d$purity_differing) == 0),
             `SCT counts, data, scale.data` = sum(is_true(d$sct_counts_identical) & is_true(d$sct_data_identical) & is_true(d$sct_scale_data_identical)),
             `SCT model` = sum(is_true(d$sct_model_identical)),
             PCA = sum(is_true(d$pca_identical)),
             `nn, snn` = sum(is_true(d$nn_graph_identical) & is_true(d$snn_graph_identical)),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tier != "T", agg$run), ]
emit("## Reference runs against the saved objects", "", "Objects whose run reproduces the saved object, per check:", "",
     md_table(agg), "")
prod <- grepl("production", vs_saved$run)
diff_rows <- vs_saved[!vs_saved$everything & prod, c("object", "run", "silhouette_differing", "purity_differing",
                                                     "sct_counts_identical", "sct_scale_data_identical", "pca_embeddings_max_abs_diff",
                                                     "snn_graph_identical", "labels_identical")]
emit(if (nrow(diff_rows)) c("Production-profile runs that differ from their saved object somewhere:", "", md_table(diff_rows))
     else sprintf("Every one of the %d production-profile runs reproduces its saved object in every file.", sum(prod)), "")
port <- vs_saved[!prod, ]
if (nrow(port)) {
  emit(sprintf(paste("Portable-profile runs (%d): labels and levels identical in %d, SCT layers and model in %d, both graphs in %d,",
                     "purities in %d; silhouette widths differ in the last bits in %d, the PCA in %d (see the floor)."),
               nrow(port), sum(is_true(port$labels_identical) & is_true(port$levels_identical)),
               sum(is_true(port$sct_counts_identical) & is_true(port$sct_data_identical) & is_true(port$sct_scale_data_identical) & is_true(port$sct_model_identical)),
               sum(is_true(port$nn_graph_identical) & is_true(port$snn_graph_identical)),
               sum(as.numeric(port$purity_differing) == 0), sum(as.numeric(port$silhouette_differing) > 0),
               sum(!is_true(port$pca_identical))), "")
}

# Run configurations
runs <- do.call(rbind, lapply(ds$object, function(o) {
  rd <- list.dirs(file.path("data/clustering", o, "runs"), recursive = FALSE, full.names = FALSE)
  rd <- rd[!grepl("\\.partial\\.", rd)]
  if (!length(rd)) return(NULL)
  data.frame(tier = ds$tier[ds$object == o], run = rd, stringsAsFactors = FALSE)
}))
rt <- as.data.frame.matrix(table(runs$run, runs$tier))
rt <- data.frame(run = rownames(rt), rt, check.names = FALSE, stringsAsFactors = FALSE)
emit("## Run configurations (objects per tier)", "", md_table(rt), "")

# Runs against runs
pairs <- cmp[!grepl("_vs_saved$", cmp$comparison), ]
emit("## Runs against runs", "")
pk <- function(prefix) pairs[startsWith(pairs$comparison, prefix), ]
same_all <- function(d) is_true(d$labels_identical) & is_true(d$levels_identical) & as.numeric(d$silhouette_differing) == 0 &
  as.numeric(d$purity_differing) == 0 & is_true(d$sct_counts_identical) & is_true(d$sct_data_identical) &
  is_true(d$sct_scale_data_identical) & is_true(d$sct_model_identical) & is_true(d$pca_identical) &
  is_true(d$nn_graph_identical) & is_true(d$snn_graph_identical)
rr <- list()
for (p in c("instr_vs_direct", "noumap_vs_umap", "v5_vs_v3", "cascadelake_vs_sapphirerapids")) {
  d <- pk(p)
  if (!nrow(d)) next
  for (tier in c("T", "P")) {
    e <- d[d$tier == tier, ]
    if (!nrow(e)) next
    rr[[length(rr) + 1]] <- data.frame(
      comparison = p, tier = tier, pairs = nrow(e), `every result file identical` = sum(same_all(e)),
      `summary result identical` = if (all(is.na(e$result_identical))) "" else sum(is_true(e$result_identical)),
      `stage arrays identical` = if (all(is.na(e$stage_arrays))) "" else paste0(sum(as.numeric(e$stage_arrays_identical), na.rm = TRUE), " / ", sum(as.numeric(e$stage_arrays), na.rm = TRUE)),
      `RNG checkpoints identical` = if (all(is.na(e$rng_checkpoints_identical))) "" else sum(is_true(e$rng_checkpoints_identical)),
      `Annoy index identical` = if (all(is.na(e$annoy_index_identical))) "" else sum(is_true(e$annoy_index_identical)),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
if (length(rr)) emit(md_table(do.call(rbind, rr)), "")

# The floor across kernels: stage arrays, then per object
fl <- pk("portable_vs_production")
fl <- fl[fl$comparison == "portable_vs_production", ]
if (nrow(fl)) {
  st_rows <- list()
  for (tier in c("T", "P")) {
    objs <- fl$object[fl$tier == tier]
    acc <- list()
    for (o in objs) {
      j <- rj(file.path("data/clustering/compare", o, "portable_vs_production.json"))
      for (a in names(j$stages$arrays)) {
        x <- j$stages$arrays[[a]]
        v <- acc[[a]]
        if (is.null(v)) v <- list(objects = 0, differing = 0, elements = 0, el_diff = 0, max_abs = 0, max_scaled = 0, len_diff = 0)
        v$objects <- v$objects + 1
        v$elements <- v$elements + if (is.null(x$elements)) 0 else x$elements
        if (!isTRUE(x$identical)) {
          v$differing <- v$differing + 1
          v$el_diff <- v$el_diff + if (is.null(x$elements_differing)) 0 else x$elements_differing
          if (!is.null(x$lengths)) v$len_diff <- v$len_diff + 1
          m <- if (!is.null(x$aligned_max_abs_diff)) x$aligned_max_abs_diff else x$max_abs_diff
          if (!is.null(m)) v$max_abs <- max(v$max_abs, m)
          if (!is.null(x$max_scaled_diff) && is.null(x$aligned_max_abs_diff)) v$max_scaled <- max(v$max_scaled, x$max_scaled_diff)
        }
        acc[[a]] <- v
      }
    }
    for (a in sort(names(acc))) {
      v <- acc[[a]]
      if (v$differing == 0) next
      st_rows[[length(st_rows) + 1]] <- data.frame(tier = tier, array = sub("\\.npy$", "", a), objects = v$objects,
        `objects differing` = v$differing,
        `elements differing` = if (v$len_diff > 0) sprintf("%s (lengths differ in %d)", fmtn(v$el_diff), v$len_diff) else fmtn(v$el_diff),
        `max abs diff` = fmtg(v$max_abs), check.names = FALSE, stringsAsFactors = FALSE)
    }
    same <- setdiff(sort(names(acc)), names(Filter(function(v) v$differing > 0, acc)))
    emit(sprintf("Tier %s, %d objects: identical across kernels in every object that dumps it: %s.", tier, length(objs),
                 paste(sub("\\.npy$", "", same), collapse = ", ")), "")
  }
  if (length(st_rows)) emit("Stage arrays that differ between the kernels. Embeddings and loadings are compared after aligning each component's sign; items and centroids are not, so their counts include whole components that changed sign (the floor statistics below align them):", "", md_table(do.call(rbind, st_rows)), "")
  per <- fl[, intersect(c("object", "tier", "cells", "labels_identical", "ari", "silhouette_differing", "purity_differing",
                          "sct_counts_identical", "sct_scale_data_identical", "pca_embeddings_max_abs_diff", "item_components_flipped",
                          "items_differing_aligned", "neighbour_lists_differing", "neighbour_distances_differing",
                          "purity_neighbor_counts_differing", "embeddings_1_30_max_abs_diff_aligned", "snn_graph_identical", "optimizer_ids_differing",
                          "rng_checkpoints_identical"), names(fl))]
  num <- function(v) suppressWarnings(as.numeric(v))
  tt <- per[per$tier == "T", ]
  if (nrow(tt)) {
    agg_t <- data.frame(object = sprintf("all %d tier-T objects", nrow(tt)), tier = "T", cells = sum(num(tt$cells)),
                        labels_identical = sum(is_true(tt$labels_identical)), ari = min(num(tt$ari)),
                        silhouette_differing = sum(num(tt$silhouette_differing)), purity_differing = sum(num(tt$purity_differing)),
                        sct_counts_identical = sum(is_true(tt$sct_counts_identical)),
                        sct_scale_data_identical = sum(is_true(tt$sct_scale_data_identical)),
                        pca_embeddings_max_abs_diff = max(num(tt$pca_embeddings_max_abs_diff)),
                        item_components_flipped = sum(num(tt$item_components_flipped)),
                        items_differing_aligned = sum(num(tt$items_differing_aligned)),
                        neighbour_lists_differing = sum(num(tt$neighbour_lists_differing)),
                        neighbour_distances_differing = sum(num(tt$neighbour_distances_differing)),
                        purity_neighbor_counts_differing = sum(num(tt$purity_neighbor_counts_differing)),
                        embeddings_1_30_max_abs_diff_aligned = max(num(tt$embeddings_1_30_max_abs_diff_aligned)),
                        snn_graph_identical = sum(is_true(tt$snn_graph_identical)),
                        optimizer_ids_differing = sum(num(tt$optimizer_ids_differing)),
                        rng_checkpoints_identical = sum(is_true(tt$rng_checkpoints_identical)), stringsAsFactors = FALSE)
    per <- rbind(agg_t[, names(per)], per[per$tier == "P", ])
  }
  emit("Per object (instrumented runs, portable against production, on Sapphire Rapids; for tier T, sums over the objects, and for the TRUE/FALSE columns the number of objects with TRUE):", "", md_table(per), "")
}
corr <- list()
for (o in ds$object) {
  j <- rj(file.path("data/clustering/compare", o, "portable_vs_production.json"))
  if (is.null(j) || isTRUE(j$sct$corrected_counts$identical)) next
  corr[[length(corr) + 1]] <- data.frame(object = o, `entries differing` = j$sct$corrected_counts$entries_differing,
    `same structure` = j$sct$corrected_counts$same_structure, check.names = FALSE)
}
if (length(corr)) emit("Corrected counts that differ between the kernels:", "", md_table(do.call(rbind, corr)), "")

# Memory and time ------------------------------------------------------------------------------------
mt <- list()
for (o in ds$object[ds$tier == "P"]) {
  for (run in c("direct-production-sapphirerapids", "direct-production-sapphirerapids-noumap")) {
    rd <- file.path("data/clustering", o, "runs", run)
    s <- rj(file.path(rd, "summary.json"))
    if (is.null(s)) next
    tv <- parse_time_v(file.path(rd, "time_v.txt"))
    ph <- s$run$phases
    step <- function(name, prev) if (is.null(ph[[name]])) NA else ph[[name]]$elapsed_s - (if (is.null(ph[[prev]])) 0 else ph[[prev]]$elapsed_s)
    prow <- panel[panel$object == o, ]
    mt[[length(mt) + 1]] <- data.frame(
      object = o, umap = if (grepl("noumap", run)) "no" else "yes", cells = fmtn(s$result$n_cells),
      `chain min` = mins(ph$chain_end$elapsed_s - ph$load$elapsed_s),
      `chain peak RSS GiB` = gib(ph$chain_end$peak_rss_kib),
      `script min` = mins(tv$elapsed_s), `script peak GiB` = gib(tv$max_rss_kib),
      `load min` = mins(ph$load$elapsed_s),
      `SCTransform min` = mins(step("SCTransform", "load")),
      `RunPCA min` = mins(step("RunPCA", "SCTransform")),
      `RunUMAP min` = if (is.null(ph$RunUMAP)) "" else mins(step("RunUMAP", "RunPCA")),
      `FindNeighbors min` = mins(ph$FindNeighbors$elapsed_s - (if (is.null(ph$RunUMAP)) ph$RunPCA$elapsed_s else ph$RunUMAP$elapsed_s)),
      `FindClusters min` = mins(step("FindClusters", "FindNeighbors")),
      `silhouette min` = mins(step("silhouette", "FindClusters")),
      `purity min` = mins(step("purity", "silhouette")),
      `workflow job` = prow$workflow_job, `workflow MaxRSS GiB` = gib(prow$workflow_max_rss_kib),
      `workflow min` = mins(sum(as.numeric(strsplit(prow$workflow_elapsed, ":")[[1]]) * c(3600, 60, 1))),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
if (length(mt)) emit("## Memory and time of R's chain (tier P, production profile, direct runs)", "", "chain: from the loaded object to the end of add_silhouette_to_metadata (VmHWM at its end); script: the whole Rscript under /usr/bin/time -v, which adds R's startup, reading the export and writing the outputs.", "", md_table(do.call(rbind, mt)), "")

# Markers ----------------------------------------------------------------------------------------------
mk <- list()
for (o in ds$object) {
  j <- rj(file.path("data/clustering", o, "markers", "markers.json"))
  if (is.null(j)) next
  rows_total <- sum(vapply(j$clusters, function(c) c$rows, numeric(1)))
  mk[[length(mk) + 1]] <- data.frame(object = o, tier = ds$tier[ds$object == o], clusters = length(j$cluster_levels),
                                     `saved cluster order` = paste(unlist(j$saved_cluster_order), collapse = " "),
                                     rows = fmtn(rows_total), identical = j$comparison$identical, reference = j$reference,
                                     `min` = mins(j$timing_s$total), check.names = FALSE, stringsAsFactors = FALSE)
}
if (length(mk)) {
  mk <- do.call(rbind, mk)
  emit("## Marker references", "", sprintf("Tier T: %d of %d identical. Tier P:", sum(mk$identical[mk$tier == "T"]), sum(mk$tier == "T")), "",
       md_table(mk[mk$tier == "P", ]), "")
}

# Optimizer --------------------------------------------------------------------------------------------
op <- list()
for (o in ds$object) {
  j <- rj(file.path("data/clustering", o, "optimizer", "optimizer.json"))
  if (is.null(j)) next
  st <- read.delim(file.path("data/clustering", o, "optimizer", "starts.tsv"), colClasses = "character")
  op[[length(op) + 1]] <- data.frame(
    object = o, tier = ds$tier[ds$object == o], nodes = fmtn(j$printed$nodes), edges = fmtn(j$printed$edges),
    Q = j$printed$max_modularity, communities = j$printed$communities,
    `distinct start modularities` = length(unique(st$modularity)), `best start` = j$best_start,
    `iterations` = paste(range(as.integer(st$iterations)), collapse = "-"),
    `= R's optimizer ids` = if (is.null(j$vs_instrumented_run)) "" else j$vs_instrumented_run$ids_identical,
    `= R's printed numbers` = if (is.null(j$vs_instrumented_run)) "" else j$vs_instrumented_run$printed_identical,
    `= saved labels` = j$saved_labels_identical,
    `= workflow log` = if (is.null(j$vs_workflow_log)) "" else all(unlist(j$vs_workflow_log[c("nodes", "edges", "max_modularity", "communities")])),
    check.names = FALSE, stringsAsFactors = FALSE)
}
if (length(op)) {
  op <- do.call(rbind, op)
  t_ok <- op$tier == "T" & op$`= R's optimizer ids` == "TRUE" & op$`= R's printed numbers` == "TRUE" & op$`= saved labels` == "TRUE"
  emit("## Instrumented optimizer", "",
       sprintf("Tier T: %d of %d objects reproduce R's optimizer ids, R's printed numbers and the saved labels (no tier-T object has singletons); in %d, the random starts end at more than one modularity. Tier P (= saved labels is FALSE where GroupSingletons merged singletons after the optimizer):",
               sum(t_ok), sum(op$tier == "T"), sum(op$tier == "T" & as.numeric(op$`distinct start modularities`) > 1)), "",
       md_table(op[op$tier == "P", ]), "")
}

writeLines(out)
