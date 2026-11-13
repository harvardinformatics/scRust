# The instrumented modularity optimizer on R's own SNN graph (docs/PLAN.md C4; milestone M8):
# Seurat 5.5.1's ModularityOptimizer.cpp in its standalone mode, built by validation/cpp/build.sh
# (data/clustering/bin/ModularityOptimizer), patched only to print each random start's
# iterations, modularity (%.17g) and partition. Reads no .rds: a login node is enough for tier T,
# validation/slurm/optimizer_clustering.sbatch runs the tier-P objects.
#
#   Rscript validation/R/optimizer_clustering.R --object <name> --out <dir>
#
# Input: the saved object's SCT_snn graph, as export_clustering.R wrote it
# (data/clustering/<object>/saved/graphs/snn.mtx.gz, entries in stored order, values %.17g). Its
# strictly lower triangle becomes the edge list main() reads, one "col<TAB>row<TAB>weight" line per
# entry with row > col, 0-based, in stored order (column by column, rows ascending): the order in
# which RunModularityClusteringCpp walks the matrix, so readInputFile() builds the network R
# builds (matrixToNetwork). The weights are copied as written (%.17g), which std::stod reads back
# exactly. readInputFile() counts the nodes as the largest index + 1; when the last cells have no
# edge, a line "n-1<TAB>n-1<TAB>0" is added, which sets the count and adds nothing (matrixToNetwork
# skips pairs that are not node1 < node2).
# Run: modularity function 1, resolution 0.8, algorithm 1 (Louvain), 10 starts, 10 iterations,
# seed 12345, output printed: FindClusters' settings in the workflow.
# Checks:
# - the final partition (the output file, clusters numbered by size) against R's optimizer result:
#   the production-profile instrumented run's K7_optimizer_ids (before GroupSingletons) where it
#   exists, and the saved seurat_clusters for objects without singletons (GroupSingletons then
#   changes nothing);
# - the printed nodes, edges, maximum modularity (%.4f) and communities against the R run's and,
#   for production objects, the workflow job's log (validation/clustering_panel.tsv);
# - the best start's partition, renumbered by size as orderClustersByNNodes does (stable sort by
#   decreasing size), against the final partition.
# Output in <dir>: starts.tsv (start, iterations, modularity as %.17g), partitions.npy (starts x
# nodes, each start's own cluster ids), ids.txt (the final partition), stdout.txt (what main()
# printed, without the SCRUST_START lines), edges.tsv.gz, optimizer.json.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args

args <- df_args(list(object = NULL, out = NULL))
obj_dir <- file.path(root, "data", "clustering", args$object)
snn <- file.path(obj_dir, "saved", "graphs", "snn.mtx.gz")
bin <- file.path(root, "data", "clustering", "bin", "ModularityOptimizer")
if (!file.exists(snn)) stop("no export of ", args$object)
if (!file.exists(bin)) stop("build the optimizer first: bash validation/cpp/build.sh")
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

saved_facts <- jsonlite::fromJSON(file.path(obj_dir, "saved", "object.json"), simplifyVector = FALSE)
gjson <- jsonlite::fromJSON(file.path(obj_dir, "saved", "graphs", "graphs.json"), simplifyVector = FALSE)
if (!isTRUE(gjson$snn$rows_sorted_within_columns)) stop("SNN rows are not sorted within columns")
n <- saved_facts$n_cells

# The edge list, from the Matrix Market text (no number is parsed)
edges <- file.path(args$out, "edges.tsv")
status <- system(sprintf(
  "zcat %s | awk 'BEGIN { OFS = \"\\t\"; dims = 0 } /^%%/ { next } dims == 0 { dims = 1; next } $1 > $2 { print $2 - 1, $1 - 1, $3; if ($1 - 1 > m) m = $1 - 1 } END { print m > \"/dev/stderr\" }' > %s 2> %s",
  shQuote(snn), shQuote(edges), shQuote(paste0(edges, ".max"))))
if (status != 0) stop("could not write the edge list")
max_node <- as.integer(readLines(paste0(edges, ".max")))
unlink(paste0(edges, ".max"))
padded <- max_node < n - 1
if (padded) cat(sprintf("%d\t%d\t0\n", n - 1L, n - 1L), file = edges, append = TRUE)
n_edges_written <- as.integer(system(sprintf("wc -l < %s", shQuote(edges)), intern = TRUE)) - padded

ids_path <- file.path(args$out, "ids.txt")
t1 <- proc.time()[["elapsed"]]
out <- system2(bin, c(shQuote(edges), shQuote(ids_path), "1", "0.8", "1", "10", "10", "12345", "1"),
               stdout = TRUE, stderr = TRUE)
run_s <- proc.time()[["elapsed"]] - t1
if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop("the optimizer failed")
starts <- grep("^SCRUST_START ", out, value = TRUE)
printed <- out[!startsWith(out, "SCRUST_START ")]
writeLines(printed, file.path(args$out, "stdout.txt"))
parts <- strsplit(starts, " ", fixed = TRUE)
st <- data.frame(start = as.integer(vapply(parts, `[`, "", 2)), iterations = as.integer(vapply(parts, `[`, "", 3)),
                 modularity = vapply(parts, `[`, "", 4), stringsAsFactors = FALSE)
writeLines(c("start\titerations\tmodularity", paste(st$start, st$iterations, st$modularity, sep = "\t")),
           file.path(args$out, "starts.tsv"))
P <- do.call(rbind, lapply(parts, function(p) as.integer(p[-(1:4)])))
if (ncol(P) != n) stop("partitions have ", ncol(P), " nodes, not ", n)
npy_write(file.path(args$out, "partitions.npy"), P)
ids <- as.integer(readLines(ids_path))
system(sprintf("gzip -f %s", shQuote(edges)))

one <- function(re) {
  hit <- grep(re, printed, value = TRUE)
  if (length(hit) != 1) stop("expected one line matching ", re)
  sub(re, "\\1", hit)
}
numbers <- list(nodes = as.integer(one("^Number of nodes: ([0-9]+)$")),
                edges = as.integer(one("^Number of edges: ([0-9]+)$")),
                max_modularity = one("^Maximum modularity in 10 random starts: ([0-9.]+)$"),
                communities = as.integer(one("^Number of communities: ([0-9]+)$")))

# The best start (strict >, so the first of equal starts), renumbered by decreasing size with a
# stable sort, as orderClustersByNNodes does
# Compared exactly as text: for values in [0, 1) written as 0.ddd..., right-padding the digits
# with zeros makes string order numeric order (R's own decimal parser is not correctly rounded).
mod <- st$modularity
if (!all(grepl("^0\\.[0-9]+$", mod))) stop("modularities outside [0, 1): ", paste(mod, collapse = " "))
padded_mod <- formatC(mod, width = max(nchar(mod)), flag = "-")
padded_mod <- gsub(" ", "0", padded_mod)
best <- 1L
for (k in seq_along(padded_mod)) if (padded_mod[k] > padded_mod[best]) best <- k  # strict >: first of equals
b <- P[best, ]
sizes <- tabulate(b + 1L, nbins = max(b) + 1L)
ord <- order(-sizes, seq_along(sizes))  # stable: equal sizes keep their order
ord <- ord[sizes[ord] > 0]
renum <- integer(length(sizes))
renum[ord] <- seq_along(ord) - 1L
best_renumbered <- renum[b + 1L]

res <- list(object = args$object, n_cells = n, edges_written = n_edges_written, padded_last_node = padded,
            printed = numbers, starts = nrow(st), best_start = best,
            best_start_renumbered_equals_output = identical(best_renumbered, ids),
            seconds = run_s)

# R's optimizer result
inst <- file.path(obj_dir, "runs", "instrumented-production-sapphirerapids", "stages", "K7_optimizer_ids.npy")
if (file.exists(inst)) {
  r_ids <- npy_read(inst)
  res$vs_instrumented_run <- list(run = "instrumented-production-sapphirerapids", ids_identical = identical(as.integer(r_ids), ids),
                                  ids_differing = sum(as.integer(r_ids) != ids))
  rs <- jsonlite::fromJSON(file.path(dirname(dirname(inst)), "summary.json"), simplifyVector = FALSE)$result$optimizer
  res$vs_instrumented_run$printed_identical <- identical(rs$nodes, numbers$nodes) && identical(rs$edges, numbers$edges) &&
    identical(rs$max_modularity, numbers$max_modularity) && identical(rs$communities, numbers$communities)
}
meta <- read.delim(file.path(obj_dir, "saved", "metadata.tsv"), colClasses = "character")
saved_labels <- as.integer(meta$seurat_clusters)
res$saved_labels_identical <- identical(saved_labels, ids)
res$saved_labels_differing <- sum(saved_labels != ids)
panel <- read.delim(file.path(root, "validation", "clustering_panel.tsv"), colClasses = "character")
prow <- panel[panel$object == args$object, ]
if (nrow(prow) == 1) {
  res$vs_workflow_log <- list(nodes = identical(as.character(numbers$nodes), prow$log_nodes),
                              edges = identical(as.character(numbers$edges), prow$log_edges),
                              max_modularity = identical(numbers$max_modularity, prow$log_max_modularity),
                              communities = identical(as.character(numbers$communities), prow$log_communities),
                              log_singletons = as.integer(prow$log_singletons))
  res$singletons <- as.integer(prow$log_singletons)
}
res$reproduces_R <- if (!is.null(res$vs_instrumented_run)) res$vs_instrumented_run$ids_identical else NA
json_write(res, file.path(args$out, "optimizer.json"))
cat(args$object, ": optimizer", numbers$nodes, "nodes", numbers$edges, "edges, Q", numbers$max_modularity,
    numbers$communities, "communities; vs R run", if (is.null(res$vs_instrumented_run)) "n/a" else res$vs_instrumented_run$ids_identical,
    "; vs saved labels", res$saved_labels_identical, "; best start", best, "renumbered = output",
    res$best_start_renumbered_equals_output, ";", round(proc.time()[["elapsed"]] - t0), "s\n")
