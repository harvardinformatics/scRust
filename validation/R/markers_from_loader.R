# K13 (docs/PLAN.md C3g, C4 "Markers", C5; milestone M14): the Seurat loader's acceptance test.
# The workflow's marker steps run on the object scrustio::load_seurat() builds from a scrust
# export, cluster by cluster, and the combined table must equal the saved <prefix>_markergenes.csv
# byte for byte, the clusters combined in the saved table's order. Reads no .rds (the object comes
# from the export's plain files); FindMarkers draws no random numbers with its defaults, and
# set.seed() is called before each cluster as each workflow job did. Run through sbatch
# (validation/slurm/m14_markers.sbatch).
#
#   Rscript validation/R/markers_from_loader.R --object <name> --export <scrust cluster --export output> \
#     --saved-markers <saved <prefix>_markergenes.csv> --out <dir>
#
# The steps are those of markers_clustering.R (M8), copied from the pipeline repo
# (scRNAseq-preprocessing, read at ea053fa; the scripts last changed in 8d34de1, 2026-10-30,
# before the test-reference and production runs):
# - workflow/scripts/find_markers.R, run by the rule find_markers once per cluster id: the same
#   packages and option, set_workflow_seed(), Idents <- "seurat_clusters", the workflow name taken
#   from the output file's name, the empty schema-correct table for clusters under 3 cells,
#   set.seed(12345) before FindMarkers(ident.1 = <cluster id>) with Seurat's defaults, rownames to
#   "genesymbol", p_val_adj <= 0.05, the cluster and workflow columns, write_csv(). The object is
#   built once and the clusters run in turn, as markers_clustering.R does. In the workflow the
#   object came from readRDS(<prefix>.rds); here from load_seurat(<export>), and the prefix (the
#   saved object's file name, which names the workflow column) from the saved table's file name.
# - workflow/scripts/combine_markers.R, run by the rule combine_markers: read_csv() of every
#   cluster's file, the non-empty ones bound with bind_rows() (all empty: the first file's empty
#   table), write_csv().
# The rule combine_markers receives the per-cluster files in the order Snakemake's
# glob_wildcards() lists the checkpoint's directory (os.walk: the filesystem's order, not the
# level order), so the clusters are combined in the order they appear in the saved table, then any
# others in level order (they add no rows when the tables agree; M8).
#
# Output in <dir>: clusters/<prefix>_markergenes_cluster<id>.csv, <prefix>_markergenes.csv (the
# combined table) and markers.json (the cluster order, rows and cells per cluster, timings, peak
# RSS, and the comparison with the saved table).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))     # set_workflow_seed
source(file.path(here, "lib", "loaders.R"))
root <- normalizePath(file.path(here, "..", ".."))

args <- df_args(list(object = NULL, export = NULL, `saved-markers` = NULL, out = NULL))
t0 <- proc.time()[["elapsed"]]
elapsed <- function() proc.time()[["elapsed"]] - t0
ld_setup(root)

# find_markers.R's setup
library("Seurat")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)
WORKFLOW_SEED <- set_workflow_seed()

prefix <- sub("_markergenes\\.csv$", "", basename(args$`saved-markers`))
cl_dir <- file.path(args$out, "clusters")
dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- scrustio::load_seurat(args$export)
load_s <- elapsed()
load_rss <- peak_rss_kib()
Idents(seurat_obj) <- "seurat_clusters"
cluster_ids <- levels(Idents(seurat_obj))
cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]

saved_tbl <- read_csv(args$`saved-markers`, show_col_types = FALSE, col_types = cols(.default = col_character()))
saved_order <- unique(saved_tbl$cluster)
if (!all(saved_order %in% cluster_ids)) stop("the saved table has clusters the object lacks: ", paste(setdiff(saved_order, cluster_ids), collapse = ", "))
order_used <- c(saved_order, setdiff(cluster_ids, saved_order))

per_cluster <- list()
for (cluster_id in order_used) {
  output_path <- file.path(cl_dir, paste0(prefix, "_markergenes_cluster", cluster_id, ".csv"))
  tc <- elapsed()
  # From find_markers.R ------------------------------------------------------------------------
  workflow_name <- file_path_sans_ext(basename(output_path))
  workflow_name <- sub("_markergenes_cluster.*$", "", workflow_name)

  n_cells_in_cluster <- sum(Idents(seurat_obj) == cluster_id)
  empty_markers <- tibble(
    genesymbol = character(),
    p_val = numeric(),
    avg_log2FC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric(),
    p_val_adj = numeric(),
    cluster = character(),
    workflow = character()
  )

  if (n_cells_in_cluster < 3) {
    message(sprintf(
      "Cluster %s has %d cell(s) (< 3); skipping marker detection and writing empty output.",
      cluster_id, n_cells_in_cluster
    ))
    write_csv(empty_markers, file = output_path)
    n_rows <- 0L
  } else {
    set.seed(WORKFLOW_SEED)
    markers <- FindMarkers(seurat_obj, ident.1 = cluster_id)
    markers_tbl <- markers %>%
      as.data.frame() %>%
      rownames_to_column(var = "genesymbol") %>%
      as_tibble()

    sig_markers <- markers_tbl %>%
      filter(p_val_adj <= 0.05) %>%
      mutate(cluster = cluster_id, workflow = workflow_name)

    write_csv(sig_markers, file = output_path)
    n_rows <- nrow(sig_markers)
  }
  # ---------------------------------------------------------------------------------------------
  per_cluster[[cluster_id]] <- list(cells = n_cells_in_cluster, rows = n_rows, seconds = elapsed() - tc)
  cat("cluster", cluster_id, ":", n_cells_in_cluster, "cells,", n_rows, "markers\n")
}
rm(seurat_obj)
invisible(gc())

# From combine_markers.R ----------------------------------------------------------------------------
output_path <- file.path(args$out, paste0(prefix, "_markergenes.csv"))
marker_paths <- file.path(cl_dir, paste0(prefix, "_markergenes_cluster", order_used, ".csv"))

marker_tables <- lapply(
  marker_paths,
  function(marker_path) read_csv(marker_path, show_col_types = FALSE)
)

non_empty <- marker_tables[vapply(marker_tables, nrow, integer(1)) > 0]
combined_markers <- if (length(non_empty) > 0) {
  bind_rows(non_empty)
} else {
  marker_tables[[1]]
}
write_csv(combined_markers, file = output_path)
# -------------------------------------------------------------------------------------------------

identical_bytes <- identical(sha256_file(output_path), sha256_file(args$`saved-markers`))
cmp <- list(identical = identical_bytes)
if (!identical_bytes) {
  mine <- read_csv(output_path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  cmp$rows <- c(saved = nrow(saved_tbl), loader = nrow(mine))
  cmp$same_columns <- identical(names(mine), names(saved_tbl))
  key_s <- paste(saved_tbl$cluster, saved_tbl$genesymbol)
  key_m <- paste(mine$cluster, mine$genesymbol)
  cmp$rows_only_saved <- sum(!key_s %in% key_m)
  cmp$rows_only_loader <- sum(!key_m %in% key_s)
  cmp$same_row_order <- identical(key_s, key_m)
  common <- intersect(key_s, key_m)
  for (col in setdiff(intersect(names(mine), names(saved_tbl)), c("genesymbol", "cluster"))) {
    a <- saved_tbl[[col]][match(common, key_s)]
    b <- mine[[col]][match(common, key_m)]
    cmp$columns[[col]] <- list(text_differing = sum(a != b))
  }
}
json_write(list(object = args$object, export = normalizePath(args$export),
                saved_markers = normalizePath(args$`saved-markers`), saved_markers_sha256 = sha256_file(args$`saved-markers`),
                combined_sha256 = sha256_file(output_path), prefix = prefix,
                cluster_levels = I(cluster_ids), saved_cluster_order = I(saved_order), cluster_order_used = I(order_used),
                rows = nrow(combined_markers), clusters = per_cluster, comparison = cmp,
                timing_s = list(load_seurat = load_s, total = elapsed()),
                peak_rss_kib = peak_rss_kib(), load_peak_rss_kib = load_rss, slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                versions = list(scrustio = as.character(packageVersion("scrustio")), Seurat = as.character(packageVersion("Seurat")),
                                SeuratObject = as.character(packageVersion("SeuratObject")), presto = as.character(packageVersion("presto")))),
           file.path(args$out, "markers.json"))
cat(args$object, ": combined table", if (identical_bytes) "IDENTICAL to" else "DIFFERS from", "the saved one;",
    length(order_used), "clusters,", nrow(combined_markers), "rows,", round(elapsed()), "s, peak RSS", peak_rss_kib(), "KiB\n")
if (!identical_bytes) quit(status = 1)
