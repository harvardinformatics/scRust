# DoubletFinder exactly as the workflow calls it (docs/PLAN.md B1), shared by
# reference_doubletfinder.R (direct run) and instrumented_doubletfinder.R. Requires lib/io.R.
#
# The call sequence below is copied from the pipeline's workflow/scripts/doubletfinder.R
# (seed 12345 from SCRNASEQ_PREPROCESS_SEED): the same packages loaded in the same order, the
# same option, the same calls with the same arguments, the same set.seed() points, and the same
# rm()/gc() between the sweep and the final run. What the workflow does after doubletFinder()
# (keeping singlets, Seurat re-processing, markers) is not part of the tool and is not run.

WORKFLOW_SEED <- 12345L

df_setup <- function() {
  suppressPackageStartupMessages({
    library("Seurat")
    library("DoubletFinder")
    library("igraph")
  })
  options(future.globals.maxSize = 16 * 1024^3)
}

# The Seurat object DoubletFinder reads. From an export directory (export_doubletfinder.R:
# matrix.mtx.gz, barcodes, features, seurat_clusters.tsv), the counts are rebuilt as the
# identical dgCMatrix and wrapped with CreateSeuratObject(); with sct = TRUE, DoubletFinder
# reads only the RNA counts layer, the meta.data row names and, for the homotypic adjustment,
# seurat_clusters. A path ending in .rds reads the workflow's object itself instead.
df_load_input <- function(input) {
  if (grepl("\\.rds$", input)) {
    seu <- readRDS(input)
    return(list(seurat = seu, source = "rds"))
  }
  counts <- mtx_read_10x(input)
  cl <- read.delim(file.path(input, "seurat_clusters.tsv"), colClasses = "character",
                   quote = "", comment.char = "")
  lev <- readLines(file.path(input, "seurat_clusters_levels.txt"))
  if (!identical(cl$barcode, colnames(counts))) stop("seurat_clusters.tsv barcodes differ from the matrix")
  seu <- CreateSeuratObject(counts = counts)
  if (!identical(LayerData(seu, assay = "RNA", layer = "counts"), counts)) {
    stop("CreateSeuratObject() changed the counts matrix")
  }
  if (!identical(rownames(seu@meta.data), colnames(counts))) stop("meta.data rows differ from the matrix columns")
  seu$seurat_clusters <- factor(cl$seurat_clusters, levels = lev)
  list(seurat = seu, source = "export")
}

# The workflow's calls. `fns` supplies the DoubletFinder functions: by default the ones on the
# search path, as in the workflow; the instrumented run passes the namespace versions it traces.
# Returns everything the outputs need; timings and peak RSS are taken after each phase.
df_run <- function(seurat, outdir, fns = NULL) {
  if (is.null(fns)) {
    fns <- list(paramSweep = paramSweep, summarizeSweep = summarizeSweep, find.pK = find.pK,
                modelHomotypic = modelHomotypic, doubletFinder = doubletFinder)
  }
  phase <- list()
  t0 <- proc.time()[["elapsed"]]
  mark <- function(name) {
    phase[[name]] <<- list(elapsed_s = proc.time()[["elapsed"]] - t0, peak_rss_kib = peak_rss_kib())
  }
  input_barcodes <- rownames(seurat@meta.data)

  # find.pK draws a base-graphics plot; the workflow left it in Rplots.pdf
  grDevices::pdf(file.path(outdir, "find_pK.pdf"))

  set.seed(WORKFLOW_SEED)
  sweep.res.list <- fns$paramSweep(seurat, PCs = 1:10, sct = TRUE)
  mark("paramSweep")
  sweep.stats <- fns$summarizeSweep(sweep.res.list, GT = FALSE)
  mark("summarizeSweep")
  bcmvn <- fns$find.pK(sweep.stats)
  mark("find.pK")
  grDevices::dev.off()
  optimal_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
  sweep_names <- names(sweep.res.list)
  # B3 outputs of the sweep, written before the workflow's rm()
  tsv_write(data.frame(pN = as.character(sweep.stats$pN), pK = as.character(sweep.stats$pK),
                       BCreal = sweep.stats$BCreal), file.path(outdir, "sweep_stats.tsv"))
  tsv_write(data.frame(ParamID = bcmvn$ParamID, pK = as.character(bcmvn$pK),
                       MeanBC = bcmvn$MeanBC, VarBC = bcmvn$VarBC, BCmetric = bcmvn$BCmetric),
            file.path(outdir, "bcmvn.tsv"))
  bc_max_row <- which.max(bcmvn$BCmetric)
  rm(sweep.res.list, sweep.stats, bcmvn)
  gc()
  mark("sweep_cleanup")

  # set value of nExp
  nExp_poi <- round(0.15*nrow(seurat@meta.data))
  homotypic.prop <- fns$modelHomotypic(seurat$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  set.seed(WORKFLOW_SEED)
  seurat <- fns$doubletFinder(seurat, PCs = 1:10, pN = 0.25,
                              pK = optimal_pk, nExp = nExp_poi.adj,
                              reuse.pANN = NULL,
                              sct=TRUE)
  mark("doubletFinder")
  rng_end <- rng_state()

  classifier_colname <- paste("DF.classifications_0.25",optimal_pk,nExp_poi.adj,sep="_")
  meta <- seurat@meta.data
  pann_col <- grep("^pANN_", colnames(meta), value = TRUE)
  class_col <- grep("^DF\\.classifications_", colnames(meta), value = TRUE)
  if (length(pann_col) != 1 || length(class_col) != 1) stop("expected one pANN and one classification column")
  if (!identical(class_col, classifier_colname)) stop("classification column ", class_col, " != ", classifier_colname)
  if (!identical(rownames(meta), input_barcodes)) stop("cell order changed")
  n <- nrow(meta)
  n_doublets <- round(n/(1 - 0.25) - n)
  list(
    cells = data.frame(barcode = rownames(meta), pANN = meta[[pann_col]],
                       classification = meta[[class_col]], stringsAsFactors = FALSE),
    result = list(
      n_cells = n,
      n_genes = nrow(LayerData(seurat, assay = "RNA", layer = "counts")),
      seed = WORKFLOW_SEED,
      PCs = 10L,
      pN = 0.25,
      sweep_pK_labels = unique(sub("^pN_[^_]+_pK_", "", sweep_names)),
      sweep_pN_labels = unique(sub("^pN_([^_]+)_pK_.*$", "\\1", sweep_names)),
      bcmvn_max_row = bc_max_row,
      pK = optimal_pk,
      pK_label = as.character(optimal_pk),
      nExp_poi = nExp_poi,
      homotypic_prop = homotypic.prop,
      nExp = nExp_poi.adj,
      n_doublets = n_doublets,
      n_cells_with_doublets = n + n_doublets,
      k = round((n + n_doublets) * optimal_pk),
      pANN_column = pann_col,
      classification_column = class_col,
      n_doublet_calls = sum(meta[[class_col]] == "Doublet"),
      rng_end_sha256 = digest::digest(paste(rng_end, collapse = ","), algo = "sha256", serialize = FALSE)
    ),
    phases = phase
  )
}

# Writes cells.tsv and summary.json (plan B3). summary.json has a deterministic "result" section
# (equal between runs that agree) and a "run" section (timings, memory, node, fingerprint).
df_write_outputs <- function(res, outdir, run_info) {
  tsv_write(res$cells, file.path(outdir, "cells.tsv"))
  json_write(list(result = res$result, run = c(run_info, list(phases = res$phases,
                                                              peak_rss_kib = peak_rss_kib()))),
             file.path(outdir, "summary.json"))
}

# Command-line arguments --name value
df_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0) stop("arguments must be --name value pairs")
  out <- defaults
  for (i in seq(1, length(a), by = 2)) {
    key <- sub("^--", "", a[i])
    if (!key %in% names(defaults)) stop("unknown argument ", a[i])
    out[[key]] <- a[i + 1]
  }
  missing <- names(out)[vapply(out, is.null, logical(1))]
  if (length(missing)) stop("missing arguments: ", paste0("--", missing, collapse = " "))
  out
}

df_run_info <- function(args, kind, input_source) {
  list(
    kind = kind,
    input = args$input,
    input_source = input_source,
    outdir = args$out,
    script = kind,
    start_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
    profile = Sys.getenv("SCRUST_PROFILE"),
    fingerprint = fingerprint()
  )
}
