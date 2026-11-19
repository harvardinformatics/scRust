set_workflow_seed <- function(seed = Sys.getenv("SCRNASEQ_PREPROCESS_SEED", "12345")) {
  seed <- suppressWarnings(as.integer(seed))
  if (length(seed) != 1 || is.na(seed)) {
    stop("SCRNASEQ_PREPROCESS_SEED must be an integer")
  }
  set.seed(seed)
  seed
}

# Guard against too-few-cells before RunPCA. PCA (npcs, default 50) and the
# downstream UMAP/neighbor graph require more cells than principal components;
# when QC/filtering removes nearly all cells (e.g. a very low-depth sample that
# fails a fixed count/feature threshold) RunPCA aborts with a cryptic
# "max(nu, nv) must be strictly less than min(nrow(A), ncol(A))" SVD error.
# Call this immediately before RunPCA to fail with an actionable message instead.
# The message embeds the stable token "[LOW_QUALITY_SAMPLE]" so the post-workflow
# quarantine step (workflow/scripts/quarantine_low_quality_samples.py) can detect
# flagged samples from their logs; keep the token in sync with that script.
require_min_cells_for_pca <- function(seurat_obj, context = "", npcs = 50L) {
  n_cells <- ncol(seurat_obj)
  if (n_cells <= npcs) {
    prefix <- if (nzchar(context)) paste0(context, ": ") else ""
    stop(sprintf(
      paste0(
        "%s[LOW_QUALITY_SAMPLE] only %d cell(s) remain - too few to compute %d principal ",
        "components (RunPCA and downstream UMAP/clustering require more cells than PCs). ",
        "This usually means upstream QC/filtering removed nearly all cells for this ",
        "sample; consider dataset-specific thresholds or excluding this sample."
      ),
      prefix, n_cells, npcs
    ), call. = FALSE)
  }
  invisible(n_cells)
}

add_silhouette_to_metadata <- function(
    seurat_obj,
    cluster_col = "seurat_clusters",
    reduction = "pca",
    dims = 1:30,
    sil_col_name = "silhouette_width",
    purity_col_name = "neighborhood_purity",
    purity_k = 50,
    purity_weighted = TRUE,
    purity_BNPARAM = NULL,
    purity_BPPARAM = NULL
) {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }

  emb <- Embeddings(seurat_obj, reduction = reduction)
  dims <- dims[dims <= ncol(emb)]
  if (length(dims) == 0) {
    stop(paste("No", reduction, "dimensions available for silhouette calculation"))
  }

  clust <- seurat_obj@meta.data[[cluster_col]]
  valid_cells <- !is.na(clust)
  sil_values <- rep(NA_real_, length(clust))
  purity_values <- rep(NA_real_, length(clust))

  if (sum(valid_cells) < 2 || length(unique(clust[valid_cells])) < 2) {
    seurat_obj@meta.data[[sil_col_name]] <- sil_values
    seurat_obj@meta.data[[purity_col_name]] <- purity_values
    return(seurat_obj)
  }

  message(sprintf(
    "Calculating approximate silhouette widths for %d cells with bluster::approxSilhouette",
    sum(valid_cells)
  ))
  sil <- bluster::approxSilhouette(
    emb[valid_cells, dims, drop = FALSE],
    clust[valid_cells]
  )
  sil_values[valid_cells] <- as.numeric(sil$width)
  seurat_obj@meta.data[[sil_col_name]] <- sil_values

  message(sprintf(
    "Calculating neighborhood purities for %d cells with bluster::neighborPurity",
    sum(valid_cells)
  ))
  bench_subphase("purity")
  purity <- bluster::neighborPurity(
    emb[valid_cells, dims, drop = FALSE],
    clust[valid_cells],
    k = purity_k,
    weighted = purity_weighted,
    BNPARAM = purity_BNPARAM,
    BPPARAM = purity_BPPARAM
  )
  purity_values[valid_cells] <- as.numeric(purity$purity)
  seurat_obj@meta.data[[purity_col_name]] <- purity_values

  seurat_obj
}
