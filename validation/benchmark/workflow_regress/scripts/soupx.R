source("workflow/scripts/bench_markers.R")
bench_phase("packages", "startup")
args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
raw <- args[2]
seurat_base <- args[3]
output <- args[4]
nclusters_output <- args[5]
cluster_ids_output <- args[6]

library("Seurat")
library("glmGamPoi")
library("SoupX")

options(future.globals.maxSize = 16 * 1024^3)

source("workflow/scripts/silhouette_utils.R")
WORKFLOW_SEED <- set_workflow_seed()

write_cluster_metadata <- function(seurat_obj, nclusters_output, cluster_ids_output) {
  cluster_ids <- levels(Idents(seurat_obj))
  if (is.null(cluster_ids) || length(cluster_ids) == 0) {
    cluster_ids <- sort(unique(as.character(Idents(seurat_obj))))
  }
  cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
  writeLines(cluster_ids, con = cluster_ids_output)
  writeLines(as.character(length(cluster_ids)), con = nclusters_output)
}

bench_phase("read", "read")
filtered_matrix <- Seurat::Read10X(filtered)
bench_phase("read_raw", "read")
raw_matrix <- Seurat::Read10X(raw)


bench_phase("read_base_rds", "read")
seurat_base <- readRDS(seurat_base)
bench_phase("soupx", "unported")
soup_channel <- SoupX::SoupChannel(tod = raw_matrix,toc=filtered_matrix,
                is10X = TRUE)
soup_channel$tod <- raw
# seurat_base carries a full SCT assay and neighbor graphs, but only its cluster labels and
# UMAP embedding are needed here. Extract those and drop the object before the memory-heavy
# autoEstCont/adjustCounts steps so its baggage is not held throughout. Results are unchanged.
soup_clusters <- as.factor(Idents(seurat_base))
soup_umap <- Seurat::Embeddings(seurat_base, "umap")
rm(seurat_base)
gc()
soup_channel <- SoupX::setClusters(soup_channel, clusters = soup_clusters)
soup_channel <- setDR(soup_channel, DR = soup_umap)
set.seed(WORKFLOW_SEED)
# autoEstCont aborts when it estimates an extremely high contamination fraction
# (> 0.8), treating it as a likely estimation failure. Across many datasets this
# hard stop kills otherwise-recoverable samples. Fall back to forceAccept = TRUE so
# the estimated fraction is used and the sample proceeds, with a clear warning.
soup_channel <- tryCatch(
  autoEstCont(soup_channel),
  error = function(e) {
    message(
      "autoEstCont failed (", conditionMessage(e),
      "); retrying with forceAccept = TRUE."
    )
    autoEstCont(soup_channel, forceAccept = TRUE)
  }
)
corrected_counts <- adjustCounts(soup_channel,roundToInt=TRUE)
seurat_soupx <- CreateSeuratObject(counts = corrected_counts)
rm(filtered_matrix, raw_matrix, soup_channel, corrected_counts)
gc()
bench_phase("percent_mt", "compute")
seurat_soupx[["percent.mt"]] <- PercentageFeatureSet(seurat_soupx, pattern = "(?i)^mt-")
bench_phase("sctransform", "compute")
seurat_soupx <- SCTransform(seurat_soupx, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("pca", "compute")
require_min_cells_for_pca(seurat_soupx, context = "soupx")
seurat_soupx <- RunPCA(seurat_soupx, seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("umap", "unported")
seurat_soupx <- RunUMAP(seurat_soupx, dims = 1:30, seed.use = WORKFLOW_SEED)
bench_phase("neighbors", "compute")
seurat_soupx <- FindNeighbors(seurat_soupx, dims = 1:30)
bench_phase("clusters", "compute")
seurat_soupx <- FindClusters(seurat_soupx, random.seed = WORKFLOW_SEED)
bench_phase("silhouette", "compute")
seurat_soupx <- add_silhouette_to_metadata(seurat_soupx)
bench_phase("write_rds", "write")
saveRDS(seurat_soupx,file=output)
bench_phase("write_cluster_ids", "write")
write_cluster_metadata(seurat_soupx, nclusters_output, cluster_ids_output)
bench_end()
