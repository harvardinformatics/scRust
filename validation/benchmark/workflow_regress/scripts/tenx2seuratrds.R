source("workflow/scripts/bench_markers.R")
bench_phase("packages", "startup")
args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("glmGamPoi")

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
filtered_loaded <- Seurat::Read10X(filtered)
seurat_obj <- CreateSeuratObject(counts = filtered_loaded)
bench_phase("percent_mt", "compute")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "(?i)^mt-")
bench_phase("sctransform", "compute")
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("pca", "compute")
require_min_cells_for_pca(seurat_obj, context = "tenx2seuratrds")
seurat_obj <- RunPCA(seurat_obj, seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("umap", "unported")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, seed.use = WORKFLOW_SEED)
bench_phase("neighbors", "compute")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
bench_phase("clusters", "compute")
seurat_obj <- FindClusters(seurat_obj, random.seed = WORKFLOW_SEED)
bench_phase("silhouette", "compute")
seurat_obj <- add_silhouette_to_metadata(seurat_obj)
bench_phase("write_rds", "write")
saveRDS(seurat_obj,file=output)
bench_phase("write_cluster_ids", "write")
write_cluster_metadata(seurat_obj, nclusters_output, cluster_ids_output)
bench_end()
