source("workflow/scripts/bench_markers.R")
bench_phase("packages", "startup")
args <- commandArgs(trailingOnly = TRUE)
raw <- args[1]
seurat_output <- args[2]
matrix_outdir <- args[3]
nclusters_output <- args[4]
cluster_ids_output <- args[5]
sample_id <- if (length(args) >= 6 && nzchar(args[6])) args[6] else basename(normalizePath(raw, mustWork = FALSE))

library("Seurat")
library("DropletUtils")
library("scater")
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

# for back conversion to matrix
library("Matrix")
library("R.utils")

# run emptyDrops
bench_phase("read", "read")
droputil_rawdata <-read10xCounts(raw)
colnames(droputil_rawdata) <- droputil_rawdata$Barcode
rownames(droputil_rawdata) <- make.unique(rowData(droputil_rawdata)$Symbol)
bench_phase("emptydrops", "unported")
set.seed(WORKFLOW_SEED)
dropletutils_out <-emptyDrops(counts(droputil_rawdata))
is.cell <-!is.na(dropletutils_out$FDR) & dropletutils_out$FDR<0.001
cell_barcodes <- colnames(droputil_rawdata)[which(is.cell)]
droputil_filtered <- droputil_rawdata[,is.cell]
droputil_filtered <- logNormCounts(droputil_filtered)
seurat_droputil_filtered <- as.Seurat(droputil_filtered, counts = "counts",data="logcounts")
seurat_droputil_filtered$Sample <- sample_id
cm <- counts(droputil_filtered)
if (!inherits(cm, "dgCMatrix")) {
    cm <- as(cm, "dgCMatrix")
}
matrix_barcodes <- colnames(droputil_filtered)
matrix_features <- rownames(droputil_filtered)
rm(droputil_rawdata, dropletutils_out, droputil_filtered, is.cell, cell_barcodes)
gc()
seurat_droputil_filtered <- RenameAssays(seurat_droputil_filtered,originalexp="RNA")
bench_phase("percent_mt", "compute")
seurat_droputil_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_droputil_filtered, pattern = "(?i)^mt-")
bench_phase("sctransform", "compute")
seurat_droputil_filtered <- SCTransform(seurat_droputil_filtered, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("pca", "compute")
require_min_cells_for_pca(seurat_droputil_filtered, context = "emptydrops")
seurat_droputil_filtered <- RunPCA(seurat_droputil_filtered, seed.use = WORKFLOW_SEED, verbose = FALSE)
bench_phase("umap", "unported")
seurat_droputil_filtered <- RunUMAP(seurat_droputil_filtered, dims = 1:30, seed.use = WORKFLOW_SEED)
bench_phase("neighbors", "compute")
seurat_droputil_filtered <- FindNeighbors(seurat_droputil_filtered, dims = 1:30)
bench_phase("clusters", "compute")
seurat_droputil_filtered <- FindClusters(seurat_droputil_filtered, random.seed = WORKFLOW_SEED)
bench_phase("silhouette", "compute")
seurat_droputil_filtered <- add_silhouette_to_metadata(seurat_droputil_filtered)
bench_phase("write_rds", "write")
saveRDS(seurat_droputil_filtered,file=seurat_output)
bench_phase("write_cluster_ids", "write")
write_cluster_metadata(seurat_droputil_filtered, nclusters_output, cluster_ids_output)


# convert emptyDrops-filtered object to 10x matrix data structure
bench_phase("write_matrix", "write")
dir.create(matrix_outdir, showWarnings = FALSE)

# 1. Write uncompressed MatrixMarket file
writeMM(cm, file.path(matrix_outdir, "matrix.mtx"))

# 2. Write barcodes
writeLines(matrix_barcodes, file.path(matrix_outdir, "barcodes.tsv"))

# 3. Write features
feature_df <- data.frame(
    id   = matrix_features,
    name = matrix_features,
    feature_type = "Gene Expression"
)
write.table(
    feature_df,
    file.path(matrix_outdir, "features.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
)

# 4. Gzip them using R's built-in compression
R.utils::gzip(file.path(matrix_outdir, "matrix.mtx"), overwrite = TRUE)
R.utils::gzip(file.path(matrix_outdir, "barcodes.tsv"), overwrite = TRUE)
R.utils::gzip(file.path(matrix_outdir, "features.tsv"), overwrite = TRUE)
bench_end()
