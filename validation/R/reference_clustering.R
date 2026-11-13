# Direct reference run of the workflow's clustering chain (docs/PLAN.md C1 and C4; milestone M8):
# SCTransform(seed.use = 12345), RunPCA(seed.use = 12345), RunUMAP(dims = 1:30, seed.use = 12345),
# FindNeighbors(dims = 1:30), FindClusters(random.seed = 12345) and the workflow's
# add_silhouette_to_metadata() (bluster's approxSilhouette and neighborPurity), with the calls and
# packages of the workflow's scripts (lib/clustering.R). Run through
# validation/slurm/reference_clustering.sbatch.
#
#   Rscript validation/R/reference_clustering.R --input <export input dir> --out <dir> \
#     --level full|standard|light --umap yes|no --assay auto|v5|v3
#
# --input is an export's input/ (export_clustering.R): the object's RNA counts as saved, rebuilt as
# the identical dgCMatrix and wrapped with CreateSeuratObject(), as a v3 Assay for the emptyDrops
# objects (input/assay.txt, with --assay auto) or a v5 Assay5 otherwise. --umap no leaves RunUMAP
# out. Outputs (lib/clustering.R cl_main and cl_write_object): the chain's results in the same
# layout as the export's saved/, so the run compares with its saved object file by file;
# cluster_ids.txt and nclusters.txt as the workflow writes them; summary.json.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "clustering.R"))

args <- df_args(list(input = NULL, out = NULL, level = NULL, umap = "yes", assay = "auto"))
cl_main(args, "validation/R/reference_clustering.R")
