# Instrumented reference run of the workflow's clustering chain (docs/PLAN.md C4; milestone M8):
# the direct run of reference_clustering.R, with the same calls and seeds, plus dumps of every
# stage into <out>/stages/ (lib/clustering_hooks.R): SCTransform as M2's S2a-S2d (and vst's
# result in `scrust sctransform`'s layout), the PCA (S3) and the RNG checkpoints, RunUMAP's RNG
# checkpoints, the float items Annoy holds, the index file from a second identical build, every
# cell's neighbours and float distances, the optimizer's result before GroupSingletons and
# GroupSingletons' result, approxSilhouette's and neighborPurity's intermediates, and correct()'s
# median depth. Run through validation/slurm/reference_clustering.sbatch.
#
#   Rscript validation/R/instrumented_clustering.R --input <export input dir> --out <dir> \
#     --level full|standard|light --umap yes|no --assay auto|v5|v3 --purity-lists auto|yes|no
#
# The hooks only read values and write files, so the final files must equal the direct run's
# byte for byte, and the "result" section of summary.json too (which includes digests of
# .Random.seed after every step); compare_clustering.R checks this for every object that has
# both runs. --purity-lists: neighborPurity's neighbour lists (auto: at level full only; the plan
# asks for them on tier T and the small tier-P object).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))  # df_args
source(file.path(here, "lib", "sct_outputs.R"))
source(file.path(here, "lib", "doubletfinder_hooks.R"))
source(file.path(here, "lib", "clustering.R"))
source(file.path(here, "lib", "clustering_hooks.R"))

args <- df_args(list(input = NULL, out = NULL, level = NULL, umap = "yes", assay = "auto",
                     `purity-lists` = "auto"))
purity_lists <- switch(args$`purity-lists`, auto = identical(args$level, "full"), yes = TRUE, no = FALSE,
                       stop("--purity-lists must be auto, yes or no"))
hooks <- list(
  install = function(seu) cli_install(file.path(args$out, "stages"), args$level, colnames(seu),
                                      rownames(LayerData(seu, assay = "RNA", layer = "counts")), purity_lists),
  finish = function(seu, umap) cli_finish(seu, umap)
)
cl_main(args, "validation/R/instrumented_clustering.R", hooks)
