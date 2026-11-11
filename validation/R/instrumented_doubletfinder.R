# Instrumented reference run of DoubletFinder 2.0.6 (docs/PLAN.md A4 and B5): the direct run of
# reference_doubletfinder.R, with the same calls and seeds, plus dumps of the internal stages
# S1-S7 into <out>/stages/. Run through validation/slurm/reference_doubletfinder.sbatch.
#
#   Rscript validation/R/instrumented_doubletfinder.R --input <export dir | .rds> --out <dir> \
#     --level full|standard|light
#
# Mechanism (lib/doubletfinder_hooks.R): trace() adds passive entry/exit hooks to the installed
# functions (DoubletFinder, Seurat's SCTransform and RunPCA methods, sctransform, glmGamPoi,
# future.apply's future_lapply, irlba, KernSmooth's bkde) that only read values and write
# files. They draw no random numbers and modify nothing, so the run's final outputs must be
# byte-identical to the direct run's: cells.tsv, sweep_stats.tsv, bcmvn.tsv and the "result"
# section of summary.json, which includes a digest of the final .Random.seed. M2 checks this on
# every object that has both runs (validation/R/compare_doubletfinder.R).
#
# Stages (index conventions in lib/doubletfinder_hooks.R):
#   S1  sweep subsample (n > 10,000), doublet index pairs per pN and for the final run,
#       .Random.seed before and after the doublet draws
#   S2a SCTransform setup: cell attributes (umi, log_umi), filtered genes, log geometric means,
#       means and variances, fitting cells and genes (draw order), the gene-sampling pool,
#       probabilities and density, RNG after cell sampling and at the fits' entry
#   S2b per-gene fits: model_pars before and after the theta floor, glmGamPoi per bin (initial
#       dispersion and beta, first and second Newton fits, overdispersion MLE, dispersion
#       trend), RNG at every bin
#   S2c regularisation: model_pars_fit, outliers, Poisson genes, the smoothing genes, bandwidth
#   S2d residual mean and variance, the top features and VariableFeatures (ordered), the PCA
#       input matrix (scale.data rows in VariableFeatures order) and, at level full, the
#       clipped residual matrix y
#   S3  PCA embeddings, loadings, sdev, irlba's d, iterations; RNG at irlba entry and exit
#   S4  per real cell, a digest of its first `ind` (sweep) or k+1 (final) order() positions;
#       the positions themselves at level full
#   S5  pANN per pK for each sweep pN, and the final pANN
#   S6  bkde() grid and density, and the evaluated density passed to bimodality_coefficient,
#       for every (pN, pK) in sweep order; BCreal and bcmvn are the direct outputs
#   S7  the final Doublet calls

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))
source(file.path(here, "lib", "doubletfinder_hooks.R"))

args <- df_args(list(input = NULL, out = NULL, level = NULL))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

df_setup()
t_load <- proc.time()[["elapsed"]]
inp <- df_load_input(args$input)
load_s <- proc.time()[["elapsed"]] - t_load
run_info <- df_run_info(args, "validation/R/instrumented_doubletfinder.R", inp$source)
run_info$load_input_s <- load_s
run_info$peak_rss_after_load_kib <- peak_rss_kib()
run_info$dump_level <- args$level

seu <- inp$seurat
rm(inp)
dfi_init(file.path(args$out, "stages"), args$level, rownames(seu@meta.data),
         rownames(LayerData(seu, assay = "RNA", layer = "counts")))
dfi_install()

fns <- list(paramSweep = DoubletFinder::paramSweep, summarizeSweep = DoubletFinder::summarizeSweep,
            find.pK = DoubletFinder::find.pK, modelHomotypic = DoubletFinder::modelHomotypic,
            doubletFinder = DoubletFinder::doubletFinder)
res <- df_run(seu, args$out, fns = fns)
rm(seu)

n_pK <- length(res$result$sweep_pK_labels)
n_pN <- length(res$result$sweep_pN_labels)
dfi_finish(expected = list(paramSweep = 1L, parallel_paramSweep = n_pN, doubletFinder = 1L,
                           summarizeSweep = 1L, SCTransform.Seurat = n_pN + 1L,
                           SCTransform.default = n_pN + 1L, vst = n_pN + 1L,
                           vst_before_gene_sampling = n_pN + 1L, get_model_pars = n_pN + 1L,
                           reg_model_pars = n_pN + 1L, RunPCA.default = n_pN + 1L,
                           irlba = n_pN + 1L, bkde = n_pN * n_pK,
                           bimodality_coefficient = n_pN * n_pK,
                           glm_gp = .dfi$counts$fit_bins))
df_write_outputs(res, args$out, run_info)
cat("pK", res$result$pK_label, "nExp", res$result$nExp, "doublets", res$result$n_doublet_calls,
    "peak RSS (KiB)", peak_rss_kib(), "\n")
