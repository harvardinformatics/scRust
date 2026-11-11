# Direct reference run of DoubletFinder 2.0.6 with the workflow's exact calls and seeds
# (docs/PLAN.md A4, B1 and B3). Run through validation/slurm/reference_doubletfinder.sbatch.
#
#   Rscript validation/R/reference_doubletfinder.R --input <export dir | object .rds> --out <dir>
#
# --input is normally a directory written by export_doubletfinder.R (10x Matrix Market counts
# plus seurat_clusters.tsv); a path ending in .rds reads the workflow's Seurat object instead.
# Outputs (plan B3):
#   cells.tsv        barcode, pANN, classification for every input cell, in input order
#   sweep_stats.tsv  pN, pK, BCreal (summarizeSweep), pN and pK as R's factor labels
#   bcmvn.tsv        ParamID, pK, MeanBC, VarBC, BCmetric (find.pK)
#   summary.json     "result": pK and its R label, nExp_poi, homotypic proportion, nExp,
#                    n_doublets, k, the column names R writes, final RNG state digest;
#                    "run": timings and peak RSS per phase, fingerprint (A5b), versions
#   find_pK.pdf      find.pK's plot (the workflow left it in Rplots.pdf)
# Doubles in TSV and JSON are written with 17 significant digits (exact; see lib/io.R).

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "doubletfinder.R"))

args <- df_args(list(input = NULL, out = NULL))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

df_setup()
t_load <- proc.time()[["elapsed"]]
inp <- df_load_input(args$input)
load_s <- proc.time()[["elapsed"]] - t_load
run_info <- df_run_info(args, "validation/R/reference_doubletfinder.R", inp$source)
run_info$load_input_s <- load_s
run_info$peak_rss_after_load_kib <- peak_rss_kib()

seu <- inp$seurat
rm(inp)
res <- df_run(seu, args$out)
rm(seu)
df_write_outputs(res, args$out, run_info)
cat("pK", res$result$pK_label, "nExp", res$result$nExp, "doublets", res$result$n_doublet_calls,
    "peak RSS (KiB)", peak_rss_kib(), "\n")
