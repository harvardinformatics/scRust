# scrustio's tests (docs/PLAN.md C3g; milestone M14). Run in the R env with the package installed
# in the repository's library (.rlib/), through sbatch or srun (validation/slurm/m14_tests.sbatch):
#   Rscript loaders/R/scrustio/tests/scrustio_tests.R
# SCRUST_ROOT (env.sh) locates the committed fixtures and .rlib/; the scrust binary (SCRUST_BIN, else
# $SCRUST_ROOT/target/release/scrust) makes the outputs the tests load: an export of a tier-T
# object, scrust sctransform --corrected-counts, scrust doubletfinder (with and without its nested
# outputs, M15) and scrust cluster without --export. Exit status 1 if a test fails.

root <- Sys.getenv("SCRUST_ROOT")
if (!nzchar(root)) stop("SCRUST_ROOT is not set (source env.sh)")
.libPaths(c(file.path(root, ".rlib"), .libPaths()))
suppressPackageStartupMessages({
  library("scrustio")
  library("Matrix")
  library("Seurat")
  library("SingleCellExperiment")
})

n_tests <- 0L
failures <- character()
check <- function(ok, what) {
  n_tests <<- n_tests + 1L
  if (!isTRUE(ok)) {
    failures <<- c(failures, what)
    cat("FAIL:", what, "\n")
  }
  invisible(ok)
}
refuses <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  check(!is.null(msg) && grepl(pattern, msg), paste0(what, if (is.null(msg)) " (no error)" else paste0(" (", msg, ")")))
}

tmp <- tempfile("scrustio-tests-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
bin <- Sys.getenv("SCRUST_BIN", file.path(root, "target", "release", "scrust"))
if (!file.exists(bin)) stop("no scrust binary at ", bin, " (cargo build --release -p scrust, or set SCRUST_BIN)")
run <- function(...) {
  out <- system2(bin, c(...), stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("scrust failed: ", paste(c(...), collapse = " "), "\n", paste(out, collapse = "\n"))
  invisible(out)
}
fx <- file.path(root, "tests", "fixtures")
counts_in <- file.path(fx, "clustering", "filtered_seurat_tenx_test", "input")
export <- file.path(tmp, "export")
run("cluster", "--counts", counts_in, "--preset", "workflow", "--export", "--out", export, "--threads", "2")
sct_out <- file.path(tmp, "sctransform")
run("sctransform", "--counts", counts_in, "--seed", "12345", "--corrected-counts", "--out", sct_out, "--threads", "2")
cl_out <- file.path(tmp, "cluster")
run("cluster", "--counts", counts_in, "--preset", "workflow", "--out", cl_out, "--threads", "2")
df_in <- file.path(fx, "doubletfinder", "soupx_tenx_test", "input")
df_out <- file.path(tmp, "doubletfinder")
run("doubletfinder", "--counts", df_in, "--clusters", file.path(df_in, "seurat_clusters.tsv"), "--out", df_out, "--threads", "2")

copy_output <- function(src, name) {
  dst <- file.path(tmp, name)
  dir.create(dst)
  invisible(file.copy(list.files(src, full.names = TRUE, all.files = FALSE, no.. = TRUE), dst, recursive = TRUE))
  dst
}
# Edits one field of the entry whose path is `path` in summary.json, as text (numbers elsewhere
# stay as written). `value` is JSON text.
tamper <- function(dir, path, field, value) {
  f <- file.path(dir, "summary.json")
  txt <- readLines(f)
  k <- grep(paste0("\"path\": \"", path, "\","), txt, fixed = TRUE)
  stopifnot(length(k) == 1)
  for (j in k:(k + 30)) {
    if (grepl(paste0("^[[:space:]]*\"", field, "\": "), txt[j])) {
      txt[j] <- sub(paste0("(\"", field, "\": ).*?(,?)$"), paste0("\\1", value, "\\2"), txt[j], perl = TRUE)
      writeLines(txt, f)
      return(invisible(TRUE))
    }
  }
  stop("no field ", field, " after ", path)
}

# 1. The slots: every wrong combination of role, scale and format is refused ---------------------
st <- slot_table()
roles <- scrustio:::.ROLES
scales <- scrustio:::.SCALES
formats <- scrustio:::.FORMATS
n_refused <- 0L
n_wrong <- 0L
for (k in seq_len(nrow(st))) {
  slot <- st$slot[k]
  for (r in roles) for (s in scales) {
    ok <- tryCatch({ scrustio:::.check_slot(slot, r, s, "f", st$format[k]); TRUE }, error = function(e) FALSE)
    want <- r == st$role[k] && s == st$scale[k]
    if (!want) {
      n_wrong <- n_wrong + 1L
      if (!ok) n_refused <- n_refused + 1L
    }
    if (ok != want) check(FALSE, paste("slot", slot, "role", r, "scale", s))
  }
  for (f in formats) {
    ok <- tryCatch({ scrustio:::.check_slot(slot, st$role[k], st$scale[k], "f", f); TRUE }, error = function(e) FALSE)
    if (f != st$format[k]) n_wrong <- n_wrong + 1L
    if (!ok && f != st$format[k]) n_refused <- n_refused + 1L
    if (ok != (f == st$format[k])) check(FALSE, paste("slot", slot, "format", f))
  }
}
check(n_refused == n_wrong, "every wrong combination refused")
cat(sprintf("slots: %d slots; %d wrong combinations of role, scale and format, %d refused\n", nrow(st), n_wrong, n_refused))
check(all(st$role[st$raw_counts] == "raw_counts") && all(st$scale[st$raw_counts] == "raw_counts"),
      "the raw-count slots take only raw counts")
check(!any(st$scale[!st$raw_counts] == "raw_counts"), "no other slot takes raw counts")

# The content rule: the comment lines
refuses(scrustio:::.check_comments("Seurat RNA counts", c("% scrust: x. NOT raw counts: y"), "f"), "NOT raw counts",
        "a file marked NOT raw counts is refused for a raw-count slot")
refuses(scrustio:::.check_comments("Seurat SCT counts", c("%metadata_json: {}"), "f"), "not so marked",
        "an unmarked file is refused for an SCT layer's slot")

# Tampered manifests, end to end (each loader must refuse)
both <- function(dir, pattern, what) {
  refuses(load_seurat(dir), pattern, paste("load_seurat:", what))
  refuses(load_sce(dir), pattern, paste("load_sce:", what))
}
d <- copy_output(export, "swap_counts")
tamper(d, "sct/corrected_counts.mtx.gz", "role", "\"raw_counts\"")
tamper(d, "sct/corrected_counts.mtx.gz", "scale", "\"raw_counts\"")
tamper(d, "counts/matrix.mtx.gz", "role", "\"sct_counts\"")
tamper(d, "counts/matrix.mtx.gz", "scale", "\"corrected_counts\"")
both(d, "NOT raw counts", "the corrected counts declared raw counts")
d <- copy_output(export, "scale_only")
tamper(d, "sct/corrected_counts.mtx.gz", "scale", "\"raw_counts\"")
both(d, "refusing to put sct/corrected_counts.mtx.gz \\(role sct_counts, scale raw_counts\\)", "the corrected counts with scale raw_counts")
d <- copy_output(export, "residuals_as_counts")
tamper(d, "counts/matrix.mtx.gz", "role", "\"stage_dump\"")
tamper(d, "sct/scale_data.npy", "role", "\"raw_counts\"")
tamper(d, "sct/scale_data.npy", "scale", "\"raw_counts\"")
both(d, "format npy", "the residuals declared raw counts")
d <- copy_output(export, "log1p_as_counts")
tamper(d, "sct/corrected_counts.mtx.gz", "role", "\"sct_data\"")
tamper(d, "sct/corrected_counts.mtx.gz", "scale", "\"log1p_corrected_counts\"")
tamper(d, "sct/log1p_corrected_counts.mtx.gz", "role", "\"sct_counts\"")
tamper(d, "sct/log1p_corrected_counts.mtx.gz", "scale", "\"corrected_counts\"")
both(d, "not log1p of the counts layer", "the log1p layer and the corrected counts swapped")

# 2. The manifest -------------------------------------------------------------------------------------
d <- copy_output(export, "schema")
txt <- readLines(file.path(d, "summary.json"))
writeLines(sub("\"schema\": \"scrust-files\"", "\"schema\": \"other-files\"", txt), file.path(d, "summary.json"))
refuses(read_manifest(d), "schema", "another schema is refused")
d <- copy_output(export, "version")
txt <- readLines(file.path(d, "summary.json"))
writeLines(sub("^(    \"version\": )1,", "\\12,", txt), file.path(d, "summary.json"))
refuses(read_manifest(d), "version 2", "another version is refused")
d <- copy_output(export, "missing")
invisible(file.remove(file.path(d, "sct", "corrected_counts.mtx.gz")))
refuses(read_manifest(d), "missing", "a listed file that is missing is refused")
d <- copy_output(export, "size")
cat("x", file = file.path(d, "pca", "stdev.tsv"), append = TRUE)
refuses(read_manifest(d), "bytes", "a listed file whose size changed is refused")
d <- copy_output(export, "sha")
f <- file.path(d, "pca", "stdev.tsv")
x <- readLines(f)
last <- substring(x[2], nchar(x[2]))
x[2] <- paste0(substring(x[2], 1, nchar(x[2]) - 1), if (last == "2") "3" else "2")
writeLines(x, f)
check(inherits(read_manifest(d), "scrust_manifest"), "a same-size change passes without the checksum check")
refuses(read_manifest(d, check_sha256 = TRUE), "SHA-256", "the checksum check catches a same-size change")
d <- copy_output(export, "shape")
tamper(d, "sct/corrected_counts.mtx.gz", "shape", "[820, 300]")
refuses(load_seurat(d), "shape gives 820 rows", "rows not matching the shape are refused")
d <- copy_output(export, "traversal")
tamper(d, "pca/stdev.tsv", "path", "\"../stdev.tsv\"")
refuses(read_manifest(d), "not a relative path", "a path outside the output is refused")
# by role, never by name: rename the corrected counts and their identifier files
d <- copy_output(export, "renamed")
invisible(file.rename(file.path(d, "sct", "corrected_counts.mtx.gz"), file.path(d, "sct", "layer_a.mtx.gz")))
invisible(file.rename(file.path(d, "sct", "sct_features.tsv.gz"), file.path(d, "sct", "rows_a.tsv.gz")))
txt <- readLines(file.path(d, "summary.json"))
txt <- gsub("sct/corrected_counts.mtx.gz", "sct/layer_a.mtx.gz", txt, fixed = TRUE)
txt <- gsub("sct/sct_features.tsv.gz", "sct/rows_a.tsv.gz", txt, fixed = TRUE)
writeLines(txt, file.path(d, "summary.json"))

# 3. Seurat -------------------------------------------------------------------------------------------
obj <- load_seurat(export)
check(identical(load_seurat(d), obj), "files are found by role: renamed files load the same object")
check(identical(class(obj[["RNA"]])[1], "Assay5") && identical(class(obj[["SCT"]])[1], "SCTAssay"), "RNA is an Assay5, SCT an SCTAssay")
check(identical(DefaultAssay(obj), "SCT"), "SCT is the default assay")
check(identical(names(obj[["SCT"]]@SCTModel.list), "model1"), "exactly one SCTModel")
lev <- readLines(file.path(export, "seurat_clusters_levels.txt"))
check(identical(levels(Idents(obj)), lev) && identical(levels(obj$seurat_clusters), lev), "identities and seurat_clusters in their levels' order")
check(identical(unname(as.character(Idents(obj))), as.character(obj$seurat_clusters)), "identities are seurat_clusters")
check(identical(levels(obj$SCT_snn_res.0.8), readLines(file.path(export, "SCT_snn_res.0.8_levels.txt"))), "res. column in its levels' order")
cnt <- LayerData(obj, assay = "RNA", layer = "counts")
m <- read_manifest(export)
check(identical(dim(cnt), as.integer(manifest_entries(m, "raw_counts")[[1]]$shape)), "RNA counts: the raw counts' shape")
check(identical(obj[["SCT"]]@data@x, log1p(obj[["SCT"]]@counts@x)), "SCT data is log1p of SCT counts")
check(identical(rownames(obj[["SCT"]]@scale.data), intersect(rownames(obj[["SCT"]]), VariableFeatures(obj[["SCT"]]))),
      "scale.data's rows in the SCT assay's gene order")
check(identical(rownames(Loadings(obj[["pca"]])), VariableFeatures(obj[["SCT"]])), "the loadings in VariableFeatures order")
check(identical(dim(Embeddings(obj[["pca"]])), c(ncol(obj), 50L)), "50 PCs")
check(all(c("SCT_nn", "SCT_snn") %in% Graphs(obj)), "both graphs")
fm <- FindMarkers(obj, ident.1 = lev[1])
check(nrow(fm) > 0 && all(rownames(fm) %in% rownames(obj[["SCT"]])), "FindMarkers runs on the SCT assay")
saved <- read.delim(gzfile(file.path(root, "tests", "fixtures", "clustering", "filtered_seurat_tenx_test", "saved", "metadata.tsv.gz")),
                    colClasses = "character", check.names = FALSE)
check(identical(as.character(obj$seurat_clusters), saved$seurat_clusters), "the labels equal the saved object's (fixture)")
f <- file.path(tmp, "object.rds")
x <- save_rds(export, f)
check(identical(readRDS(f), obj) && identical(x, obj), "save_rds: the .rds reloads into the loader's object")
w <- tryCatch({ scrustio:::.check_versions(c(Seurat = "0.0.1")); NULL }, warning = function(w) conditionMessage(w))
check(!is.null(w) && grepl("untested", w), "another Seurat version warns")

# 3b. The SCTModel's cells_step1: every cell up to ncells; above it, drawn again from the
# manifest's seed and ncells and checked against the step-1 genes the export lists; the caller's
# random-number state left as it was --------------------------------------------------------------
check(all(obj[["SCT"]]@SCTModel.list[[1]]@cell.attributes$cells_step1), "cells_step1: every cell when there are at most ncells")
drawn <- file.path(tmp, "export_ncells100")
run("cluster", "--counts", counts_in, "--preset", "workflow", "--ncells", "100", "--export", "--dump-stages", "light",
    "--out", drawn, "--threads", "2")
rust_step1 <- scrustio:::.read_npy(file.path(drawn, "stages", "S2a_cells_step1.npy")) + 1L
# loads drawn's Seurat object after setup(); the cells_step1 it rebuilds, and whether the caller's
# generator kinds and .Random.seed (or its absence) are as before, with no warning
rng_kept <- function(setup) {
  setup()
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  before <- if (had) get(".Random.seed", envir = env) else NULL
  kind <- RNGkind()
  warned <- FALSE
  o <- withCallingHandlers(load_seurat(drawn), warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  has <- exists(".Random.seed", envir = env, inherits = FALSE)
  kept <- identical(RNGkind(), kind) && identical(has, had) && (!had || identical(get(".Random.seed", envir = env), before))
  list(step1 = o[["SCT"]]@SCTModel.list[[1]]@cell.attributes$cells_step1, kept = kept && !warned)
}
r1 <- rng_kept(function() { suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection")); set.seed(1) })
check(sum(r1$step1) == 100 && identical(which(r1$step1), sort(as.integer(rust_step1))),
      "cells_step1: the 100 cells scrust drew (its own record, stages/S2a_cells_step1.npy)")
check(r1$kept, "cells_step1: the caller's generator state is restored (R's default kinds)")
r2 <- rng_kept(function() { suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rounding")); set.seed(2) })
check(identical(r2$step1, r1$step1), "cells_step1: the same draw whatever the caller's generator kinds")
check(r2$kept, "cells_step1: the caller's generator kinds and state are restored (L'Ecuyer-CMRG, Box-Muller, Rounding)")
r3 <- rng_kept(function() {
  suppressWarnings(RNGkind("default", "default", "default"))
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
})
check(identical(r3$step1, r1$step1) && r3$kept, "cells_step1: no .Random.seed before the load, none after")
suppressWarnings(RNGkind("default", "default", "default"))
check(identical(as.data.frame(colData(altExp(load_sce(drawn), "SCT")))$cells_step1, r1$step1),
      "cells_step1: the same in the SingleCellExperiment's altExp SCT")
check(identical(load_seurat(drawn)[["SCT"]]@SCTModel.list[[1]]@cell.attributes[, c("umi", "log_umi")],
                obj[["SCT"]]@SCTModel.list[[1]]@cell.attributes[, c("umi", "log_umi")]),
      "cells_step1: the other cell attributes as without the draw")
# the self-check: two of the listed step-1 genes swapped (the file's size unchanged)
d <- copy_output(drawn, "step1_swapped")
f <- file.path(d, "sct", "model_pars.tsv")
x <- readLines(f)
k <- which(!startsWith(x, "#"))[2:3]
x[k] <- x[rev(k)]
writeLines(x, f)
refuses(load_seurat(d), "cells_step1 cannot be rebuilt", "cells_step1: a draw that does not give the listed step-1 genes is refused")
refuses(load_sce(d), "cells_step1 cannot be rebuilt", "cells_step1: the SingleCellExperiment refuses it too")
# another seed, in both places the manifest records it: another draw, which the check catches
tamper_line <- function(dir, anchors, field, value) {
  f <- file.path(dir, "summary.json")
  txt <- readLines(f)
  k <- 0L
  for (a in anchors) k <- k + match(TRUE, grepl(a, txt[(k + 1):length(txt)], fixed = TRUE))
  j <- k + match(TRUE, grepl(paste0("\"", field, "\": "), txt[(k + 1):length(txt)], fixed = TRUE))
  txt[j] <- sub(paste0("(\"", field, "\": )[^,]*"), paste0("\\1", value), txt[j])
  writeLines(txt, f)
}
d <- copy_output(drawn, "step1_seed")
tamper_line(d, c("\"files\": {", "\"provenance\": {", "\"sctransform\": {"), "seed_use", "12346")
tamper_line(d, c("\"files\": {", "\"provenance\": {", "\"seeds\": {"), "SCTransform seed.use", "12346")
refuses(load_seurat(d), "cells_step1 cannot be rebuilt", "cells_step1: a draw from another seed is caught by the check")
d <- copy_output(drawn, "step1_seeds_differ")
tamper_line(d, c("\"files\": {", "\"provenance\": {", "\"seeds\": {"), "SCTransform seed.use", "12346")
refuses(load_seurat(d), "seeds differ", "cells_step1: the manifest's two records of the seed must agree")
d <- copy_output(drawn, "step1_n")
tamper_line(d, c("\"result\": {", "\"sctransform\": {"), "n_cells_step1", "101")
refuses(load_seurat(d), "n_cells_step1 \\(101\\) is not min\\(ncells, cells\\) = 100", "cells_step1: the summary's n_cells_step1 must be min(ncells, cells)")
# The check runs only under the versions it mirrors (`tested`, injectable as .check_versions()'s
# is); with others the cells are drawn without it, with one warning, and nothing is refused.
check(is.null(scrustio:::.step1_check_skip()), "cells_step1: this session has the versions the check mirrors (R, sctransform, Linux)")
why <- scrustio:::.step1_check_skip(c(R = "4.3.0", sctransform = "0.4.3", os = "linux"))
check(!is.null(why) && grepl(paste("this session has R", getRversion()), why, fixed = TRUE), "cells_step1: another R skips the check, saying so")
other <- c(R = as.character(getRversion()), sctransform = "0.0.1", os = "linux")
m_d <- read_manifest(drawn)
feat_d <- scrustio:::.counts_features(m_d)
cells_d <- scrustio:::.cells_of(m_d)
cnt_d <- read_for_slot(m_d, "raw_counts", "Seurat RNA counts")
rownames(cnt_d) <- feat_d$seurat
# .cells_step1() under `other` versions, from a known generator state: its result, its warnings and
# whether the state is as before
skip_run <- function(m, tabs) {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(4)
  before <- get(".Random.seed", envir = globalenv())
  kind <- RNGkind()
  warns <- character()
  x <- withCallingHandlers(scrustio:::.cells_step1(m, cells_d, cnt_d, tabs, tested = other), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(step1 = x, warns = warns, kept = identical(get(".Random.seed", envir = globalenv()), before) && identical(RNGkind(), kind))
}
k1 <- skip_run(m_d, scrustio:::.sct_tables(m_d, feat_d))
check(identical(k1$step1, r1$step1), "cells_step1 under other versions: the same cells, drawn without the check")
check(length(k1$warns) == 1 && grepl("without its check.*this session has sctransform", k1$warns),
      "cells_step1 under other versions: one warning, saying why")
check(k1$kept, "cells_step1 under other versions: the caller's generator state restored")
m_s <- read_manifest(file.path(tmp, "step1_swapped"))
k2 <- skip_run(m_s, scrustio:::.sct_tables(m_s, scrustio:::.counts_features(m_s)))
check(identical(k2$step1, r1$step1) && length(k2$warns) == 1,
      "cells_step1 under other versions: an output the check would refuse loads, with one warning")
suppressWarnings(RNGkind("default", "default", "default"))

# 4. SingleCellExperiment ----------------------------------------------------------------------------
sce <- load_sce(export)
check(identical(assay(sce, "counts"), cnt), "SCE counts: the raw counts")
check(identical(altExpNames(sce), c("SCT", "SCT_scale_data")), "SCE altExps SCT and SCT_scale_data")
check(identical(assayNames(altExp(sce, "SCT")), c("corrected_counts", "logcounts")), "SCT altExp's assays")
check(identical(assay(altExp(sce, "SCT"), "corrected_counts"), obj[["SCT"]]@counts), "corrected_counts equal Seurat's SCT counts")
check(identical(rownames(assay(altExp(sce, "SCT_scale_data"), "scaledata")), VariableFeatures(obj[["SCT"]])),
      "SCT_scale_data in VariableFeatures order")
e <- reducedDim(sce, "PCA")
check(identical(attr(e, "rotation"), Loadings(obj[["pca"]])), "reducedDim PCA's rotation is the loadings")
check(identical(as(colPair(sce, "SCT_snn", asSparse = TRUE), "CsparseMatrix")@x, obj[["SCT_snn"]]@x), "colPair SCT_snn's values")
check(identical(levels(sce$seurat_clusters), lev), "colData factors in their levels' order")
f <- file.path(tmp, "sce.rds")
x <- save_rds(export, f, object = "sce")
check(identical(readRDS(f), sce) && identical(x, sce), "save_rds(object = \"sce\"): the .rds reloads into the loader's object")

# 5. Other commands' outputs --------------------------------------------------------------------------
so <- load_seurat(sct_out)
check(identical(SeuratObject::Assays(so), "SCT") && identical(DefaultAssay(so), "SCT"), "scrust sctransform: an SCT assay alone")
check(identical(so[["SCT"]]@counts, obj[["SCT"]]@counts) && identical(so[["SCT"]]@SCTModel.list, obj[["SCT"]]@SCTModel.list),
      "scrust sctransform: the same SCT layers and model as the export")
check(all(c("nCount_RNA", "nFeature_RNA") %in% colnames(so[[]])), "scrust sctransform: the counts' per-cell columns")
ss <- load_sce(sct_out)
check(nrow(ss) == 0 && identical(assay(altExp(ss, "SCT"), "logcounts"), obj[["SCT"]]@data), "scrust sctransform as an SCE")
dm <- read_manifest(df_out)
res <- dm$summary$result
dfo <- load_seurat(df_out)
check(identical(SeuratObject::Assays(dfo), "RNA"), "scrust doubletfinder: the RNA assay")
singlets <- read.delim(gzfile(file.path(df_out, "singlets", "barcodes.tsv.gz")), header = FALSE)$V1
check(identical(colnames(dfo), singlets), "scrust doubletfinder: the singlets are the cells")
check(all(c(res$pANN_column, res$classification_column) %in% colnames(dfo[[]])), "DoubletFinder's column names")
check(all(dfo[[res$classification_column, drop = TRUE]] == "Singlet"), "the singlets' calls")
dfs <- load_sce(df_out)
check(identical(colnames(dfs), singlets) && res$pANN_column %in% colnames(colData(dfs)), "scrust doubletfinder as an SCE")
refuses(load_seurat(cl_out), "neither raw counts nor", "scrust cluster without --export holds no assay data")
# add_scrust: DoubletFinder's calls onto an object of all the input cells; clusters onto a bare object
raw <- Read10X(df_in)
base <- CreateSeuratObject(raw)
with_calls <- add_scrust(base, df_out)
check(sum(with_calls[[res$classification_column, drop = TRUE]] == "Doublet") == res$n_doublet_calls, "add_scrust: the doublets' calls")
bare <- CreateSeuratObject(Read10X(counts_in))
with_cl <- suppressMessages(add_scrust(bare, cl_out))
check(identical(with_cl$seurat_clusters, obj$seurat_clusters), "add_scrust: scrust cluster's labels with their levels")

# 6. Nested outputs (M15) ----------------------------------------------------------------------------
# scrust doubletfinder clustering its own input (clusters/) and re-clustering the singlets with the
# export (singlets_clusters/): the DoubletFinder output loads as before, by its own files' roles; each nested
# output loads on its own
dfn <- file.path(tmp, "doubletfinder_nested")
run("doubletfinder", "--counts", df_in, "--out", dfn, "--cluster-stats", "silhouette,purity", "--recluster",
    "--recluster-export", "--threads", "2")
nm <- read_manifest(dfn)
roles <- vapply(nm$entries, `[[`, "", "role")
check(sum(roles == "nested_output") == 2 && sum(roles == "raw_counts") == 1,
      "nested outputs: listed by their summary.json; the parent's one raw_counts")
dfno <- load_seurat(dfn)
check(identical(colnames(dfno), colnames(dfo)) && identical(SeuratObject::Assays(dfno), "RNA") &&
        identical(dfno[[res$pANN_column, drop = TRUE]], dfo[[res$pANN_column, drop = TRUE]]),
      "nested outputs: the DoubletFinder output loads as before")
rc <- load_seurat(file.path(dfn, "singlets_clusters"))
check(identical(colnames(rc), singlets) && identical(DefaultAssay(rc), "SCT") && "pca" %in% names(rc@reductions) &&
        is.factor(rc$seurat_clusters), "nested outputs: the re-clustering loads on its own, with the export")
with_hook <- suppressMessages(add_scrust(base, file.path(dfn, "clusters")))
check(identical(levels(with_hook$seurat_clusters), readLines(file.path(df_in, "seurat_clusters_levels.txt"))) &&
        all(c("silhouette_width", "neighborhood_purity") %in% colnames(with_hook[[]])),
      "nested outputs: clustering on the fly's labels, levels and scores, by add_scrust")

cat(sprintf("scrustio %s: %d tests, %d failed\n", as.character(packageVersion("scrustio")), n_tests, length(failures)))
if (length(failures)) quit(status = 1)
