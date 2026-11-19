#!/usr/bin/env Rscript
# Milestone M15b (docs/PLAN.md C3, "vars.to.regress and percent.mt"): the code path of
# SCTransform's second regression and of PercentageFeatureSet() as the pinned env runs them, and
# reference vectors for the Rust port (rcompat::linpack, scio's percent.mt).
#
# 1. The code path. The installed functions are compared with the source tarballs that
#    third_party/fetch_sources.sh fetched (their bodies parsed from the tarball's files must be
#    identical() to the installed ones): Seurat 5.5.1's SCTransform.default, .Assay, .Seurat,
#    .StdAssay, FetchResiduals.SCTAssay, FetchResidualSCTModel, ScaleData.default,
#    RegressOutMatrix, FastRowScale, PercentageFeatureSet; R 4.5.3's lm, lm.fit, qr.resid; and
#    how R's grep (TRE, perl = FALSE) treats "(?i)^mt-".
# 2. tests/fixtures/rcompat/linpack_qr.tsv and linpack_resid.tsv: the QR that
#    RegressOutMatrix() takes once, lm(GENE ~ <vars>, qr = TRUE)$qr (LINPACK dqrdc2 through
#    Cdqrls, tol = 1e-7), and qr.resid() of several rows (dqrrsd -> dqrsl), with
#    RegressOutMatrix()'s own result for the same rows. dqrdc2 and dqrsl call BLAS ddot, daxpy,
#    dnrm2 and dscal: run with OPENBLAS_CORETYPE=Haswell, the portable profile (env.sh's default).
# 3. tests/fixtures/clustering/percent_mt/: a small 10x directory with tricky feature names
#    (case, underscores and pipes that Seurat turns into dashes, duplicates that Read10X makes
#    unique, a mitochondrial gene with no counts, non-integer counts in a second matrix) and
#    PercentageFeatureSet()'s results, by the pattern and by a list of features.
#
# Reads no .rds. Runs in the production R env with the hermetic flags, through srun or sbatch:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     "$R_ENV/bin/Rscript" validation/R/regress_vectors.R <repo root>
#
# Doubles are written as their IEEE-754 bit patterns (16 hex digits), lists comma-separated.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
fx <- file.path(root, "tests", "fixtures")
stopifnot(R.version$major == "4", R.version$minor == "5.3")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
stopifnot(as.character(packageVersion("Seurat")) == "5.5.1")
cat("OpenBLAS kernel:", Sys.getenv("OPENBLAS_CORETYPE"), " threads:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
cat("La_library:", La_library(), "\n")
cat("BLAS (extSoftVersion):", extSoftVersion()[["BLAS"]], "\n")

hex <- function(x) {
  x <- as.double(x)
  if (!length(x)) return(character())
  b <- writeBin(x, raw(), size = 8, endian = "big")
  apply(matrix(as.character(b), nrow = 8), 2, paste, collapse = "")
}
lst <- function(v) if (length(v)) paste(v, collapse = ",") else "-"
hexl <- function(x) lst(hex(x))

# 1. The code path ----------------------------------------------------------------------------

src <- file.path(root, "third_party", "src")
from_source <- function(file, name) {
  e <- new.env()
  exprs <- parse(file = file, keep.source = FALSE)
  for (ex in exprs) {
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) && identical(ex[[2]], as.name(name))) {
      return(eval(ex[[3]], e))
    }
  }
  stop(name, " not found in ", file)
}
same_fun <- function(installed, file, name) {
  f <- from_source(file, name)
  ok <- identical(body(installed), body(f)) && identical(formals(installed), formals(f))
  cat(sprintf("%-28s %s (%s)\n", name, if (ok) "identical to the source" else "DIFFERS", basename(file)))
  ok
}
seu_r <- file.path(src, "Seurat-5.5.1", "R")
checks <- c(
  same_fun(Seurat:::SCTransform.default, file.path(seu_r, "preprocessing.R"), "SCTransform.default"),
  same_fun(Seurat:::SCTransform.Assay, file.path(seu_r, "preprocessing.R"), "SCTransform.Assay"),
  same_fun(Seurat:::SCTransform.Seurat, file.path(seu_r, "preprocessing.R"), "SCTransform.Seurat"),
  same_fun(Seurat:::SCTransform.StdAssay, file.path(seu_r, "preprocessing5.R"), "SCTransform.StdAssay"),
  same_fun(Seurat:::FetchResiduals.SCTAssay, file.path(seu_r, "preprocessing5.R"), "FetchResiduals.SCTAssay"),
  same_fun(Seurat:::FetchResidualSCTModel, file.path(seu_r, "preprocessing5.R"), "FetchResidualSCTModel"),
  same_fun(Seurat:::ScaleData.default, file.path(seu_r, "preprocessing.R"), "ScaleData.default"),
  same_fun(Seurat:::RegressOutMatrix, file.path(seu_r, "preprocessing.R"), "RegressOutMatrix"),
  same_fun(Seurat::FastRowScale, file.path(seu_r, "utilities.R"), "FastRowScale"),
  same_fun(Seurat::PercentageFeatureSet, file.path(seu_r, "utilities.R"), "PercentageFeatureSet"),
  same_fun(stats::lm, file.path(src, "R-4.5.3", "src", "library", "stats", "R", "lm.R"), "lm"),
  same_fun(stats::lm.fit, file.path(src, "R-4.5.3", "src", "library", "stats", "R", "lm.R"), "lm.fit"),
  same_fun(base::qr.resid, file.path(src, "R-4.5.3", "src", "library", "base", "R", "qr.R"), "qr.resid")
)
if (!all(checks)) stop("an installed function differs from its source")
cat("future plan workers:", future::nbrOfWorkers(), "\n")

# grep's TRE engine (perl = FALSE) and the inline flag
names_tre <- c("mt-Nd1", "MT-ND1", "Mt-Co1", "mT-Cytb", "mt-", "MT-", "mt", "MT", "mt.1", "xmt-a",
               " mt-a", "mt -a", "mt‐a", "ｍt-a", "m-t-a", "MT--", "mt-mt-", "Mtor", "MTX-1",
               "mt_a", "mt|a")
hit <- grepl("(?i)^mt-", names_tre)
cat("grep(\"(?i)^mt-\", perl = FALSE):\n")
for (i in seq_along(names_tre)) cat(sprintf("  %-10s %s\n", encodeString(names_tre[i], quote = "\""), hit[i]))
want <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
          FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
if (!identical(hit, want)) stop("TRE's matches are not the ones expected")
# the same with perl = TRUE (PCRE): for the record
stopifnot(identical(grepl("(?i)^mt-", names_tre, perl = TRUE), want))

# 2. The QR and the residuals ------------------------------------------------------------------

set.seed(20261119)
cases <- list()
add <- function(name, latent) cases[[name]] <<- as.data.frame(latent)
add("full_n5", list(v = rnorm(5)))
add("percent_n40", list(percent.mt = round(runif(40) * 100, 3)))
add("count_score_n151", list(nCount_RNA = round(runif(151) * 20000) + 500, score = runif(151) * 5))
add("zero_n151", list(percent.mt = rep(0, 151)))
add("constant_n151", list(v = rep(7.25, 151)))
v151 <- rnorm(151)
add("collinear_n151", list(v = v151, w = 2 * v151))
pz <- ifelse(runif(613) < 0.3, 0, rexp(613) * 3)
add("percent_count_n613", list(percent.mt = pz, nCount_RNA = round(rexp(613) * 4000) + 200))
add("square_n2", list(v = c(1, 3)))
add("square_n3", list(v = rnorm(3), w = rnorm(3)))
add("wide_n2", list(v = rnorm(2), w = rnorm(2)))
add("outlier_n151", list(v = c(rnorm(150), 1e6)))
add("block_n64", list(v = runif(64)))
add("odd_n17", list(v = runif(17), w = rnorm(17)))
add("nearly_n151", list(v = rep(1, 151) + c(1e-9, rep(0, 150))))

qr_rows <- character()
resid_rows <- character()
for (name in names(cases)) {
  L <- cases[[name]]
  n <- nrow(L)
  Y <- rbind(round(rnorm(n) * 2, 6), pmin(pmax(rnorm(n) * 3, -3), 3), rep(0.5, n))
  rownames(Y) <- paste0("g", seq_len(nrow(Y)))
  colnames(Y) <- paste0("c", seq_len(n))
  rownames(L) <- colnames(Y)
  # as RegressOutMatrix takes it (Seurat 5.5.1, preprocessing.R): the formula, then lm's qr once
  fmla <- as.formula(paste("GENE ~", paste(colnames(L), collapse = "+")))
  regression.mat <- cbind(L, Y[1, ])
  colnames(regression.mat) <- c(colnames(L), "GENE")
  qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
  X <- model.matrix(fmla, data = regression.mat)
  qr_rows <- c(qr_rows, paste(name, n, ncol(X), paste(colnames(L), collapse = ","), hexl(X),
                              hexl(qr$qr), hexl(qr$qraux), lst(qr$pivot), qr$rank, hexl(qr$tol), sep = "\t"))
  seurat <- Seurat:::RegressOutMatrix(data.expr = Y, latent.data = L, features.regress = rownames(Y),
                                      model.use = "linear", use.umi = FALSE, verbose = FALSE)
  for (k in seq_len(nrow(Y))) {
    r <- qr.resid(qr = qr, y = Y[k, ])
    if (!identical(unname(r), unname(seurat[k, ]))) stop(name, " row ", k, ": qr.resid and RegressOutMatrix differ")
    resid_rows <- c(resid_rows, paste(name, k, hexl(Y[k, ]), hexl(r), sep = "\t"))
  }
  cat(sprintf("%-20s n %4d, p %d, rank %d, pivot %s\n", name, n, ncol(X), qr$rank, lst(qr$pivot)))
}
dir.create(file.path(fx, "rcompat"), showWarnings = FALSE, recursive = TRUE)
gz_write <- function(lines, path) {
  writeLines(lines, path)
  if (system2("gzip", c("-n", "-9", "-f", shQuote(path))) != 0) stop("gzip ", path)
}
gz_write(c("# case\tn\tp\tvariables\tx (model.matrix, column-major)\tqr\tqraux\tpivot\trank\ttol", qr_rows),
         file.path(fx, "rcompat", "linpack_qr.tsv"))
gz_write(c("# case\trow\ty\tqr.resid(qr, y) (= RegressOutMatrix's row)", resid_rows),
         file.path(fx, "rcompat", "linpack_resid.tsv"))

# 3. percent.mt -------------------------------------------------------------------------------

pdir <- file.path(fx, "clustering", "percent_mt")
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
feat <- data.frame(
  id = sprintf("ID%02d", 1:14),
  name = c("mt-Co1", "MT-ND2", "Mt_Atp6", "mT|Cytb", "Mtor", "xmt-a", "mt", "mt", "Gene_1", "Gene_1",
           "mt-Zero", "Actb", "mtx-1", "MT-"),
  type = "Gene Expression", stringsAsFactors = FALSE)
cells <- sprintf("CELL%d-1", 1:7)
set.seed(7)
m <- matrix(rpois(14 * 7, 3), 14, 7)
m[11, ] <- 0  # mt-Zero: no counts
m[, 7] <- 0   # a cell with no mitochondrial counts
m[c(1, 12), 7] <- c(0, 9)
m[1, 6] <- 0
counts_int <- as(m, "CsparseMatrix")
mreal <- m * 1
mreal[mreal > 0] <- mreal[mreal > 0] + c(0.1, 0.2, 0.7, 0.3)[seq_len(sum(mreal > 0)) %% 4 + 1]
counts_real <- as(mreal, "CsparseMatrix")
write_10x <- function(cm, d) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  writeMM(cm, file.path(d, "matrix.mtx"))
  writeLines(cells, file.path(d, "barcodes.tsv"))
  write.table(feat, file.path(d, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  # Read10X reads Cell Ranger v3's gzipped names
  for (f in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
    if (system2("gzip", c("-n", "-9", "-f", shQuote(file.path(d, f)))) != 0) stop("gzip ", f)
  }
}
ids_file <- file.path(pdir, "mt_genes.txt")
writeLines(c("# an ID list: ID03 and ID04 (renamed by Seurat), ID11 (no counts), ID02 twice,",
             "# two IDs not among the features, and a blank line", "ID03", "ID04", "", "ID11", "ID02",
             "ID02", "NOPE1", "NOPE2"), ids_file)
pm_rows <- character()
for (kind in c("integer", "real")) {
  d <- file.path(pdir, kind)
  write_10x(if (kind == "integer") counts_int else counts_real, d)
  seu <- suppressWarnings(CreateSeuratObject(counts = Read10X(d)))
  rn <- rownames(seu)
  by_pattern <- PercentageFeatureSet(seu, pattern = "(?i)^mt-")
  # the list: the features whose id (column 1) is listed, in the counts' row order
  ids <- readLines(ids_file)
  ids <- ids[nzchar(ids) & !startsWith(ids, "#")]
  listed <- rn[feat$id %in% ids]
  by_list <- suppressWarnings(PercentageFeatureSet(seu, features = listed))
  pm_rows <- c(pm_rows,
               paste(kind, "pattern", lst(grep("(?i)^mt-", rn, value = TRUE)), hexl(unname(by_pattern)), sep = "\t"),
               paste(kind, "list", lst(listed), hexl(unname(by_list)), sep = "\t"),
               paste(kind, "nCount_RNA", "-", hexl(seu$nCount_RNA), sep = "\t"))
  if (kind == "integer") writeLines(rn, file.path(pdir, "seurat_names.txt"))
}
# make.unique(), as Read10X() applies it
mu <- list(c("a", "a", "b", "a", "a.1", "b"), c("a", "a.1", "a"), c("x", "x", "x.1", "x.1", "x"),
           c("g", "g.2", "g", "g", "g.1"), c("mt", "mt", "mt.1", "MT", "mt"))
writeLines(c("# names\tmake.unique(names)", vapply(mu, function(v) paste(lst(v), lst(make.unique(v)), sep = "\t"), "")),
           file.path(pdir, "make_unique.tsv"))
writeLines(c("# counts\tdefinition\tfeatures (Seurat's names)\tpercent.mt per cell", pm_rows),
           file.path(pdir, "expected.tsv"))
cat("Seurat's names:", lst(readLines(file.path(pdir, "seurat_names.txt"))), "\n")
cat("done\n")
