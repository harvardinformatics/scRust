#!/usr/bin/env Rscript
# Reference vectors for the pca crate (docs/PLAN.md B3b, milestone M5): the installed irlba 2.4.1
# and Seurat 5.5.1 on small synthetic inputs, for the branches the reference data never reach.
#
# - irlba.tsv: irlba(A, nv, ...) after set.seed(seed): d, u, v, iter, mprod and .Random.seed
#   afterwards, for inputs that take the tiny-problem svd(), nv = 1 (where R drops the restart's
#   Bsvd$v[, 1:k] to a vector), the work cap (work >= min(dim)), rank-deficient matrices (the
#   near-invariant-subspace branches, which draw more rnorm() vectors), clustered spectra (the
#   working-dimension expansion) and maxit reached (no convergence);
# - rowvar.tsv: Seurat's RowVar (Eigen) on rows of every length up to 9 and longer ones, with
#   offsets that make the summation order visible, and constant rows.
# Every input is a recipe Rust rebuilds exactly: rnorm() from a seed, R's %*% (OpenBLAS), and
# exact diagonals. irlba's products are OpenBLAS calls: the vectors hold for the kernel they were
# made with, the portable profile's Haswell.
#
# Reads no .rds, so it may run on a login node. From the repository root, with the production R
# env, the hermetic flags and the portable OpenBLAS profile:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 \
#     "$R_ENV/bin/Rscript" validation/R/pca_vectors.R tests/fixtures/pca
#
# Output: tab-separated files; lines starting with "#" describe the columns. Doubles are written
# as their IEEE-754 bit patterns (16 hex digits), lists comma-separated, "-" for an empty list.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path("tests", "fixtures", "pca")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

stopifnot(R.version$major == "4", R.version$minor == "5.3")
suppressPackageStartupMessages({
    library(irlba)
    library(Seurat)
})
stopifnot(packageVersion("irlba") == "2.4.1", packageVersion("Seurat") == "5.5.1")
stopifnot(Sys.getenv("OPENBLAS_CORETYPE") == "Haswell", Sys.getenv("OPENBLAS_NUM_THREADS") == "1")

hex <- function(x) {
    x <- as.double(x)
    if (!length(x)) return(character())
    b <- writeBin(x, raw(), size = 8, endian = "big")
    apply(matrix(as.character(b), nrow = 8), 2, paste, collapse = "")
}
lst <- function(v) if (length(v)) paste(v, collapse = ",") else "-"
hexl <- function(x) lst(hex(x))

write_tsv <- function(name, header, rows) {
    path <- file.path(out, name)
    writeLines(c(paste0("# ", header), unlist(rows)), path)
    cat(sprintf("%-12s %4d lines\n", name, length(unlist(rows))))
}

# ---------------------------------------------------------------------------------------------
# Inputs. recipe: "rnorm" (matrix(rnorm(m * n), m, n) after set.seed(data_seed)), "lowrank:r"
# (X %*% Y, X = matrix(rnorm(m * r), m, r) then Y = matrix(rnorm(r * n), r, n)), "diag:a:b:c"
# (an m x n matrix of zeros whose first a diagonal entries are 2 and next b are 1, plus c times
# matrix(rnorm(m * n), m, n): D + c * N), "cluster:k:gap:rest" (zeros but a diagonal of k values
# 1, 1 - gap, 1 - 2 gap, ..., then rest values 0.5 (1 - j / rest), j = 0, 1, ...: a tight
# cluster of the largest singular values).

make_input <- function(recipe, m, n, data_seed) {
    set.seed(data_seed)
    parts <- strsplit(recipe, ":", fixed = TRUE)[[1]]
    switch(parts[1],
        rnorm = matrix(rnorm(m * n), m, n),
        lowrank = {
            r <- as.integer(parts[2])
            X <- matrix(rnorm(m * r), m, r)
            Y <- matrix(rnorm(r * n), r, n)
            X %*% Y
        },
        cluster = {
            k <- as.integer(parts[2]); gap <- as.numeric(parts[3]); rest <- as.integer(parts[4])
            vals <- c(1 - (0:(k - 1)) * gap, 0.5 * (1 - (0:(rest - 1)) / rest))
            D <- matrix(0, m, n)
            for (i in seq_along(vals)) D[i, i] <- vals[i]
            D
        },
        diag = {
            a <- as.integer(parts[2]); b <- as.integer(parts[3]); c <- as.numeric(parts[4])
            D <- matrix(0, m, n)
            for (i in seq_len(a)) D[i, i] <- 2
            for (i in a + seq_len(b)) D[i, i] <- 1
            N <- matrix(rnorm(m * n), m, n)
            if (c == 0) D else D + c * N
        },
        stop("unknown recipe ", recipe))
}

# ---------------------------------------------------------------------------------------------
# irlba(A, nv, work, maxit) after set.seed(seed)
cases <- list(
    list("random", "rnorm", 60, 40, 1, 5, NA, 1000, 42),
    list("random_nv1", "rnorm", 60, 40, 1, 1, NA, 1000, 42),
    list("wide", "rnorm", 30, 70, 2, 4, NA, 1000, 7),
    list("tiny", "rnorm", 12, 5, 3, 2, NA, 1000, 42),
    list("tiny_rows", "rnorm", 5, 9, 3, 3, NA, 1000, 42),
    list("work_cap", "rnorm", 40, 20, 4, 10, NA, 1000, 42),
    list("work_small", "rnorm", 50, 45, 5, 3, 6, 1000, 11),
    list("lowrank3", "lowrank:3", 50, 30, 6, 2, NA, 1000, 42),
    list("lowrank1", "lowrank:1", 40, 25, 7, 1, NA, 1000, 5),
    list("lowrank6_nv5", "lowrank:6", 80, 60, 8, 5, NA, 1000, 42),
    list("clustered", "diag:15:15:0", 40, 30, 9, 5, NA, 1000, 42),
    list("clustered_noisy", "diag:12:10:1e-9", 60, 40, 10, 6, NA, 1000, 3),
    list("cluster_a", "cluster:10:1e-9:30", 60, 45, 14, 5, NA, 1000, 42),
    list("cluster_b", "cluster:8:1e-10:40", 70, 50, 15, 3, NA, 1000, 42),
    list("cluster_c", "cluster:12:2e-9:28", 60, 45, 16, 10, NA, 1000, 42),
    list("cluster_d", "cluster:6:1e-11:30", 50, 40, 17, 4, NA, 1000, 9),
    list("maxit1", "rnorm", 80, 60, 11, 10, NA, 1, 42),
    list("maxit2", "rnorm", 80, 60, 12, 10, NA, 2, 42),
    list("seurat_like", "rnorm", 300, 120, 13, 10, NA, 1000, 42)
)
rows <- lapply(cases, function(cs) {
    names(cs) <- c("case", "recipe", "m", "n", "data_seed", "nv", "work", "maxit", "seed")
    A <- make_input(cs$recipe, cs$m, cs$n, cs$data_seed)
    set.seed(cs$seed)
    warned <- character()
    r <- withCallingHandlers(
        if (is.na(cs$work)) irlba(A, nv = cs$nv, maxit = cs$maxit) else irlba(A, nv = cs$nv, work = cs$work, maxit = cs$maxit),
        warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
    converged <- !any(grepl("did not converge", warned))
    paste(cs$case, cs$recipe, cs$m, cs$n, cs$data_seed, cs$nv, if (is.na(cs$work)) "-" else cs$work, cs$maxit,
          cs$seed, r$iter, r$mprod, as.integer(converged), hexl(r$d), hexl(r$u), hexl(r$v),
          lst(.Random.seed), sep = "\t")
})
write_tsv("irlba.tsv", paste("case", "recipe", "m", "n", "data_seed", "nv", "work", "maxit", "seed", "iter",
                             "mprod", "converged", "d", "u", "v", "random_seed_after", sep = "\t"), rows)

# ---------------------------------------------------------------------------------------------
# Seurat's RowVar on rnorm rows (offset + scale * rnorm, row-major in the recipe's order: the
# matrix is matrix(values, nrow, ncol, byrow = TRUE)) and constant rows.
rv <- list(
    list("rnorm", 1, 3, 1, 1, 0), list("rnorm", 2, 3, 2, 1, 0), list("rnorm", 3, 3, 3, 1, 0),
    list("rnorm", 4, 3, 4, 1, 1e8), list("rnorm", 5, 3, 5, 1, 1e8), list("rnorm", 6, 3, 6, 1, 1e8),
    list("rnorm", 7, 3, 7, 1, 1e8), list("rnorm", 8, 3, 8, 1, 1e8), list("rnorm", 9, 3, 9, 1, 1e8),
    list("rnorm", 10, 2, 316, 1, 1e8), list("rnorm", 11, 2, 4144, 3, 0.5), list("rnorm", 12, 4, 1001, 1e-3, 10),
    list("const", 0.1, 2, 316), list("const", 1 / 3, 1, 17), list("const", 2^-10, 1, 316),
    list("const", 0, 1, 50), list("const", -7.25, 1, 9)
)
rows <- lapply(rv, function(cs) {
    if (cs[[1]] == "rnorm") {
        set.seed(cs[[2]])
        x <- matrix(cs[[6]] + cs[[5]] * rnorm(cs[[3]] * cs[[4]]), cs[[3]], cs[[4]], byrow = TRUE)
        recipe <- paste("rnorm", cs[[2]], hex(cs[[5]]), hex(cs[[6]]), sep = ":")
    } else {
        x <- matrix(cs[[2]], cs[[3]], cs[[4]])
        recipe <- paste("const", hex(cs[[2]]), sep = ":")
    }
    v <- Seurat:::RowVar(x)
    paste(recipe, nrow(x), ncol(x), hexl(v), hex(sum(v)), sep = "\t")
})
write_tsv("rowvar.tsv", paste("recipe", "nrow", "ncol", "rowvar", "sum", sep = "\t"), rows)
