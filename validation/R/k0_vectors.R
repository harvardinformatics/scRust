# R's part of the K0 vectors (docs/PLAN.md C2 and C5; milestone M8): the draw GroupSingletons makes
# when a singleton is equally connected to several clusters. Seurat 5.5.1's GroupSingletons calls
# set.seed(1) before each candidate cluster's connectivity, then, after the loop,
#   closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
# so every tie is broken by the first draw after set.seed(1), from a character vector of the tied
# clusters' names in candidate order. The result depends only on how many clusters tie.
#
# Reads no .rds and uses no BLAS, so it may run on a login node, with the hermetic flags:
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     "$R_ENV/bin/Rscript" validation/R/k0_vectors.R <out dir>
# Writes groupsingletons_sample.tsv: n (the number of tied names), the 1-based position sample()
# returns, the name it returns (names are "c1".."cn"), and .Random.seed[2] afterwards: the
# Mersenne-Twister's position, which counts the uniforms the draw consumed (set.seed(1) leaves it
# at 624; M1's vectors check whole states).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: k0_vectors.R <out dir>")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))

ns <- as.integer(c(1:200, 255, 256, 257, 500, 1000, 1023, 1024, 1025, 4096, 65535, 65536, 65537, 100000))
rows <- vapply(ns, function(n) {
  x <- paste0("c", seq_len(n))
  set.seed(1)
  picked <- sample(x = x, 1)
  pos <- get(".Random.seed", envir = globalenv())[2]
  paste(n, match(picked, x), picked, pos, sep = "\t")
}, character(1))
set.seed(1)
stopifnot(get(".Random.seed", envir = globalenv())[2] == 624L)
writeLines(c("# n\tposition (1-based)\tname\t.Random.seed[2] after set.seed(1); sample(x, 1) (624 after set.seed(1) alone)",
             rows), file.path(out, "groupsingletons_sample.tsv"))
cat("wrote", length(rows), "rows to", file.path(out, "groupsingletons_sample.tsv"), "\n")
