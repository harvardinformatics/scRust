#!/usr/bin/env Rscript
# Reference vectors for the rcompat crate (docs/PLAN.md A5 and B6, milestones M1 and M3; stage S0
# in validation/tolerances.toml): R's random-number stream, sampling, rnorm/qnorm, round, ^, seq,
# order, sum/mean/var/sd and number formatting (M1); as.numeric() of strings, fft, dnorm and
# approxfun (M3); quantile, median, mad, IQR, density (bw = "nrd"), nlminb, gamma functions,
# matrixStats, dnbinom (mu form) and optimize (M4); over many seeds,
# sizes and edge cases.
#
# Reads no .rds, so it may run on a login node. Run from the repository root with the production
# R env and the hermetic flags:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
#
# Output: one tab-separated file per family; lines starting with "#" describe the columns. Doubles
# are written as their IEEE-754 bit patterns (16 hex digits), lists comma-separated, "-" for an
# empty list. Results do not depend on BLAS, except nlminb.tsv (PORT's inner products and norms
# are BLAS ddot and dnrm2): run with OPENBLAS_CORETYPE=Haswell, the portable profile, as env.sh
# sets it. The files are the same on any x86-64 node with this R build and glibc (2.28 on Cannon).

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path("tests", "fixtures", "rcompat")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

stopifnot(R.version$major == "4", R.version$minor == "5.3")
stopifnot(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
stopifnot(getOption("scipen") == 0)

# IEEE-754 bit patterns, 16 hex digits each
hex <- function(x) {
    x <- as.double(x)
    if (!length(x)) return(character())
    b <- writeBin(x, raw(), size = 8, endian = "big")
    apply(matrix(as.character(b), nrow = 8), 2, paste, collapse = "")
}
unhex <- function(h) vapply(h, function(s)
    readBin(as.raw(strtoi(substring(s, seq(1, 15, 2), seq(2, 16, 2)), 16L)), "double",
            size = 8, endian = "big"), 0, USE.NAMES = FALSE)
lst <- function(v) if (length(v)) paste(v, collapse = ",") else "-"
hexl <- function(x) lst(hex(x))
# integers that may exceed .Machine$integer.max, without scientific notation
intl <- function(x) lst(formatC(as.double(x), format = "f", digits = 0))

write_tsv <- function(name, header, rows) {
    path <- file.path(out, name)
    writeLines(c(paste0("# ", header), unlist(rows)), path)
    cat(sprintf("%-16s %6d lines\n", name, length(unlist(rows))))
}

# ---------------------------------------------------------------------------------------------
# set.seed: the full .Random.seed

seeds_full <- c(0L, 1L, 42L, 12345L, 1448145L, -1L, -2147483647L, 2147483647L)
write_tsv("random_seed.tsv", "seed\t.Random.seed (626 integers)",
          lapply(seeds_full, function(s) {
              set.seed(s)
              paste(s, lst(.Random.seed), sep = "\t")
          }))

# .Random.seed checkpoints part-way through the stream (k uniforms drawn), around the point where
# the 624-word state is regenerated, and the three uniforms that follow

ckpt_cases <- list(c(1L, 1L), c(42L, 623L), c(42L, 624L), c(42L, 625L), c(12345L, 1000L),
                   c(1448145L, 1248L), c(-7L, 5000L))
write_tsv("random_seed_after.tsv",
          "seed\tk\t.Random.seed after set.seed(seed); runif(k)\tthe next three uniforms",
          lapply(ckpt_cases, function(cs) {
              set.seed(cs[1]); runif(cs[2])
              state <- .Random.seed
              paste(cs[1], cs[2], lst(state), hexl(runif(3)), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# runif streams: 10 values for many seeds; 2000 (three state regenerations) for a few

unif_cases <- c(lapply(c(0:20, 42L, 123L, 12345L, 1448145L, -1L, -100L, 987654321L,
                         2147483647L, -2147483647L), function(s) list(s, 10L)),
                lapply(c(1L, 42L, 12345L), function(s) list(s, 2000L)))
write_tsv("unif.tsv", "seed\tn\trunif(n) after set.seed(seed)",
          lapply(unif_cases, function(cs) {
              set.seed(cs[[1]])
              paste(cs[[1]], cs[[2]], hexl(runif(cs[[2]])), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# sample.int: uniform and weighted, with and without replacement. After each draw, the next three
# uniforms check that exactly R's number of uniforms was consumed. Long results are checked by
# their head, tail and a positional checksum.
#
# Weight specifications (Rust rebuilds the same doubles):
#   int:a:b:c       ((1:n * a) %% b) + c
#   unif:s          set.seed(s); runif(n)
#   mix:s:k         set.seed(s); c(runif(k), rep(1e-6, n - k))   (only k entries have n*p > 0.1)

make_prob <- function(spec, n) {
    if (spec == "-") return(NULL)
    f <- strsplit(spec, ":", fixed = TRUE)[[1]]
    switch(f[1],
        int = ((seq_len(n) * as.numeric(f[2])) %% as.numeric(f[3])) + as.numeric(f[4]),
        unif = { set.seed(as.integer(f[2])); runif(n) },
        mix = { set.seed(as.integer(f[2])); k <- as.integer(f[3]); c(runif(k), rep(1e-6, n - k)) },
        stop("unknown spec ", spec))
}

LONG <- 2000
ck1 <- function(r) sum((as.double(r) * seq_along(r)) %% 1000003)

sample_cases <- list()
add_sample <- function(seed, n, size, replace, prob = "-")
    sample_cases[[length(sample_cases) + 1]] <<- list(seed, n, size, replace, prob)

# uniform, with replacement: R_unif_index across bit widths (2^15, 2^16 and their neighbours)
i <- 0
for (n in c(1, 2, 3, 10, 255, 256, 257, 32767, 32768, 32769, 65535, 65536, 65537, 1e5,
            2^24 + 1, 2147483647, 2147483648, 3e9, 2^40 + 3, 4.5e15)) {
    i <- i + 1
    add_sample(1000L + i, n, if (n == 1e5) 1000 else 50, TRUE)
}
add_sample(12345L, 11831, 3944, TRUE)           # doublet parents: n.real = 11831, pN = 0.25
# uniform, without replacement (partial Fisher-Yates)
for (cs in list(c(1, 1), c(2, 2), c(5, 5), c(10, 3), c(100, 100), c(1000, 999), c(1000, 1),
                c(32769, 100), c(65537, 1000), c(1e5, 1e5), c(0, 0), c(1e7, 10),
                c(1.2e7, 6000001))) {
    i <- i + 1
    add_sample(2000L + i, cs[1], cs[2], FALSE)
}
add_sample(12345L, 11831, 10000, FALSE)         # paramSweep's subsample
add_sample(1448145L, 11831, 5000, FALSE)        # SCTransform's fitting cells
# uniform, without replacement, n > 1e7 and size <= n/2: the hash-based sample2
for (cs in list(c(1e7 + 1, 5), c(2e7, 1000), c(2147483647, 100), c(3e9, 100), c(2e7, 2e5))) {
    i <- i + 1
    add_sample(3000L + i, cs[1], cs[2], FALSE)
}
# weighted, with replacement: inversion (<= 200 entries with n*p > 0.1) or Walker's alias method
for (cs in list(list(1, 5, "int:1:7:1"), list(2, 20, "int:1:3:1"), list(5, 30, "int:3:7:1"),
                list(50, 100, "unif:7"), list(200, 300, "int:1:1:1"), list(1000, 500, "mix:9:150"),
                list(201, 500, "int:1:1:1"), list(250, 1000, "unif:11"), list(1000, 2000, "int:7:13:0"),
                list(5000, 3000, "unif:12"), list(20000, 100, "unif:13"), list(10, 0, "unif:3"))) {
    i <- i + 1
    add_sample(4000L + i, cs[[1]], cs[[2]], TRUE, cs[[3]])
}
# weighted, without replacement (ProbSampleNoReplace with revsort)
for (cs in list(list(1, 1, "int:1:7:1"), list(2, 2, "int:1:3:1"), list(10, 3, "int:3:7:1"),
                list(10, 10, "int:3:7:1"), list(100, 50, "int:7:13:0"), list(1000, 200, "unif:14"),
                list(15000, 2000, "unif:15"), list(5, 0, "unif:16"))) {
    i <- i + 1
    add_sample(5000L + i, cs[[1]], cs[[2]], FALSE, cs[[3]])
}

# many seeds: small draws of every kind
sizes_n <- c(1, 2, 3, 7, 100, 1000, 32768, 32769, 65536, 1e5, 2^31)
for (k in 1:300) {
    n <- sizes_n[(k %% length(sizes_n)) + 1]
    replace <- k %% 2 == 0
    size <- if (replace) (k %% 9) + 1 else min(n, (k %% 9) + 1)
    prob <- if (k %% 3 == 0 && n <= 1e5) c("int:7:13:1", sprintf("unif:%d", 100 + k))[(k %% 2) + 1] else "-"
    add_sample(7919L * k, n, size, replace, prob)
}

write_tsv("sample.tsv",
          paste("seed\tn\tsize\treplace\tprob\tlength\tchecksum\thead (all if length <= 2000)",
                "\ttail (last 50 if longer)\tnext three uniforms (results are R's, 1-based)"),
          lapply(sample_cases, function(cs) {
              prob <- make_prob(cs[[5]], cs[[2]])
              set.seed(cs[[1]])
              r <- if (is.null(prob)) sample.int(cs[[2]], cs[[3]], cs[[4]])
                   else sample.int(cs[[2]], cs[[3]], cs[[4]], prob)
              after <- runif(3)
              long <- length(r) > LONG
              paste(cs[[1]], intl(cs[[2]]), intl(cs[[3]]), as.integer(cs[[4]]), cs[[5]],
                    length(r), if (long) intl(ck1(r)) else "-",
                    intl(if (long) head(r, 50) else r), if (long) intl(tail(r, 50)) else "-",
                    hexl(after), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# rnorm by inversion, and qnorm

rnorm_cases <- c(list(list(42L, 3000L, 0, 1), list(1L, 10L, 0, 1), list(12345L, 100L, 5, 2),
                      list(7L, 5L, 0, 0), list(8L, 5L, Inf, 1), list(9L, 200L, -3, 1e-3)),
                 lapply(1:50, function(k) list(104729L * k, 5L, 0, 1)))
write_tsv("rnorm.tsv", "seed\tn\tmean\tsd\trnorm(n, mean, sd)\tnext three uniforms",
          lapply(rnorm_cases, function(cs) {
              set.seed(cs[[1]])
              x <- rnorm(cs[[2]], cs[[3]], cs[[4]])
              paste(cs[[1]], cs[[2]], hex(cs[[3]]), hex(cs[[4]]), hexl(x), hexl(runif(3)), sep = "\t")
          }))

set.seed(20261110)
p_fixed <- c(0, 1, 0.5, 1e-300, 1e-100, 1e-50, 1e-20, 1e-16, 1e-12, exp(-25), 1e-10, 1e-5, 0.001,
             0.01, 0.05, 0.074, 0.075, 0.0751, 0.25, 0.425, 0.575, 0.9, 0.924, 0.925, 0.926, 0.95,
             0.99, 0.999999, 1 - 1e-12, 1 - 2^-53, 2^-1074, .Machine$double.xmin, -0.1, 1.1, NaN)
p_rand <- c(runif(1000), 10^runif(200, -300, 0))
lp <- -c(0, 1e-300, 1e-20, 1e-5, 0.1, 0.5, 0.69, 0.7, 1, 2, 5, 10, 24, 25, 26, 100, 700, 728, 730,
         1000, 3000, 5000, 12000, 20000, 1e5, 7e5, 1e6, 1e8, 1e9, 1.3e9, 1e10, 1e12, 1e15, 1e17,
         4e17, 4.1e17, 5e17, 1e20, 1e100, 1e300, Inf)
lp_rand <- -10^runif(100, -3, 20)
qcase <- function(p, lower, logp, mu = 0, sigma = 1)
    paste(hex(p), as.integer(lower), as.integer(logp), hex(mu), hex(sigma),
          hex(suppressWarnings(qnorm(p, mu, sigma, lower, logp))), sep = "\t")
write_tsv("qnorm.tsv", "p\tlower.tail\tlog.p\tmean\tsd\tqnorm(p, mean, sd, lower.tail, log.p)",
          c(lapply(c(p_fixed, p_rand), qcase, TRUE, FALSE),
            lapply(c(p_fixed, p_rand[1:100]), qcase, FALSE, FALSE),
            lapply(c(lp, lp_rand), qcase, TRUE, TRUE),
            lapply(c(lp, lp_rand), qcase, FALSE, TRUE),
            lapply(p_rand[1:50], qcase, TRUE, FALSE, 3, 2),
            list(qcase(0.3, TRUE, FALSE, 1, 0), qcase(0.3, TRUE, FALSE, 1, -1))))

# ---------------------------------------------------------------------------------------------
# round(x) and round(x, digits); x ^ y

set.seed(1)
x_halves <- c((-6:6) + 0.5, 4503599627370495.5, 2^52, 2^53 + 2, 2.675, 1.005, 0.125, 0.135, 0.285,
              1.0005, 2.345, 5.015, 0.0005, 0.15, 0.25, 0.35, 123.455, -2.675, 1e300, 5e-324, -0,
              1e-320, NaN, NA, Inf, -Inf)
digits_all <- c(0, -2, -1, 1, 2, 3, 5, 10, 15, 16, 20, 300, 309, 310, 400, -308, -309, -400, 1.4,
                2.5, -0.5)
x_rand <- c(runif(300, -1000, 1000), 10^runif(100, -10, 10), round(runif(100, 0, 1e4)) / 1000 + 5e-4)
round_rows <- c(
    unlist(lapply(x_halves, function(x)
        lapply(digits_all, function(d) paste(hex(x), hex(d), hex(round(x, d)), sep = "\t")))),
    unlist(lapply(x_rand, function(x) {
        d <- sample(c(0, 0, 1, 2, 3, 4, -1, 6), 1)
        paste(hex(x), hex(d), hex(round(x, d)), sep = "\t")
    })))
write_tsv("round.tsv", "x\tdigits\tround(x, digits)", round_rows)

set.seed(2)
pow_x <- c(0, -0, 1, -1, 2, -2, 0.5, 10, Inf, -Inf, NaN, NA)
pow_y <- c(0, -0, 1, -1, 2, 3, 0.5, -0.5, 2.5, 1/3, Inf, -Inf, NaN, NA, 1e300, -1e300)
pow_rows <- c(
    unlist(lapply(pow_x, function(x) lapply(pow_y, function(y)
        paste(hex(x), hex(y), hex(x^y), sep = "\t")))),
    mapply(function(x, y) paste(hex(x), hex(y), hex(x^y), sep = "\t"),
           c(runif(500, 0, 10), -runif(100, 0, 3), 10 + 0 * 1:100, runif(100)),
           c(runif(500, -20, 20), round(runif(100, -10, 10)), runif(100, -320, 308), rep(c(2, 3, 0.5, 1/3, 4), 20))))
write_tsv("pow.tsv", "x\ty\tx^y", pow_rows)

# ---------------------------------------------------------------------------------------------
# seq(from, to, by) and seq(from, to, length.out)

seq_err <- function(expr) tryCatch(hexl(expr), error = function(e) "error")
set.seed(3)
by_cases <- list(c(0.01, 0.3, 0.01), c(0.05, 0.3, 0.05), c(0, 1, 0.1), c(1, 0, -0.1),
                 c(0.1, 0.7, 0.2), c(1e6, 1e6 + 1, 0.1), c(-5, 5, 0.25), c(1, 2, 3), c(1, 1, 0),
                 c(0, 0, 1), c(5, 5, 1), c(1, 1 + 1e-15, 1e-16), c(0, 100, 1), c(0, 1, 1/3),
                 c(0.3, 0.1, -0.1), c(-1e308, 1e308, 1e307), c(0, 1, -0.1), c(0, 1, 0),
                 c(0, 1, 1e-12), c(0, -1, -0.07))
for (k in 1:30) {
    from <- runif(1, -10, 10); by <- runif(1, 0.001, 1) * sample(c(-1, 1), 1)
    by_cases[[length(by_cases) + 1]] <- c(from, from + by * runif(1, 1, 50), by)
}
write_tsv("seq_by.tsv", "from\tto\tby\tseq(from, to, by) (\"error\" if R stops)",
          lapply(by_cases, function(cs)
              paste(hex(cs[1]), hex(cs[2]), hex(cs[3]), seq_err(seq(cs[1], cs[2], by = cs[3])), sep = "\t")))

len_cases <- list(c(0, 1, 401), c(-3.2, 7.9, 401), c(1, 1, 5), c(0, 1, 1), c(0, 1, 2), c(0, 1, 0),
                  c(0, 1, 2.5), c(-1e308, 1e308, 5), c(5, -5, 11), c(0.0123, 0.9876, 300))
for (k in 1:20) {
    a <- runif(1, -100, 100)
    len_cases[[length(len_cases) + 1]] <- c(a, a + runif(1, -50, 50), sample(3:1000, 1))
}
write_tsv("seq_len.tsv", "from\tto\tlength.out\tseq(from, to, length.out = length.out)",
          lapply(len_cases, function(cs)
              paste(hex(cs[1]), hex(cs[2]), hex(cs[3]), seq_err(seq(cs[1], cs[2], length.out = cs[3])), sep = "\t")))

# ---------------------------------------------------------------------------------------------
# order(x) and order(x, decreasing = TRUE)

set.seed(4)
order_vecs <- list(
    c(2, 1, 2, NaN, -0, 0, 1, NA, Inf, -Inf, 3, 1),
    5,
    numeric(0),
    round(runif(150) * 10) / 10,                         # < 200 elements: insertion sort
    round(runif(3000) * 50) / 7,                         # radix sort, many ties
    ((1:2500 * 7) %% 31) / 30,                           # pANN-like multiples of 1/k
    rev(sort(runif(500))),                               # strictly decreasing
    sort(runif(500)),                                    # strictly increasing
    rep(1.5, 300),
    c(runif(200), rep(NA, 50), rep(NaN, 50), rep(c(0, -0), 50), rep(0.5, 50)),
    c(-Inf, Inf, 1e308, -1e308, 5e-324, -5e-324, 0, -0))
write_tsv("order.tsv", "decreasing\tx\torder(x, decreasing = decreasing) (1-based)",
          unlist(lapply(order_vecs, function(x) lapply(c(FALSE, TRUE), function(dec)
              paste(as.integer(dec), hexl(x), lst(order(x, decreasing = dec)), sep = "\t")))))

# ---------------------------------------------------------------------------------------------
# sum, mean, var, sd (each also with na.rm = TRUE)

set.seed(5)
summary_vecs <- list(
    c(1, rep(1e-16, 10)), c(1e308, 1e308, -1e308), c(1e308, 1e308), c(1.7e308, 1.7e308, -1.7e308, 1e308),
    c(-Inf, Inf), c(1, Inf), c(NA, 1), c(NaN, 1), c(NaN, NA), c(NA, NaN), c(NaN, NaN, 2),
    numeric(0), 5, c(3, 3, 3), c(1e16, 1, -1e16), c(0.1, 0.2, 0.3), c(5e-324, 5e-324, -5e-324),
    1e9 + runif(1000), runif(1000) * 1e6, ((1:500 * 7) %% 31) / 30, rep(1 / 3, 1000), (1:1000) / 7,
    10^runif(200, -300, 300), rnorm(777, 1e5, 1e-3), c(runif(100), NA, runif(10)),
    c(2, 4, 4, 4, 5, 5, 7, 9))
summ <- function(x) hex(c(sum(x), mean(x), var(x), sd(x),
                          sum(x, na.rm = TRUE), mean(x, na.rm = TRUE), var(x, na.rm = TRUE),
                          sd(x, na.rm = TRUE)))
write_tsv("summary.tsv",
          paste("x\tsum\tmean\tvar\tsd\tsum(na.rm)\tmean(na.rm)\tvar(na.rm)\tsd(na.rm);",
                "x is a hex list or unif:seed:n:a:b (set.seed(seed); runif(n, a, b))"),
          c(lapply(summary_vecs, function(x) paste(c(hexl(x), suppressWarnings(summ(x))), collapse = "\t")),
            lapply(list(c(11L, 1e5, 0, 1), c(12L, 1e5, 1e9, 1e9 + 1), c(13L, 54321, -1e-3, 1e-3)),
                   function(cs) {
                       set.seed(cs[1]); x <- runif(cs[2], cs[3], cs[4])
                       paste(c(sprintf("unif:%d:%d:%s:%s", cs[1], cs[2], hex(cs[3]), hex(cs[4])), summ(x)),
                             collapse = "\t")
                   })))

# ---------------------------------------------------------------------------------------------
# as.character (as paste writes numbers)

set.seed(6)
pk_grid <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
fmt_named <- c(pk_grid, seq(0.05, 0.3, by = 0.05), 0, -0, 1, -1, 0.1, 1/3, 2/3, 1e5, 1e-5, 123456,
               1234567, 1e15, 1e16, 1e-15, 100000, 99999, 0.0001, 0.00012, 1e21, 1e22, 1.5e300,
               5e-324, .Machine$double.xmax, .Machine$double.xmin, NaN, NA, Inf, -Inf, 0.1 + 0.2,
               99999.99999999999, 9.999999999999998, 9.99999999999999, 99999999999999.99,
               999999999999999.9, 1e-300, 2^60, 123456789.123456789, -2.5e-7, 1234.5678, 1e-13,
               9.9e-14, 1e42, 9.99e41, 1e41, 0.1234, 100, 1e-4, 1e-3 * 3, 314159.26535, 2^-1074 * 3,
               0.95, 1/0.95 - 1, 11831 / 0.95 - 11831)
n_r <- 1500
fmt_rand <- c(
    (runif(n_r) * 9 + 1) * 10^round(runif(n_r, -320, 308)),   # all magnitudes
    runif(500, -1e6, 1e6),
    round(runif(500) * 1e4) / 10^sample(0:8, 500, replace = TRUE),   # short decimals
    round(runif(200) * 1e15) / 1e15,
    (round(runif(200) * 1e14) + 0.5) / 10^sample(10:20, 200, replace = TRUE))
# Numbers whose string depends on glibc's powl(10, k) being one ulp off for some k: here the
# format of |x| < 1e-13 or >= 1e42 scales by powl, and a correctly rounded 10^k would give
# different digits. Found by searching near rounding ties in rcompat for those k.
fmt_powl <- unhex(c(
    "03516acc7153d1d0", "0674209494806ef0", "0674a010b87438e8", "07b73e9039e221d6",
    "07b94d05fa201b13", "08509e2bc672c583", "0854760a121a614f", "09c53699e3b2cb9d",
    "09c40fde14d6fcca", "0a95e8c6def0e6f3", "0a957b585fe9c6e3", "0d471609891587a9",
    "0d53500210b38383", "0e511921a38bbe41", "0e5aecf4cf1a3322", "0ef2923e5d136f20",
    "0ef48ccf5cd0ee33", "0fcbef1eef19d231", "0fc5abdaf7e2fda5", "12126ed1927cd53e",
    "1216a97ffc8c71b6", "1495019516233497", "1494738f91792c7d", "16d71f20dfd4bbdd",
    "16d91a6909dd0b79", "1a2cbd20944efc6d", "1a2fab8cdc26523a", "1a998fad8d07303a",
    "1a997010a4578dbe", "1c1172e6d2620d91", "1c0a6ea4d68b2fb6", "1c733d5fdb9fbfe6",
    "1c746ba310692f7e", "1d4f6d5bb69c6f31", "1d4a11747bde7989", "1eecc6a3fc64df2f",
    "1ef75bfc92ae1a95", "266f37763972609d", "266b2fbde5a18155", "2ca231f36162b4d2",
    "2ca232106570563d", "2e1f5ea1fa854ed8", "2e1e7c92da1dcd03", "327c8d679185f625",
    "327b15e2f25bcb72", "35cf35fe040ff332", "35cf3ed614fb182e", "36389205f5284aa9",
    "3634b1c10bc93ebd", "3b29bb671f93e678", "3b30a8c649aace3c", "4bd11c209a152438",
    "4bc8e71bbb628f10", "516aeeec5325a9ed", "516c27dde385b14c", "520596aae626fba6",
    "520c58d140b68a75", "569a160792c26f2b", "569a4d505b249b4e", "5878299b952a7488",
    "587cfaa11f02c5c9", "5e0e1c0bde1f7584", "5e15e3028da04098", "5f20ac4dc3c80da1",
    "5f15e41b3a4716a6", "6021ffb5ed37a780", "6022070b09cbbb15", "60949ddf3ce4423a",
    "60916ccae7c38fac", "6267c6cd7b556c67", "626f48673ca83b4a", "63b24540a4d72022",
    "63ab9800b8e2a756", "64f2445fec0c7ca0", "64ec424d0f0674c3", "676b26a2a8819ada",
    "67707f6fdbc00e29", "694b04b97361c24c", "694b58c8cb593b94", "6ac0e300c7d522ff",
    "6abd70731672e741", "6c630df4c6d2a1c1", "6c62cd2aef8b8384", "6fbd48cfb28b861f",
    "6fb97378e1a44ebb", "70f625c548b78a06", "70f6030050603e79", "7629b96f2d630660",
    "762b0d91d0c5f217", "76f6502d17618cbb", "77011d369c6a28a7", "77cca322e6d3c2e3",
    "77d4e11af8cdbe0f", "7840a283d8158d58", "784150c7d64cc112", "78a36b8dd56e56e2",
    "78a2b5b6fb85e899", "7941d827d2c78dc2", "7943dfd2640ce3ac", "7bf5451b12586241",
    "7d64ffb2f4c06bcd", "7d64906100b58ba7", "7fe867cf1ffb547e", "7fe4fcc4b5984a40"
))
write_tsv("format.tsv", "x\tas.character(x)",
          lapply(c(fmt_named, fmt_rand, fmt_powl), function(x) paste(hex(x), as.character(x), sep = "\t")))

# ---------------------------------------------------------------------------------------------
# as.numeric() of strings: R's own decimal parser R_strtod (milestone M3). It is not correctly
# rounded: digits accumulate in a long double, a power of ten is built in long double, and the
# quotient is rounded to double. DoubletFinder parses pK and pN back from paste()d strings
# (summarizeSweep, find.pK's as.numeric(as.character(...))). Columns: the string, as.numeric() of
# it, and as.character() of that value. Strings contain no tab or newline.

set.seed(7)
pn_grid <- seq(0.05, 0.3, by = 0.05)
strtod_named <- c(
    as.character(pk_grid), as.character(pn_grid),             # paste() of paramSweep's grids
    sprintf("%.17g", pk_grid), sprintf("%.17g", pn_grid),
    as.character(seq(0.0001, 0.5, by = 0.0001)),              # every 4-decimal label below 0.5
    "0", "-0", "+0", "00", "0.", ".0", ".5", "5.", "-.5", "+.25e1", "000123.4500", "1e5", "1E+05",
    "1e-05", "5e-04", "1e-4", "0.0005", "0.1", "0.2", "0.3", "0.7", "0.29", "0.95", "1.05", "3.14159",
    " 0.25", "0.25 ", "  -1.5  ", "1e", "1e+", "e5", "abc", "1.2.3", "--1", "+-1", "", " ", "0x",
    "0x1p-2", "0x1.8p1", "0X1A", "0xA.8", "0x.8", "0x1P+10", "0x1p", "0xg", "0x1.fffffffffffffp1023",
    "0x1p-1074", "0x1.8p-1074", "0x10000000000000000001", "0x1.000000000000000001p0",
    "Inf", "inf", "-Inf", "infinity", "-INFINITY", "infinit", "NaN", "nan", "-NaN", "NA", "na",
    "1.7976931348623157e308", "1.7976931348623158e308", "1.797693134862316e308", "1e308", "9e307",
    "1e309", "-1e400", "1e99999", "1e-99999", "0e99999", "0e-99999", "1e-310", "123e-400",
    "4.9406564584124654e-324", "2.4703282292062328e-324", "2.4703282292062327e-324", "1e-323",
    "2.2250738585072014e-308", "2.2250738585072011e-308", "1e-300", "1e-307", "1e-308",
    paste0("0.", strrep("0", 320), "123"), paste0("1", strrep("0", 30)), paste0(strrep("9", 40), ".5"),
    "0.1000000000000000055511151231257827021181583404541015625",
    "9007199254740993", "9007199254740992.5", "18446744073709551617", "123456789012345678901234567890")
n_s <- 2500
rand_mag <- function(n) (runif(n) * 9 + 1) * 10^round(runif(n, -330, 310))
strtod_rand <- c(
    sprintf("%.*g", sample(1:20, n_s, replace = TRUE), rand_mag(n_s)),  # 1-20 significant digits
    sprintf("%.17g", runif(1000)),          # round-trip text, as the reference files write it
    sprintf("%.15g", runif(500)),           # as.character's precision
    sprintf("%.*g", sample(16:19, 1500, replace = TRUE), runif(1500) * 10^sample(-5:5, 1500, replace = TRUE)),
    sprintf("%.*f", sample(0:12, 1000, replace = TRUE), runif(1000) * 10^sample(-3:6, 1000, replace = TRUE)),
    sprintf("%.40f", runif(50)),
    sprintf("%.*e", sample(0:25, 500, replace = TRUE), rand_mag(500)),
    as.character(round(runif(500), sample(1:6, 500, replace = TRUE))))
# Strings R does not round correctly (the long-double power of ten is inexact for exponents
# beyond 27, and the quotient is rounded twice): found by parsing 1.6 million random strings
# (set.seed(20261111); "%.17g" of runif(n) and of runif(n) * 10^k for k in -30:30, "%.15g" to
# "%.20g" of all magnitudes, "%.16e" to "%.22e" of runif(n); n = 4e5 each) and keeping those whose
# value differs from the correctly rounded double. Within [0, 1] none did.
strtod_misrounded <- c(
    "5.25528057198971e+210", "1.71470329165459e+249", "2.242252947296947e+270",
    "6.860132599016651e-270", "4.510290635284036e-246", "3.502626012312248e+195",
    "4.92278569051996e+168", "7.34419400524348e-289", "5.62239264859818e-124",
    "4.65777662117034e+283", "2.085569987306371e+188", "5.118861177237704e+293",
    "2.086596223525703e+194", "5.65985701070167e+254", "4.698357041692361e+263",
    "7.87398234778084e-176", "5.019062082748861e+276", "7.473196343984456e-295",
    "5.522084330441431e-273", "8.70618498325348e-274", "3.83974621864036e+282", "23361.110619735",
    "5.44503967859782e+298", "1.72871372662485e+120", "5.97657346026972e+65",
    "7.76994571927935e-186", "3.5953846455086e-176", "7.37715537007898e-121",
    "2.99026781436987e-258", "1.50350060709752e-307", "7.6716720117256e-152",
    "7.465276727918535e-283", "7.83056729123928e-259", "6.55742067098618e-298",
    "1.22114257398061e-297", "3.362189100123942e+127", "6.45266169845127e-304",
    "7.61927206208929e-272", "1.74368128064089e-255", "4.46175417094491e+21",
    "6.75499872490764e+202", "5.980723627610132e+184", "9.66649661841802e+228",
    "1.730707140173763e+172", "7.373004509136082e+248", "3.769152639899403e+213",
    "4.19709940254688e+284", "6.86438280856237e-143", "8.17883213260211e-143",
    "2.70772597659379e-244", "1.583063934929669e-154", "2.526549922069535e-253",
    "4.07147940271534e-219", "5.79826105781831e-211", "1.76951495371759e-282",
    "2.641541251679882e-252", "5.82878880831413e+269", "4.5610746210441e+241",
    "4.18201254680753e+173", "8.59632390853949e-266", "3.197703327983618e+44",
    "4.92872675624676e-304", "7.392714625690132e-202", "3.676960405427962e-98",
    "5.132704433053732e-90")
strtod_strings <- c(strtod_named, strtod_rand, strtod_misrounded)
stopifnot(!any(grepl("[\t\n\r]", strtod_strings)))
write_tsv("strtod.tsv", "string\tas.numeric(string)\tas.character(as.numeric(string))",
          lapply(strtod_strings, function(s) {
              v <- suppressWarnings(as.numeric(s))
              paste(s, hex(v), as.character(v), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# fft(z) and fft(z, inverse = TRUE) (Singleton's mixed-radix algorithm; M3: KernSmooth::bkde
# convolves through it at lengths 512 and 1024). Input: set.seed(n); complex(real = runif(n),
# imaginary = runif(n) - 0.5), which Rust regenerates. Output: real and imaginary parts
# interleaved.

fft_n <- c(1:32, 36, 45, 48, 49, 60, 64, 81, 97, 100, 105, 121, 125, 128, 169, 210, 243, 256, 289,
           343, 360, 401, 512, 1024, 1155)
write_tsv("fft.tsv", "n\tinverse\tfft(z, inverse) as re,im,re,im,... (z from set.seed(n))",
          unlist(lapply(fft_n, function(n) {
              set.seed(n)
              z <- complex(real = runif(n), imaginary = runif(n) - 0.5)
              vapply(c(FALSE, TRUE), function(inv) {
                  f <- fft(z, inverse = inv)
                  paste(n, as.integer(inv), hexl(as.vector(rbind(Re(f), Im(f)))), sep = "\t")
              }, "")
          })))

# ---------------------------------------------------------------------------------------------
# dnorm(x, mean, sd, log) (bkde's normal kernel weights; M3)

set.seed(8)
dn_x <- c(seq(0, 4, by = 1 / 64), runif(300, 0, 5), runif(100, 4.9, 40), -runif(50, 0, 6),
          0, -0, 5, 4.999999999999999, 38.5, 38.6, 38.58601, 40, 1e-300, 1e300, Inf, -Inf, NaN, NA)
dn_cases <- c(lapply(dn_x, function(x) list(x, 0, 1, FALSE)),
              lapply(c(0.5, 3, 7, 50), function(x) list(x, 0, 1, TRUE)),
              list(list(1, 2, 0.5, FALSE), list(1, 1, 0, FALSE), list(1, 2, 0, FALSE),
                   list(1, 0, -1, FALSE), list(1, 0, Inf, FALSE), list(Inf, Inf, 1, FALSE),
                   list(0.3, -1.25, 2.5, FALSE), list(0.3, -1.25, 2.5, TRUE)))
write_tsv("dnorm.tsv", "x\tmean\tsd\tlog\tdnorm(x, mean, sd, log)",
          lapply(dn_cases, function(cs) paste(hex(cs[[1]]), hex(cs[[2]]), hex(cs[[3]]), as.integer(cs[[4]]),
                                               hex(suppressWarnings(dnorm(cs[[1]], cs[[2]], cs[[3]], cs[[4]]))),
                                               sep = "\t")))

# ---------------------------------------------------------------------------------------------
# approxfun(x, y, method, rule, f)(xout): the grid of bkde (seq(a, b, length.out = 401)) and
# irregular, unsorted and tied abscissae (M3). x, y and xout are written out.

set.seed(9)
ap_cases <- list()
add_ap <- function(x, y, xout, method = "linear", rule = 1, f = 0)
    ap_cases[[length(ap_cases) + 1]] <<- list(x, y, xout, method, rule, f)
g <- seq(-0.123, 1.456, length.out = 401)
add_ap(g, dnorm(g, 0.6, 0.2), seq(min(g) + 0.5, max(g) - 0.3, length.out = 777))
add_ap(g, runif(401), c(g[c(1, 2, 200, 401)], seq(min(g), max(g), length.out = 300), min(g) - 1, max(g) + 1, NA, NaN))
xi <- sort(runif(50, -3, 3))
add_ap(xi, runif(50), runif(200, -3.5, 3.5))
add_ap(xi, runif(50), runif(200, -3.5, 3.5), rule = 2)
add_ap(xi, runif(50), runif(200, -3.5, 3.5), method = "constant", f = 0.25)
add_ap(xi, runif(50), runif(200, -3.5, 3.5), method = "constant", f = 1, rule = 2)
add_ap(c(3, 1, 2, 2, 5, 3, 3), c(1, 2, 3, 4, 5, 6, 7.5), seq(0.5, 5.5, by = 0.25))
add_ap(c(1, NA, 2, 4), c(1, 5, NaN, 2), seq(0, 5, by = 0.5))
write_tsv("approx.tsv", "x\ty\txout\tmethod\trule\tf\tapproxfun(x, y, method, rule = rule, f = f)(xout)",
          lapply(ap_cases, function(cs) {
              fun <- suppressWarnings(approxfun(cs[[1]], cs[[2]], method = cs[[4]], rule = cs[[5]], f = cs[[6]]))
              paste(hexl(cs[[1]]), hexl(cs[[2]]), hexl(cs[[3]]), cs[[4]], cs[[5]], hex(cs[[6]]),
                    hexl(fun(cs[[3]])), sep = "\t")
          }))

# =============================================================================================
# M4: the stats functions SCTransform v2 reaches (sctransform 0.4.3, glmGamPoi 1.22.0)
# =============================================================================================

# ---------------------------------------------------------------------------------------------
# quantile(x, c(0.25, 0.75), type = 7), median, mad, IQR: odd and even lengths, ties, negative
# values; the vectors are written out.

set.seed(41)
os_cases <- list(c(3, 1, 2), c(4, 1, 3, 2), c(1, 1, 1, 1), c(-2.5, 7, 0, 0, 0.5),
                 runif(1), runif(2), runif(101, -5, 5), round(rnorm(500), 1), rexp(1000),
                 log10(1 + rpois(777, 3)), c(rep(0, 50), runif(51)))
write_tsv("orderstat.tsv", "x\tquantile(x, c(0.25, 0.75), type = 7)\tmedian(x)\tmad(x)\tIQR(x)",
          lapply(os_cases, function(x)
              paste(hexl(x), hexl(quantile(x, c(0.25, 0.75), names = FALSE)), hex(median(x)),
                    hex(mad(x)), hex(IQR(x)), sep = "\t")))

# ---------------------------------------------------------------------------------------------
# density(x, bw = "nrd", adjust = 1)$x, $y and $bw (512 points; R >= 4.4 coordinates), as
# sctransform weights its gene sampling: log10 geometric means and other shapes.

set.seed(43)
de_cases <- list(log10(1 + rexp(800, 2)), rnorm(300), c(runif(200), runif(30, 5, 6)),
                 log10(1 + rpois(1200, 0.3) + runif(1200)), c(1, 2), seq(0, 1, length.out = 7))
write_tsv("density.tsv", "x\tbw.nrd(x)\tdensity(x, bw = \"nrd\")$x\t$y",
          lapply(de_cases, function(x) {
              d <- density(x, bw = "nrd", adjust = 1)
              stopifnot(identical(d$bw, bw.nrd(x)))
              paste(hexl(x), hex(d$bw), hexl(d$x), hexl(d$y), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# nlminb(start, objective, gradient, hessian, lower, upper, control = list(iter.max = m)): R's
# PORT library (DRMNHB; DRMNGB without the Hessian, glmGamPoi's fallback), as glmGamPoi calls
# it, on analytic test functions in one and two
# dimensions whose arithmetic Rust repeats operation for operation (always with finite bounds:
# without them R calls DRMNH, which glmGamPoi never reaches): the result (par, objective,
# convergence, iterations, evaluations, the PORT return code iv[1]) and the trace of every point
# at which the objective ("f") or the gradient and Hessian ("g") were evaluated.
#   quad:a:c        a * ((x - c) * (x - c))
#   quart:a:c       (x - c)^4 + a * (x - c)^2, written as products
#   explin:b        exp(x) - b * x
#   flat            exp(-x) (no minimum inside large bounds)
#   rosen           100 * (x2 - x1 * x1)^2 + (1 - x1)^2 (two dimensions)
#   nbll:s:m        a negative binomial profile in log(theta): s * log(1 + m * exp(x)) - ...

nl_fun <- function(spec) {
    f <- strsplit(spec, ":", fixed = TRUE)[[1]]
    num <- as.numeric(f[-1])
    switch(f[1],
        quad = list(
            fn = function(x) num[1] * ((x - num[2]) * (x - num[2])),
            gr = function(x) 2 * num[1] * (x - num[2]),
            he = function(x) matrix(2 * num[1], 1, 1)),
        quart = list(
            fn = function(x) { d <- x - num[2]; d2 <- d * d; d2 * d2 + num[1] * d2 },
            gr = function(x) { d <- x - num[2]; 4 * (d * d * d) + 2 * num[1] * d },
            he = function(x) { d <- x - num[2]; matrix(12 * (d * d) + 2 * num[1], 1, 1) }),
        explin = list(
            fn = function(x) exp(x) - num[1] * x,
            gr = function(x) exp(x) - num[1],
            he = function(x) matrix(exp(x), 1, 1)),
        flat = list(
            fn = function(x) exp(-x),
            gr = function(x) -exp(-x),
            he = function(x) matrix(exp(-x), 1, 1)),
        rosen = list(
            fn = function(x) { a <- x[2] - x[1] * x[1]; b <- 1 - x[1]; 100 * (a * a) + b * b },
            gr = function(x) { a <- x[2] - x[1] * x[1]; b <- 1 - x[1];
                               c(-400 * x[1] * a - 2 * b, 200 * a) },
            he = function(x) { a <- x[2] - x[1] * x[1];
                               h11 <- 1200 * (x[1] * x[1]) - 400 * x[2] + 2
                               matrix(c(h11, -400 * x[1], -400 * x[1], 200), 2, 2) }),
        nbll = list(
            # theta = exp(x): -( s * x - (s + m) * log(1 + exp(x)) ) style concave profile
            fn = function(x) (num[1] + num[2]) * log(1 + exp(x)) - num[1] * x,
            gr = function(x) (num[1] + num[2]) * (exp(x) / (1 + exp(x))) - num[1],
            he = function(x) { e <- exp(x); matrix((num[1] + num[2]) * (e / ((1 + e) * (1 + e))), 1, 1) }),
        stop("unknown function ", spec))
}

nl_cases <- list(
    list("quad:1:3", 0, -10, 10, 200), list("quad:1:3", 0, -10, 1, 200),
    list("quad:2.5:-1.25", 4, -1e300, 1e300, 200), list("quad:1e-3:7", -5, -20, 20, 200),
    list("quart:0:2", -1, -10, 10, 200), list("quart:0.01:0.5", 3, -5, 5, 200),
    list("quart:1e-8:1", 10, -36.841361487904734, 36.841361487904734, 200),
    list("explin:2", 0, -10, 10, 200), list("explin:1e-6", 0, -36.841361487904734, 36.841361487904734, 200),
    list("explin:1e8", 0, -36.841361487904734, 36.841361487904734, 200),
    list("flat", 0, -36.841361487904734, 36.841361487904734, 200),
    list("flat", 0, -36.841361487904734, 36.841361487904734, 3),
    list("explin:3", 5, -10, 10, 1), list("explin:3", 5, -10, 10, 2),
    list("nbll:5:20", log(0.5), -36.841361487904734, 36.841361487904734, 200),
    list("nbll:0.5:300", 1, -36.841361487904734, 36.841361487904734, 200),
    list("nbll:50:2", -3, -36.841361487904734, 36.841361487904734, 200),
    list("rosen", c(-1.2, 1), c(-50, -50), c(50, 50), 200),
    list("rosen", c(-1.2, 1), c(-2, -2), c(0.5, 2), 200),
    list("rosen", c(2, 2), c(1.5, -1), c(3, 3), 200),
    list("rosen", c(0, 0), c(-1, -1), c(1, 1), 5))
# the same cases without the Hessian (DRMNGB: glmGamPoi's fallback)
nl_cases <- c(lapply(nl_cases, function(cs) c(cs, TRUE)), lapply(nl_cases, function(cs) c(cs, FALSE)))
write_tsv("nlminb.tsv",
          paste("function\tstart\tlower\tupper\titer.max\thessian\tpar\tobjective\tconvergence",
                "iterations\tevaluations (f,g)\tiv[1] code\ttrace kinds\ttrace points", sep = "\t"),
          lapply(nl_cases, function(cs) {
              fun <- nl_fun(cs[[1]])
              kinds <- character(); pts <- numeric()
              res <- nlminb(cs[[2]],
                            objective = function(x) { kinds <<- c(kinds, "f"); pts <<- c(pts, x); fun$fn(x) },
                            gradient = function(x) { kinds <<- c(kinds, "g"); pts <<- c(pts, x); fun$gr(x) },
                            hessian = if (cs[[6]]) function(x) fun$he(x),
                            lower = cs[[3]], upper = cs[[4]], control = list(iter.max = cs[[5]]))
              code <- match(res$message, vapply(c(3:10, 14:16, 63, 65), function(i)
                  stats:::port_msg(i), ""))
              code <- c(3:10, 14:16, 63, 65)[code]
              if (is.na(code)) code <- -1L
              paste(cs[[1]], hexl(cs[[2]]), hexl(cs[[3]]), hexl(cs[[4]]), cs[[5]], as.integer(cs[[6]]), hexl(res$par),
                    hex(res$objective), res$convergence, res$iterations,
                    paste(res$evaluations, collapse = ","), code, lst(kinds), hexl(pts), sep = "\t")
          }))

# ---------------------------------------------------------------------------------------------
# lgamma (R's lgammafn), digamma and trigamma (dpsifn), at the arguments glmGamPoi's
# likelihood reaches (count + 1/theta for theta from 1e-16 to 1e16, counts and half counts) and
# across the branches of the implementations (tiny, <= 10, > 10, huge, negative non-integers).

set.seed(47)
gm_x <- c(1e-310, 1e-300, 1e-20, 1e-8, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 9.999, 10, 10.5, 11, 50.25,
          1e3, 4934720, 4934721, 1e10, 1e17, 1.1e17, 1e300,
          -0.5, -1.5, -2.25, -9.5, -10.5, -11.25, -100.5,
          exp(runif(150, log(1e-16), log(1e16))),
          outer(c(0, 0.5, 1, 1.5, 2, 3, 5, 7.5, 12, 40, 155.5),
                exp(c(-36.8, -20, -9, -3, -1, 0, 1, 3, 9, 20, 36.8)), "+"),
          runif(100, 0, 30))
write_tsv("gamma.tsv", "x\tlgamma(x)\tdigamma(x)\ttrigamma(x)",
          lapply(gm_x, function(x) paste(hex(x), hex(lgamma(x)), hex(digamma(x)),
                                         hex(trigamma(x)), sep = "\t")))

# ---------------------------------------------------------------------------------------------
# matrixStats 1.5.0 as glmGamPoi reaches it through DelayedMatrixStats/MatrixGenerics (checked
# identical here): rowMeans2, colMeans2 and rowVars of matrices (refine = TRUE; the installed
# package accumulates in double, as HAVE_LONG_DOUBLE is not visible to packages), and
# weightedMedian(x, w, interpolate = FALSE) with ties in x (R_qsort_I is not stable), as
# loc_median_fit calls it with dnorm weights, and weightedMedian(x, w) (interpolate = TRUE). Matrices are written column-major.

suppressPackageStartupMessages(library(matrixStats))
set.seed(53)
ms_mats <- list(matrix(rpois(60, 2) / 2, 6, 10), matrix(runif(35) * 1e3, 5, 7),
                matrix(exp(rnorm(400, 0, 3)), 20, 20), matrix(c(0.1, 0.1, 0.1, 0.2, 0.3, 0.1), 2, 3),
                matrix(rep(1 / 3, 64), 4, 16))
write_tsv("matrixstats.tsv", "nrow\tncol\tx\trowMeans2\tcolMeans2\trowVars",
          lapply(ms_mats, function(m) {
              stopifnot(identical(DelayedMatrixStats::rowMeans2(m), matrixStats::rowMeans2(m)),
                        identical(DelayedMatrixStats::colMeans2(m), matrixStats::colMeans2(m)),
                        identical(DelayedMatrixStats::rowVars(m), matrixStats::rowVars(m)))
              paste(nrow(m), ncol(m), hexl(m), hexl(rowMeans2(m)), hexl(colMeans2(m)),
                    hexl(rowVars(m)), sep = "\t")
          }))
wm_cases <- list()
for (h in c(1, 2, 5, 10, 50)) {
    w <- dnorm(seq(-3, 3, length.out = 2 * h + 1))
    wm_cases[[length(wm_cases) + 1]] <- list(runif(2 * h + 1), w)
    wm_cases[[length(wm_cases) + 1]] <- list(sample(c(0, 0, 0, 1e16, runif(2 * h + 1)))[seq_len(2 * h + 1)], w)
    wm_cases[[length(wm_cases) + 1]] <- list(round(runif(2 * h + 1) * 3) / 4, w)
}
wm_cases[[length(wm_cases) + 1]] <- list(c(1, 2, 3), c(1, 1, 1))
wm_cases[[length(wm_cases) + 1]] <- list(c(1, 2, 3, 4), c(1, 1, 1, 1))
wm_cases[[length(wm_cases) + 1]] <- list(rep(0.25, 31), dnorm(seq(-3, 3, length.out = 31)))
# loc_median_fit's fallback when its window is wider than the data: weightedMedian(y, w =
# dnorm(seq(-3, 3, length.out = length(y)))) with the default interpolate = TRUE
for (n in c(2, 3, 4, 7, 10, 25)) {
    w <- dnorm(seq(-3, 3, length.out = n))
    wm_cases[[length(wm_cases) + 1]] <- list(runif(n), w)
    wm_cases[[length(wm_cases) + 1]] <- list(round(runif(n) * 3) / 4, w)
}
write_tsv("weighted_median.tsv", "x\tw\tweightedMedian(x, w, interpolate = FALSE)\tweightedMedian(x, w)",
          lapply(wm_cases, function(cs) paste(hexl(cs[[1]]), hexl(cs[[2]]),
                 hex(weightedMedian(cs[[1]], w = cs[[2]], interpolate = FALSE)),
                 hex(weightedMedian(cs[[1]], w = cs[[2]])), sep = "\t")))

# ---------------------------------------------------------------------------------------------
# dnbinom(x, size, mu = mu, log = TRUE) (nmath dnbinom_mu, with dbinom_raw, stirlerr, bd0 and
# lgamma1p), as glmGamPoi's optimize() fallback evaluates it: sizes on both sides of every
# stirlerr() branch and series cut-off, x = 0, x < 1e-10 * size, non-integer x (-Inf), mu from
# 1e-12 to 1e8. One line per (size, mu) pair.

set.seed(61)
nb_size <- c(1e-6, 0.006, 0.1, 0.3, 0.49, 0.5, 0.51, 0.75, 0.99, 1, 1.5, 2, 2.25, 3.7, 5, 5.25,
             5.3, 5.5, 6, 6.05, 6.3, 6.6, 6.9, 7.3, 8, 8.9, 10.25, 12.3, 12.5, 12.8, 13.3, 15, 15.5,
             20, 23.5, 24, 26.9, 27, 50, 86, 100, 205, 1000, 6180, 1e4, 15.7e6, 2e7, 1e11, 1e13,
             1 / c(166.61835542451206, 2.869887632263776, 4.69942116835039),
             runif(10, 0, 30), exp(rnorm(10, 0, 4)))
nb_x <- c(0, 1, 2, 3, 4, 7, 10, 13, 30, 100, 1234, 1e5, 0.5, 2.5, 3 + 1e-12, 1e-300)
nb_mu <- c(1e-12, 1e-4, 0.03, 0.5, 1, 2.7, 10, 333, 1e5, 1e8)
nb_rows <- list()
for (size in nb_size) for (mu in nb_mu)
    nb_rows[[length(nb_rows) + 1]] <- paste(hex(size), hex(mu), hexl(nb_x),
        hexl(suppressWarnings(dnbinom(nb_x, size = size, mu = mu, log = TRUE))), sep = "\t")
write_tsv("dnbinom.tsv", "size\tmu\tx\tdnbinom(x, size, mu = mu, log = TRUE)", nb_rows)

# ---------------------------------------------------------------------------------------------
# optimize() (Brent_fmin; fcn1 turns -Inf into -DBL_MAX, NaN and Inf into DBL_MAX) on functions
# whose every value Rust computes the same way: each line has the points optimize() evaluated,
# in order (without the objective at the result that optimize() adds), and the result.

opt_fns <- list(
    list("sq", function(x) (x - 1/3)^2, -3, 5),
    list("abs", function(x) abs(x - 0.7), 0, 2),
    list("floor", function(x) floor(4 * x) / 4, -1, 1),
    list("inf_right", function(x) if (x > 0.3) Inf else (x - 0.5)^2, -2, 2),
    list("sin", function(x) sin(3 * x) + x / 5, -4, 4),
    list("neg_inf", function(x) -Inf, -30, 30),
    list("xexp", function(x) exp(-x) * x, 0, 10),
    list("log", function(x) log(x), -1, 2))
opt_rows <- list()
for (fn in opt_fns) for (maximum in c(FALSE, TRUE)) {
    xs <- numeric(0)
    g <- function(x) { xs <<- c(xs, x); fn[[2]](x) }
    r <- suppressWarnings(optimize(g, lower = fn[[3]], upper = fn[[4]], maximum = maximum))
    res <- if (maximum) r$maximum else r$minimum
    stopifnot(identical(xs[length(xs)], res))
    opt_rows[[length(opt_rows) + 1]] <- paste(fn[[1]], as.integer(maximum), hex(fn[[3]]), hex(fn[[4]]),
                                              hexl(xs[-length(xs)]), hex(res), sep = "\t")
}
write_tsv("optimize.tsv", "function\tmaximum\tlower\tupper\tevaluated at\tresult", opt_rows)

cat("written to", normalizePath(out), "\n")
