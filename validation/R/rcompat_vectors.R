#!/usr/bin/env Rscript
# Reference vectors for the rcompat crate (docs/PLAN.md A5 and B6, milestone M1; stage S0 in
# validation/tolerances.toml): R's random-number stream, sampling, rnorm/qnorm, round, ^, seq,
# order, sum/mean/var/sd and number formatting, over many seeds, sizes and edge cases.
#
# Reads no .rds, so it may run on a login node. Run from the repository root with the production
# R env and the hermetic flags:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
#     "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
#
# Output: one tab-separated file per family; lines starting with "#" describe the columns. Doubles
# are written as their IEEE-754 bit patterns (16 hex digits), lists comma-separated, "-" for an
# empty list. Results do not depend on BLAS, so the files are the same on any x86-64 node with
# this R build and glibc (2.28 on Cannon).

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

cat("written to", normalizePath(out), "\n")
