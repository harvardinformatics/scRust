# Markdown tables for docs/validation/doubletfinder_m5.md from the M5 outputs: `scrust validate
# pca` on every instrumented reference run (data/doubletfinder/<object>/rust/<label>/<run>/:
# validate.json, sources.txt, time_v.txt), the runs' own stage scalars (stages.json, index.tsv),
# and the benchmark (data/doubletfinder/<object>/rust/m5-benchmark/final/). Reads only plain files.
#
#   Rscript validation/R/report_doubletfinder_m5.R <label> > <tables.md>

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
setwd(normalizePath(file.path(here, "..", "..")))
args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[[1]] else stop("usage: report_doubletfinder_m5.R <label>")

objects <- c("soupx_tenx_test", "soupx_emptydrops_test", "decontx_raw_tenx_test", "decontx_raw_emptydrops_test",
             "cellbender_fromraw_test", "soupx_tenx_SRR18256555", "soupx_tenx_10xmousebrain",
             "soupx_emptydrops_4500871lane", "soupx_emptydrops_sscovelli", "cellbender_fromraw_ghirsutum2")
js <- function(p) if (file.exists(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
time_v <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p)
  g <- function(pat) sub(".*: ", "", grep(pat, l, value = TRUE, fixed = TRUE))
  pt <- as.numeric(strsplit(g("Elapsed (wall clock) time"), ":")[[1]])
  secs <- if (length(pt) == 3) pt[1] * 3600 + pt[2] * 60 + pt[3] else pt[1] * 60 + pt[2]
  list(rss_kib = as.numeric(g("Maximum resident set size")), wall_s = secs,
       user_s = as.numeric(g("User time (seconds)")), sys_s = as.numeric(g("System time (seconds)")))
}
fmt_s <- function(s) if (s >= 120) sprintf("%d:%05.2f", as.integer(s %/% 60), s %% 60) else sprintf("%.1f s", s)
fmt_gib <- function(kib) sprintf("%.1f", kib / 1024^2)

# Per run and context ------------------------------------------------------------------------
tot <- list()
for (obj in objects) {
  runs <- sort(list.files(file.path("data/doubletfinder", obj, "rust", label), pattern = "^instrumented-"))
  runs <- runs[!grepl("partial", runs)]
  for (run in runs) {
    d <- file.path("data/doubletfinder", obj, "rust", label, run)
    v <- js(file.path(d, "validate.json"))
    st <- js(file.path("data/doubletfinder", obj, "runs", run, "stages", "stages.json"))
    idx <- read.delim(file.path("data/doubletfinder", obj, "runs", run, "stages", "index.tsv"), colClasses = "character")
    src <- read.delim(file.path(d, "sources.txt"), header = FALSE, colClasses = "character", col.names = c("ctx", "src"))
    tv <- time_v(file.path(d, "time_v.txt"))
    ctxs <- idx$ctx[idx$file == "S3_embeddings.npy"]
    cat(sprintf("\n#### %s, `%s`\n\nValidation: %s wall, %s GiB peak RSS (%s threads); %s.\n\n", obj, run, fmt_s(tv$wall_s),
                fmt_gib(tv$rss_kib), 8, if (all(src$src == "R")) "teacher-forced on R's PCA input" else
                  "teacher-forced on Rust's S2d (R dumped no PCA input)"))
    rows <- lapply(ctxs, function(ctx) {
      sc <- st$scalars[[ctx]]
      cells <- as.integer(idx$nrow[idx$file == "S3_embeddings.npy" & idx$ctx == ctx])
      cc <- Filter(function(c) startsWith(c$name, paste0(ctx, ":")), v$checks)
      agg <- function(stage) {
        x <- Filter(function(c) c$stage == stage, cc)
        ok <- sum(vapply(x, function(c) c$status == "identical", logical(1)))
        n <- sum(vapply(x, function(c) c$n, numeric(1)))
        nd <- sum(vapply(x, function(c) c$n_differ, numeric(1)))
        tot[[stage]] <<- (if (is.null(tot[[stage]])) c(0, 0, 0, 0) else tot[[stage]]) + c(ok, length(x), n, nd)
        rng <- Filter(function(c) grepl("RNG", c$name), x)
        tot[[paste0(stage, "_rng")]] <<- (if (is.null(tot[[paste0(stage, "_rng")]])) c(0, 0) else tot[[paste0(stage, "_rng")]]) +
          c(sum(vapply(rng, function(c) c$status == "identical", logical(1))), length(rng))
        sprintf("%d/%d, %s", ok, length(x), num(n))
      }
      data.frame(context = ctx, cells = num(cells), features = num(sc$pca_n_features), iter = sc$irlba_iter,
                 mprod = sc$irlba_mprod, S3 = agg("S3"), `S3-chain` = agg("S3-chain"), check.names = FALSE)
    })
    cat(md_table(do.call(rbind, rows)), sep = "\n")
  }
}
cat("\n### Totals\n\n")
cat(md_table(data.frame(stage = c("S3", "S3-chain"),
  `checks identical` = vapply(c("S3", "S3-chain"), function(s) sprintf("%s of %s", num(tot[[s]][1]), num(tot[[s]][2])), ""),
  `elements compared` = vapply(c("S3", "S3-chain"), function(s) num(tot[[s]][3]), ""),
  `elements differing` = vapply(c("S3", "S3-chain"), function(s) num(tot[[s]][4]), ""),
  `RNG checkpoints identical` = vapply(c("S3", "S3-chain"), function(s) sprintf("%s of %s", num(tot[[paste0(s, "_rng")]][1]), num(tot[[paste0(s, "_rng")]][2])), ""),
  check.names = FALSE, row.names = NULL)), sep = "\n")

# Benchmark -----------------------------------------------------------------------------------
cat("\n### Benchmark\n\n")
brow <- list()
for (obj in c("soupx_emptydrops_sscovelli", "cellbender_fromraw_ghirsutum2")) {
  b <- file.path("data/doubletfinder", obj, "rust", "m5-benchmark", "final")
  if (!dir.exists(b)) next
  for (side in sort(list.files(b, pattern = "^(r|rust_t[0-9]+)$"))) {
    tv <- time_v(file.path(b, side, "time_v.txt"))
    if (side == "r") {
      t <- js(file.path(b, side, "out", "timing.json"))
      pca_s <- t$runpca_elapsed_s; load_rss <- t$peak_rss_after_input_kib; threads <- 1
      shape <- sprintf("%s x %s", num(t$n_cells), num(t$n_features)); it <- sprintf("%d, %d", t$irlba_iter, t$irlba_mprod)
    } else {
      s <- js(file.path(b, side, "out", "summary.json"))
      ph <- s$run$phases
      pca_s <- ph$pca$elapsed_s - ph$load$elapsed_s; load_rss <- ph$load$peak_rss_kib; threads <- s$run$threads
      shape <- sprintf("%s x %s", num(s$result$n_cells), num(s$result$n_features))
      it <- sprintf("%d, %d", s$result$irlba_iter, s$result$irlba_mprod)
    }
    brow[[length(brow) + 1]] <- data.frame(`object (final run)` = obj, `cells x features` = shape, side = if (side == "r") "R" else "Rust",
      threads = threads, `irlba iter, products` = it, `RunPCA (s)` = sprintf("%.1f", pca_s), `process wall (s)` = sprintf("%.1f", tv$wall_s),
      `process CPU (s)` = sprintf("%.1f", tv$user_s + tv$sys_s), `peak RSS (GiB)` = fmt_gib(tv$rss_kib),
      `RSS after loading (GiB)` = fmt_gib(load_rss), check.names = FALSE)
  }
}
if (length(brow)) cat(md_table(do.call(rbind, brow)), sep = "\n")
