# Tables for docs/validation/doubletfinder_m7.md (milestone M7, the full panel) from the M7 outputs.
# For each object of validation/doubletfinder_panel.tsv (its workflow job: node, CPU type, time
# and peak RSS from sacct):
# - the exports (validation/slurm/m7_export.sbatch): data/doubletfinder/<object>/input/export.json,
#   saved/saved.json, saved/singlets/export.json;
# - the Rust run and its comparison with the saved workflow output
#   (validation/slurm/m7_panel.sbatch), data/doubletfinder/<object>/rust/<label>/production/:
#   vs_saved.json, run/summary.json, fingerprint.json, node_features.txt, time_v.txt.
# Reads only plain files. Writes the markdown tables to stdout; with --tsv, also one row per
# object as TSV; with --same-as, compares each object's run with the run of another label (made
# with another binary) file by file: cells.tsv, sweep_stats.tsv, bcmvn.tsv and the three files of
# singlets/ byte for byte, summary.json's result section field by field.
#
#   Rscript validation/R/report_doubletfinder_m7.R <label> [--tsv <objects.tsv>] [--same-as <label>] > <tables.md>
#
# An object passes when all nine checks are made and identical: pANN_column,
# classification_column, pK_label, nExp, the doublet set, the kept cells' pANN, and the
# singlets' barcodes, feature lines and matrix.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
setwd(normalizePath(file.path(here, "..", "..")))
args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: report_doubletfinder_m7.R <label> [--tsv <objects.tsv>] [--same-as <label>]"
if (length(args) < 1 || startsWith(args[[1]], "--")) stop(usage)
label <- args[[1]]
opt <- function(name) { k <- match(name, args); if (is.na(k)) NA else if (k < length(args)) args[[k + 1]] else stop(usage) }
tsv_out <- opt("--tsv")
same_as <- opt("--same-as")

panel <- read.delim("validation/doubletfinder_panel.tsv", colClasses = "character", quote = "")
expected <- c("pANN_column", "classification_column", "pK_label", "nExp", "doublet set", "singlet pANN",
              "singlets barcodes", "singlets features", "singlets matrix")
arms <- c("soupx_tenx", "soupx_emptydrops", "decontx_raw_tenx", "decontx_raw_emptydrops", "cellbender_fromraw")

js <- function(p) if (file.exists(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else NULL
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
out <- function(lines) cat(lines, sep = "\n")
num <- function(x) ifelse(is.na(x), "", format(x, big.mark = ",", scientific = FALSE, trim = TRUE))
secs <- function(txt) {
  pt <- as.numeric(strsplit(txt, ":")[[1]])
  if (length(pt) == 3) pt[1] * 3600 + pt[2] * 60 + pt[3] else pt[1] * 60 + pt[2]
}
time_v <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p)
  g <- function(pat) sub(".*: ", "", grep(pat, l, value = TRUE, fixed = TRUE))
  list(rss_kib = as.numeric(g("Maximum resident set size")), wall_s = secs(g("Elapsed (wall clock) time")),
       cpu_s = as.numeric(g("User time (seconds)")) + as.numeric(g("System time (seconds)")))
}
fmt_s <- function(s) {
  if (is.null(s) || is.na(s)) return("")
  if (s >= 3600) sprintf("%d:%02d:%02d", as.integer(s %/% 3600), as.integer((s %% 3600) %/% 60), as.integer(round(s %% 60)))
  else if (s >= 120) sprintf("%d:%04.1f", as.integer(s %/% 60), s %% 60) else sprintf("%.1f s", s)
}
fmt_gib <- function(kib) sprintf("%.2f", kib / 1024^2)
cpu_short <- function(m) sub("^Intel\\(R\\) Xeon\\(R\\) (Platinum|Gold|Silver) ", "", m)

# One record per object --------------------------------------------------------------------------
recs <- list()
for (k in seq_len(nrow(panel))) {
  p <- panel[k, ]
  obj <- p$object
  base <- file.path("data/doubletfinder", obj)
  d <- file.path(base, "rust", label, "production")
  ein <- js(file.path(base, "input", "export.json"))
  esv <- js(file.path(base, "saved", "saved.json"))
  esg <- js(file.path(base, "saved", "singlets", "export.json"))
  v <- js(file.path(d, "vs_saved.json"))
  sm <- js(file.path(d, "run", "summary.json"))
  fp <- js(file.path(d, "fingerprint.json"))
  tv <- time_v(file.path(d, "time_v.txt"))
  nf <- if (file.exists(file.path(d, "node_features.txt"))) strsplit(readLines(file.path(d, "node_features.txt"), n = 1), " ")[[1]] else c(NA, NA)
  ch <- if (is.null(v)) list() else v$checks
  names(ch) <- vapply(ch, function(c) c$name, "")
  st <- vapply(expected, function(n) if (is.null(ch[[n]])) "missing" else ch[[n]]$status, "")
  n_of <- function(n) if (is.null(ch[[n]])) NA else as.numeric(ch[[n]]$n)
  made <- st %in% c("identical", "differ")
  rust_cpu_type <- if (is.na(nf[2])) NA else {
    f <- strsplit(nf[2], ",")[[1]]
    f <- f[f %in% c("sapphirerapids", "emeraldrapids", "icelake", "cascadelake", "skylake", "genoa", "milan", "rome")]
    if (length(f)) f[1] else NA
  }
  recs[[k]] <- list(
    object = obj, tier = p$tier, arm = p$arm, sample = p$sample,
    cells = if (is.null(ein)) NA else ein$n_cells, genes = if (is.null(ein)) NA else ein$n_features,
    run_cells = if (is.null(sm)) NA else sm$result$n_cells,
    singlets = if (is.null(esv)) NA else esv$n_kept, doublets = if (is.null(esv)) NA else esv$n_doublets,
    pK_label = if (is.null(esv)) NA else esv$pK_label, nExp = if (is.null(esv)) NA else esv$nExp,
    run_pK_label = if (is.null(sm)) NA else sm$result$pK_label, run_nExp = if (is.null(sm)) NA else sm$result$nExp,
    k = if (is.null(sm)) NA else sm$result$k, pANN_column = if (is.null(esv)) NA else esv$pANN_column,
    singlet_entries = if (is.null(esg)) NA else esg$nnz, singlet_explicit_zeros = if (is.null(esg)) NA else esg$explicit_zeros,
    saved_equals_input_subset = if (is.null(esg)) NA else isTRUE(esg$equals_input_subset),
    status = st, n = vapply(expected, n_of, numeric(1)),
    checks = sum(made), checks_identical = sum(st == "identical"),
    values_compared = sum(vapply(ch, function(c) if (c$status %in% c("identical", "differ")) as.numeric(c$n) else 0, numeric(1))),
    values_differing = sum(vapply(ch, function(c) if (c$status %in% c("identical", "differ")) as.numeric(c$n_differ) else 0, numeric(1))),
    pass = !is.null(v) && all(st == "identical"),
    workflow_job = p$workflow_job, workflow_run = p$workflow_run, workflow_node = p$workflow_node,
    workflow_cpu = p$workflow_cpu, workflow_wall_s = secs(p$workflow_elapsed),
    workflow_max_rss_kib = as.numeric(p$workflow_max_rss_kib),
    rust_job = if (is.null(sm)) NA else sm$run$slurm_job_id, rust_node = nf[1], rust_cpu_type = rust_cpu_type,
    rust_cpu_model = if (is.null(fp)) NA else fp$cpu_model, openblas_core = if (is.null(fp)) NA else fp$openblas_corename,
    openblas_threads = if (is.null(fp)) NA else fp$openblas_num_threads,
    threads = if (is.null(sm)) NA else sm$run$threads,
    wall_s = if (is.null(tv)) NA else tv$wall_s, cpu_s = if (is.null(tv)) NA else tv$cpu_s,
    peak_rss_kib = if (is.null(tv)) NA else tv$rss_kib)
}
ok <- vapply(recs, function(r) r$pass, logical(1))
done <- vapply(recs, function(r) !is.na(r$wall_s), logical(1))
# one field over all objects: character if any object's is, else numeric (logicals count as 0/1)
fld <- function(name) {
  x <- lapply(recs, function(r) r[[name]])
  na <- vapply(x, function(v) is.null(v) || (length(v) == 1 && is.na(v)), logical(1))
  if (any(vapply(x, is.character, logical(1)))) {
    vapply(seq_along(x), function(i) if (na[i]) NA_character_ else as.character(x[[i]]), "")
  } else {
    vapply(seq_along(x), function(i) if (na[i]) NA_real_ else as.numeric(x[[i]]), numeric(1))
  }
}

# 1. Per object ------------------------------------------------------------------------------------
# A check that is not identical shows its status: differs, missing (not made), n/a.
bad <- function(s) sprintf("**%s**", if (s == "differ") "differs" else s)
cell <- function(r, name) {
  s <- r$status[[name]]
  if (s == "identical") sprintf("= (%s)", num(r$n[[name]])) else bad(s)
}
first_bad <- function(r, names) { s <- r$status[names]; if (all(s == "identical")) NA else s[s != "identical"][[1]] }
rows <- lapply(recs, function(r) {
  same_cpu <- !is.na(r$rust_cpu_type) && r$rust_cpu_type == r$workflow_cpu
  pk_bad <- first_bad(r, "pK_label")
  nexp_bad <- first_bad(r, c("nExp", "pANN_column", "classification_column"))
  sgl_bad <- first_bad(r, c("singlets barcodes", "singlets features", "singlets matrix"))
  data.frame(
    object = r$object, cells = num(r$cells), singlets = num(r$singlets),
    `doublet set` = cell(r, "doublet set"), `singlet pANN` = cell(r, "singlet pANN"),
    pK = if (is.na(pk_bad)) r$pK_label else bad(pk_bad),
    nExp = if (is.na(nexp_bad)) num(r$nExp) else bad(nexp_bad),
    `singlets' counts` = if (is.na(sgl_bad)) sprintf("= (%s)", num(r$n[["singlets matrix"]])) else bad(sgl_bad),
    `workflow job: node, CPU type` = sprintf("%s, %s", r$workflow_node, r$workflow_cpu),
    `Rust: node, CPU, kernel` = if (is.na(r$rust_node)) "" else sprintf("%s, %s%s, %s", r$rust_node, cpu_short(r$rust_cpu_model),
                                                                     if (same_cpu) "" else " **(other CPU type)**", r$openblas_core),
    wall = fmt_s(r$wall_s), `peak RSS GiB` = if (is.na(r$peak_rss_kib)) "" else fmt_gib(r$peak_rss_kib),
    pass = if (r$pass) "yes" else "**no**", check.names = FALSE)
})
out(c("### Per object", ""))
out(md_table(do.call(rbind, rows)))

# 2. Totals ----------------------------------------------------------------------------------------
tot_checks <- sum(fld("checks"))
check_tot <- lapply(expected, function(n) {
  s <- vapply(recs, function(r) r$status[[n]], "")
  data.frame(check = n, objects = num(sum(s %in% c("identical", "differ"))), identical = num(sum(s == "identical")),
             `values compared` = num(sum(vapply(recs, function(r) if (r$status[[n]] %in% c("identical", "differ")) r$n[[n]] else 0, numeric(1)))),
             check.names = FALSE)
})
out(c("", "### Totals", ""))
out(md_table(data.frame(
  objects = num(length(recs)), compared = num(sum(done)), passed = num(sum(ok)),
  `checks identical` = sprintf("%s of %s", num(sum(fld("checks_identical"))), num(tot_checks)),
  `values compared` = num(sum(fld("values_compared"))), `values differing` = num(sum(fld("values_differing"))),
  `not made or missing` = num(length(recs) * length(expected) - tot_checks), check.names = FALSE)))
out(c("", "### By check", ""))
out(md_table(do.call(rbind, check_tot)))

# 3. By arm ------------------------------------------------------------------------------------------
rng <- function(x) if (all(is.na(x))) "" else sprintf("%s to %s", num(min(x, na.rm = TRUE)), num(max(x, na.rm = TRUE)))
arm_rows <- lapply(c(arms, "all"), function(a) {
  sel <- if (a == "all") rep(TRUE, length(recs)) else fld("arm") == a
  tot <- function(name) num(sum(fld(name)[sel], na.rm = TRUE))
  data.frame(arm = a, objects = num(sum(sel)), passed = num(sum(ok[sel])), cells = rng(fld("cells")[sel]),
             `cells, total` = tot("cells"), `doublets, total` = tot("doublets"), `singlets, total` = tot("singlets"),
             `singlets' stored entries` = tot("singlet_entries"), check.names = FALSE)
})
out(c("", "### By arm", ""))
out(md_table(do.call(rbind, arm_rows)))

# 4. CPU types and kernels -----------------------------------------------------------------------------
wf_cpu <- table(fld("workflow_cpu"))
ru_cpu <- table(fld("rust_cpu_type"), useNA = "ifany")
cores <- table(fld("openblas_core"), useNA = "ifany")
models <- table(fld("rust_cpu_model"), useNA = "ifany")
same <- sum(fld("rust_cpu_type") == fld("workflow_cpu"), na.rm = TRUE)
out(c("", "### CPU types and OpenBLAS kernels", ""))
out(md_table(data.frame(
  side = c("workflow jobs (sacct node, SLURM features)", "Rust runs (node's SLURM features; scrust fingerprint)"),
  objects = c(num(length(recs)), num(sum(done))),
  nodes = c(num(length(unique(fld("workflow_node")))), num(length(unique(na.omit(fld("rust_node")))))),
  `CPU type` = c(paste(sprintf("%s (%s)", names(wf_cpu), wf_cpu), collapse = ", "),
                 paste(sprintf("%s (%s)", names(ru_cpu), ru_cpu), collapse = ", ")),
  `CPU model` = c("", paste(sprintf("%s (%s)", cpu_short(names(models)), models), collapse = ", ")),
  `OpenBLAS kernel` = c("not recorded; OPENBLAS_CORETYPE unset (production profile)",
                        paste(sprintf("%s (%s)", names(cores), cores), collapse = ", ")),
  check.names = FALSE)))
out(c("", sprintf("Rust runs on the CPU type of the object's workflow job: %s of %s.", num(same), num(sum(done)))))

# 5. Cost ----------------------------------------------------------------------------------------------
w <- fld("wall_s"); m <- fld("peak_rss_kib"); cpu <- fld("cpu_s")
ww <- fld("workflow_wall_s"); wm <- fld("workflow_max_rss_kib")
imax <- function(x) which.max(replace(x, is.na(x), -Inf))
out(c("", "### Cost", ""))
out(md_table(data.frame(
  side = c(sprintf("Rust (scrust doubletfinder, %s threads, singlets/ included)", paste(unique(na.omit(fld("threads"))), collapse = "/")),
           "the workflow's jobs (DoubletFinder plus doubletfinder.R's Seurat processing; 1 CPU)"),
  `wall, total` = c(fmt_s(sum(w, na.rm = TRUE)), fmt_s(sum(ww))),
  `wall, median` = c(fmt_s(median(w, na.rm = TRUE)), fmt_s(median(ww))),
  `wall, largest` = c(sprintf("%s (%s)", fmt_s(max(w, na.rm = TRUE)), recs[[imax(w)]]$object), sprintf("%s (%s)", fmt_s(max(ww)), recs[[imax(ww)]]$object)),
  `CPU time, total` = c(fmt_s(sum(cpu, na.rm = TRUE)), ""),
  `peak RSS GiB, median` = c(fmt_gib(median(m, na.rm = TRUE)), fmt_gib(median(wm))),
  `peak RSS GiB, largest` = c(sprintf("%s (%s)", fmt_gib(max(m, na.rm = TRUE)), recs[[imax(m)]]$object), sprintf("%s (%s)", fmt_gib(max(wm)), recs[[imax(wm)]]$object)),
  check.names = FALSE)))

# 6. The exports -----------------------------------------------------------------------------------------
eq <- fld("saved_equals_input_subset")
out(c("", sprintf("Saved RNA counts equal to the input's restricted to the kept cells, slot for slot: %s of %s objects; explicit zeros among the singlets' stored entries: %s.",
                  num(sum(eq, na.rm = TRUE)), num(sum(!is.na(eq))), num(sum(fld("singlet_explicit_zeros"), na.rm = TRUE)))))

# 7. The same objects run with another binary (--same-as) ---------------------------------------------
same_files <- c("cells.tsv", "sweep_stats.tsv", "bcmvn.tsv", "singlets/matrix.mtx.gz", "singlets/barcodes.tsv.gz",
                "singlets/features.tsv.gz")
if (!is.na(same_as)) {
  sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA
  for (k in seq_along(recs)) {
    a <- file.path("data/doubletfinder", recs[[k]]$object, "rust", label, "production", "run")
    b <- file.path("data/doubletfinder", recs[[k]]$object, "rust", same_as, "production", "run")
    if (!dir.exists(a) || !dir.exists(b)) { recs[[k]]$same_as <- NA; next }
    files_same <- vapply(same_files, function(f) { x <- sha(file.path(a, f)); !is.na(x) && identical(x, sha(file.path(b, f))) }, logical(1))
    ra <- js(file.path(a, "summary.json"))$result
    rb <- js(file.path(b, "summary.json"))$result
    recs[[k]]$same_files <- sum(files_same)
    recs[[k]]$same_result <- !is.null(ra) && identical(ra, rb)
    recs[[k]]$same_as <- all(files_same) && recs[[k]]$same_result
  }
  sa <- lapply(recs, function(r) r$same_as)
  both <- !vapply(sa, function(x) is.null(x) || is.na(x), logical(1))
  out(c("", sprintf("### Against the runs of label %s", same_as), ""))
  out(md_table(data.frame(
    `objects run with both binaries` = num(sum(both)),
    `files identical` = sprintf("%s of %s", num(sum(vapply(recs[both], function(r) r$same_files, numeric(1)))), num(sum(both) * length(same_files))),
    `result sections identical` = sprintf("%s of %s", num(sum(vapply(recs[both], function(r) r$same_result, logical(1)))), num(sum(both))),
    `objects identical` = sprintf("%s of %s", num(sum(vapply(recs[both], function(r) r$same_as, logical(1)))), num(sum(both))),
    check.names = FALSE)))
  differ <- vapply(recs[both], function(r) if (isTRUE(r$same_as)) NA_character_ else r$object, "")
  if (any(!is.na(differ))) out(c("", paste("Not identical:", paste(na.omit(differ), collapse = ", "))))
}

# One row per object, as TSV ----------------------------------------------------------------------------
if (!is.na(tsv_out)) {
  cols <- c("object", "tier", "arm", "sample", "cells", "genes", "singlets", "doublets", "pK_label", "nExp", "k",
            "pANN_column", "singlet_entries", "singlet_explicit_zeros", "saved_equals_input_subset", "checks",
            "checks_identical", "values_compared", "values_differing", "pass", "workflow_job", "workflow_run",
            "workflow_node", "workflow_cpu", "workflow_wall_s", "workflow_max_rss_kib", "rust_job", "rust_node",
            "rust_cpu_type", "rust_cpu_model", "openblas_core", "openblas_threads", "threads", "wall_s", "cpu_s",
            "peak_rss_kib")
  df <- as.data.frame(lapply(setNames(cols, cols), function(cn) vapply(recs, function(r) {
    x <- r[[cn]]
    if (is.null(x) || (length(x) == 1 && is.na(x))) "" else if (is.logical(x)) tolower(as.character(x)) else format(x, scientific = FALSE, trim = TRUE, digits = 15)
  }, "")), check.names = FALSE)
  # the status of each check, by name
  for (n in expected) df[[paste0("check: ", n)]] <- vapply(recs, function(r) r$status[[n]], "")
  if (!is.na(same_as)) {
    df[[paste0("identical to ", same_as)]] <- vapply(recs, function(r) if (is.null(r$same_as) || is.na(r$same_as)) "" else tolower(r$same_as), "")
  }
  write.table(df, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
