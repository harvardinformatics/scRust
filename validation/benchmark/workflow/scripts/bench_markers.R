# Phase markers for the benchmark copies of the preprocessing workflow's R scripts
# (scRust data/benchmarks; not part of the workflow). Sourced first by each copy:
#   source("workflow/scripts/bench_markers.R")
# and called between the scripts' unchanged statements:
#   bench_phase("<phase>", "<kind>")   ends the current phase and starts the next
#   bench_subphase("<phase>")          the same, keeping the current phase's kind
#   bench_end()                        ends the last phase
# kind: startup (R and the packages), read, compute (the steps the Rust port does), recluster
# (doubletfinder.R's clustering of the singlets: ported by `scrust cluster`, not by
# `scrust doubletfinder`), unported (steps no port does: SoupX, emptyDrops, RunUMAP), write.
#
# Each mark appends one row per event to $BENCH_MARKERS (tab-separated): the phase's end
# (its CPU time so far, its peak RSS = VmHWM since the phase began, the RSS at its end, the
# /proc/self/io counters) and the next phase's start (after VmHWM is reset to the current RSS
# by writing 5 to /proc/self/clear_refs, so the next VmHWM read is that phase's own peak).
# CPU times are proc.time(): getrusage of this process (all its threads) and of its waited-for
# children. Without $BENCH_MARKERS the functions do nothing.

.bench <- new.env(parent = emptyenv())
.bench$file <- Sys.getenv("BENCH_MARKERS", "")
.bench$phase <- "startup"
.bench$kind <- "startup"
.bench$seq <- 0L

.bench_counters <- function() {
  st <- readLines("/proc/self/status")
  kib <- function(key) {
    l <- st[startsWith(st, paste0(key, ":"))]
    if (length(l) == 0) return(NA_real_)
    as.numeric(strsplit(trimws(sub("^[^:]*:", "", l[1])), "[[:space:]]+")[[1]][1])
  }
  io <- readLines("/proc/self/io")
  iov <- as.numeric(sub("^[^:]*:[[:space:]]*", "", io))
  names(iov) <- sub(":.*", "", io)
  pt <- proc.time()
  list(
    epoch_s = sprintf("%.6f", as.numeric(Sys.time())),
    elapsed_s = sprintf("%.3f", pt[["elapsed"]]),
    user_s = sprintf("%.3f", pt[["user.self"]]),
    sys_s = sprintf("%.3f", pt[["sys.self"]]),
    child_user_s = sprintf("%.3f", pt[["user.child"]]),
    child_sys_s = sprintf("%.3f", pt[["sys.child"]]),
    vmrss_kib = sprintf("%.0f", kib("VmRSS")), vmhwm_kib = sprintf("%.0f", kib("VmHWM")),
    rssanon_kib = sprintf("%.0f", kib("RssAnon")), rssfile_kib = sprintf("%.0f", kib("RssFile")),
    rchar = sprintf("%.0f", iov[["rchar"]]), wchar = sprintf("%.0f", iov[["wchar"]]),
    syscr = sprintf("%.0f", iov[["syscr"]]), syscw = sprintf("%.0f", iov[["syscw"]]),
    read_bytes = sprintf("%.0f", iov[["read_bytes"]]), write_bytes = sprintf("%.0f", iov[["write_bytes"]])
  )
}

.bench_row <- function(event, phase, kind, reset_ok = NA) {
  c0 <- .bench_counters()
  .bench$seq <- .bench$seq + 1L
  header <- !file.exists(.bench$file)
  # every value is a string already (sprintf): no padding, no scientific notation
  row <- c(seq = as.character(.bench$seq), event = event, phase = phase, kind = kind, unlist(c0),
           hwm_reset = if (is.na(reset_ok)) "" else as.character(reset_ok))
  con <- file(.bench$file, open = "a")
  on.exit(close(con))
  if (header) writeLines(paste(names(row), collapse = "\t"), con)
  writeLines(paste(row, collapse = "\t"), con)
  invisible(NULL)
}

.bench_reset_hwm <- function() {
  ok <- tryCatch({
    cat("5\n", file = "/proc/self/clear_refs")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  ok
}

bench_phase <- function(phase, kind) {
  if (!nzchar(.bench$file)) return(invisible(NULL))
  .bench_row("end", .bench$phase, .bench$kind)
  ok <- .bench_reset_hwm()
  .bench$phase <- phase
  .bench$kind <- kind
  .bench_row("start", phase, kind, ok)
  invisible(NULL)
}

# a phase of the same kind as the current one (used inside the shared helper silhouette_utils.R,
# whose add_silhouette_to_metadata() runs in the chain and in DoubletFinder's re-clustering)
bench_subphase <- function(phase) bench_phase(phase, .bench$kind)

bench_end <- function() {
  if (!nzchar(.bench$file)) return(invisible(NULL))
  .bench_row("end", .bench$phase, .bench$kind)
  .bench$phase <- "(after the script)"
  invisible(NULL)
}
