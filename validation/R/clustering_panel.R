# Writes the clustering crate's full panel table (milestone M8; docs/PLAN.md C4):
# validation/clustering_panel.tsv, one row per Seurat object the production run clustered
# (840: 30 samples x 28 objects), with its pipeline paths and the SLURM job that wrote it.
#
#   Rscript validation/R/clustering_panel.R <out.tsv> <snakemake log> [<snakemake log> ...]
#
# For M8, the two Snakemake runs of the workflow's rerun #2 (8 and 9 November 2026):
#   .snakemake/log/2026-11-08T214133.605174.snakemake.log .snakemake/log/2026-11-09T220250.045638.snakemake.log
#
# It reads the pipeline repo's Snakemake logs, samplesheet, config and job logs (read-only), the
# saved objects' modification times (no .rds is opened), and queries sacct and scontrol, so a
# login node is enough. Run it with the hermetic flags (validation/R/README.md).
#
# - The jobs: each Snakemake log's "Job <n> has been submitted with SLURM jobid <id> (log:
#   .../rule_<rule>/<wildcards>/<id>.log)", for the 13 rules that run the clustering chain
#   (docs/PLAN.md C1). An object submitted more than once (a retry, or again in a later run) gets
#   the job that wrote its saved .rds: its last completed job. Snakemake touches each output once
#   it sees the job finish, so the file's modification time must come after that job's end, by at
#   most 30 minutes; every other completed job of the object must have ended before it started.
#   The others are listed in superseded_jobs (job:state).
# - sacct: the allocation's node, partition, state, exit code, start, end and elapsed time; the
#   MaxRSS of step .0 (the srun that runs the R script), in KiB.
# - scontrol: the node's CPU-type feature (sapphirerapids, icelake, cascadelake, ...).
# - The job's log (results/logs/...): the modularity optimizer's numbers as Seurat prints them
#   (nodes, edges, maximum modularity of the 10 starts with 4 decimals, communities) and
#   GroupSingletons' message (singletons, final clusters; 0 and the communities when absent).
#   A log must hold exactly one optimizer block.
#
# Columns: object, tier (P for the seven tier-P objects of validation/clustering_datasets.tsv, else
# F), sample, kind, rule, assay (v3 for the emptyDrops base objects, which hold a v3 Assay, else
# v5), input (the job's data inputs, ';'-separated: pipeline-relative, or absolute for Cell
# Ranger's matrices), saved, log, markers (the combined marker table), workflow_job, workflow_run
# (the Snakemake log's time stamp), superseded_jobs, workflow_node, workflow_cpu,
# workflow_partition, workflow_start, workflow_end, workflow_elapsed, workflow_max_rss_kib,
# saved_mtime, log_nodes, log_edges, log_max_modularity, log_communities, log_singletons,
# log_final_clusters.
# Rows: the samplesheet's samples in order (those the config excludes are skipped), each with its
# 28 objects in a fixed order (KINDS below).

PIPELINE <- "/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: clustering_panel.R <out.tsv> <snakemake log> [<snakemake log> ...]")
out_path <- args[1]
logs <- args[-1]

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "lib", "io.R"))

# The samples -----------------------------------------------------------------------------------
ss <- read.delim(file.path(PIPELINE, "samplesheet.tsv"), colClasses = "character")
cfg <- readLines(file.path(PIPELINE, "config", "config.yaml"))
excl_line <- grep("^excluded_samples:", cfg, value = TRUE)
excluded <- regmatches(excl_line, gregexpr('"[^"]*"', excl_line))[[1]]
excluded <- gsub('"', "", excluded)
samples <- ss$sampleid[!ss$sampleid %in% excluded]
tenx_dir <- setNames(ss$tenx_datadir, ss$sampleid)

# The objects: kind, rule, wildcard directory, object name and paths ---------------------------
ARMS <- c("soupx_tenx", "soupx_emptydrops", "decontx_raw_tenx", "decontx_raw_emptydrops", "cellbender_fromraw")
arm_rds <- function(arm, s) {
  if (arm == "cellbender_fromraw") return(sprintf("results/cellbender_fromraw/seurat_cellbender_fromraw_%s.rds", s))
  d <- if (startsWith(arm, "soupx")) "soupx" else "decontx_raw"
  sprintf("results/%s/seurat_%s_%s.rds", d, arm, s)
}
objects_of <- function(s) {
  fil <- file.path(tenx_dir[[s]], "filtered_feature_bc_matrix")
  raw <- file.path(tenx_dir[[s]], "raw_feature_bc_matrix")
  ed_mat <- sprintf("results/emptydrops/%s_emptydrops_filtered_matrix", s)
  base_tenx <- sprintf("results/seurat_filtered/filtered_seurat_tenx_%s.rds", s)
  base_ed <- sprintf("results/emptydrops/filtered_seurat_emptydrops_%s.rds", s)
  rows <- list(
    list("base_tenx", "tenx2seuratrds", s, sprintf("filtered_seurat_tenx_%s", s), fil, base_tenx,
         sprintf("results/logs/tenx2seuratrds/tenx2seuratrds_%s.log", s)),
    list("base_emptydrops", "emptydrops", s, sprintf("filtered_seurat_emptydrops_%s", s), raw, base_ed,
         sprintf("results/logs/emptydrops/emptydrops_%s.log", s)),
    list("base_cellbendercells", "cellbendercells2seurat", s, sprintf("filtered_seurat_cellbendercells_%s", s),
         sprintf("results/cellbendercells/%s_cellbendercells_filtered_matrix", s),
         sprintf("results/cellbendercells/filtered_seurat_cellbendercells_%s.rds", s),
         sprintf("results/logs/cellbendercells2seurat/cellbendercells2seurat_%s.log", s)),
    list("soupx", "soupx", s, sprintf("soupx_tenx_%s", s), paste(fil, raw, base_tenx, sep = ";"),
         arm_rds("soupx_tenx", s), sprintf("results/logs/soupx/soupx_tenx_%s.log", s)),
    list("soupx", "soupx_emptydrops", s, sprintf("soupx_emptydrops_%s", s), paste(ed_mat, raw, base_ed, sep = ";"),
         arm_rds("soupx_emptydrops", s), sprintf("results/logs/soupx/soupx_emptydrops_%s.log", s)),
    list("decontx", "decontx_raw", sprintf("tenx_%s", s), sprintf("decontx_raw_tenx_%s", s),
         paste(fil, raw, base_tenx, sep = ";"), arm_rds("decontx_raw_tenx", s),
         sprintf("results/logs/decontx/decontx_raw_tenx_%s.log", s)),
    list("decontx", "decontx_raw", sprintf("emptydrops_%s", s), sprintf("decontx_raw_emptydrops_%s", s),
         paste(ed_mat, raw, base_ed, sep = ";"), arm_rds("decontx_raw_emptydrops", s),
         sprintf("results/logs/decontx/decontx_raw_emptydrops_%s.log", s)),
    list("cellbender_fromraw", "cellbender2seurat", s, sprintf("cellbender_fromraw_%s", s),
         sprintf("results/cellbender/cellbender_%s_filtered.h5", s), arm_rds("cellbender_fromraw", s),
         sprintf("results/logs/cellbender2seurat/cellbender2seurat_%s.log", s))
  )
  for (tool in c("doubletfinder", "scdblfinder")) {
    for (arm in ARMS) {
      cb <- arm == "cellbender_fromraw"
      rows[[length(rows) + 1]] <- list(
        tool, if (cb) paste0(tool, "_cellbender") else tool, if (cb) s else paste0(arm, "_", s),
        sprintf("%s_%s_%s", tool, arm, s), arm_rds(arm, s),
        sprintf("results/%s/seurat_%s_%s_%s.rds", tool, tool, arm, s),
        sprintf("results/logs/%s/%s_%s_%s.log", tool, tool, arm, s))
    }
  }
  for (tool in c("doubletfinder", "scdblfinder")) {
    for (arm in ARMS) {
      cb <- arm == "cellbender_fromraw"
      rows[[length(rows) + 1]] <- list(
        "posthoc", if (cb) "posthocfilter_threshold_cellbender" else "posthocfilter_threshold",
        if (cb) paste0(tool, "_", s) else paste0(tool, "_", arm, "_", s),
        sprintf("posthocfilt_threshold_%s_%s_%s", tool, arm, s),
        sprintf("results/%s/seurat_%s_%s_%s.rds", tool, tool, arm, s),
        sprintf("results/posthocfilter/seurat_posthocfilt_threshold_%s_%s_%s.rds", tool, arm, s),
        sprintf("results/logs/posthocfilter/posthocfilter_threshold_%s_%s_%s.log", tool, arm, s))
    }
  }
  df <- as.data.frame(do.call(rbind, lapply(rows, unlist)), stringsAsFactors = FALSE)
  names(df) <- c("kind", "rule", "wildcards", "object", "input", "saved", "log")
  df$sample <- s
  df
}
obj <- do.call(rbind, lapply(samples, objects_of))
if (anyDuplicated(obj$object)) stop("duplicated object names")
obj$assay <- ifelse(obj$kind == "base_emptydrops", "v3", "v5")
obj$markers <- sub("\\.rds$", "_markergenes.csv", obj$saved)
obj$key <- paste(obj$rule, obj$wildcards)

# The jobs, from the Snakemake logs --------------------------------------------------------------
jobs <- do.call(rbind, lapply(logs, function(log) {
  txt <- readLines(file.path(PIPELINE, log), warn = FALSE)
  m <- regmatches(txt, regexpr("has been submitted with SLURM jobid [0-9]+ \\(log: [^)]*\\)", txt))
  rule <- sub(".*/rule_([a-z0-9_]+)/[^/]+/[0-9]+\\.log\\)$", "\\1", m)
  wc <- sub(".*/rule_[a-z0-9_]+/([^/]+)/[0-9]+\\.log\\)$", "\\1", m)
  id <- sub("^has been submitted with SLURM jobid ([0-9]+) .*", "\\1", m)
  data.frame(key = paste(rule, wc), job = id, run = basename(sub("\\.snakemake\\.log$", "", log)),
             stringsAsFactors = FALSE)
}))
jobs <- jobs[jobs$key %in% obj$key, ]
missing <- setdiff(obj$key, jobs$key)
if (length(missing)) stop("objects without a job: ", paste(head(missing), collapse = ", "))

# sacct and scontrol ----------------------------------------------------------------------------
sacct <- function(ids, fmt, steps = FALSE) {
  out <- character()
  for (chunk in split(ids, ceiling(seq_along(ids) / 200))) {
    a <- c("-j", paste(chunk, collapse = ","), "-n", "-P", "--units=K", paste0("--format=", fmt))
    if (!steps) a <- c(a, "-X")
    out <- c(out, system2("sacct", a, stdout = TRUE))
  }
  m <- do.call(rbind, strsplit(out, "|", fixed = TRUE))
  as.data.frame(m, stringsAsFactors = FALSE)
}
alloc <- sacct(unique(jobs$job), "JobID,NodeList,Partition,State,ExitCode,Start,End,Elapsed")
names(alloc) <- c("job", "node", "partition", "state", "exit", "start", "end", "elapsed")
steps <- sacct(unique(jobs$job), "JobID,MaxRSS", steps = TRUE)
names(steps) <- c("step", "maxrss")
steps <- steps[grepl("\\.0$", steps$step), ]
rss <- setNames(sub("K$", "", steps$maxrss), sub("\\.0$", "", steps$step))
jobs <- merge(jobs, alloc, by = "job", all.x = TRUE, sort = FALSE)
if (anyNA(jobs$state)) stop("jobs unknown to sacct: ", paste(jobs$job[is.na(jobs$state)], collapse = ", "))
nodes <- unique(jobs$node[jobs$node != "None assigned"])
cpu_of <- vapply(nodes, function(n) {
  f <- system2("scontrol", c("show", "node", n), stdout = TRUE)
  feats <- sub(".*AvailableFeatures=([^ ]*).*", "\\1", grep("AvailableFeatures=", f, value = TRUE))
  hit <- intersect(strsplit(feats, ",")[[1]], c("sapphirerapids", "emeraldrapids", "icelake",
                                                 "cascadelake", "skylake", "broadwell", "haswell",
                                                 "genoa", "milan", "rome", "turin"))
  if (length(hit)) hit[1] else "unknown"
}, character(1))

# The job that wrote each saved object ---------------------------------------------------------
to_time <- function(x) as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%S", tz = "America/New_York")
obj$saved_mtime_posix <- file.mtime(file.path(PIPELINE, obj$saved))
if (anyNA(obj$saved_mtime_posix)) stop("saved objects missing: ", paste(obj$saved[is.na(obj$saved_mtime_posix)], collapse = ", "))
# Snakemake's main process touches each output once it sees the job finish (mtime = ctime, tens of
# seconds after sacct's end), so the modification time follows the writing job's end instead of
# falling inside it. The writer is the last job that completed; the file's time must come after
# its end, within MAX_TOUCH_S, and every other completed job of the object must have ended before
# the writer started.
MAX_TOUCH_S <- 1800
pick <- lapply(seq_len(nrow(obj)), function(i) {
  j <- jobs[jobs$key == obj$key[i], ]
  mt <- obj$saved_mtime_posix[i]
  done <- j$state == "COMPLETED" & j$exit == "0:0"
  if (!any(done)) stop(obj$object[i], ": no completed job")
  w <- j[done, ][which.max(to_time(j$end[done])), ]
  lag <- as.numeric(difftime(mt, to_time(w$end), units = "secs"))
  if (lag < 0 || lag > MAX_TOUCH_S) {
    stop(obj$object[i], ": the saved file's time ", format(mt), " is ", lag, " s after job ", w$job,
         " ended (", w$end, ")")
  }
  others <- j[j$job != w$job, ]
  od <- others$state == "COMPLETED"
  if (any(od) && any(to_time(others$end[od]) > to_time(w$start))) {
    stop(obj$object[i], ": another completed job overlaps the writer ", w$job)
  }
  list(w = w, lag = lag,
       superseded = if (nrow(others)) paste0(others$job, ":", others$state, collapse = ",") else "-")
})
lags <- vapply(pick, `[[`, 0, "lag")
cat("saved file time after the writer's end: median", median(lags), "s, max", max(lags), "s\n")
w <- do.call(rbind, lapply(pick, `[[`, "w"))
obj$workflow_job <- w$job
obj$workflow_run <- w$run
obj$superseded_jobs <- vapply(pick, `[[`, "", "superseded")
obj$workflow_node <- w$node
obj$workflow_cpu <- unname(cpu_of[w$node])
obj$workflow_partition <- w$partition
obj$workflow_start <- w$start
obj$workflow_end <- w$end
obj$workflow_elapsed <- w$elapsed
obj$workflow_max_rss_kib <- unname(rss[w$job])
if (anyNA(obj$workflow_max_rss_kib)) stop("no step .0 MaxRSS for some jobs")
obj$saved_mtime <- format(obj$saved_mtime_posix, "%Y-%m-%dT%H:%M:%S")

# The optimizer's numbers from each job's log ----------------------------------------------------
one <- function(txt, re, what, path) {
  hit <- grep(re, txt, value = TRUE)
  if (length(hit) != 1) stop(path, ": ", length(hit), " lines with ", what)
  sub(re, "\\1", hit)
}
lognum <- do.call(rbind, lapply(obj$log, function(path) {
  txt <- readLines(file.path(PIPELINE, path), warn = FALSE)
  if (sum(grepl("^Modularity Optimizer version 1\\.3\\.0", txt)) != 1) stop(path, ": not exactly one optimizer block")
  nodes <- one(txt, "^Number of nodes: ([0-9]+)$", "the nodes", path)
  edges <- one(txt, "^Number of edges: ([0-9]+)$", "the edges", path)
  q <- one(txt, "^Maximum modularity in 10 random starts: ([0-9.]+)$", "the modularity", path)
  comm <- one(txt, "^Number of communities: ([0-9]+)$", "the communities", path)
  sg <- grep("^[0-9]+ singletons identified\\. [0-9]+ final clusters\\.$", txt, value = TRUE)
  if (length(sg) > 1) stop(path, ": more than one singleton message")
  c(nodes, edges, q, comm,
    if (length(sg)) sub("^([0-9]+) singletons.*", "\\1", sg) else "0",
    if (length(sg)) sub(".*identified\\. ([0-9]+) final clusters\\.$", "\\1", sg) else comm)
}))
obj$log_nodes <- lognum[, 1]
obj$log_edges <- lognum[, 2]
obj$log_max_modularity <- lognum[, 3]
obj$log_communities <- lognum[, 4]
obj$log_singletons <- lognum[, 5]
obj$log_final_clusters <- lognum[, 6]

# Tier P (validation/clustering_datasets.tsv) ------------------------------------------------------
ds_path <- file.path(dirname(here), "clustering_datasets.tsv")
tierP <- if (file.exists(ds_path)) {
  ds <- read.delim(ds_path, colClasses = "character")
  ds$object[ds$tier == "P"]
} else character()
obj$tier <- ifelse(obj$object %in% tierP, "P", "F")

cols <- c("object", "tier", "sample", "kind", "rule", "assay", "input", "saved", "log", "markers",
          "workflow_job", "workflow_run", "superseded_jobs", "workflow_node", "workflow_cpu",
          "workflow_partition", "workflow_start", "workflow_end", "workflow_elapsed",
          "workflow_max_rss_kib", "saved_mtime", "log_nodes", "log_edges", "log_max_modularity",
          "log_communities", "log_singletons", "log_final_clusters")
tsv_write(obj[, cols], out_path)
cat("wrote", nrow(obj), "objects to", out_path, "\n")
cat("CPU types:", paste(names(table(obj$workflow_cpu)), table(obj$workflow_cpu), collapse = ", "), "\n")
cat("objects with superseded jobs:", sum(obj$superseded_jobs != "-"), "\n")
cat("objects with singletons:", sum(obj$log_singletons != "0"), "max", max(as.integer(obj$log_singletons)), "\n")
