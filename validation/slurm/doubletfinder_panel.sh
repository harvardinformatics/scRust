#!/bin/bash
# Writes the full panel's table (milestone M7; docs/PLAN.md A4 and B4) to stdout: one row per
# production DoubletFinder object of the workflow, with its pipeline paths and the SLURM job that
# wrote its saved output. Reads the pipeline repo's Snakemake logs and samplesheet (read-only), and
# queries sacct and scontrol; reads no .rds, so a login node is enough. From the repository's top
# directory:
#   bash validation/slurm/doubletfinder_panel.sh <snakemake log> [<snakemake log> ...] > validation/doubletfinder_panel.tsv
# For M7, the two runs of the workflow's rerun #2:
#   .snakemake/log/2026-11-08T214133.605174.snakemake.log .snakemake/log/2026-11-09T220250.045638.snakemake.log
#
# - The jobs: each log's "Job <n> has been submitted with SLURM jobid <id> (log: .../rule_doubletfinder
#   [_cellbender]/<wildcards>/<id>.log)". An object submitted more than once (a retry, or a later
#   run) is an error: the table needs the one job that wrote the saved output.
# - Columns: object, tier (P if the object is in validation/doubletfinder_datasets.tsv, else F), arm,
#   sample, input, saved, log (pipeline-relative), workflow_job, workflow_run (the log's time stamp),
#   workflow_node (sacct NodeList), workflow_cpu (the node's CPU-type feature, from scontrol),
#   workflow_elapsed (sacct Elapsed), workflow_max_rss_kib (sacct MaxRSS of step .0, the srun that runs
#   doubletfinder.R, in KiB).
# - Rows: the samplesheet's samples in order (samples without DoubletFinder jobs, such as those the
#   config excludes, are skipped), then the arms soupx_tenx, soupx_emptydrops, decontx_raw_tenx,
#   decontx_raw_emptydrops, cellbender_fromraw.
# It stops if a job did not complete with exit code 0.

set -euo pipefail
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
[ $# -ge 1 ] || { echo "usage: doubletfinder_panel.sh <snakemake log> [<snakemake log> ...]" >&2; exit 2; }
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# object, job, run
for log in "$@"; do
    run=$(basename "$log" .snakemake.log)
    grep -o "has been submitted with SLURM jobid [0-9]* (log: [^)]*/rule_doubletfinder[a-z_]*/[^/]*/[0-9]*\.log)" "$PIPELINE/$log" |
        sed -E "s#.*SLURM jobid ([0-9]+) \(log: .*/(rule_doubletfinder[a-z_]*)/([^/]+)/[0-9]+\.log\)#\2\t\3\t\1\t$run#"
done | awk -F'\t' 'BEGIN { OFS = "\t" } { print ($1 == "rule_doubletfinder_cellbender" ? "cellbender_fromraw_" $2 : $2), $3, $4 }' > "$tmp/jobs.tsv"
dup=$(cut -f1 "$tmp/jobs.tsv" | sort | uniq -d)
[ -z "$dup" ] || { echo "objects with more than one DoubletFinder job: $dup" >&2; exit 1; }

# sacct: the allocation (node, state, exit code, elapsed) and step .0 (MaxRSS in KiB)
ids=$(cut -f2 "$tmp/jobs.tsv" | paste -sd,)
sacct -j "$ids" -X -n -P --format=JobID,NodeList,State,ExitCode,Elapsed > "$tmp/alloc.txt"
sacct -j "$ids" -n -P --units=K --format=JobID,MaxRSS | awk -F'|' '$1 ~ /\.0$/ { sub(/\.0$/, "", $1); sub(/K$/, "", $2); print $1 "|" $2 }' > "$tmp/rss.txt"
bad=$(awk -F'|' '$3 != "COMPLETED" || $4 != "0:0"' "$tmp/alloc.txt")
[ -z "$bad" ] || { echo "jobs that did not complete with exit code 0: $bad" >&2; exit 1; }
for node in $(cut -d'|' -f2 "$tmp/alloc.txt" | sort -u); do
    feats=$(scontrol show node "$node" | sed -n 's/.*AvailableFeatures=\([^ ]*\).*/\1/p')
    cpu=$(tr ',' '\n' <<< "$feats" | grep -xE 'sapphirerapids|emeraldrapids|icelake|cascadelake|skylake|genoa|milan|rome' | head -n 1 || true)
    echo "$node|${cpu:-unknown}"
done > "$tmp/cpu.txt"

tail -n +2 "$PIPELINE/samplesheet.tsv" | cut -f1 > "$tmp/samples.txt"
awk -F'\t' 'NR > 1 { print $1 }' validation/doubletfinder_datasets.tsv > "$tmp/datasets.txt"
awk -F'\t' -v OFS='\t' '
    FILENAME ~ /cpu.txt$/      { split($0, f, "|"); cpu[f[1]] = f[2]; next }
    FILENAME ~ /alloc.txt$/    { split($0, f, "|"); node[f[1]] = f[2]; el[f[1]] = f[5]; next }
    FILENAME ~ /rss.txt$/      { split($0, f, "|"); rss[f[1]] = f[2]; next }
    FILENAME ~ /datasets.txt$/ { ds[$1] = 1; next }
    FILENAME ~ /jobs.tsv$/     { job[$1] = $2; run[$1] = $3; next }
    FILENAME ~ /samples.txt$/  { order[++ns] = $1; next }
    END {
        print "object", "tier", "arm", "sample", "input", "saved", "log", "workflow_job", "workflow_run",
              "workflow_node", "workflow_cpu", "workflow_elapsed", "workflow_max_rss_kib"
        na = split("soupx_tenx soupx_emptydrops decontx_raw_tenx decontx_raw_emptydrops cellbender_fromraw", arms, " ")
        for (s = 1; s <= ns; s++) for (a = 1; a <= na; a++) {
            o = arms[a] "_" order[s]
            if (!(o in job)) continue
            arm = arms[a]
            if (arm == "cellbender_fromraw") input = "results/cellbender_fromraw/seurat_cellbender_fromraw_" order[s] ".rds"
            else { d = (arm ~ /^soupx/) ? "soupx" : "decontx_raw"; input = "results/" d "/seurat_" o ".rds" }
            j = job[o]
            print o, (o in ds) ? "P" : "F", arm, order[s], input, "results/doubletfinder/seurat_doubletfinder_" o ".rds",
                  "results/logs/doubletfinder/doubletfinder_" o ".log", j, run[o], node[j], cpu[node[j]], el[j], rss[j]
            n++
        }
        if (n != length(job)) { print "jobs of objects not in the samplesheet x arms: " length(job) - n > "/dev/stderr"; exit 1 }
    }' "$tmp/cpu.txt" "$tmp/alloc.txt" "$tmp/rss.txt" "$tmp/datasets.txt" "$tmp/jobs.tsv" "$tmp/samples.txt"
