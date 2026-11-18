#!/bin/bash
# Benchmark (scRust data/benchmarks): the workflow's own jobs for the benchmarked objects, from
# sacct: the job's step .0 (Snakemake's job step, which runs the rule's Rscript) and the job.
#   bash sacct_workflow.sh > sacct_workflow.tsv
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
printf 'tool\tobject\tjob\tstep\tnode\tpartition\talloc_cpus\treq_mem\telapsed\ttotal_cpu\tuser_cpu\tsystem_cpu\tmax_rss\tstart\tend\tstate\n'
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["tier"] == "P" {print "cluster", $c["object"], $c["workflow_job"]}' \
    "$S/validation/clustering_panel.tsv" > /tmp/bench_jobs.$$
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["tier"] == "P" {print "doubletfinder", $c["object"], $c["workflow_job"]}' \
    "$S/validation/doubletfinder_datasets.tsv" >> /tmp/bench_jobs.$$
while read -r tool obj job; do
    case $obj in soupx_tenx_SRR18256555|soupx_tenx_10xmousebrain|filtered_seurat_tenx_10xmousebrain|soupx_emptydrops_4500871lane|soupx_emptydrops_sscovelli|filtered_seurat_emptydrops_sscovelli|cellbender_fromraw_ghirsutum2) ;; *) continue ;; esac
    sacct -j "$job" -P -n --format=JobID,NodeList,Partition,AllocCPUS,ReqMem,Elapsed,TotalCPU,UserCPU,SystemCPU,MaxRSS,Start,End,State |
        awk -F'|' -v t="$tool" -v o="$obj" -v j="$job" 'BEGIN {OFS = "\t"} ($1 "") == (j "") || ($1 "") == (j ".0") {print t, o, j, (($1 "") == (j "") ? "job" : "step.0"), $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}'
done < /tmp/bench_jobs.$$
rm -f /tmp/bench_jobs.$$
