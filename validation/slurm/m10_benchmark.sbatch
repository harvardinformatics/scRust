#!/bin/bash
# Milestone M10's benchmark (docs/PLAN.md C5, "Performance"; C6): FindClusters on the two
# 79,409-cell objects, R against Rust, on one node. From the repository's top directory, after
# `cargo build --release -p scrust`:
#   sbatch validation/slurm/m10_benchmark.sbatch soupx_emptydrops_sscovelli filtered_seurat_emptydrops_sscovelli
# For each object, on its saved SCT_snn graph (data/clustering/<object>/saved/graphs/):
#   R     FindClusters(obj, random.seed = 12345) alone (cluster_references.R --mode bench): its time,
#         and VmHWM before and after; the whole Rscript under /usr/bin/time -v
#   Rust  scrust cluster --snn at the workflow's settings, one thread; and with the five
#         resolutions of the `resolutions` setting (0.2, 0.4, 0.8, 1.2, 2.0) at one thread and at
#         $SLURM_CPUS_PER_TASK threads, whose files must be byte-identical (summary.json: its result
#         section). Each under /usr/bin/time -v
# Outputs in data/clustering/<object>/rust/m10/bench/: r.json, r_time_v.txt, and per Rust run its
# output directory and <run>.time_v.txt; threads.txt. R's time for five resolutions is the
# `resolutions` reference's (cluster_refs/resolutions/summary.json, run.seconds), made on another
# node. SCRUST_BIN names the binary (default target/release/scrust).
#
#SBATCH --job-name=scrust-m10-bench
#SBATCH -A informatics
#SBATCH --partition=test
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary $bin sha256 $(sha256sum "$bin" | cut -d' ' -f1)"
threads=${SLURM_CPUS_PER_TASK:-1}
status=0

for obj in "$@"; do
    base=data/clustering/$obj
    out=$base/rust/m10/bench
    [ -e "$out" ] && { echo "skip $obj: $out exists"; continue; }
    tmp=$out.partial.${SLURM_JOB_ID:-local}
    rm -rf "$tmp"
    mkdir -p "$tmp"
    /usr/bin/time -v -o "$tmp/r_time_v.txt" \
        "$R_ENV/bin/Rscript" validation/R/cluster_references.R --mode bench --object "$obj" --out "$tmp/r.json"
    g=$base/saved/graphs
    /usr/bin/time -v -o "$tmp/workflow_t1.time_v.txt" \
        "$bin" cluster --snn "$g" --out "$tmp/workflow_t1" --seed 12345 --stats none --threads 1 2> "$tmp/workflow_t1.stderr.txt"
    for t in 1 "$threads"; do
        /usr/bin/time -v -o "$tmp/resolutions_t$t.time_v.txt" \
            "$bin" cluster --snn "$g" --out "$tmp/resolutions_t$t" --seed 12345 --resolution 0.2,0.4,0.8,1.2,2.0 \
            --stats none --threads "$t" 2> "$tmp/resolutions_t$t.stderr.txt"
    done
    same=1
    for f in "$tmp"/resolutions_t1/*; do
        n=$(basename "$f")
        [ "$n" = summary.json ] && continue
        cmp -s "$f" "$tmp/resolutions_t$threads/$n" || { echo "$n differs"; same=0; }
    done
    a=$(sed -n '/^  "result": {/,/^  },$/p' "$tmp/resolutions_t1/summary.json" | sha256sum)
    b=$(sed -n '/^  "result": {/,/^  },$/p' "$tmp/resolutions_t$threads/summary.json" | sha256sum)
    [ "$a" = "$b" ] || { echo "summary.json result differs"; same=0; }
    if [ "$same" = 1 ]; then
        echo "identical at 1 and $threads threads: every file, and summary.json's result section" > "$tmp/threads.txt"
    else
        echo "DIFFERS between 1 and $threads threads" > "$tmp/threads.txt"
        status=1
    fi
    cat "$tmp/threads.txt"
    for r in workflow_t1 resolutions_t1 "resolutions_t$threads"; do
        printf '%s %s: %s, peak %s KiB\n' "$obj" "$r" \
            "$(sed -n 's/^\tElapsed (wall clock) time (h:mm:ss or m:ss): //p' "$tmp/$r.time_v.txt")" \
            "$(sed -n 's/^\tMaximum resident set size (kbytes): //p' "$tmp/$r.time_v.txt")"
    done
    mv "$tmp" "$out"
done
echo "status $status"
exit $status
