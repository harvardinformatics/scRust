#!/bin/bash
# Milestone M9 (docs/PLAN.md C5, "Threads" and "Performance"): FindNeighbors on the 79,409-cell
# objects, R against Rust, on one node. From the repository's top directory, after
# `cargo build --release -p scrust`:
#   sbatch validation/slurm/m9_benchmark.sbatch <object> [<object> ...]
# For each object, on its production instrumented run's embeddings (runs/<RUN>/pca, columns 1:30;
# RUN defaults to instrumented-production-sapphirerapids), into data/clustering/<object>/rust/m9-bench/:
#   r/            R's FindNeighbors alone, the workflow's call (neighbors_references.R --mode bench):
#                 findneighbors.json (its time, VmHWM before and after), time_v.txt
#   t1/, t<N>/    scrust neighbors with 1 and N threads (N = SLURM_CPUS_PER_TASK), no index file:
#                 the outputs, summary.json, time_v.txt
#   a1/, a<N>/    the same with --annoy yes, compared file by file (threads.txt): the outputs,
#                 annoy.ann included, must not depend on the thread count
# BLAS is not involved; the CPU type is recorded in the log and in each summary.json.
#
#SBATCH --job-name=scrust-m9-bench
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary $bin sha256 $(sha256sum "$bin" | cut -d' ' -f1)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
run=${RUN:-instrumented-production-sapphirerapids}
nt=${SLURM_CPUS_PER_TASK:-1}
status=0

for obj in "$@"; do
    base=data/clustering/$obj
    pca=$base/runs/$run/pca
    out=$base/rust/m9-bench
    [ -e "$out" ] && { echo "$out exists" >&2; exit 1; }
    mkdir -p "$out/r"
    /usr/bin/time -v -o "$out/r/time_v.txt" "$R_ENV/bin/Rscript" validation/R/neighbors_references.R \
        --mode bench --runs "$base/runs/$run" --out "$out/r/findneighbors.json"
    for t in 1 "$nt"; do
        /usr/bin/time -v -o "$out/time_v.t$t.txt" "$bin" neighbors --pca "$pca" --out "$out/t$t" --dims 1:30 --threads "$t"
        mv "$out/time_v.t$t.txt" "$out/t$t/time_v.txt"
        /usr/bin/time -v -o "$out/time_v.a$t.txt" "$bin" neighbors --pca "$pca" --out "$out/a$t" --dims 1:30 --annoy yes --threads "$t"
        mv "$out/time_v.a$t.txt" "$out/a$t/time_v.txt"
    done
    : > "$out/threads.txt"
    for f in knn_index.tsv knn_distance.tsv nn.mtx.gz snn.mtx.gz annoy.ann; do
        if cmp -s "$out/a1/$f" "$out/a$nt/$f"; then echo "$f identical at 1 and $nt threads"; else echo "$f DIFFERS at 1 and $nt threads"; fi >> "$out/threads.txt"
    done
    for f in knn_index.tsv knn_distance.tsv nn.mtx.gz snn.mtx.gz; do
        if cmp -s "$out/t1/$f" "$out/a1/$f"; then echo "$f identical with and without the index file"; else echo "$f DIFFERS with the index file"; fi >> "$out/threads.txt"
    done
    echo "annoy.ann sha256 $(sha256sum "$out/a1/annoy.ann" | cut -d' ' -f1), R's $(sha256sum "$base/runs/$run/stages/annoy.ann" | cut -d' ' -f1)" >> "$out/threads.txt"
    cat "$out/threads.txt"
    if grep -q DIFFERS "$out/threads.txt"; then status=1; fi
    # the index files duplicate R's (compared above); keep their digests only
    rm "$out/a1/annoy.ann" "$out/a$nt/annoy.ann"
done
echo "status $status"
exit $status
