#!/bin/bash
# Milestone M3 benchmark (docs/PLAN.md B5 "Performance", B6): wall time and peak RSS of the Rust
# DoubletFinder core on R's own embeddings, at 1 thread and at N threads, on
# soupx_emptydrops_sscovelli (79,409 cells; 105,879 with the artificial doublets of the final
# run), where R's fields::rdist matrix alone is 89.7 GB. From the repository's top directory,
# after `cargo build --release -p scrust`:
#   sbatch [--cpus-per-task=N] validation/slurm/m3_benchmark.sbatch [<object>] [<run>]
# Defaults: soupx_emptydrops_sscovelli, instrumented-portable-sapphirerapids. Each measurement
# runs under /usr/bin/time -v, into data/doubletfinder/<object>/rust/m3-benchmark/<job id>/:
#   knn_t<T>/   scrust knn: the final run's neighbour search alone (79,409 queries among
#               105,879 cells, k + 1 = 20,118 positions; doublet counts, no sort)
#   core_t<T>/  scrust doubletfinder --from-stages: the whole core (six sweep iterations, the
#               sweep summary and pK, nExp, the final run, the outputs), no stage dumps
# with T = 1 and T = N (SLURM_CPUS_PER_TASK). Every core run's cells.tsv, sweep_stats.tsv and
# bcmvn.tsv are checked against the R run's.
#
#SBATCH --job-name=scrust-m3-bench
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
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

obj=${1:-soupx_emptydrops_sscovelli}
run=${2:-instrumented-portable-sapphirerapids}
bin=target/release/scrust
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
base=data/doubletfinder/$obj
ref=$base/runs/$run
stages=$ref/stages
out=$base/rust/m3-benchmark/${SLURM_JOB_ID:-local}
[ -e "$out" ] && { echo "$out exists" >&2; exit 1; }
mkdir -p "$out"
n=${SLURM_CPUS_PER_TASK:-1}

# the final run's sizes, from the R run's summary
n_cells=$(sed -n 's/^    "n_cells": \([0-9]*\),$/\1/p' "$ref/summary.json")
k=$(sed -n 's/^    "k": \([0-9]*\),$/\1/p' "$ref/summary.json")
echo "$obj: $n_cells real cells, k = $k"

status=0
for t in 1 "$n"; do
    d=$out/knn_t$t
    mkdir -p "$d"
    echo "== knn, $t thread(s)"
    /usr/bin/time -v -o "$d/time_v.txt" "$bin" knn --from-stages "$stages" --context final \
        --queries "$n_cells" --neighbours $((k + 1)) --output count --out "$d/counts.npy" --threads "$t" 2> "$d/knn.log"
    grep -E "Elapsed|Maximum resident" "$d/time_v.txt"

    d=$out/core_t$t
    echo "== core, $t thread(s)"
    /usr/bin/time -v -o "$out/core_t$t.time_v.txt" "$bin" doubletfinder --counts "$base/input" \
        --clusters "$base/input/seurat_clusters.tsv" --from-stages "$stages" --out "$d" --threads "$t"
    mv "$out/core_t$t.time_v.txt" "$d/time_v.txt"
    grep -E "Elapsed|Maximum resident" "$d/time_v.txt"
    for f in cells.tsv sweep_stats.tsv bcmvn.tsv; do
        if cmp -s "$ref/$f" "$d/$f"; then echo "$f identical to R's"; else echo "$f DIFFERS from R's"; status=1; fi
    done
done
exit $status
