#!/bin/bash
# Milestone M11 (docs/PLAN.md A5b mechanism 3): the R-vs-R floor across OpenBLAS kernels, measured
# with the statistics kernel-gated grading uses (scio::grade, `scrust validate floor`), on the
# instrumented R runs M2 and M8 made in both profiles on Sapphire Rapids:
#   - clustering: every object with instrumented-production-sapphirerapids and
#     instrumented-portable-sapphirerapids runs (27 tier-T, 7 tier-P);
#   - DoubletFinder: every object with both (5 tier-T, 3 tier-P).
# No R is run: the measurement reads the runs' stage dumps. From the repository's top directory,
# after `cargo build --release -p scrust`:
#   sbatch validation/slurm/m11_floor.sbatch
# Outputs (data/, listed in the manifest): data/clustering/floor_kernel_m11/{pairs,blocks}.tsv and
# summary.md; data/doubletfinder/compare/floor_kernel_m11/ likewise. summary.md is the table
# behind [clustering.floor.kernel] and [doubletfinder.floor.kernel] (validation/tolerances.toml).
#
#SBATCH --job-name=scrust-m11-floor
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=0-03:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary $bin sha256 $(sha256sum "$bin" | cut -d' ' -f1)"

# <tool data dir> <out dir> <object list>
measure() {
    local data=$1 out=$2
    shift 2
    mkdir -p "$out"
    {
        printf 'label\ta\tb\n'
        for o in "$@"; do
            local a=$data/$o/runs/instrumented-production-sapphirerapids/stages
            local b=$data/$o/runs/instrumented-portable-sapphirerapids/stages
            if [ -f "$a/index.tsv" ] && [ -f "$b/index.tsv" ]; then printf '%s\t%s\t%s\n' "$o" "$a" "$b"; fi
        done
    } > "$out/pairs.tsv"
    echo "$(($(wc -l < "$out/pairs.tsv") - 1)) pairs into $out"
    /usr/bin/time -v "$bin" validate floor --pairs "$out/pairs.tsv" --out "$out/blocks.tsv" --summary "$out/summary.md" 2> "$out/time_v.txt"
    cat "$out/summary.md"
}

clustering_objects=$(tail -n +2 validation/clustering_datasets.tsv | cut -f1)
measure data/clustering data/clustering/floor_kernel_m11 $clustering_objects
df_objects=$(ls data/doubletfinder | grep -v -e '^compare$' -e '^logs$')
measure data/doubletfinder data/doubletfinder/compare/floor_kernel_m11 $df_objects
