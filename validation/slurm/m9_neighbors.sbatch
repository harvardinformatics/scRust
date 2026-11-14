#!/bin/bash
# Milestone M9 (docs/PLAN.md C5 K3-K6, C6): `scrust neighbors` on R's own embeddings (teacher
# forcing, plan A5b) for every R FindNeighbors run of the given objects, each compared with that
# run by `scrust validate neighbors`. From the repository's top directory, after
# `cargo build --release -p scrust`:
#   sbatch [--cpus-per-task=N] validation/slurm/m9_neighbors.sbatch <object> [<object> ...]
# The R runs, each with its own embeddings (runs/<run>/pca/embeddings.npy, columns 1:30):
#   data/clustering/<object>/runs/instrumented-*/           annoy.metric euclidean (M8's chain)
#   data/clustering/<object>/neighbors/<run>/<metric>/      cosine and manhattan, where
#                                                           neighbors_references.sbatch made them
# Outputs, data/clustering/<object>/rust/m9/<run>/<metric>/:
#   knn_index.tsv, knn_distance.tsv, nn.mtx.gz, snn.mtx.gz, summary.json   Rust's outputs
#   annoy.ann.sha256   the index's sha256 (the index itself is deleted once it has been found
#                      byte-identical to R's, which then holds the same bytes)
#   time_v.txt         /usr/bin/time -v of the Rust run
#   validate.json, validate.md   scrust validate neighbors against the R run
# EUCLIDEAN_ONLY=1 skips the cosine and manhattan references (e.g. while they are being made).
# An existing output directory is never overwritten: its run is skipped. Work in progress goes
# to <dir>.partial.<job id>, renamed when complete. Exit status 1 if any comparison differs.
#
# No BLAS is involved, so any CPU and several threads are fine; the results do not depend on
# either.
#
#SBATCH --job-name=scrust-m9-nn
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
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

bin=target/release/scrust
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
threads=${SLURM_CPUS_PER_TASK:-1}
status=0

# one Rust run against one R run: <object> <run> <metric> <R dir>
one() {
    local obj=$1 run=$2 metric=$3 ref=$4
    local base=data/clustering/$obj
    local out=$base/rust/m9/$run/$metric
    if [ -e "$out" ]; then
        echo "skip $out: exists"
        return 0
    fi
    local tmp=$out.partial.${SLURM_JOB_ID:-local}
    rm -rf "$tmp"
    mkdir -p "$(dirname "$out")"
    /usr/bin/time -v -o "$tmp.time_v.txt" \
        "$bin" neighbors --pca "$base/runs/$run/pca" --out "$tmp" --dims 1:30 --metric "$metric" \
        --annoy yes --threads "$threads"
    mv "$tmp.time_v.txt" "$tmp/time_v.txt"
    local ok=0
    "$bin" validate neighbors --reference "$ref" --candidate "$tmp" \
        --json "$tmp/validate.json" --markdown "$tmp/validate.md" || ok=1
    sha256sum "$tmp/annoy.ann" | sed "s|$tmp/||" > "$tmp/annoy.ann.sha256"
    if [ "$ok" -eq 0 ] && grep -q '"pass": true' "$tmp/validate.json"; then
        rm "$tmp/annoy.ann"
    else
        echo "DIFFERS: $obj $run $metric (index kept in $out)"
        status=1
    fi
    mv "$tmp" "$out"
}

for obj in "$@"; do
    base=data/clustering/$obj
    [ -d "$base/runs" ] || { echo "no runs for $obj" >&2; status=1; continue; }
    for r in "$base"/runs/instrumented-*/; do
        r=${r%/}
        run=$(basename "$r")
        case $run in *.partial.*) continue ;; esac
        one "$obj" "$run" euclidean "$r"
    done
    if [ -d "$base/neighbors" ] && [ "${EUCLIDEAN_ONLY:-0}" != 1 ]; then
        for d in "$base"/neighbors/*/*/; do
            d=${d%/}
            metric=$(basename "$d")
            run=$(basename "$(dirname "$d")")
            case $metric in *.partial.*) continue ;; esac
            one "$obj" "$run" "$metric" "$d"
        done
    fi
done
echo "status $status"
exit $status
