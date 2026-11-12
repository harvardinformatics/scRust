#!/bin/bash
# Milestone M5 benchmark (docs/PLAN.md B5 "Performance", B6): RunPCA (npcs = 10, DoubletFinder's)
# on the scale.data of one DoubletFinder context, the same input for both sides, each under
# /usr/bin/time -v:
#   rust: scrust pca --sct <dir> --npcs 10 at 1 thread and at N threads (SLURM_CPUS_PER_TASK),
#         the two outputs compared byte for byte;
#   r:    validation/R/runpca_reference.R (Seurat 5.5.1's RunPCA() through DoubletFinder's path,
#         timed alone).
# The input is an SCTransform output directory in `scrust sctransform`'s layout; by default R's
# own SCTransform of the object's final run from the M4 benchmark,
# data/doubletfinder/<object>/rust/m4-benchmark/final/r/out (byte-identical to Rust's). From the
# repository's top directory, after `cargo build --release -p scrust` (pass a pinned copy of
# the binary as SCRUST_BIN, as for m5_pca.sbatch):
#   SCRUST_BIN=<copy> sbatch --cpus-per-task=32 --mem=32G validation/slurm/m5_benchmark.sbatch rust <object> [<sct dir>]
#   sbatch --cpus-per-task=1 --mem=48G validation/slurm/m5_benchmark.sbatch r <object> [<sct dir>]
# OpenBLAS runs the portable profile (OPENBLAS_CORETYPE=Haswell), the profile of the instrumented
# runs and of the M4 benchmark that made the input, so the Rust and R outputs compare byte for
# byte (compare them with cmp once both jobs are done). Outputs, never overwritten:
#   data/doubletfinder/<object>/rust/<BENCH_LABEL, default m5-benchmark>/final/{rust_t1,rust_t<N>,r}/
#     the scrust pca (or R) outputs in out/, time_v.txt; R adds out/timing.json (RunPCA alone)
#
#SBATCH --job-name=scrust-m5-bench
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
export OPENBLAS_CORETYPE=Haswell
source ./env.sh

side=${1:?usage: m5_benchmark.sbatch rust|r <object> [<sct dir>]}
obj=${2:?usage: m5_benchmark.sbatch rust|r <object> [<sct dir>]}
base=data/doubletfinder/$obj
sct=${3:-$base/rust/m4-benchmark/final/r/out}
[ -f "$sct/scale_data.npy" ] || { echo "no $sct/scale_data.npy" >&2; exit 1; }
out=$base/rust/${BENCH_LABEL:-m5-benchmark}/final
mkdir -p "$out"
n=${SLURM_CPUS_PER_TASK:-1}

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown) side=$side $obj sct=$sct"

status=0
case $side in
    rust)
        bin=${SCRUST_BIN:-target/release/scrust}
        [ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
        echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
        for t in 1 "$n"; do
            d=$out/rust_t$t
            [ -e "$d" ] && { echo "$d exists; skipped"; continue; }
            tmp=$d.partial.${SLURM_JOB_ID:-local}
            mkdir -p "$tmp"
            echo "== scrust pca, $t thread(s)"
            "$bin" fingerprint > "$tmp/fingerprint.json"
            /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" pca --sct "$sct" --out "$tmp/out" --npcs 10 --threads "$t" || status=1
            mv "$tmp" "$d"
            sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$d/time_v.txt"
        done
        if [ -d "$out/rust_t1/out" ] && [ -d "$out/rust_t$n/out" ]; then
            for f in embeddings.tsv loadings.tsv stdev.tsv random_seed.txt; do
                cmp "$out/rust_t1/out/$f" "$out/rust_t$n/out/$f" && echo "1 and $n threads: $f identical" || status=1
            done
        fi
        ;;
    r)
        R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
        d=$out/r
        [ -e "$d" ] && { echo "$d exists" >&2; exit 1; }
        tmp=$d.partial.${SLURM_JOB_ID:-local}
        mkdir -p "$tmp"
        echo "== R RunPCA"
        /usr/bin/time -v -o "$tmp/time_v.txt" env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB \
            PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" validation/R/runpca_reference.R "$sct" "$tmp/out" 10 || status=1
        mv "$tmp" "$d"
        sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$d/time_v.txt"
        ;;
    *) echo "side must be rust or r" >&2; exit 2 ;;
esac
if [ -d "$out/r/out" ]; then
    for t in "$out"/rust_t*/out; do
        [ -d "$t" ] || continue
        for f in embeddings.tsv loadings.tsv stdev.tsv random_seed.txt; do
            cmp "$out/r/out/$f" "$t/$f" && echo "R and ${t#"$out"/}: $f identical" || status=1
        done
    done
fi
exit $status
