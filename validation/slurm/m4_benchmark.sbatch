#!/bin/bash
# Milestone M4 benchmark (docs/PLAN.md B5 "Performance", B6): SCTransform on one context of a
# DoubletFinder run (by default the final run: all cells and its artificial doublets), the same
# input for both sides, each under /usr/bin/time -v:
#   rust: scrust sctransform --from-stages ... --context ... at 1 thread and at N threads
#         (SLURM_CPUS_PER_TASK), each output compared byte for byte with the other;
#   r:    validation/R/sctransform_benchmark.R (Seurat 5.5.1's SCTransform(), timed alone).
# From the repository's top directory, after `cargo build --release -p scrust` (copy the binary
# and pass it as SCRUST_BIN, as for m4_sctransform.sbatch):
#   sbatch --cpus-per-task=32 --mem=64G  validation/slurm/m4_benchmark.sbatch rust <object> [<run>] [<context>]
#   sbatch --cpus-per-task=1  --mem=250G validation/slurm/m4_benchmark.sbatch r    <object> [<run>] [<context>]
# Defaults: run instrumented-portable-sapphirerapids, context final. OpenBLAS runs under the
# run's profile (portable: OPENBLAS_CORETYPE=Haswell), so the Rust and R outputs compare byte for
# byte (compare them with cmp once both jobs are done). Outputs, never overwritten:
#   data/doubletfinder/<object>/rust/m4-benchmark/<context>/{rust_t1,rust_t<N>,r}/
#     the scrust sctransform (or R) outputs, time_v.txt; R adds timing.json (SCTransform alone)
#
#SBATCH --job-name=scrust-m4-bench
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

side=${1:?usage: m4_benchmark.sbatch rust|r <object> [<run>] [<context>]}
obj=${2:?usage: m4_benchmark.sbatch rust|r <object> [<run>] [<context>]}
run=${3:-instrumented-portable-sapphirerapids}
ctx=${4:-final}
base=data/doubletfinder/$obj
stages=$base/runs/$run/stages
[ -d "$stages" ] || { echo "no $stages" >&2; exit 1; }
case $run in
    instrumented-portable-*) export OPENBLAS_CORETYPE=Haswell ;;
    instrumented-production-*) unset OPENBLAS_CORETYPE ;;
    *) echo "unknown profile in $run" >&2; exit 1 ;;
esac
out=$base/rust/m4-benchmark/$ctx
mkdir -p "$out"
n=${SLURM_CPUS_PER_TASK:-1}

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown) side=$side $obj $run $ctx"

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
            echo "== scrust sctransform, $t thread(s)"
            "$bin" fingerprint > "$tmp/fingerprint.json"
            /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" sctransform --counts "$base/input" \
                --from-stages "$stages" --context "$ctx" --out "$tmp/out" --threads "$t" || status=1
            mv "$tmp" "$d"
            sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$d/time_v.txt"
        done
        if [ -d "$out/rust_t1/out" ] && [ -d "$out/rust_t$n/out" ]; then
            for f in cell_attr.tsv gene_attr.tsv model_pars.tsv model_pars_fit.tsv variable_features.tsv scale_data.npy; do
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
        echo "== R SCTransform"
        /usr/bin/time -v -o "$tmp/time_v.txt" env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB \
            PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" validation/R/sctransform_benchmark.R \
            "$base/input" "$stages" "$ctx" "$tmp/out" || status=1
        mv "$tmp" "$d"
        sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$d/time_v.txt"
        ;;
    *) echo "side must be rust or r" >&2; exit 2 ;;
esac
exit $status
