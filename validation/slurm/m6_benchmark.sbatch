#!/bin/bash
# Milestone M6 benchmark (docs/PLAN.md B5, "Performance"): the whole `scrust doubletfinder` on an
# object's exported input (the counts and the cluster labels, the same files the direct R run
# read), at each given thread count, under /usr/bin/time -v, in the production profile:
# OPENBLAS_CORETYPE unset on a Sapphire Rapids node, as the direct R run of M2
# (runs/direct-production-sapphirerapids, whose time_v.txt holds R's wall time and peak RSS).
# Each run's cells.tsv, sweep_stats.tsv, bcmvn.tsv and summary.json result are compared with that
# R run's (`scrust validate doubletfinder`), and the runs with each other byte for byte (from M6b
# on, singlets/ too). From the
# repository's top directory, with a pinned copy of the binary (as for m6_doubletfinder.sbatch):
#   SCRUST_BIN=<copy> sbatch --cpus-per-task=32 --mem=32G validation/slurm/m6_benchmark.sbatch <label> <object> <threads> [<threads> ...]
# Output, never overwritten: data/doubletfinder/<object>/rust/<label>/t<threads>/
#   run/                    the outputs (no stage dumps)
#   time_v.txt, stderr.txt  /usr/bin/time -v, the run's messages
#   fingerprint.json        scrust fingerprint (the OpenBLAS kernel must be SapphireRapids)
#   vs_direct.json, .md     against the direct R run
# The job's other CPUs stay idle while a run uses fewer threads, so each run has its cores.
#
#SBATCH --job-name=scrust-m6-bench
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
# production profile: OpenBLAS detects the kernel (env.sh keeps it unset for "auto")
export OPENBLAS_CORETYPE=auto
source ./env.sh

[ $# -ge 3 ] || { echo "usage: m6_benchmark.sbatch <label> <object> <threads> [<threads> ...]" >&2; exit 2; }
label=$1
obj=$2
shift 2
base=data/doubletfinder/$obj
ref=$base/runs/direct-production-sapphirerapids
[ -f "$ref/cells.tsv" ] || { echo "no $ref" >&2; exit 1; }
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown) $obj"
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
echo "R (direct production run): $(sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: //p' "$ref/time_v.txt" | paste -sd' ')"

status=0
done_runs=()
for t in "$@"; do
    d=$base/rust/$label/t$t
    if [ -e "$d" ]; then
        echo "$d exists; skipped"
        done_runs+=("$d")
        continue
    fi
    tmp=$d.partial.${SLURM_JOB_ID:-local}
    mkdir -p "$tmp"
    "$bin" fingerprint > "$tmp/fingerprint.json"
    core=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
    if [ "$core" != SapphireRapids ]; then
        echo "OpenBLAS reports core '$core', expected SapphireRapids" >&2
        rm -rf "$tmp"
        exit 1
    fi
    echo "== scrust doubletfinder, $t thread(s)"
    if ! /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" doubletfinder --counts "$base/input" \
            --clusters "$base/input/seurat_clusters.tsv" --out "$tmp/run" --threads "$t" 2> "$tmp/stderr.txt"; then
        cat "$tmp/stderr.txt" >&2
        status=1
        continue
    fi
    cat "$tmp/stderr.txt"
    sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$tmp/time_v.txt"
    "$bin" validate doubletfinder --reference "$ref" --candidate "$tmp/run" \
        --json "$tmp/vs_direct.json" --markdown "$tmp/vs_direct.md" || status=1
    mv "$tmp" "$d"
    done_runs+=("$d")
done
# the runs against each other, byte for byte
first=${done_runs[0]:-}
for d in "${done_runs[@]:1}"; do
    for f in cells.tsv sweep_stats.tsv bcmvn.tsv singlets/matrix.mtx.gz singlets/barcodes.tsv.gz singlets/features.tsv.gz; do
        [ -e "$first/run/$f" ] || continue   # singlets/ from M6b on
        cmp "$first/run/$f" "$d/run/$f" && echo "${first##*/} and ${d##*/}: $f identical" || status=1
    done
done
exit $status
