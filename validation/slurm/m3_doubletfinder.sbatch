#!/bin/bash
# Milestone M3 (docs/PLAN.md B6): the Rust DoubletFinder core on R's own embeddings (teacher
# forcing, plan A5b) for every instrumented reference run of the given objects, each compared
# with that same run. From the repository's top directory, after
# `cargo build --release -p scrust`:
#   sbatch [--cpus-per-task=N] validation/slurm/m3_doubletfinder.sbatch <object> [<object> ...]
# For each data/doubletfinder/<object>/runs/instrumented-<profile>-<arch>/:
#   data/doubletfinder/<object>/rust/m3/instrumented-<profile>-<arch>/
#     cells.tsv, sweep_stats.tsv, bcmvn.tsv, summary.json   Rust's outputs (reference layout)
#     stages/          Rust's S4 digests (and order positions where R dumped them), S5, S6, S7
#     time_v.txt       /usr/bin/time -v of the Rust run
#     validate.json, validate.md   scrust validate doubletfinder against the R run
# An existing output directory is never overwritten: its run is skipped. Work in progress goes
# to <dir>.partial.<job id>, renamed when complete. Exit status 1 if any comparison differs.
#
# The core uses no BLAS, so several threads are fine; the results do not depend on the count.
#
#SBATCH --job-name=scrust-m3-df
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
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

bin=target/release/scrust
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
threads=${SLURM_CPUS_PER_TASK:-1}
status=0

for obj in "$@"; do
    base=data/doubletfinder/$obj
    [ -d "$base/runs" ] || { echo "no runs for $obj" >&2; status=1; continue; }
    for r in "$base"/runs/instrumented-*/; do
        r=${r%/}
        run=$(basename "$r")
        case $run in *.partial.*) continue ;; esac
        out=$base/rust/m3/$run
        if [ -e "$out" ]; then
            echo "skip $out: exists"
            continue
        fi
        # full stage dumps (S4 order positions too) where R dumped the positions
        level=light
        if compgen -G "$r/stages/S4_positions.*.npy" > /dev/null; then level=full; fi
        tmp=$out.partial.${SLURM_JOB_ID:-local}
        mkdir -p "$(dirname "$out")"
        echo "== $obj $run (dump level $level, $threads threads)"
        /usr/bin/time -v -o "$tmp.time_v.txt" "$bin" doubletfinder \
            --counts "$base/input" --clusters "$base/input/seurat_clusters.tsv" \
            --from-stages "$r/stages" --out "$tmp" --dump-stages "$level" --threads "$threads"
        mv "$tmp.time_v.txt" "$tmp/time_v.txt"
        if ! "$bin" validate doubletfinder --reference "$r" --candidate "$tmp" \
                --json "$tmp/validate.json" --markdown "$tmp/validate.md" --threads "$threads"; then
            status=1
        fi
        mv "$tmp" "$out"
    done
done
exit $status
