#!/bin/bash
# Milestone M5 (docs/PLAN.md B6): Rust's RunPCA (irlba) against every instrumented DoubletFinder
# reference run of the given objects, stage S3 in all 7 contexts (`scrust validate pca`):
# teacher-forced on R's PCA input (for the light runs, which did not dump it, Rust's S2d on R's
# S2c parameters) and after Rust's own SCTransform (S3-chain). From the repository's top
# directory, after `cargo build --release -p scrust` and copying the binary to a file of its own:
#   SCRUST_BIN=<copy> sbatch [--cpus-per-task=N] [--mem=...] validation/slurm/m5_pca.sbatch <label> <object> [<object> ...]
# For each data/doubletfinder/<object>/runs/instrumented-<profile>-<arch>/:
#   data/doubletfinder/<object>/rust/<label>/instrumented-<profile>-<arch>/
#     validate.json, validate.md   the comparison (every check, per stage)
#     sources.txt                  per context, where the teacher-forced PCA input came from
#     fingerprint.json             scrust fingerprint under the run's profile
#     time_v.txt                   /usr/bin/time -v of the validation
# An existing output directory is never overwritten: its run is skipped. Work in progress goes
# to <dir>.partial.<job id>, renamed when complete. Exit status 1 if any comparison differs.
#
# SCRUST_BIN (environment) names the scrust binary to run (default target/release/scrust): a
# pinned copy, so that later builds do not change it. PCA_MODE (environment) selects the checks:
# both (default), teacher or chained, passed to --mode.
#
# BLAS profile per run (A5b): each validation runs under the profile of the R run it checks:
# portable runs with OPENBLAS_CORETYPE=Haswell, production runs with the variable unset, so that
# OpenBLAS detects SapphireRapids (this job must run on a `sapphire` node). The fingerprint is
# checked before each run. OpenBLAS runs one thread; the products are split across the job's
# CPUs only where that keeps R's bits (crates/pca/src/products.rs), so the CPU count changes
# nothing but the time.
#
#SBATCH --job-name=scrust-m5-pca
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS PCA_MODE=${PCA_MODE:-both}" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -ge 2 ] || { echo "usage: m5_pca.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
threads=${SLURM_CPUS_PER_TASK:-1}
mode=${PCA_MODE:-both}
status=0

for obj in "$@"; do
    base=data/doubletfinder/$obj
    [ -d "$base/runs" ] || { echo "no runs for $obj" >&2; status=1; continue; }
    for r in "$base"/runs/instrumented-*/; do
        r=${r%/}
        run=$(basename "$r")
        case $run in *.partial.*) continue ;; esac
        out=$base/rust/$label/$run
        if [ -e "$out" ]; then
            echo "skip $out: exists"
            continue
        fi
        case $run in
            instrumented-portable-*) core=Haswell; envset=(env OPENBLAS_CORETYPE=Haswell) ;;
            instrumented-production-*) core=SapphireRapids; envset=(env -u OPENBLAS_CORETYPE) ;;
            *) echo "unknown profile in $run" >&2; status=1; continue ;;
        esac
        tmp=$out.partial.${SLURM_JOB_ID:-local}
        mkdir -p "$tmp"
        "${envset[@]}" "$bin" fingerprint > "$tmp/fingerprint.json"
        got=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
        if [ "$got" != "$core" ]; then
            echo "$run: OpenBLAS reports core '$got', expected $core; skipped" >&2
            status=1
            rm -rf "$tmp"
            continue
        fi
        echo "== $obj $run ($core)"
        if ! "${envset[@]}" /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" validate pca \
                --input "$base/input" --stages "$r/stages" --mode "$mode" \
                --json "$tmp/validate.json" --markdown "$tmp/validate.md" --threads "$threads" \
                2> "$tmp/stderr.txt"; then
            status=1
        fi
        cat "$tmp/stderr.txt" >&2
        sed -n 's/: teacher-forced PCA input from /\t/p' "$tmp/stderr.txt" > "$tmp/sources.txt"
        rm -f "$tmp/stderr.txt"
        mv "$tmp" "$out"
    done
done
exit $status
