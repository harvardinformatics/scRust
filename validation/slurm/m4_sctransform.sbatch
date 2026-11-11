#!/bin/bash
# Milestone M4 (docs/PLAN.md B6): Rust's SCTransform against every instrumented DoubletFinder
# reference run of the given objects, stage S2 in all 7 contexts, each sub-stage teacher-forced
# on R's own inputs (`scrust validate sctransform`). From the repository's top directory, after
# `cargo build --release -p scrust`:
#   sbatch [--mem=...] validation/slurm/m4_sctransform.sbatch <label> <object> [<object> ...]
# For each data/doubletfinder/<object>/runs/instrumented-<profile>-<arch>/:
#   data/doubletfinder/<object>/rust/<label>/instrumented-<profile>-<arch>/
#     validate.json, validate.md   the comparison (every check, per stage)
#     fingerprint.json             scrust fingerprint under the run's profile
#     time_v.txt                   /usr/bin/time -v of the validation
# An existing output directory is never overwritten: its run is skipped. Work in progress goes
# to <dir>.partial.<job id>, renamed when complete. Exit status 1 if any comparison differs.
#
# BLAS profile per run (A5b): each validation runs under the profile of the R run it checks:
# portable runs with OPENBLAS_CORETYPE=Haswell, production runs with the variable unset, so that
# OpenBLAS detects SapphireRapids (this job must run on a `sapphire` node). The fingerprint is
# checked before each run. One CPU: the residuals' products are OpenBLAS calls.
#
#SBATCH --job-name=scrust-m4-sct
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-06:00:00
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
[ $# -ge 2 ] || { echo "usage: m4_sctransform.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
threads=${SLURM_CPUS_PER_TASK:-1}
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
        if ! "${envset[@]}" /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" validate sctransform \
                --input "$base/input" --stages "$r/stages" \
                --json "$tmp/validate.json" --markdown "$tmp/validate.md" --threads "$threads"; then
            status=1
        fi
        mv "$tmp" "$out"
    done
done
exit $status
