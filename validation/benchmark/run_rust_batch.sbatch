#!/bin/bash
# Benchmark (scRust data/benchmarks): several of one object's Rust runs in one job, one after
# another (run_rust.sbatch, run with bash inside this allocation: the same settings, checks and
# outputs; each run evicts its inputs from the page cache first). One job per object instead of
# a chain of jobs, because every 33-CPU job waits in the queue.
#   sbatch --time=<t> --job-name=bench-rust-batch run_rust_batch.sbatch <object> <tool>:<rep> ...
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a[32608,24102,24105-24112,24201-24212,26105-26112,26201-26212]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=33
#SBATCH --mem=48G
#SBATCH --output=data/benchmarks/logs/%x-%j.log

set -uo pipefail
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
obj=$1
shift
status=0
for spec in "$@"; do
    tool=${spec%%:*} rep=${spec#*:}
    echo "===== $tool $obj rep$rep ====="
    bash "$S/validation/benchmark/run_rust.sbatch" "$tool" "$obj" "$rep" || status=1
done
exit $status
