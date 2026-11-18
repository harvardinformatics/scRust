#!/bin/bash
# Benchmark (scRust data/benchmarks): submits the Rust runs (run_rust.sbatch), repeats
# <first>..<last> of <tool> for each object, as one chain per object (each job starts after the
# previous one ends; [after <job>] puts the chain after that job), so no two of these runs read
# the same input files at once. Records the job ids in logs/submitted.tsv.
#   bash submit_rust.sh <tool> <first rep> <last rep> <object>[:<after job>] ...
set -euo pipefail
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
tool=$1 first=$2 last=$3
shift 3
cd "$S"
for spec in "$@"; do
    obj=${spec%%:*}
    dep=
    [ "$spec" != "$obj" ] && dep=${spec#*:}
    for rep in $(seq "$first" "$last"); do
        j=$(sbatch --parsable --job-name="bench-rust-$tool" ${dep:+--dependency=afterany:$dep} \
            validation/benchmark/run_rust.sbatch "$tool" "$obj" "$rep")
        printf '%s\trust\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$tool" "$obj" "$rep" "$j" "${dep:--}" | tee -a "$B/logs/submitted.tsv"
        dep=$j
    done
done
