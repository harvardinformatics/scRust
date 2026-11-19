#!/bin/bash
# Benchmark (scRust data/benchmarks): submits the R runs of r_runs.tsv, repeats <first>..<last>,
# each (tool, object) as a chain (rep k+1 starts after rep k ends: repeats of one object never
# share a node's page cache or Lustre stream). Records the job ids in logs/submitted.tsv.
#   bash submit_r.sh <first rep> <last rep> [<tool> [<object> ...]]
set -euo pipefail
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
first=$1 last=$2
shift 2
want_tool=${1:-}
[ $# -gt 0 ] && shift
want_objs=("$@")
cd "$S"
tail -n +2 validation/benchmark/r_runs.tsv | while IFS=$'\t' read -r tool obj rule script mem tlimit inputs args; do
    [ -n "$want_tool" ] && [ "$tool" != "$want_tool" ] && continue
    if [ ${#want_objs[@]} -gt 0 ]; then
        hit=0; for o in "${want_objs[@]}"; do [ "$o" = "$obj" ] && hit=1; done
        [ $hit = 1 ] || continue
    fi
    dep=
    for rep in $(seq "$first" "$last"); do
        j=$(sbatch --parsable --mem="$mem" --time="$tlimit" --job-name="bench-r-$tool" \
            ${dep:+--dependency=afterany:$dep} validation/benchmark/run_r.sbatch "$tool" "$obj" "$rep")
        printf '%s\tr\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$tool" "$obj" "$rep" "$j" "${dep:--}" | tee -a "$B/logs/submitted.tsv"
        dep=$j
    done
done
