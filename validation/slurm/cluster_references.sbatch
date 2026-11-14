#!/bin/bash
# Milestone M10 (docs/PLAN.md C4, C5 K7 and C7): R references for Seurat's FindClusters that M8
# did not make, by validation/R/cluster_references.R in the production R env, each with the
# instrumented optimizer's per-start results (data/clustering/bin/ModularityOptimizer, built by
# validation/cpp/build.sh). From the repository's top directory, after
# `cargo build --release -p scrust` (for `scrust fingerprint`):
#   sbatch validation/slurm/cluster_references.sbatch settings <object> [<object> ...]
#   sbatch validation/slurm/cluster_references.sbatch edge
#   sbatch validation/slurm/cluster_references.sbatch bench <object> [<object> ...]
# settings: the settings of cluster_references.R (the workflow's call with one argument changed;
#   SETTINGS=a,b,... picks some) on the saved SCT_snn graph, into
#   data/clustering/<object>/cluster_refs/<setting>/. The tier (validation/clustering_datasets.tsv)
#   decides which settings run.
# edge: the constructed edge cases, into data/clustering/edge/ (cases.tsv, <case>/input/,
#   <case>/<setting>/).
# bench: FindClusters alone at the workflow's settings on the saved graph, timed, into
#   data/clustering/<object>/cluster_refs/bench.json and bench_time_v.txt.
# No .rds is read and no BLAS is involved, so the profile and CPU type do not matter (env.sh's
# portable default is kept). Existing outputs are never overwritten: the object is skipped.
#
#SBATCH --job-name=scrust-m10-ref
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
[ -x data/clustering/bin/ModularityOptimizer ] || { echo "build the optimizer first: bash validation/cpp/build.sh" >&2; exit 1; }
mode=${1:?settings, edge or bench}
shift
job=${SLURM_JOB_ID:-local}
status=0
case $mode in
    settings)
        for obj in "$@"; do
            tier=$(awk -F'\t' -v o="$obj" '$1 == o { print $2 }' validation/clustering_datasets.tsv)
            [ -n "$tier" ] || { echo "unknown object $obj" >&2; status=1; continue; }
            out=data/clustering/$obj/cluster_refs
            if [ -e "$out/done" ]; then
                echo "skip $obj: $out complete"
                continue
            fi
            [ -e "$out" ] && [ -z "${SETTINGS:-}" ] && { echo "skip $obj: $out exists but is incomplete; remove it" >&2; status=1; continue; }
            mkdir -p "$out"
            if /usr/bin/time -v -o "$out/time_v.$job.txt" \
                "$R_ENV/bin/Rscript" validation/R/cluster_references.R --mode settings --object "$obj" \
                --tier "$tier" --out "$out" --settings "${SETTINGS:-all}"; then
                [ -z "${SETTINGS:-}" ] && date -Iseconds > "$out/done"
            else
                echo "FAILED: $obj" >&2
                status=1
            fi
        done
        ;;
    edge)
        out=data/clustering/edge
        [ -e "$out" ] && { echo "$out exists: remove it first to regenerate" >&2; exit 1; }
        tmp=$out.partial.$job
        /usr/bin/time -v -o "data/clustering/edge_time_v.$job.txt" \
            "$R_ENV/bin/Rscript" validation/R/cluster_references.R --mode edge --out "$tmp"
        mv "$tmp" "$out"
        ;;
    bench)
        for obj in "$@"; do
            out=data/clustering/$obj/cluster_refs
            mkdir -p "$out"
            [ -e "$out/bench.json" ] && { echo "skip $obj: bench.json exists"; continue; }
            /usr/bin/time -v -o "$out/bench_time_v.txt" \
                "$R_ENV/bin/Rscript" validation/R/cluster_references.R --mode bench --object "$obj" --out "$out/bench.json"
        done
        ;;
    *) echo "unknown mode $mode" >&2; exit 2 ;;
esac
echo "status $status"
exit $status
