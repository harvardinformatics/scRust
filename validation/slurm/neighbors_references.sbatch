#!/bin/bash
# Milestone M9 (docs/PLAN.md C4, C5 K3-K6 and C7): R references for Seurat's FindNeighbors that
# M8 did not make, by validation/R/neighbors_references.R in the production R env. From the
# repository's top directory, after `cargo build --release -p scrust` (for `scrust fingerprint`):
#   sbatch validation/slurm/neighbors_references.sbatch edge
#   sbatch validation/slurm/neighbors_references.sbatch real <object>/<run> [<object>/<run> ...]
# edge: the constructed edge cases, all three metrics, into data/neighbors/edge/ (cases.tsv,
#   expected.tsv, <case>/input/, <case>/<metric>/).
# real: annoy.metric = "cosine" and "manhattan" (set METRICS to change) from each instrumented
#   run's own embeddings, data/clustering/<object>/runs/<run>/, into
#   data/clustering/<object>/neighbors/<run>/<metric>/.
# No .rds is read and no BLAS is involved, so the profile and CPU type do not matter (env.sh's
# portable default is kept). Existing outputs are never overwritten: the job stops instead.
#
#SBATCH --job-name=scrust-m9-nnref
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=0-03:00:00
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
mkdir -p data/neighbors
mode=${1:?edge or real}
shift
case $mode in
    edge)
        out=data/neighbors/edge
        [ -e "$out" ] && { echo "$out exists: remove it first to regenerate" >&2; exit 1; }
        tmp=$out.partial.${SLURM_JOB_ID:-local}
        /usr/bin/time -v -o "data/neighbors/edge_time_v.${SLURM_JOB_ID:-local}.txt" \
            "$R_ENV/bin/Rscript" validation/R/neighbors_references.R --mode edge --out "$tmp"
        mv "$tmp" "$out"
        ;;
    real)
        metrics=${METRICS:-cosine,manhattan}
        runs=()
        for ob_run in "$@"; do
            r=data/clustering/${ob_run%%/*}/runs/${ob_run#*/}
            [ -d "$r/pca" ] || { echo "no run $r" >&2; exit 1; }
            for m in ${metrics//,/ }; do
                o=data/clustering/${ob_run%%/*}/neighbors/${ob_run#*/}/$m
                [ -e "$o" ] && { echo "$o exists" >&2; exit 1; }
            done
            runs+=("$r")
        done
        list=$(IFS=,; echo "${runs[*]}")
        /usr/bin/time -v -o "data/neighbors/real_time_v.${SLURM_JOB_ID:-local}.txt" \
            "$R_ENV/bin/Rscript" validation/R/neighbors_references.R --mode real --runs "$list" --metrics "$metrics"
        ;;
    *) echo "unknown mode $mode" >&2; exit 2 ;;
esac
echo "done"
