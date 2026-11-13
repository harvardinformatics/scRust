#!/bin/bash
# Exports DoubletFinder's inputs and the compared fields of the workflow's saved outputs to plain
# files (validation/R/export_doubletfinder.R; docs/PLAN.md A4 and B4). Reads .rds files, so it
# runs only as a batch job. Objects and their pipeline paths: validation/doubletfinder_datasets.tsv
# (tiers T and P), else validation/doubletfinder_panel.tsv (the full panel, M7).
# Submit from the repository's top directory:
#   sbatch [--mem=..] [--time=..] validation/slurm/export_doubletfinder.sbatch <object> [<object> ...]
# Output: data/doubletfinder/<object>/{input,saved}/ and export_run.json; an existing export is
# never overwritten.
#
#SBATCH --job-name=scrust-df-export
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
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

PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
[ -x target/release/scrust ] || { echo "build scrust first: cargo build --release -p scrust" >&2; exit 2; }

# lookup <object> <column>: the column of the object's row, by header name, from the datasets
# table, else the full panel's
lookup() {
    local t
    for t in validation/doubletfinder_datasets.tsv validation/doubletfinder_panel.tsv; do
        awk -F'\t' -v o="$1" -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
            $1 == o && (c in h) { print $h[c]; found = 1; exit } END { exit !found }' "$t" && return 0
    done
    return 1
}

for obj in "$@"; do
    input=$(lookup "$obj" input) && saved=$(lookup "$obj" saved) || { echo "unknown object $obj" >&2; exit 2; }
    out=data/doubletfinder/$obj
    if [ -e "$out/input" ] || [ -e "$out/saved" ]; then
        echo "$out already has an export; not overwriting" >&2
        exit 1
    fi
    partial=$out/export.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"   # this job's own incomplete output, if it was requeued
    mkdir -p "$partial"
    /usr/bin/time -v -o "$partial/time_v.txt" \
        "$R_ENV/bin/Rscript" validation/R/export_doubletfinder.R --object "$obj" \
        --input "$PIPELINE/$input" --saved "$PIPELINE/$saved" --out "$partial"
    mv "$partial/input" "$partial/saved" "$out/"
    mv "$partial/export_run.json" "$out/export_run.json"
    mv "$partial/time_v.txt" "$out/export_time_v.txt"
    rmdir "$partial"
done
