#!/bin/bash
# Exports the singlets' counts of the workflow's saved DoubletFinder outputs (milestone M6b;
# validation/R/export_doubletfinder_singlets.R) next to milestone M2's export of the same
# objects, which must exist. Reads .rds files, so it runs only as a batch job. Objects and their
# pipeline paths: validation/doubletfinder_datasets.tsv. From the repository's top directory:
#   sbatch [--mem=..] [--time=..] validation/slurm/export_doubletfinder_singlets.sbatch <object> [<object> ...]
# Output: data/doubletfinder/<object>/saved/singlets/ (matrix.mtx.gz, barcodes.tsv.gz,
# features.tsv.gz, export.json), export_singlets_run.json and export_singlets_time_v.txt. An
# existing export is never overwritten.
#
#SBATCH --job-name=scrust-df-export-singlets
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
status=0

for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/doubletfinder_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; status=2; continue; }
    saved=$(cut -f6 <<< "$row")
    base=data/doubletfinder/$obj
    if [ ! -f "$base/saved/saved.json" ] || [ ! -f "$base/input/export.json" ]; then
        echo "$base: no M2 export (validation/slurm/export_doubletfinder.sbatch)" >&2
        status=1
        continue
    fi
    if [ -e "$base/saved/singlets" ]; then
        echo "$base/saved/singlets already exists; not overwriting" >&2
        status=1
        continue
    fi
    partial=$base/export_singlets.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"   # this job's own incomplete output, if it was requeued
    mkdir -p "$partial"
    echo "== $obj"
    if ! /usr/bin/time -v -o "$partial/time_v.txt" \
            "$R_ENV/bin/Rscript" validation/R/export_doubletfinder_singlets.R --object "$obj" \
            --saved "$PIPELINE/$saved" --export "$base" --out "$partial"; then
        echo "$obj: export failed; partial output left in $partial" >&2
        status=1
        continue
    fi
    mv "$partial/singlets" "$base/saved/singlets"
    mv "$partial/export_singlets_run.json" "$base/export_singlets_run.json"
    mv "$partial/time_v.txt" "$base/export_singlets_time_v.txt"
    rmdir "$partial"
done
exit $status
