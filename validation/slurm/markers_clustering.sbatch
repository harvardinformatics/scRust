#!/bin/bash
# The marker references (validation/R/markers_clustering.R; docs/PLAN.md C4, milestone M8): the
# workflow's find_markers.R and combine_markers.R steps on each saved object, compared byte for
# byte with its saved <prefix>_markergenes.csv. Reads .rds files.
# Submit from the repository's top directory:
#   sbatch [--mem=..] [--time=..] validation/slurm/markers_clustering.sbatch <object> [<object> ...]
# Objects and pipeline paths are in validation/clustering_datasets.tsv. Output:
# data/clustering/<object>/markers/ (markers.json, the combined table, reference_markergenes.csv,
# clusters/), written to a .partial directory first; an existing one is never overwritten.
# FindMarkers uses no BLAS, so any CPU type will do (-p shared,sapphire starts sooner). Peak RSS
# at most 9.6 GiB (cellbender_fromraw_ghirsutum2) and 4 minutes per object: 4 GB covers tier T,
# 16 GB every P object.
#
#SBATCH --job-name=scrust-cl-markers
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; exit 2; }
    saved=$(cut -f6 <<< "$row")
    out=data/clustering/$obj/markers
    if [ -e "$out" ]; then
        echo "$out exists; not overwriting" >&2
        continue
    fi
    partial=$out.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"
    echo "== markers of $obj"
    /usr/bin/time -v -o "$partial.time_v.txt" \
        "$R_ENV/bin/Rscript" validation/R/markers_clustering.R --object "$obj" --rds "$PIPELINE/$saved" \
        --saved-markers "$PIPELINE/${saved%.rds}_markergenes.csv" --out "$partial"
    mv "$partial.time_v.txt" "$partial/time_v.txt"
    mv "$partial" "$out"
done
