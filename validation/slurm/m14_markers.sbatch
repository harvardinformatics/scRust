#!/bin/bash
# Milestone M14 (docs/PLAN.md C4 "Markers", C5's K13), as an sbatch array: the workflow's marker
# steps on the Seurat loader's object (validation/R/markers_from_loader.R), one object per array
# task, built from the export of an earlier validation run (label <label>, profile <profile>):
# data/clustering/<object>/rust/<label>/<profile>/run. The combined table must equal the saved
# <prefix>_markergenes.csv byte for byte. Reads no .rds. From the repository's top directory:
#   sbatch --array=0-<n-1> [--mem=..] validation/slurm/m14_markers.sbatch <label> <profile> <object> [<object> ...]
# Task i takes the i-th object. Output: data/clustering/<object>/loaders/<label>/<profile>/markers/
# (markers.json, the combined table, clusters/), written to a .partial directory first.
# FindMarkers uses no BLAS. Memory: 4 GB covers tier T; the largest P objects need 24 GB
# (the loaded object and presto's working copies).
#
#SBATCH --job-name=scrust-m14-markers
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 3 ] || { echo "usage: m14_markers.sbatch <label> <profile> <object> [<object> ...]" >&2; exit 2; }
label=$1
profile=$2
shift 2
objects=("$@")
i=${SLURM_ARRAY_TASK_ID:-0}
[ "$i" -lt "${#objects[@]}" ] || { echo "array task $i: only ${#objects[@]} objects" >&2; exit 2; }
obj=${objects[$i]}
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
[ -n "$row" ] || { echo "unknown object $obj" >&2; exit 2; }
saved=$(cut -f6 <<< "$row")
export_dir=data/clustering/$obj/rust/$label/$profile/run
out=data/clustering/$obj/loaders/$label/$profile/markers
partial=$out.partial.${SLURM_JOB_ID:-$$}
mkdir -p "$(dirname "$out")"
rm -rf "$partial"
echo "== K13: markers of $obj from the Seurat loader ($label, $profile)"
status=0
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 /usr/bin/time -v -o "$partial.time_v.txt" \
    "$R_ENV/bin/Rscript" validation/R/markers_from_loader.R --object "$obj" --export "$export_dir" \
    --saved-markers "$PIPELINE/${saved%.rds}_markergenes.csv" --out "$partial" || status=$?
if [ -d "$partial" ]; then
    mv "$partial.time_v.txt" "$partial/time_v.txt"
    rm -rf "$out"
    mv "$partial" "$out"
else
    echo "$obj: no output (status $status)" >&2
    status=1
fi
exit $status
