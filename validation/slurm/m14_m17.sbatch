#!/bin/bash
# Milestone M14: the R loader's comparison as milestone M17 will run it (docs/PLAN.md C4 "Full
# panel"), timed: one Rust export against its saved .rds through the Seurat loader, without the
# file and round-trip checks (validation/R/loader_compare.R --loaders seurat --files no
# --round-trip no), under /usr/bin/time -v, to size M17's sbatch arrays. For each <object> of
# validation/clustering_datasets.tsv, the export data/clustering/<object>/rust/<label>/<profile>/run.
# Reads .rds files. From the repository's top directory:
#   sbatch [--mem=..] validation/slurm/m14_m17.sbatch <label> <profile> <object> [<object> ...]
# Output: data/clustering/<object>/loaders/<label>/<profile>/m17/ (k11.tsv, compare.json, log.txt,
# time_v.txt). The objects run one after another, so each one's peak is its own.
#
#SBATCH --job-name=scrust-m14-m17
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 3 ] || { echo "usage: m14_m17.sbatch <label> <profile> <object> [<object> ...]" >&2; exit 2; }
label=$1
profile=$2
shift 2
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
status=0
for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; status=1; continue; }
    saved=$PIPELINE/$(cut -f6 <<< "$row")
    out=data/clustering/$obj/loaders/$label/$profile/m17
    rm -rf "$out"
    mkdir -p "$out"
    echo "== $obj"
    if ! env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 /usr/bin/time -v -o "$out/time_v.txt" \
            "$R_ENV/bin/Rscript" validation/R/loader_compare.R --object "$obj" --export "data/clustering/$obj/rust/$label/$profile/run" \
            --rds "$saved" --out "$out" --loaders seurat --files no --round-trip no > "$out/log.txt" 2>&1; then
        echo "$obj: comparison FAILED" >&2
        status=1
    fi
    tail -n 1 "$out/log.txt"
    grep -E "Elapsed|Maximum resident" "$out/time_v.txt"
done
exit $status
