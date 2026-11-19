#!/bin/bash
# Milestone M15b: percent.mt as the production run saved it (validation/R/export_percent_mt.R), from
# the saved workflow objects: tasks 1-34 are the 34 clustering reference objects
# (validation/clustering_datasets.tsv: tier T and P; their counts in data/clustering/<object>/input),
# tasks 35-184 the 150 saved DoubletFinder objects, the singlets re-clustered by doubletfinder.R
# (validation/clustering_panel.tsv's doubletfinder rows; their counts in M7's
# data/doubletfinder/<object>/saved/singlets). Reads .rds files; not a new R run of the workflow.
# Submit from the repository's top directory:
#   sbatch --array=1-184%40 validation/slurm/m15b_saved.sbatch
# Output: data/clustering/m15b/saved/<object>/percent_mt.tsv and percent_mt.json.
#
#SBATCH --job-name=scrust-m15b-saved
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
i=${SLURM_ARRAY_TASK_ID:?submit as an array job}
if [ "$i" -le 34 ]; then
    read -r obj saved < <(awk -F'\t' -v i="$((i + 1))" 'NR == i { print $1, $6 }' validation/clustering_datasets.tsv)
    counts=data/clustering/$obj/input
else
    read -r obj saved < <(awk -F'\t' -v i="$((i - 34))" 'NR == 1 { for (c = 1; c <= NF; c++) h[$c] = c; next }
        $h["kind"] == "doubletfinder" { n++; if (n == i) print $h["object"], $h["saved"] }' validation/clustering_panel.tsv)
    counts=data/doubletfinder/${obj#doubletfinder_}/saved/singlets
fi
[ -n "${obj:-}" ] || { echo "no object for task $i" >&2; exit 2; }
out=data/clustering/m15b/saved/$obj
if [ -e "$out/percent_mt.json" ]; then
    echo "skip $out: exists"
    exit 0
fi
partial=$out.partial.${SLURM_JOB_ID:-$$}
rm -rf "$partial"
mkdir -p "$partial"
echo "== task $i: $obj from $saved (counts: $counts)"
/usr/bin/time -v -o "$partial/time_v.txt" "$R_ENV/bin/Rscript" validation/R/export_percent_mt.R \
    --rds "$PIPELINE/$saved" --out "$partial" --counts-from "$counts"
rm -rf "$out"
mv "$partial" "$out"
