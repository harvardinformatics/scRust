#!/bin/bash
# Milestone M7 (docs/PLAN.md B4 and B6): the full panel's one-time exports. One array task per
# object of validation/doubletfinder_panel.tsv (task i: the table's data row i). For its object
# the task runs, in this job:
#   1. export_doubletfinder.sbatch (validation/R/export_doubletfinder.R): the input object's
#      counts, cluster labels and their level order, and the saved output's compared fields
#      (kept barcodes and pANN, the doublet set, the column names), unless
#      data/doubletfinder/<object>/ already holds that export;
#   2. export_doubletfinder_singlets.sbatch (validation/R/export_doubletfinder_singlets.R): the
#      saved output's RNA counts, the singlets', into saved/singlets/, unless they are there.
# Both read .rds files and need no particular CPU type, so `shared` (Cascade Lake) serves as well
# as `sapphire`. 24 GB covers every object (the largest input .rds is 2.3 GB; the tier-P
# exports peaked at 5.8 GB). From the repository's top directory, after
# `cargo build --release -p scrust` (the exports record `scrust fingerprint`):
#   sbatch --array=1-150 [--output=<dir>/%x-%A_%a.sbatch.log] validation/slurm/m7_export.sbatch
# The two wrappers run with bash inside this job, where their #SBATCH lines are comments; they
# never overwrite an export. Exit status: theirs.
#
#SBATCH --job-name=scrust-m7-export
#SBATCH -A informatics
#SBATCH --partition=shared,sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

panel=validation/doubletfinder_panel.tsv
i=${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch --array=1-150 validation/slurm/m7_export.sbatch}
obj=$(awk -F'\t' -v i="$i" 'NR == i + 1 { print $1 }' "$panel")
[ -n "$obj" ] || { echo "$panel has no data row $i" >&2; exit 2; }
base=data/doubletfinder/$obj
echo "== task $i: $obj"

if [ -f "$base/input/export.json" ] && [ -f "$base/saved/saved.json" ]; then
    echo "$base: input and saved fields already exported"
else
    bash validation/slurm/export_doubletfinder.sbatch "$obj"
fi
if [ -f "$base/saved/singlets/export.json" ]; then
    echo "$base/saved/singlets: already exported"
else
    bash validation/slurm/export_doubletfinder_singlets.sbatch "$obj"
fi
