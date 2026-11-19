#!/bin/bash
# Milestone M15b: how the saved DoubletFinder objects' percent.mt was made
# (validation/R/percent_mt_sct.R): on each of the 150 saved DoubletFinder inputs
# (validation/doubletfinder_panel.tsv, column input), PercentageFeatureSet() on the SCT assay and on
# the RNA assay, at the singlets, against the saved DoubletFinder object's percent.mt (exported by
# validation/slurm/m15b_saved.sbatch). Reads .rds files; not a new R run of the workflow.
#   sbatch --array=1-150%40 validation/slurm/m15b_percent_mt_sct.sbatch [<out name> <scrust label>]
# Output: data/clustering/m15b/<out name, default percent_mt_sct>/<object>.json. With a scrust
# label, also against scrust's percent.mt of the same singlets
# (data/clustering/m15b/saved_rust/<label>/doubletfinder_<object>/run/singlets_clusters/metadata.tsv).
#
#SBATCH --job-name=scrust-m15b-sct
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=0-00:40:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
i=${SLURM_ARRAY_TASK_ID:?submit as an array job}
read -r obj input < <(awk -F'\t' -v i="$((i + 1))" 'NR == 1 { for (c = 1; c <= NF; c++) h[$c] = c; next }
    NR == i { print $h["object"], $h["input"] }' validation/doubletfinder_panel.tsv)
[ -n "${obj:-}" ] || { echo "no object $i" >&2; exit 2; }
name=${1:-percent_mt_sct}
rust=()
[ -n "${2:-}" ] && rust=(--rust "data/clustering/m15b/saved_rust/$2/doubletfinder_$obj/run/singlets_clusters/metadata.tsv")
out=data/clustering/m15b/$name/$obj.json
[ -e "$out" ] && { echo "skip $out"; exit 0; }
mkdir -p "data/clustering/m15b/$name"
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB "$R_ENV/bin/Rscript" validation/R/percent_mt_sct.R --rds "$PIPELINE/$input" \
    --saved "data/clustering/m15b/saved/doubletfinder_$obj/percent_mt.tsv" --out "$out.partial" "${rust[@]}"
mv "$out.partial" "$out"
