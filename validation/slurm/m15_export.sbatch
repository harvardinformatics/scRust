#!/bin/bash
# Milestone M15 (docs/PLAN.md C5's K12 and the singlets' re-clustering): one-time exports of the
# 150 saved production DoubletFinder objects' clustering, the workflow's re-clustering of the
# singlets (doubletfinder.R after DoubletFinder): their metadata (seurat_clusters,
# SCT_snn_res.0.8, silhouette_width, neighborhood_purity, nCount_SCT, nFeature_SCT and the RNA
# columns) and the factors' levels, by validation/R/export_clustering.R at level "labels". Their
# chain input, the singlets' RNA counts, was exported by M7 (data/doubletfinder/<object>/saved/
# singlets/); input/ links to it, and the export checks the saved object's counts against it, slot
# for slot. M17's full panel reuses these exports. Reads .rds files.
#
# One array task per object: task i is the i-th row of validation/clustering_panel.tsv whose kind
# is doubletfinder (150 rows, in the table's order). Submit from the repository's top directory:
#   sbatch --array=1-150%50 [--output=<dir>/%x-%A_%a.sbatch.log] validation/slurm/m15_export.sbatch
# Output: data/clustering/<object>/{input,saved}/ and export_run.json (<object> as the panel names
# it, doubletfinder_<arm>_<sample>), written to data/clustering/<object>.partial.<job>/ first and
# moved into place when complete; an existing export is skipped. No particular CPU type is needed.
#
#SBATCH --job-name=scrust-m15-export
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
panel=validation/clustering_panel.tsv
i=${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch --array=1-150 validation/slurm/m15_export.sbatch}
# the i-th doubletfinder row: the object, its saved .rds and its marker table, by header name
read -r obj saved markers < <(awk -F'\t' -v i="$i" 'NR == 1 { for (c = 1; c <= NF; c++) h[$c] = c; next }
    $h["kind"] == "doubletfinder" { n++; if (n == i) print $h["object"], $h["saved"], $h["markers"] }' "$panel") || true
[ -n "${obj:-}" ] || { echo "$panel has no doubletfinder row $i" >&2; exit 2; }
# the DoubletFinder object's own name (M7's data/doubletfinder/<name>), without the prefix
df=${obj#doubletfinder_}
counts_from=data/doubletfinder/$df/saved/singlets
for f in matrix.mtx.gz barcodes.tsv.gz features.tsv.gz; do
    [ -f "$counts_from/$f" ] || { echo "$obj: no $counts_from/$f (M7's export of the saved singlets)" >&2; exit 1; }
done
out=data/clustering/$obj
if [ -e "$out/saved/object.json" ]; then
    echo "skip $out: exists"
    exit 0
fi
partial=data/clustering/$obj.partial.${SLURM_JOB_ID:-$$}
rm -rf "$partial"
mkdir -p "$partial"
echo "== task $i: export $obj (labels) from $saved"
/usr/bin/time -v -o "$partial/time_v.txt" "$R_ENV/bin/Rscript" validation/R/export_clustering.R --object "$obj" \
    --rds "$PIPELINE/$saved" --markers "$PIPELINE/$markers" --out "$partial" --level labels --counts-from "$counts_from"
sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$partial/time_v.txt"
# input/'s links are relative, and the partial directory sits at the same depth as the final one
mkdir -p "$out"
for d in input saved; do
    [ ! -e "$out/$d" ] || { echo "$out/$d exists" >&2; exit 1; }
    mv "$partial/$d" "$out/$d"
done
mv "$partial/export_run.json" "$out/export_run.json"
mv "$partial/time_v.txt" "$out/export_time_v.txt"
rmdir "$partial"
for f in matrix.mtx.gz barcodes.tsv.gz features.tsv.gz; do
    [ -r "$out/input/$f" ] || { echo "$out/input/$f does not resolve" >&2; exit 1; }
done
echo "done $obj"
