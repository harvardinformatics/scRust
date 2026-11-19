#!/bin/bash
# Milestone M15b: scrust's percent.mt by the workflow's pattern, without regression, against the
# saved workflow objects' (validation/R/export_percent_mt.R; validation/slurm/m15b_saved.sbatch).
# Tasks 1-34: the 34 clustering reference objects, scrust cluster --counts
# data/clustering/<object>/input --preset workflow --percent-mt --stats none (percent.mt is
# computed from the counts before SCTransform, as the workflow's scripts do). Tasks 35-184: the 150
# DoubletFinder objects, whose saved objects are the singlets re-clustered by doubletfinder.R (which
# computes percent.mt on them before its SCTransform): scrust doubletfinder --counts <input>
# --clusters <input's labels> --recluster --recluster-params <percent_mt = true>, in the production
# profile on Sapphire Rapids (the singlets are DoubletFinder's, M7's criteria). Then
# validation/python/m15b_percent_mt.py compares each percent.mt column with the saved one.
#   SCRUST_BIN=<copy> sbatch --array=1-184%40 validation/slurm/m15b_saved_rust.sbatch <label>
# Output, never overwritten: data/clustering/m15b/saved_rust/<label>/<object>/: run/ (for the
# DoubletFinder objects, run/singlets_clusters/ holds the re-clustering), time_v.txt, stderr.txt,
# fingerprint.json, compare.json.
#
#SBATCH --job-name=scrust-m15b-pm
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
root=$PWD
label=${1:?usage: m15b_saved_rust.sbatch <label>}
export OPENBLAS_CORETYPE=auto
source ./env.sh
bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' kernel=$core binary=$bin"
[ "$core" = SapphireRapids ] || { echo "the production profile needs SapphireRapids, got '$core'" >&2; exit 3; }
threads=${SLURM_CPUS_PER_TASK:-1}
i=${SLURM_ARRAY_TASK_ID:?submit as an array job}
if [ "$i" -le 34 ]; then
    obj=$(awk -F'\t' -v i="$((i + 1))" 'NR == i { print $1 }' validation/clustering_datasets.tsv)
else
    obj=$(awk -F'\t' -v i="$((i - 34))" 'NR == 1 { for (c = 1; c <= NF; c++) h[$c] = c; next }
        $h["kind"] == "doubletfinder" { n++; if (n == i) print $h["object"] }' validation/clustering_panel.tsv)
fi
[ -n "$obj" ] || { echo "no object for task $i" >&2; exit 2; }
out=data/clustering/m15b/saved_rust/$label/$obj
[ -e "$out" ] && { echo "skip $out: exists"; exit 0; }
tmp=$out.partial.${SLURM_JOB_ID:-$$}
rm -rf "$tmp"
mkdir -p "$tmp"
"$bin" fingerprint > "$tmp/fingerprint.json"
if [ "$i" -le 34 ]; then
    args=(cluster --counts "data/clustering/$obj/input" --preset workflow --percent-mt --stats none --out "$tmp/run" --threads "$threads")
    meta=$tmp/run/metadata.tsv
else
    df=${obj#doubletfinder_}
    printf 'preset = "workflow"\n[sctransform]\npercent_mt = true\n' > "$tmp/percent_mt.toml"
    args=(doubletfinder --counts "data/doubletfinder/$df/input" --clusters "data/doubletfinder/$df/input/seurat_clusters.tsv"
          --out "$tmp/run" --recluster --recluster-params "$tmp/percent_mt.toml" --threads "$threads")
    meta=$tmp/run/singlets_clusters/metadata.tsv
fi
echo "== $obj: scrust ${args[*]}"
/usr/bin/time -v -o "$tmp/time_v.txt" "$bin" "${args[@]}" 2> "$tmp/stderr.txt" || { cat "$tmp/stderr.txt" >&2; exit 1; }
tail -n 2 "$tmp/stderr.txt"
# the large intermediates are not needed: only metadata.tsv and the summaries are compared
rm -rf "$tmp/run/pca" "$tmp/run/neighbors" "$tmp/run/singlets_clusters/pca" "$tmp/run/singlets_clusters/neighbors"
status=0
PYTHONNOUSERSITE=1 "$root/.envs/loaders-py/bin/python" -I validation/python/m15b_percent_mt.py one \
    "$meta" "data/clustering/m15b/saved/$obj" --json "$tmp/compare.json" || status=1
mv "$tmp" "$out"
exit $status
