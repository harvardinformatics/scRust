#!/bin/bash
# One-time exports of clustered objects (validation/R/export_clustering.R; docs/PLAN.md C4,
# milestone M8): each object's chain input (its RNA counts, or a link to an earlier M2 or M6b
# export of the same counts) and what the chain left in it (labels and levels, silhouette and
# purity, the SCT assay and model, the PCA, both graphs), at the object's level, plus a copy of its
# combined marker table. Reads .rds files.
# Submit from the repository's top directory:
#   sbatch [--mem=..] [--time=..] validation/slurm/export_clustering.sbatch <object> [<object> ...]
# Objects, levels, pipeline paths and reused exports are in validation/clustering_datasets.tsv.
# Output: data/clustering/<object>/{input,saved}/ and export_run.json, written to
# data/clustering/<object>.partial.<job>/ first and moved into place when complete; an existing
# export is never overwritten. Exports need no particular CPU type (-p shared,sapphire starts
# sooner); peak RSS at most 11.9 GiB: 4 GB covers tier T, 24 GB every P object.
#
#SBATCH --job-name=scrust-cl-export
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-02:00:00
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
    level=$(cut -f4 <<< "$row")
    saved=$(cut -f6 <<< "$row")
    counts_from=$(cut -f7 <<< "$row")
    out=data/clustering/$obj
    if [ -e "$out/saved/object.json" ]; then
        echo "$out exists; not overwriting" >&2
        continue
    fi
    partial=data/clustering/$obj.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"
    mkdir -p "$partial"
    echo "== export $obj ($level) from $saved"
    "$R_ENV/bin/Rscript" validation/R/export_clustering.R --object "$obj" --rds "$PIPELINE/$saved" \
        --markers "$PIPELINE/${saved%.rds}_markergenes.csv" --out "$partial" --level "$level" \
        --counts-from "$counts_from"
    # input/'s links are relative, and the partial directory sits at the same depth as the final one
    mkdir -p "$out"
    for d in input saved; do
        [ ! -e "$out/$d" ] || { echo "$out/$d exists" >&2; exit 1; }
        mv "$partial/$d" "$out/$d"
    done
    mv "$partial/export_run.json" "$out/export_run.json"
    rmdir "$partial"
    for f in matrix.mtx.gz barcodes.tsv.gz features.tsv.gz; do
        [ -r "$out/input/$f" ] || { echo "$out/input/$f does not resolve" >&2; exit 1; }
    done
done
