#!/bin/bash
# The portable profile's SCT scale.data against the saved object's (validation/R/floor_scale_data.R;
# M8's R-vs-R floor) for objects whose reference runs kept scale.data only as a digest (dump level
# light). Reads .rds files. Portable profile (OPENBLAS_CORETYPE=Haswell, env.sh's default) on
# Sapphire Rapids, like the portable reference runs it stands for.
# Submit from the repository's top directory:
#   sbatch [--partition=..] [--mem=..] validation/slurm/floor_scale_data.sbatch <object> [<object> ...]
# Output: data/clustering/<object>/floor/scale_data_portable.json; an existing one is not
# overwritten. Memory as the reference runs: 96 GB for the 79,409-cell objects, 160 GB for
# cellbender_fromraw_ghirsutum2. M8 also ran it on soupx_tenx_SRR18256555 (8 GB), whose runs kept
# scale.data in full, as a cross-check: its numbers equal the element-by-element comparison of
# the saved and run scale_data.npy.
#
#SBATCH --job-name=scrust-cl-floor
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
export SCRUST_PROFILE=portable
fp=$(target/release/scrust fingerprint)
core=$(sed -n 's/^  "openblas_corename": "\(.*\)",$/\1/p' <<< "$fp")
[ "$core" = Haswell ] || { echo "the portable profile needs OpenBLAS core Haswell, got '$core'" >&2; exit 3; }
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' core=$core"

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; exit 2; }
    out=data/clustering/$obj/floor/scale_data_portable.json
    if [ -e "$out" ]; then
        echo "$out exists; not overwriting" >&2
        continue
    fi
    mkdir -p "data/clustering/$obj/floor"
    /usr/bin/time -v -o "$out.time_v.partial" \
        "$R_ENV/bin/Rscript" validation/R/floor_scale_data.R --object "$obj" --rds "$PIPELINE/$(cut -f6 <<< "$row")" \
        --input "data/clustering/$obj/input" --out "$out.partial"
    mv "$out.time_v.partial" "data/clustering/$obj/floor/time_v.txt"
    mv "$out.partial" "$out"
done
