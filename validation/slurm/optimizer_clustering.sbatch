#!/bin/bash
# The instrumented modularity optimizer (validation/R/optimizer_clustering.R; docs/PLAN.md C4,
# milestone M8) on each object's saved SNN graph: Seurat's ModularityOptimizer.cpp standalone,
# built by validation/cpp/build.sh into data/clustering/bin/, compared with R's optimizer result
# (the production-profile instrumented run), the saved labels and the workflow job's log.
# Reads no .rds; tier T can run it with bash on a login node (a second each).
# From the repository's top directory, after the object's export and production instrumented run:
#   sbatch validation/slurm/optimizer_clustering.sbatch <object> [<object> ...]
# Output: data/clustering/<object>/optimizer/; an existing one is never overwritten. The optimizer
# uses no BLAS and its own generator, so any CPU type will do.
#
#SBATCH --job-name=scrust-cl-optimizer
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh
fi
cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
[ -x data/clustering/bin/ModularityOptimizer ] || { echo "build it first: bash validation/cpp/build.sh" >&2; exit 2; }
for obj in "$@"; do
    out=data/clustering/$obj/optimizer
    if [ -e "$out" ]; then
        echo "$out exists; not overwriting" >&2
        continue
    fi
    env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 LC_ALL=C \
        "$R_ENV/bin/Rscript" validation/R/optimizer_clustering.R --object "$obj" --out "$out.partial"
    mv "$out.partial" "$out"
done
