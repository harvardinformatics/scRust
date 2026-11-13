#!/bin/bash
# M8's comparisons (validation/R/compare_clustering.R; docs/PLAN.md C4): for each object, every
# reference run in data/clustering/<object>/runs/ against the saved object (its export's saved/),
# and the pairs of runs that answer M8's questions:
#   instr_vs_direct_<profile>-<arch>    the hooks change nothing (instrumented against direct)
#   portable_vs_production              the R-vs-R floor across OpenBLAS kernels on one CPU type
#                                       (instrumented runs, Haswell against SapphireRapids)
#   cascadelake_vs_sapphirerapids_<mode>  the floor across CPU types with one kernel (portable)
#   noumap_vs_umap                      RunUMAP changes nothing downstream (direct production runs)
#   v5_vs_v3                            the emptyDrops objects' v3 assay path against the v5 path
#                                       (instrumented production runs)
# into data/clustering/compare/<object>/<name>.json. Reads only plain files: small objects can run it
# with bash on a login node; the large ones through sbatch (32 GB).
# From the repository's top directory:
#   sbatch validation/slurm/compare_clustering.sbatch <object> [<object> ...]
#
#SBATCH --job-name=scrust-cl-compare
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh
fi
cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
rs() { env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" "$@"; }

for obj in "$@"; do
    base=data/clustering/$obj
    out=data/clustering/compare/$obj
    mkdir -p "$out"
    [ -d "$base/saved" ] || { echo "no export of $obj" >&2; exit 2; }
    pair() {  # name a b
        [ -d "$2" ] && [ -d "$3" ] || return 0
        rs validation/R/compare_clustering.R "$obj" "$2" "$3" "$out/$1.json"
    }
    for run in "$base"/runs/*; do
        [ -d "$run" ] || continue
        case $run in *.partial.*) continue ;; esac
        pair "$(basename "$run")_vs_saved" "$base/saved" "$run"
    done
    for pa in portable-sapphirerapids production-sapphirerapids portable-cascadelake; do
        pair "instr_vs_direct_$pa" "$base/runs/direct-$pa" "$base/runs/instrumented-$pa"
    done
    pair portable_vs_production "$base/runs/instrumented-production-sapphirerapids" "$base/runs/instrumented-portable-sapphirerapids"
    pair portable_vs_production_direct "$base/runs/direct-production-sapphirerapids" "$base/runs/direct-portable-sapphirerapids"
    for mode in direct instrumented; do
        pair "cascadelake_vs_sapphirerapids_$mode" "$base/runs/$mode-portable-sapphirerapids" "$base/runs/$mode-portable-cascadelake"
    done
    pair noumap_vs_umap "$base/runs/direct-production-sapphirerapids" "$base/runs/direct-production-sapphirerapids-noumap"
    pair v5_vs_v3 "$base/runs/instrumented-production-sapphirerapids" "$base/runs/instrumented-production-sapphirerapids-v5"
done
