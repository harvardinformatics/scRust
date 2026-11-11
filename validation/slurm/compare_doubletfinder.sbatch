#!/bin/bash
# M2 checks for DoubletFinder's reference runs (validation/R/compare_doubletfinder.R; docs/PLAN.md
# A5b and B6). Reads only plain files, so small objects may be checked directly with bash; large
# stage files are better checked as a batch job. From the repository's top directory:
#   sbatch [--mem=..] validation/slurm/compare_doubletfinder.sbatch <object> [<object> ...]
#   bash validation/slurm/compare_doubletfinder.sbatch <object> ...
# For each object, under data/doubletfinder/compare/<object>/:
#   crosscheck_<run>.json         every run against the saved workflow output (the core M2 gate)
#   instr_vs_direct_<profile>-<arch>.*   the instrumented run against the direct run (identical)
#   <mode>_<arch>_vs_<primary>.*  the portable profile on another CPU type against Sapphire Rapids
#                                 (identical files: the portable references are CPU-independent)
#   portable_vs_production.*      Haswell against SapphireRapids kernels on the same CPU type: the
#                                 R-vs-R floor per stage
#
#SBATCH --job-name=scrust-df-compare
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
cmp() { "$R_ENV/bin/Rscript" validation/R/compare_doubletfinder.R "$@"; }
primary=sapphirerapids

for obj in "$@"; do
    base=data/doubletfinder/$obj
    out=data/doubletfinder/compare/$obj
    mkdir -p "$out"
    runs=$(cd "$base/runs" && ls -d */ 2>/dev/null | tr -d / | grep -v '\.partial\.' || true)
    for r in $runs; do
        cmp crosscheck "$base/runs/$r" "$base" "$out/crosscheck_$r.json"
    done
    for r in $runs; do
        case $r in instrumented-*)
            d=direct-${r#instrumented-}
            [ -d "$base/runs/$d" ] && cmp runs "$base/runs/$d" "$base/runs/$r" "$out/instr_vs_direct_${r#instrumented-}"
            ;;
        esac
    done
    for r in $runs; do
        case $r in *-portable-$primary) ;; *-portable-*)
            mode=${r%%-*}
            p=$mode-portable-$primary
            [ -d "$base/runs/$p" ] && cmp runs "$base/runs/$p" "$base/runs/$r" "$out/${mode}_${r##*-}_vs_$primary"
            ;;
        esac
    done
    a=$base/runs/instrumented-portable-$primary
    b=$base/runs/instrumented-production-$primary
    if [ -d "$a" ] && [ -d "$b" ]; then
        cmp runs "$a" "$b" "$out/portable_vs_production"
    elif [ -d "$a" ] && [ -d "$base/runs/direct-production-$primary" ]; then
        cmp runs "$a" "$base/runs/direct-production-$primary" "$out/portable_vs_production"
    fi
done
