#!/bin/bash
# Template for scRust's sbatch wrappers (docs/PLAN.md A2, A4 and A5b). Copy it to
# validation/slurm/<job>.sbatch, set the job name and resources, and put the work at the end.
# Submit from the repository's top directory, or export SCRUST_ROOT=<checkout> first:
#   sbatch validation/slurm/<job>.sbatch
# The log is <job name>-<job id>.sbatch.log in the submit directory (gitignored).
#
# Keep in every wrapper:
# - the node guard and --exclude. Some nodes have lost their connection to one holylfs05 storage
#   target (OST 51), and any read of a file stored there hangs unkillably; the guard exits with
#   code 75 instead, and --exclude skips the nodes known to be affected;
# - env.sh: the project env, OpenBLAS pinned to one kernel and one thread for R and Rust alike,
#   and no user-level R, Perl or Python libraries.
#
# Reference profiles (A5b): "portable" is env.sh's default (OPENBLAS_CORETYPE=Haswell) and runs
# on any node. "production" lets OpenBLAS detect the kernel, as the production run did, on
# Sapphire Rapids nodes: set OPENBLAS_CORETYPE=auto below, and keep --partition=sapphire (or use
# --constraint=sapphirerapids on another partition).
#
#SBATCH --job-name=scrust-template
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
# Production profile: uncomment.
# export OPENBLAS_CORETYPE=auto
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

# The work ------------------------------------------------------------------------------------
# Rust, for example:
#   cargo run --release -p scrust -- doubletfinder ...
# An original tool, in the R reference env (envs/r-reference/README.md); env.sh has already
# removed user-level R, Perl and Python libraries:
#   R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
#   "$R_ENV/bin/Rscript" validation/R/<script>.R ...
