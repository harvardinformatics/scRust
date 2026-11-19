#!/bin/bash
# Milestone M15b: the code path of SCTransform's second regression and PercentageFeatureSet() in
# the pinned env, and the reference vectors for the Rust port (validation/R/regress_vectors.R):
# tests/fixtures/rcompat/linpack_qr.tsv, linpack_resid.tsv and tests/fixtures/clustering/percent_mt/.
# Portable profile (OPENBLAS_CORETYPE=Haswell, env.sh's default): dqrdc2 and dqrsl call BLAS.
# Submit from the repository's top directory:
#   sbatch validation/slurm/m15b_vectors.sbatch
#
#SBATCH --job-name=scrust-m15b-vec
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-00:30:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
[ "$OPENBLAS_CORETYPE" = Haswell ] || { echo "the vectors are made under the Haswell kernel" >&2; exit 2; }

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
mkdir -p data/clustering/logs/m15b
"$R_ENV/bin/Rscript" validation/R/regress_vectors.R . 2>&1 | tee data/clustering/logs/m15b/regress_vectors.txt
