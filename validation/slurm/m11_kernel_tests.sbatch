#!/bin/bash
# Milestone M11 (docs/PLAN.md A5b; C6's M11 row): `cargo test --workspace` under another OpenBLAS
# kernel than the portable fixtures' (Haswell). The tests that compare a BLAS-dependent result
# with a portable fixture must stop with a message naming both kernels; every other test
# (teacher-forced and BLAS-free, or Rust against Rust) must pass. From the repository's top
# directory, with the tests built (`cargo test --workspace --release --no-run`):
#   sbatch [--partition=test] validation/slurm/m11_kernel_tests.sbatch [auto|<kernel>]
# auto (the default) lets OpenBLAS detect the CPU's kernel: SapphireRapids on a Sapphire Rapids
# node. Output: data/clustering/m11_kernel_tests/<kernel>/cargo_test.txt and results.tsv (each
# test, its outcome, and the kernels its failure named).
#
#SBATCH --job-name=scrust-m11-kern
#SBATCH -A informatics
#SBATCH --partition=test
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -uo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
export OPENBLAS_CORETYPE=${1:-auto}
source ./env.sh
core=$(target/release/scrust fingerprint | sed -n 's/.*"openblas_corename": "\([A-Za-z0-9]*\)".*/\1/p')
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' kernel in use: $core commit=$(git rev-parse --short HEAD)"
out=data/clustering/m11_kernel_tests/$core
mkdir -p "$out"
cargo test --workspace --release --no-fail-fast -- --test-threads=4 > "$out/cargo_test.txt" 2>&1
echo "cargo test exit status $?"
# each test's outcome, and the kernels a kernel check's message names
{
    printf 'test\toutcome\tmessage\n'
    grep -E '^test .* \.\.\. (ok|FAILED)$' "$out/cargo_test.txt" | sed -E 's/^test (.*) \.\.\. (ok|FAILED)$/\1\t\2\t/'
} > "$out/results.tsv"
grep -E "was made under kernel" "$out/cargo_test.txt" | sed -E 's/^.*panicked at [^:]*:[0-9]+:[0-9]+:$//' | sort | uniq -c | head -20
grep -c $'\tok\t' "$out/results.tsv"; grep -c $'\tFAILED\t' "$out/results.tsv"
