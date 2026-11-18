#!/bin/bash
# Milestone M15: `cargo test --workspace` in the default (debug) profile and in --release, under
# env.sh's OpenBLAS kernel (Haswell, the portable fixtures' kernel), as every branch needs before
# it is pushed. From the repository's top directory:
#   sbatch validation/slurm/m15_cargo_tests.sbatch
# Output: data/clustering/m15/cargo_tests/<job id>/{debug,release}.txt and results.tsv (each
# test's outcome in each profile), and the warnings count of each build.
#
#SBATCH --job-name=scrust-m15-cargo
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0-03:00:00
#SBATCH --output=%x-%j.sbatch.log

set -uo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
out=data/clustering/m15/cargo_tests/${SLURM_JOB_ID:-$$}
mkdir -p "$out"
core=$(target/release/scrust fingerprint 2>/dev/null | sed -n 's/.*"openblas_corename": "\([A-Za-z0-9]*\)".*/\1/p')
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} kernel=$core commit=$(git rev-parse --short HEAD) $(rustc --version)"
status=0
for profile in debug release; do
    flag=$([ "$profile" = release ] && echo --release || true)
    cargo test --workspace $flag -- --test-threads=4 > "$out/$profile.txt" 2>&1 || status=1
    passed=$(grep -E '^test result:' "$out/$profile.txt" | awk '{p += $4; f += $6; i += $8} END {print p, f, i}')
    warnings=$(grep -c '^warning' "$out/$profile.txt" || true)
    echo "$profile: passed failed ignored = $passed; warnings $warnings; suites $(grep -c '^test result:' "$out/$profile.txt")"
done
{
    printf 'profile\ttest\toutcome\n'
    for profile in debug release; do
        grep -E '^test .* \.\.\. (ok|FAILED|ignored)$' "$out/$profile.txt" | sed -E "s/^test (.*) \.\.\. (ok|FAILED|ignored)$/$profile\t\1\t\2/"
    done
} > "$out/results.tsv"
exit $status
