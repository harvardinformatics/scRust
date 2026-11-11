#!/bin/bash
# Direct or instrumented reference runs of DoubletFinder (validation/R/reference_doubletfinder.R,
# validation/R/instrumented_doubletfinder.R; docs/PLAN.md A4, A5b and B5) on exported inputs
# (data/doubletfinder/<object>/input, made by export_doubletfinder.sbatch).
# Submit from the repository's top directory:
#   sbatch [--partition=..] [--constraint=..] [--mem=..] [--time=..] \
#     validation/slurm/reference_doubletfinder.sbatch <portable|production> <direct|instrumented> <object> [<object> ...]
# Profiles (A5b):
#   portable    OPENBLAS_CORETYPE=Haswell (env.sh's default), any node with AVX2
#   production  OPENBLAS_CORETYPE unset, so OpenBLAS detects the CPU as the production jobs did;
#               Sapphire Rapids nodes only (--partition=sapphire, this file's default). Setting
#               OPENBLAS_CORETYPE=SapphireRapids or =auto literally selects Cooperlake instead.
# The job checks the kernel OpenBLAS reports (scrust fingerprint) and stops if it is not the
# profile's: Haswell, or SapphireRapids.
# Instrumented runs dump the level given in validation/doubletfinder_datasets.tsv, or $DF_LEVEL.
# Output: data/doubletfinder/<object>/runs/<mode>-<profile>-<cpu arch>/, written to a .partial
# directory first and renamed when complete; an existing run directory is never overwritten.
# /usr/bin/time -v goes to time_v.txt in the run directory.
#
#SBATCH --job-name=scrust-df-ref
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
profile=$1
mode=$2
shift 2
case $profile in
    portable) expected_core=Haswell ;;
    production) export OPENBLAS_CORETYPE=auto; expected_core=SapphireRapids ;;
    *) echo "profile must be portable or production" >&2; exit 2 ;;
esac
case $mode in
    direct) script=validation/R/reference_doubletfinder.R ;;
    instrumented) script=validation/R/instrumented_doubletfinder.R ;;
    *) echo "mode must be direct or instrumented" >&2; exit 2 ;;
esac
source ./env.sh
export SCRUST_PROFILE=$profile

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
[ -x target/release/scrust ] || { echo "build scrust first: cargo build --release -p scrust" >&2; exit 2; }
fp=$(target/release/scrust fingerprint)
echo "$fp"
core=$(sed -n 's/^  "openblas_corename": "\(.*\)",$/\1/p' <<< "$fp")
if [ "$core" != "$expected_core" ]; then
    echo "profile $profile needs OpenBLAS core $expected_core, got '$core'; stopping" >&2
    exit 3
fi
arch=$(scontrol show node "${SLURMD_NODENAME:-$(hostname -s)}" 2>/dev/null \
    | grep -o 'AvailableFeatures=[^ ]*' | cut -d= -f2 | tr ',' '\n' \
    | grep -E -x 'sapphirerapids|emeraldrapids|icelake|cascadelake|skylake|broadwell|haswell|genoa|milan|rome|turin' \
    | head -n 1)
arch=${arch:-unknown}

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/doubletfinder_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; exit 2; }
    input=data/doubletfinder/$obj/input
    [ -e "$input/matrix.mtx.gz" ] || { echo "no export in $input" >&2; exit 2; }
    out=data/doubletfinder/$obj/runs/$mode-$profile-$arch
    if [ -e "$out" ]; then
        echo "$out exists; not overwriting" >&2
        exit 1
    fi
    partial=$out.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"   # this job's own incomplete output, if it was requeued
    mkdir -p "$partial"
    echo "$fp" > "$partial/fingerprint_wrapper.json"
    extra=()
    if [ "$mode" = instrumented ]; then
        extra=(--level "${DF_LEVEL:-$(cut -f4 <<< "$row")}")
    fi
    echo "== $obj $mode $profile ($arch) -> $out"
    /usr/bin/time -v -o "$partial/time_v.txt" \
        "$R_ENV/bin/Rscript" "$script" --input "$input" --out "$partial" "${extra[@]}"
    mv "$partial" "$out"
done
