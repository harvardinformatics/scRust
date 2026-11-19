#!/bin/bash
# Milestone M15b: R reference runs of the workflow's clustering chain with SCTransform's
# vars.to.regress (validation/R/reference_clustering.R with --vars-to-regress, --percent-mt,
# --mt-genes, --metadata; validation/R/lib/clustering.R cl_covariates), one array task per row of
# validation/m15b_cases.tsv (validation/python/m15b_cases.py writes it and the cases' inputs).
# Approved for M15b: the chain with vars.to.regress, in the production R env, by sbatch.
# Submit from the repository's top directory, the rows as array indices (1 = the first case):
#   sbatch --array=1-108 validation/slurm/m15b_references.sbatch production
#   sbatch --array=109,110 --mem=40G validation/slurm/m15b_references.sbatch production
# Profiles (A5b): production (OPENBLAS_CORETYPE unset: OpenBLAS detects SapphireRapids, which the
# job checks, on the sapphire partition or test) or portable (Haswell).
# Output: data/clustering/m15b/runs/<object>/<case>/direct-<profile>-<cpu arch>/, written to a
# .partial directory first and renamed when complete; never overwritten. /usr/bin/time -v goes to
# time_v.txt in the run directory, the wrapper's fingerprint to fingerprint_wrapper.json.
#
#SBATCH --job-name=scrust-m15b-ref
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
profile=$1
case $profile in
    portable) ;;
    production) export OPENBLAS_CORETYPE=auto ;;
    *) echo "profile must be portable or production" >&2; exit 2 ;;
esac
source ./env.sh
export SCRUST_PROFILE=$profile

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no scrust binary at $bin (for the fingerprint)" >&2; exit 2; }
fp=$("$bin" fingerprint)
core=$(sed -n 's/^  "openblas_corename": "\(.*\)",$/\1/p' <<< "$fp")
arch=$(scontrol show node "${SLURMD_NODENAME:-$(hostname -s)}" 2>/dev/null \
    | grep -o 'AvailableFeatures=[^ ]*' | cut -d= -f2 | tr ',' '\n' \
    | grep -E -x 'sapphirerapids|emeraldrapids|icelake|cascadelake|skylake|broadwell|haswell|genoa|milan|rome|turin' \
    | head -n 1)
arch=${arch:-unknown}
if [ "$profile" = portable ] && [ "$core" != Haswell ]; then
    echo "profile portable needs OpenBLAS core Haswell, got '$core'; stopping" >&2
    exit 3
fi
if [ "$profile" = production ] && { [ "$arch" != sapphirerapids ] || [ "$core" != SapphireRapids ]; }; then
    echo "production needs a sapphirerapids node with core SapphireRapids, got $arch/'$core'; stopping" >&2
    exit 3
fi

row=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n' validation/m15b_cases.tsv)
[ -n "$row" ] || { echo "no case $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
# (cut, not read: read collapses the empty fields' tabs)
field() { cut -f"$1" <<< "$row"; }
obj=$(field 1); tier=$(field 2); input=$(field 3); assay=$(field 4); level=$(field 5); cs=$(field 6)
vars=$(field 7); pm=$(field 8); mt=$(field 9); md=$(field 10)
out=data/clustering/m15b/runs/$obj/$cs/direct-$profile-$arch
if [ -e "$out" ]; then
    echo "$out exists; not overwriting" >&2
    exit 0
fi
partial=$out.partial.${SLURM_JOB_ID:-$$}
rm -rf "$partial"
mkdir -p "$partial"
echo "$fp" > "$partial/fingerprint_wrapper.json"
extra=(--vars-to-regress "$vars" --percent-mt "$pm")
[ -z "$mt" ] || extra+=(--mt-genes "$mt")
[ -z "$md" ] || extra+=(--metadata "$md")
echo "== $obj case $cs ($tier, $assay, level $level): vars.to.regress $vars, percent.mt $pm, mt_genes '${mt}', metadata '${md}' -> $out"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
/usr/bin/time -v -o "$partial/time_v.txt" \
    "$R_ENV/bin/Rscript" validation/R/reference_clustering.R --input "$input" --out "$partial" --level "$level" \
    --umap yes --assay "$assay" "${extra[@]}"
mv "$partial" "$out"
