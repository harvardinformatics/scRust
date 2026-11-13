#!/bin/bash
# Direct or instrumented reference runs of the workflow's clustering chain
# (validation/R/reference_clustering.R, validation/R/instrumented_clustering.R; docs/PLAN.md A5b and
# C4, milestone M8) on exported inputs (data/clustering/<object>/input, made by
# export_clustering.sbatch).
# Submit from the repository's top directory:
#   sbatch [--partition=..] [--constraint=..] [--mem=..] [--time=..] \
#     validation/slurm/reference_clustering.sbatch <portable|production> <direct|instrumented> <object> [<object> ...]
# Profiles (A5b):
#   portable    OPENBLAS_CORETYPE=Haswell (env.sh's default), any node with AVX2
#   production  OPENBLAS_CORETYPE unset, so OpenBLAS detects the CPU as the workflow's jobs did, on
#               the CPU type of the object's workflow job (validation/clustering_datasets.tsv: all
#               Sapphire Rapids, so --partition=sapphire, this file's default, or test, or
#               --constraint=sapphirerapids elsewhere). The job stops unless the node has that CPU
#               feature and OpenBLAS reports its core (SapphireRapids). Setting
#               OPENBLAS_CORETYPE=SapphireRapids or =auto literally selects Cooperlake instead.
# Variants, through the environment:
#   CL_UMAP=no            leave RunUMAP out (run directory suffix -noumap)
#   CL_ASSAY=v5|v3        force the RNA assay's class instead of the object's own (input/assay.txt;
#                         suffix -v5 or -v3): the emptyDrops objects' v3 path against the v5 path
#   CL_LEVEL=...          dump level instead of the table's (full, standard, light)
#   CL_PURITY_LISTS=yes|no  neighborPurity's neighbour lists in instrumented runs (default: at level
#                         full and for the small tier-P object)
# Output: data/clustering/<object>/runs/<mode>-<profile>-<cpu arch><suffix>/, written to a .partial
# directory first and renamed when complete; an existing run directory is never overwritten.
# /usr/bin/time -v goes to time_v.txt in the run directory, the wrapper's fingerprint to
# fingerprint_wrapper.json.
#
#SBATCH --job-name=scrust-cl-ref
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
profile=$1
mode=$2
shift 2
case $profile in
    portable) ;;
    production) export OPENBLAS_CORETYPE=auto ;;
    *) echo "profile must be portable or production" >&2; exit 2 ;;
esac
case $mode in
    direct) script=validation/R/reference_clustering.R ;;
    instrumented) script=validation/R/instrumented_clustering.R ;;
    *) echo "mode must be direct or instrumented" >&2; exit 2 ;;
esac
umap=${CL_UMAP:-yes}
case $umap in yes|no) ;; *) echo "CL_UMAP must be yes or no" >&2; exit 2 ;; esac
source ./env.sh
export SCRUST_PROFILE=$profile

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
[ -x target/release/scrust ] || { echo "build scrust first: cargo build --release -p scrust" >&2; exit 2; }
fp=$(target/release/scrust fingerprint)
echo "$fp"
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

R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
for obj in "$@"; do
    row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
    [ -n "$row" ] || { echo "unknown object $obj" >&2; exit 2; }
    role=$(cut -f3 <<< "$row")
    level=${CL_LEVEL:-$(cut -f4 <<< "$row")}
    want_cpu=$(cut -f9 <<< "$row")
    if [ "$profile" = production ]; then
        case $want_cpu in
            sapphirerapids) want_core=SapphireRapids ;;
            *) echo "no OpenBLAS core known for CPU type $want_cpu" >&2; exit 2 ;;
        esac
        if [ "$arch" != "$want_cpu" ] || [ "$core" != "$want_core" ]; then
            echo "$obj: production needs a $want_cpu node with core $want_core, got $arch/'$core'; stopping" >&2
            exit 3
        fi
    fi
    input=data/clustering/$obj/input
    [ -r "$input/matrix.mtx.gz" ] || { echo "no export in $input" >&2; exit 2; }
    suffix=""
    [ "$umap" = yes ] || suffix="$suffix-noumap"
    assay=auto
    if [ -n "${CL_ASSAY:-}" ]; then
        assay=$CL_ASSAY
        suffix="$suffix-$CL_ASSAY"
    fi
    [ -z "${CL_LEVEL:-}" ] || suffix="$suffix-$CL_LEVEL"
    out=data/clustering/$obj/runs/$mode-$profile-$arch$suffix
    if [ -e "$out" ]; then
        echo "$out exists; not overwriting" >&2
        continue
    fi
    partial=$out.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$partial"
    mkdir -p "$partial"
    echo "$fp" > "$partial/fingerprint_wrapper.json"
    extra=()
    if [ "$mode" = instrumented ]; then
        lists=${CL_PURITY_LISTS:-auto}
        [ -n "${CL_PURITY_LISTS:-}" ] || [ "$role" != small ] || lists=yes
        extra=(--purity-lists "$lists")
    fi
    echo "== $obj $mode $profile ($arch) level $level umap $umap assay $assay -> $out"
    /usr/bin/time -v -o "$partial/time_v.txt" \
        "$R_ENV/bin/Rscript" "$script" --input "$input" --out "$partial" --level "$level" \
        --umap "$umap" --assay "$assay" "${extra[@]}"
    mv "$partial" "$out"
done
