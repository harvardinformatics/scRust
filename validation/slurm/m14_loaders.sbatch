#!/bin/bash
# Milestone M14 (docs/PLAN.md C3g, C4 "Loader tests", C5's K11): each loader's object against the
# saved object, for scrust outputs made by an earlier validation run (label <label>, profile
# <profile>). For each <object>:
#   - a clustering object (validation/clustering_datasets.tsv): its export
#     data/clustering/<object>/rust/<label>/<profile>/run (scrust cluster --export):
#       R: validation/R/loader_compare.R: the Seurat and SingleCellExperiment loaders against the
#          saved .rds, read directly; the files loaded exactly; the .rds round trips;
#       Python: validation/python/loader_compare.py: the AnnData loader against M8's export of the
#          saved object (data/clustering/<object>/saved, input); the files; the .h5ad round trip;
#       the 10x readers: Seurat's Read10X (validation/R/guard_read10x.R) and scanpy's read_10x_mtx
#          (validation/python/guard_read10x.py, every prefix) must fail on every directory of the
#          export but counts/, which they must read as the input;
#   - df:<input> (validation/doubletfinder_datasets.tsv): the scrust doubletfinder output
#     data/doubletfinder/<input>/rust/<label>/<profile>/run against the saved DoubletFinder object
#     (R: the .rds; Python: M6b's and M7's export of it, data/doubletfinder/<input>/saved).
# Reads .rds files. From the repository's top directory, after installing scrustio into .rlib/
# (R CMD INSTALL --library=.rlib loaders/R/scrustio) and creating .envs/loaders-py:
#   sbatch [--mem=..] validation/slurm/m14_loaders.sbatch <label> <profile> <object> [<object> ...]
# Output: data/clustering/<object>/loaders/<label>/<profile>/{r,python,guard}/ (df:
# data/doubletfinder/<input>/loaders/<label>/<profile>/{r,python}/), each step's time_v.txt.
# No step uses BLAS, so any CPU type will do. Memory: 4 GB covers tier T, 96 GB the largest P
# objects (both loaders, the saved object and the round trips at once).
#
#SBATCH --job-name=scrust-m14-loaders
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0-12:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 3 ] || { echo "usage: m14_loaders.sbatch <label> <profile> <object> [<object> ...]" >&2; exit 2; }
label=$1
profile=$2
shift 2
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PIPELINE=/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing
PY=${LOADERS_PY:-.envs/loaders-py/bin/python}
[ -x "$PY" ] || { echo "no $PY: create it from envs/loaders-py.yml" >&2; exit 1; }
[ -d .rlib/scrustio ] || { echo "scrustio is not installed in .rlib" >&2; exit 1; }
# Caches of matplotlib, numba and the XDG locations inside the job's own directory, never $HOME
cache=data/clustering/m14/cache/${SLURM_JOB_ID:-$$}
mkdir -p "$cache"
export MPLCONFIGDIR=$PWD/$cache/mpl NUMBA_CACHE_DIR=$PWD/$cache/numba XDG_CACHE_HOME=$PWD/$cache/xdg XDG_CONFIG_HOME=$PWD/$cache/xdgconfig
renv() { env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 "$@"; }
penv() { env -u PYTHONPATH -u PYTHONHOME PYTHONNOUSERSITE=1 "$@"; }
status=0

for o in "$@"; do
    if [[ $o == df:* ]]; then
        obj=${o#df:}
        row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/doubletfinder_datasets.tsv)
        [ -n "$row" ] || { echo "unknown DoubletFinder object $obj" >&2; status=1; continue; }
        saved=$PIPELINE/$(cut -f6 <<< "$row")
        export_dir=data/doubletfinder/$obj/rust/$label/$profile/run
        out=data/doubletfinder/$obj/loaders/$label/$profile
        kind=doubletfinder
        name=doubletfinder_$obj
    else
        obj=$o
        row=$(awk -F'\t' -v o="$obj" '$1 == o' validation/clustering_datasets.tsv)
        [ -n "$row" ] || { echo "unknown clustering object $obj" >&2; status=1; continue; }
        saved=$PIPELINE/$(cut -f6 <<< "$row")
        export_dir=data/clustering/$obj/rust/$label/$profile/run
        out=data/clustering/$obj/loaders/$label/$profile
        kind=cluster
        name=$obj
    fi
    [ -f "$export_dir/summary.json" ] || { echo "$o: no output in $export_dir" >&2; status=1; continue; }
    mkdir -p "$out"
    tmp=$out/tmp.${SLURM_JOB_ID:-$$}
    echo "== $o ($kind, $profile): R loaders"
    rm -rf "$out/r"
    mkdir -p "$out/r"
    if ! renv /usr/bin/time -v -o "$out/r/time_v.txt" "$R_ENV/bin/Rscript" validation/R/loader_compare.R --object "$name" \
            --export "$export_dir" --rds "$saved" --out "$out/r" --kind "$kind" --tmp "$tmp" > "$out/r/log.txt" 2>&1; then
        echo "$o: R loaders FAILED" >&2
        status=1
    fi
    tail -n 30 "$out/r/log.txt" | grep -E "components|FAIL|Error" || true
    echo "== $o ($kind, $profile): Python loader"
    rm -rf "$out/python"
    mkdir -p "$out/python"
    if [ "$kind" = cluster ]; then
        pargs=(--saved "data/clustering/$obj/saved" --input "data/clustering/$obj/input")
    else
        pargs=(--saved "data/doubletfinder/$obj/saved" --kind doubletfinder)
    fi
    if ! penv /usr/bin/time -v -o "$out/python/time_v.txt" "$PY" -I -W ignore validation/python/loader_compare.py \
            --object "$name" --export "$export_dir" "${pargs[@]}" --out "$out/python" --tmp "$tmp" > "$out/python/log.txt" 2>&1; then
        echo "$o: Python loader FAILED" >&2
        status=1
    fi
    tail -n 30 "$out/python/log.txt" | grep -E "components|Error" || true
    rm -rf "$tmp"
    if [ "$kind" = cluster ]; then
        echo "== $o ($profile): the 10x readers on the export"
        mkdir -p "$out/guard"
        if renv "$R_ENV/bin/Rscript" validation/R/guard_read10x.R "$out/guard/read10x_seurat.tsv" "$export_dir" "data/clustering/$obj/input" \
                > "$out/guard/read10x_seurat.log" 2>&1; then
            tail -n 1 "$out/guard/read10x_seurat.log"
        else
            echo "$o: Read10X NOT as required" >&2
            status=1
        fi
        if penv "$PY" -I validation/python/guard_read10x.py "$out/guard/read10x_scanpy.tsv" "$export_dir" "data/clustering/$obj/input" \
                > "$out/guard/read10x_scanpy.log" 2>&1; then
            tail -n 1 "$out/guard/read10x_scanpy.log"
        else
            echo "$o: read_10x_mtx NOT as required" >&2
            status=1
        fi
    fi
done
rm -rf "$cache"
exit $status
