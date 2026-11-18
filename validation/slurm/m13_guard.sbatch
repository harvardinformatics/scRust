#!/bin/bash
# Milestone M13 (docs/PLAN.md C3f, C5's K11, C6's M13 row): the export's guard against being read
# as raw counts, tried with the readers people use. For each object, on the export of its M13
# validation run (data/clustering/<object>/rust/<label>/<profile>/run, `scrust cluster --export`)
# and on a `scrust sctransform --corrected-counts` output made here:
#   - Seurat 5.5.1's Read10X() in the production R env (validation/R/guard_read10x.R): it must
#     fail on every directory but the export's counts/, which it must read as the input;
#   - scanpy's read_10x_mtx() in the hermetic env envs/loaders-py.yml (.envs/loaders-py; with
#     PYTHONNOUSERSITE=1 and python -I), without prefix= and with every prefix that could name a
#     file there (validation/python/guard_read10x.py): the same.
# R reads plain files only (no .rds). From the repository's top directory:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m13_guard.sbatch <label> <profile> <object> [<object> ...]
# Output, data/clustering/<object>/rust/<label>/guard/: sctransform_corrected/ (the second output),
# read10x_seurat.tsv, read10x_scanpy.tsv, and their logs. Exit status 1 if either reader behaves
# otherwise anywhere.
#
#SBATCH --job-name=scrust-m13-guard
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 3 ] || { echo "usage: m13_guard.sbatch <label> <profile> <object> [<object> ...]" >&2; exit 2; }
label=$1
profile=$2
shift 2
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PY=.envs/loaders-py/bin/python
[ -x "$PY" ] || { echo "no $PY: create it from envs/loaders-py.yml" >&2; exit 1; }
threads=${SLURM_CPUS_PER_TASK:-1}
status=0

for o in "$@"; do
    d=data/clustering/$o
    export_dir=$d/rust/$label/$profile/run
    [ -f "$export_dir/README.txt" ] || { echo "$o: no export in $export_dir" >&2; status=1; continue; }
    out=$d/rust/$label/guard
    mkdir -p "$out"
    rm -rf "$out/sctransform_corrected"
    "$bin" sctransform --counts "$d/input" --seed 12345 --corrected-counts --out "$out/sctransform_corrected" --threads "$threads" 2> "$out/sctransform.log"
    echo "== $o: Seurat::Read10X"
    if env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
            "$R_ENV/bin/Rscript" validation/R/guard_read10x.R "$out/read10x_seurat.tsv" "$export_dir" "$d/input" "$out/sctransform_corrected" \
            > "$out/read10x_seurat.log" 2>&1; then
        tail -1 "$out/read10x_seurat.log"
    else
        echo "$o: Read10X NOT as required" >&2
        cat "$out/read10x_seurat.log" >&2
        status=1
    fi
    echo "== $o: scanpy.read_10x_mtx"
    if env -u PYTHONPATH -u PYTHONHOME PYTHONNOUSERSITE=1 "$PY" -I validation/python/guard_read10x.py \
            "$out/read10x_scanpy.tsv" "$export_dir" "$d/input" "$out/sctransform_corrected" > "$out/read10x_scanpy.log" 2>&1; then
        tail -1 "$out/read10x_scanpy.log"
    else
        echo "$o: read_10x_mtx NOT as required" >&2
        tail -30 "$out/read10x_scanpy.log" >&2
        status=1
    fi
done
exit $status
