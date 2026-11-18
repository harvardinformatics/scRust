#!/bin/bash
# Milestone M14: the loaders' own tests, in their hermetic envs (docs/PLAN.md A2, C3g):
#   - installs loaders/R/scrustio into the repository's library .rlib/ (never into the R env) and
#     runs its tests in the production R env (loaders/R/scrustio/tests/scrustio_tests.R);
#   - runs the Python package's tests (pytest) in .envs/loaders-py, with PYTHONNOUSERSITE=1 and
#     python -I; the package is put on the path by its pyproject.toml, not installed.
# Both make their test outputs with the scrust binary (SCRUST_BIN, else target/release/scrust)
# from the committed fixtures. From the repository's top directory:
#   sbatch validation/slurm/m14_tests.sbatch
# Output: data/clustering/m14/tests/<job id>/ (r.log, python.log, install.log).
#
#SBATCH --job-name=scrust-m14-tests
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
PY=${LOADERS_PY:-$PWD/.envs/loaders-py/bin/python}
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
export SCRUST_BIN=$PWD/$bin
echo "binary=$SCRUST_BIN sha256=$(sha256sum "$SCRUST_BIN" | cut -c1-16) $("$SCRUST_BIN" --version)"
out=data/clustering/m14/tests/${SLURM_JOB_ID:-$$}
mkdir -p "$out" .rlib
cache=$PWD/$out/cache
export MPLCONFIGDIR=$cache/mpl NUMBA_CACHE_DIR=$cache/numba XDG_CACHE_HOME=$cache/xdg XDG_CONFIG_HOME=$cache/xdgconfig
status=0

echo "== R: install scrustio into .rlib and run its tests"
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
    "$R_ENV/bin/R" CMD INSTALL --no-docs --no-multiarch --library=.rlib loaders/R/scrustio > "$out/install.log" 2>&1 \
    || { echo "R CMD INSTALL failed" >&2; cat "$out/install.log" >&2; exit 1; }
if env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 TMPDIR="$PWD/$out" \
        "$R_ENV/bin/Rscript" loaders/R/scrustio/tests/scrustio_tests.R > "$out/r.log" 2>&1; then
    tail -n 2 "$out/r.log"
else
    echo "R tests FAILED" >&2
    tail -n 40 "$out/r.log" >&2
    status=1
fi

echo "== Python: pytest"
if (cd loaders/python/scrustio && env -u PYTHONPATH -u PYTHONHOME PYTHONNOUSERSITE=1 \
        "$PY" -I -m pytest -p no:cacheprovider -q --basetemp="$PWD/../../../$out/pytest" -W ignore::DeprecationWarning) \
        > "$out/python.log" 2>&1; then
    tail -n 3 "$out/python.log"
else
    echo "Python tests FAILED" >&2
    tail -n 60 "$out/python.log" >&2
    status=1
fi
rm -rf "$cache" "$out/pytest"
exit $status
