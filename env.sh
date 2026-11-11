# shellcheck shell=bash
# scRust development environment (docs/PLAN.md A2 and A5b). Source it, from any directory, in
# an interactive shell or at the top of an sbatch script:
#   source /path/to/scRust/env.sh
#
# It:
# - activates the project-local conda env .envs/rust-dev by path (create it first: README.md);
#   on Cannon, conda comes from `module load python`, loaded here if needed;
# - sets CARGO_HOME to .cargo-home/ in the repository, so cargo writes nothing to the home
#   directory (the conda rust package has no rustup);
# - pins OpenBLAS's kernel and threads, for R and Rust alike. OPENBLAS_CORETYPE defaults to
#   Haswell, the "portable" reference profile; set it before sourcing to choose another kernel.
#   OPENBLAS_CORETYPE=auto unsets it, so OpenBLAS detects the CPU: the "production" profile,
#   run on Sapphire Rapids nodes. OPENBLAS_NUM_THREADS is always 1;
# - keeps user-level libraries out of R, Perl and Python: unsets R_LIBS_USER, R_LIBS_SITE, R_LIBS
#   and PERL5LIB, and sets PYTHONNOUSERSITE=1.

if [ -z "${BASH_VERSION:-}" ]; then
    echo "env.sh: this file needs bash" >&2
    return 1 2>/dev/null || exit 1
fi
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "env.sh: source this file instead of running it: source $0" >&2
    exit 1
fi

_scrust_env_setup() {
    local root env
    root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || return 1
    env="$root/.envs/rust-dev"
    if [ ! -d "$env/conda-meta" ]; then
        echo "env.sh: no conda env at $env; create it first (see README.md)" >&2
        return 1
    fi

    # conda's shell function: from the python module, or from the conda executable's hook
    # (in batch jobs the module may be loaded already while its shell functions are not defined)
    if [ "$(type -t conda)" != function ]; then
        if [ -z "${CONDA_EXE:-}" ] && type module >/dev/null 2>&1; then
            module load python || return 1
        fi
        if [ "$(type -t conda)" != function ]; then
            eval "$("${CONDA_EXE:-conda}" shell.bash hook)" || return 1
        fi
    fi
    conda activate "$env" || return 1

    export SCRUST_ROOT="$root"
    export CARGO_HOME="$root/.cargo-home"

    case "${OPENBLAS_CORETYPE:-}" in
        "") export OPENBLAS_CORETYPE=Haswell ;;
        auto) unset OPENBLAS_CORETYPE ;;
        *) export OPENBLAS_CORETYPE ;;
    esac
    export OPENBLAS_NUM_THREADS=1

    export PYTHONNOUSERSITE=1
    unset R_LIBS_USER R_LIBS_SITE R_LIBS PERL5LIB
}

# conda's and Lmod's shell functions do not support `set -u`: relax it while they run.
_scrust_opts=$-
set +u
if _scrust_env_setup; then _scrust_status=0; else _scrust_status=1; fi
case $_scrust_opts in *u*) set -u ;; esac
unset -f _scrust_env_setup
unset _scrust_opts
if [ "$_scrust_status" -ne 0 ]; then
    unset _scrust_status
    return 1
fi
unset _scrust_status
