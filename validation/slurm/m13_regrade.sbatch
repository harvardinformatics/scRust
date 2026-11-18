#!/bin/bash
# Milestone M13: earlier reports re-grade unchanged. M12's own Rust runs (label m12-4231635: every
# tier-T and tier-P object, both profiles) compared again with M13's `scrust validate cluster`,
# against the same R runs with the same options as m12_cluster.sbatch (the profile's R runs with
# the instrumented optimizer, the saved object and, on tier P, the panel; the other profile's
# instrumented run, the mismatched kernels). No Rust or R run is made. M12's runs made no
# corrected counts, so the new K10 checks are not applicable to them, except K10 teacher-forced,
# which needs only R's run and the counts. Then validation/slurm/m13_unchanged.sh compares each
# new report with M12's, check by check. From the repository's top directory:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m13_regrade.sbatch <label> <object> [<object> ...]
# Output: data/clustering/<object>/rust/m12-4231635/<profile>/<label>/vs_*.json, .md, .log.
#
#SBATCH --job-name=scrust-m13-regrade
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto}"
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -ge 2 ] || { echo "usage: m13_regrade.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
old=m12-4231635
status=0
for o in "$@"; do
    d=data/clustering/$o
    tier=$(awk -F'\t' -v o="$o" '$1 == o {print $2}' validation/clustering_datasets.tsv)
    for profile in portable production; do
        m=$d/rust/$old/$profile
        [ -d "$m/run" ] || continue
        out=$m/$label
        mkdir -p "$out"
        other=$([ "$profile" = portable ] && echo production || echo portable)
        panel=()
        [ "$tier" = P ] && panel=(--panel validation/clustering_panel.tsv --object "$o")
        for f in "$m"/vs_*.json; do
            name=$(basename "$f" .json)
            case $name in
                vs_other_*) r=${name#vs_other_}; opts=(--optimizer "$d/optimizer") ;;
                *) r=${name#vs_}; opts=(--optimizer "$d/optimizer" --saved "$d/saved" "${panel[@]}") ;;
            esac
            if "$bin" validate cluster --reference "$d/runs/$r" --candidate "$m/run" "${opts[@]}" \
                    --json "$out/$name.json" --markdown "$out/$name.md" 2> "$out/$name.log"; then :; else
                case $name in vs_other_*) ;; *) status=1; echo "$o $profile $name: FAIL" >&2 ;; esac
            fi
        done
        echo "$o $profile: $(ls "$out"/vs_*.json | wc -l) reports"
    done
done
exit $status
