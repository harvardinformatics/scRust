#!/bin/bash
#SBATCH --job-name=scrust-m13-noexport
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log
# Milestone M13: without --export, scrust cluster writes M12's files. M12's validation run of each
# object in the portable profile (label m12-4231635, validation/slurm/m12_cluster.sbatch: the counts,
# the workflow's preset, the stage dumps at the level of M8's instrumented runs, --keep-stages on
# tier T, 8 threads) is made again with M13's binary and the same options, and every file is
# compared with M12's byte for byte (summary.json by its result section: its run section records
# the run, and M13 adds the files section, the manifest). The new run's manifest must describe its
# directory (scrust validate manifest). From the repository's top directory, after building and
# copying the binary:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m13_noexport.sbatch <label> <object> [<object> ...]
# Output: data/clustering/<object>/rust/<label>-noexport/portable/ (run/, identical.txt,
# manifest.json) and one row per object in data/clustering/compare/m13_noexport_<label>.tsv.
# Exit status 1 if any file differs.

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh
export OPENBLAS_CORETYPE=Haswell

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -ge 2 ] || { echo "usage: m13_noexport.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
old_label=m12-4231635
threads=${SLURM_CPUS_PER_TASK:-8}
status=0
summary=data/clustering/compare/m13_noexport_$label.tsv
mkdir -p "$(dirname "$summary")"
printf 'object\ttier\tprofile\tcompared_with\tfiles\tidentical\tdiffering\tresult_section\tmanifest\n' > "$summary"
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
[ "$core" = Haswell ] || { echo "OpenBLAS reports core '$core', expected Haswell" >&2; exit 1; }

result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }

for o in "$@"; do
    d=data/clustering/$o
    old=$d/rust/$old_label/portable/run
    [ -d "$old" ] || { echo "$o: no $old" >&2; status=1; continue; }
    row=$(awk -F'\t' -v o="$o" '$1 == o' validation/clustering_datasets.tsv)
    tier=$(echo "$row" | cut -f2)
    level=$(echo "$row" | cut -f4)
    role=$(echo "$row" | cut -f3)
    lists=auto
    [ "$role" = small ] && lists=yes
    keep=()
    [ "$tier" = T ] && keep=(--keep-stages)
    out=$d/rust/$label-noexport/portable
    if [ -e "$out" ]; then echo "skip $out: exists"; continue; fi
    work=$out.partial.${SLURM_JOB_ID:-local}
    rm -rf "$work"
    mkdir -p "$work"
    echo "== $o (tier $tier, dumps $level, portable, $threads threads, no --export)"
    if ! "$bin" cluster --counts "$d/input" --preset workflow --out "$work/run" \
            --dump-stages "$level" --purity-lists "$lists" "${keep[@]}" --threads "$threads" 2> "$work/stderr.txt"; then
        tail -5 "$work/stderr.txt" >&2
        echo "$o: scrust cluster failed" >&2
        status=1
        continue
    fi
    n=0
    dif=0
    : > "$work/identical.txt"
    while IFS= read -r f; do
        [ "$f" = summary.json ] && continue
        n=$((n + 1))
        if [ -f "$old/$f" ] && cmp -s "$work/run/$f" "$old/$f"; then
            echo "identical	$f" >> "$work/identical.txt"
        else
            echo "differs	$f" >> "$work/identical.txt"
            dif=$((dif + 1))
        fi
    done < <(cd "$work/run" && find . -type f | sed 's|^\./||' | sort)
    while IFS= read -r f; do
        [ "$f" = summary.json ] && continue
        [ -f "$work/run/$f" ] || { echo "missing	$f" >> "$work/identical.txt"; dif=$((dif + 1)); n=$((n + 1)); }
    done < <(cd "$old" && find . -type f | sed 's|^\./||' | sort)
    rs=identical
    cmp -s <(result_section "$work/run/summary.json") <(result_section "$old/summary.json") || { rs=differs; status=1; }
    man=pass
    "$bin" validate manifest --dir "$work/run" --json "$work/manifest.json" 2> "$work/manifest.log" || { man=FAIL; status=1; }
    printf '%s\t%s\tportable\t%s\t%s\t%s\t%s\t%s\t%s\n' "$o" "$tier" "$old_label" "$n" "$((n - dif))" "$dif" "$rs" "$man" >> "$summary"
    echo "   $n files: $dif differing; result section $rs; manifest $man"
    [ "$dif" = 0 ] || { grep -v '^identical' "$work/identical.txt" | head >&2; status=1; }
    mv "$work" "$out"
done
exit $status
