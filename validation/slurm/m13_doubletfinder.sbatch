#!/bin/bash
# Milestone M13: Step 1 unchanged. M13 changes code `scrust doubletfinder` runs (SCTransform's
# residual pass, which correct() now shares; the output manifest), so the command is run again
# with M13's binary on the objects M11 re-graded, with the options of the run it is compared
# with, and every file it writes is compared with that run's byte for byte:
#   m6: the ten objects of M6 (tier T and P), both profiles, with the stage dumps at the level of
#       the profile's instrumented R run (M6's options): against M6b's runs (label m6b-13592a6:
#       cells.tsv, sweep_stats.tsv, bcmvn.tsv, singlets/, stages/);
#   m7: M7's objects, production profile, no dumps (M7's options): against M7's runs (label
#       m7-b8c8a8c: cells.tsv, sweep_stats.tsv, bcmvn.tsv, singlets/).
# summary.json is compared by its result section (its run section records the run, and M13 adds
# the files section, the manifest). The new run's manifest must describe its directory (scrust
# validate manifest). From the repository's top directory, after building and copying the binary:
#   SCRUST_BIN=<copy> sbatch [--partition=...] validation/slurm/m13_doubletfinder.sbatch <label> m6|m7 <object> [...]
# Production runs need a Sapphire Rapids node (OpenBLAS detects SapphireRapids there).
# Output: data/doubletfinder/<object>/rust/<label>/<profile>/: run/, identical.txt (one line per
# file: identical or differs), manifest.json, time_v.txt, stderr.txt; and one row per run in
# data/doubletfinder/compare/m13_step1_<label>.tsv. Exit status 1 if any file differs.
#
#SBATCH --job-name=scrust-m13-df
#SBATCH -A informatics
#SBATCH --partition=test
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -ge 3 ] || { echo "usage: m13_doubletfinder.sbatch <label> m6|m7 <object> [<object> ...]" >&2; exit 2; }
label=$1
set_=$2
shift 2
threads=${SLURM_CPUS_PER_TASK:-1}
status=0
summary=data/doubletfinder/compare/m13_step1_$label.tsv
mkdir -p "$(dirname "$summary")"
[ -f "$summary" ] || printf 'object\tprofile\tcompared_with\tfiles\tidentical\tdiffering\tresult_section\tmanifest\twall\tpeak_rss_kib\n' > "$summary"

result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }

compare_runs() {  # new run dir, old run dir, output file; prints differing count
    local new=$1 old=$2 out=$3 n=0 d=0
    : > "$out"
    while IFS= read -r f; do
        [ "$f" = summary.json ] && continue
        n=$((n + 1))
        if [ -f "$old/$f" ] && cmp -s "$new/$f" "$old/$f"; then
            echo "identical	$f" >> "$out"
        else
            echo "differs	$f" >> "$out"
            d=$((d + 1))
        fi
    done < <(cd "$new" && find . -type f | sed 's|^\./||' | sort)
    # files the old run has that the new one lacks
    while IFS= read -r f; do
        [ "$f" = summary.json ] && continue
        [ -f "$new/$f" ] || { echo "missing	$f" >> "$out"; d=$((d + 1)); n=$((n + 1)); }
    done < <(cd "$old" && find . -type f | sed 's|^\./||' | sort)
    echo "$n $d"
}

for obj in "$@"; do
    base=data/doubletfinder/$obj
    case $set_ in
        m6) profiles="portable production"; old_label=m6b-13592a6 ;;
        m7) profiles=production; old_label=m7-b8c8a8c ;;
        *) echo "set: m6 or m7, not $set_" >&2; exit 2 ;;
    esac
    for profile in $profiles; do
        old=$base/rust/$old_label/$profile/run
        [ -d "$old" ] || continue
        case $profile in
            portable) core=Haswell; envset=(env OPENBLAS_CORETYPE=Haswell) ;;
            production) core=SapphireRapids; envset=(env -u OPENBLAS_CORETYPE) ;;
        esac
        out=$base/rust/$label/$profile
        if [ -e "$out" ]; then echo "skip $out: exists"; continue; fi
        tmp=$out.partial.${SLURM_JOB_ID:-local}
        rm -rf "$tmp"
        mkdir -p "$tmp"
        "${envset[@]}" "$bin" fingerprint > "$tmp/fingerprint.json"
        got=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
        if [ "$got" != "$core" ]; then
            echo "$obj $profile: OpenBLAS reports core '$got', expected $core; skipped" >&2
            status=1
            rm -rf "$tmp"
            continue
        fi
        # the old run's options: its dump level (M6), none (M7)
        dumps=()
        level=$(sed -n 's/.*"dump_stages": *"\([a-z]*\)".*/\1/p' "$old/summary.json" | head -n 1)
        [ -n "$level" ] && [ "$level" != none ] && dumps=(--dump-stages "$level")
        echo "== $obj $profile ($core), ${dumps[*]:-no dumps}, $threads threads, against $old_label"
        if ! "${envset[@]}" /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" doubletfinder \
                --counts "$base/input" --clusters "$base/input/seurat_clusters.tsv" \
                --out "$tmp/run" "${dumps[@]}" --threads "$threads" 2> "$tmp/stderr.txt"; then
            cat "$tmp/stderr.txt" >&2
            echo "$obj $profile: scrust doubletfinder failed" >&2
            status=1
            continue
        fi
        read -r n d < <(compare_runs "$tmp/run" "$old" "$tmp/identical.txt")
        idn=$((n - d))
        rs=identical
        cmp -s <(result_section "$tmp/run/summary.json") <(result_section "$old/summary.json") || { rs=differs; d=$((d + 1)); }
        man=pass
        "$bin" validate manifest --dir "$tmp/run" --json "$tmp/manifest.json" 2> "$tmp/manifest.log" || { man=FAIL; status=1; }
        wall=$(sed -n 's/^\s*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$tmp/time_v.txt")
        rss=$(sed -n 's/^\s*Maximum resident set size (kbytes): //p' "$tmp/time_v.txt")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$obj" "$profile" "$old_label" "$n" "$idn" "$d" "$rs" "$man" "$wall" "$rss" >> "$summary"
        echo "   $n files: $d differing; result section $rs; manifest $man; $wall, $rss KiB"
        [ "$d" = 0 ] || { grep -v '^identical' "$tmp/identical.txt" >&2; status=1; }
        mv "$tmp" "$out"
    done
done
exit $status
