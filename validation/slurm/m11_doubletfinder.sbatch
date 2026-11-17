#!/bin/bash
# Milestone M11 (docs/PLAN.md A5b mechanism 3): `scrust validate doubletfinder` with kernel-gated
# grading, on Step 1's own Rust runs. No Rust or R run is made: the validator reads M6's runs
# (label m6-bc1b48e: every tier-T and tier-P object, both profiles, with stage dumps) and M7's
# (label m7-b8c8a8c: the full panel, production profile) and the R references. From the
# repository's top directory, after `cargo build --release -p scrust` and copying the binary:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m11_doubletfinder.sbatch <label> [<M7 object> ...]
#
# 1. Same kernels: every comparison M6 made (each R run of the run's profile, and the saved
#    output), again with this binary, into data/doubletfinder/<object>/rust/m6-bc1b48e/<profile>/
#    <label>/vs_<R run>.json; and M7's comparison with the saved output for the M7 objects given,
#    into .../rust/m7-b8c8a8c/production/<label>/vs_saved.json. Each new report must give the
#    same verdict as the old one, check by check (stage, name, status): verdicts.tsv.
# 2. Mismatched kernels: each M6 Rust run against the other profile's instrumented R run on
#    Sapphire Rapids (portable Rust, Haswell, against R's SapphireRapids run; production Rust
#    against R's Haswell run): tolerance mode, every BLAS-dependent difference within the floor
#    (vs_other_<R run>.json).
# Exit status 1 if a verdict changed, or a same-kernel comparison failed.
#
#SBATCH --job-name=scrust-m11-df
#SBATCH -A informatics
#SBATCH --partition=shared
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
[ $# -ge 1 ] || { echo "usage: m11_doubletfinder.sbatch <label> [<M7 object> ...]" >&2; exit 2; }
label=$1
shift
threads=${SLURM_CPUS_PER_TASK:-1}
status=0
summary=data/doubletfinder/compare/m11_verdicts_$label.tsv
mkdir -p "$(dirname "$summary")"
printf 'object\trun\tprofile\tcomparison\told_pass\tnew_pass\tmode\tchecks\tsame_status\tchanged_status\twithin_floor\n' > "$summary"

# the (stage, name, status) of every check of a report, one per line
triples() { awk '/"stage": /{s=$0} /"name": /{n=$0} /"status": /{print s "|" n "|" $0}' "$1"; }
jfield() { grep -m1 "\"$2\": " "$1" | sed 's/.*": "\{0,1\}\([^",]*\)"\{0,1\},\{0,1\}$/\1/'; }

# <object> <run label> <profile> <comparison name> <old json> <new json>
verdict() {
    local o=$1 run=$2 profile=$3 name=$4 old=$5 new=$6
    local same changed n within
    n=$(triples "$new" | wc -l)
    within=$(grep -c '"status": "within floor"' "$new" || true)
    if [ -f "$old" ]; then
        same=$(comm -12 <(triples "$old" | grep -v '"OpenBLAS kernel"' | sort) <(triples "$new" | grep -v '"OpenBLAS kernel"' | sort) | wc -l)
        changed=$(comm -23 <(triples "$old" | grep -v '"OpenBLAS kernel"' | sort) <(triples "$new" | grep -v '"OpenBLAS kernel"' | sort) | wc -l)
        old_pass=$(jfield "$old" pass)
    else
        same=-; changed=-; old_pass=-
    fi
    local mode
    mode=$(awk '/"grading": \{/{g=1} g && /"mode": /{gsub(/[",]/,""); print $2; exit}' "$new")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$o" "$run" "$profile" "$name" "$old_pass" "$(jfield "$new" pass)" "${mode:--}" "$n" "$same" "$changed" "$within" >> "$summary"
    if [ "$changed" != - ] && [ "$changed" != 0 ]; then
        echo "  $o $profile $name: $changed checks changed status" >&2
        status=1
    fi
}

for o in $(ls -d data/doubletfinder/*/rust/m6-bc1b48e 2>/dev/null | cut -d/ -f3); do
    d=data/doubletfinder/$o
    for profile in portable production; do
        m6=$d/rust/m6-bc1b48e/$profile
        [ -d "$m6/run" ] || continue
        out=$m6/$label
        mkdir -p "$out"
        for ref in "$d"/runs/*-"$profile"-*; do
            r=$(basename "$ref")
            rs=()
            [ -f "$ref/stages/index.tsv" ] && rs=(--reference-stages "$ref/stages")
            if "$bin" validate doubletfinder --reference "$ref" "${rs[@]}" --candidate "$m6/run" --json "$out/vs_$r.json" --markdown "$out/vs_$r.md" --threads "$threads" 2> "$out/vs_$r.log"; then :; else status=1; echo "$o $profile vs $r: FAIL" >&2; fi
            verdict "$o" m6-bc1b48e "$profile" "vs_$r" "$m6/vs_$r.json" "$out/vs_$r.json"
        done
        if "$bin" validate doubletfinder --saved "$d/saved" --candidate "$m6/run" --json "$out/vs_saved.json" --markdown "$out/vs_saved.md" 2> "$out/vs_saved.log"; then :; else status=1; echo "$o $profile vs saved: FAIL" >&2; fi
        verdict "$o" m6-bc1b48e "$profile" vs_saved "$m6/vs_saved.json" "$out/vs_saved.json"
        other=$([ "$profile" = portable ] && echo production || echo portable)
        ref=$d/runs/instrumented-$other-sapphirerapids
        if [ -d "$ref" ]; then
            r=$(basename "$ref")
            "$bin" validate doubletfinder --reference "$ref" --candidate "$m6/run" --json "$out/vs_other_$r.json" --markdown "$out/vs_other_$r.md" --threads "$threads" 2> "$out/vs_other_$r.log" || true
            verdict "$o" m6-bc1b48e "$profile" "vs_other_$r" - "$out/vs_other_$r.json"
        fi
        echo "$o $profile: $(grep -c . "$summary") rows so far"
    done
done

for o in "$@"; do
    d=data/doubletfinder/$o
    m7=$d/rust/m7-b8c8a8c/production
    [ -d "$m7/run" ] || { echo "$o: no M7 run" >&2; status=1; continue; }
    out=$m7/$label
    mkdir -p "$out"
    if "$bin" validate doubletfinder --saved "$d/saved" --candidate "$m7/run" --json "$out/vs_saved.json" --markdown "$out/vs_saved.md" 2> "$out/vs_saved.log"; then :; else status=1; echo "$o vs saved: FAIL" >&2; fi
    verdict "$o" m7-b8c8a8c production vs_saved "$m7/vs_saved.json" "$out/vs_saved.json"
done
column -t -s $'\t' "$summary" | head -80
exit $status
