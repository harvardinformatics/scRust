#!/bin/bash
# Milestone M11 (docs/PLAN.md A5b mechanism 3): the instrumented production-profile R runs of
# DoubletFinder's two large objects (soupx_emptydrops_sscovelli, cellbender_fromraw_ghirsutum2;
# made with reference_doubletfinder.sbatch on Sapphire Rapids, OpenBLAS kernel detected), checked
# as references before the R-vs-R floor is measured on them (m11_floor.sbatch):
#   1. against the object's direct production run, as M2 checked its instrumented runs
#      (compare_doubletfinder.sbatch, instr_vs_direct): cells.tsv, sweep_stats.tsv and bcmvn.tsv
#      byte-identical, the result section of summary.json identical (it holds a digest of the
#      final .Random.seed), and the run's own hook checks clean (stages.json: problems empty);
#   2. M6's production Rust run (label m6-bc1b48e) against it with `scrust validate doubletfinder`,
#      stage by stage: exact mode (both runs record OpenBLAS SapphireRapids), every check passing.
# No R is run. From the repository's top directory, after `cargo build --release -p scrust`:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m11_df_references.sbatch <label> <object> [<object> ...]
# Outputs: data/doubletfinder/compare/<object>/instr_vs_direct_production-sapphirerapids.tsv, and
# data/doubletfinder/<object>/rust/m6-bc1b48e/production/<label>/
# vs_instrumented-production-sapphirerapids.{json,md,log}. Exit status 1 if any check fails.
#
#SBATCH --job-name=scrust-m11-dfref
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-01:00:00
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
[ $# -ge 2 ] || { echo "usage: m11_df_references.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
threads=${SLURM_CPUS_PER_TASK:-1}
status=0
result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }

for o in "$@"; do
    d=data/doubletfinder/$o
    ref=$d/runs/instrumented-production-sapphirerapids
    direct=$d/runs/direct-production-sapphirerapids
    [ -d "$ref" ] && [ -d "$direct" ] || { echo "$o: no instrumented or direct production run" >&2; status=1; continue; }

    # 1. the instrumented run against the direct run
    out=data/doubletfinder/compare/$o
    mkdir -p "$out"
    t=$out/instr_vs_direct_production-sapphirerapids.tsv
    printf 'check\tidentical\n' > "$t"
    for f in cells.tsv sweep_stats.tsv bcmvn.tsv; do
        if cmp -s "$direct/$f" "$ref/$f"; then printf '%s\ttrue\n' "$f"; else printf '%s\tfalse\n' "$f"; fi >> "$t"
    done
    if [ -n "$(result_section "$ref/summary.json")" ] && cmp -s <(result_section "$direct/summary.json") <(result_section "$ref/summary.json"); then
        printf 'summary.json result section\ttrue\n' >> "$t"
    else
        printf 'summary.json result section\tfalse\n' >> "$t"
    fi
    if grep -q '"problems": \[\]' "$ref/stages/stages.json"; then
        printf 'stages.json problems empty\ttrue\n' >> "$t"
    else
        printf 'stages.json problems empty\tfalse\n' >> "$t"
    fi
    for k in direct ref; do
        printf 'OpenBLAS kernel (%s run)\t%s\n' "$k" "$(grep -m1 '"openblas_corename"' "${!k}/summary.json" | sed 's/.*: "\(.*\)",*/\1/')" >> "$t"
    done
    echo "== $o: instrumented production run against the direct production run"
    cat "$t"
    grep -q $'\tfalse$' "$t" && { echo "$o: the instrumented run differs from the direct run" >&2; status=1; }

    # 2. M6's production Rust run against it, stage by stage
    m6=$d/rust/m6-bc1b48e/production
    [ -d "$m6/run" ] || { echo "$o: no M6 production run" >&2; status=1; continue; }
    v=$m6/$label
    mkdir -p "$v"
    r=$(basename "$ref")
    if "$bin" validate doubletfinder --reference "$ref" --reference-stages "$ref/stages" --candidate "$m6/run" \
        --json "$v/vs_$r.json" --markdown "$v/vs_$r.md" --threads "$threads" 2> "$v/vs_$r.log"; then :; else status=1; echo "$o vs $r: FAIL" >&2; fi
    mode=$(awk '/"grading": \{/{g=1} g && /"mode": /{gsub(/[",]/,""); print $2; exit}' "$v/vs_$r.json")
    echo "== $o: M6's production Rust run against it: pass $(grep -m1 '"pass": ' "$v/vs_$r.json" | sed 's/.*"pass": \([a-z]*\).*/\1/'), mode $mode"
    grep -o '"status": "[^"]*"' "$v/vs_$r.json" | sort | uniq -c
    [ "$mode" = exact ] || { echo "$o: not graded in exact mode" >&2; status=1; }
    sed -n '/^| stage | checks identical/,/^$/p' "$v/vs_$r.md"
done
exit $status
