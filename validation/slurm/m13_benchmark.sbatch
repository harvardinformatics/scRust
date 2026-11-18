#!/bin/bash
# Milestone M13 (docs/PLAN.md C5, Performance and Threads): M12's benchmark (m12_benchmark.sbatch)
# with the export: `scrust cluster --counts --preset workflow --export`, the whole chain with both
# scores and the export (counts/, sct/ with the corrected counts and their log1p, the README and
# the manifest), on the seven tier-P objects, at one thread and at many, in the production profile
# on a Sapphire Rapids node (Intel Xeon Platinum 8480CL, the CPU of M8's R runs), against R's
# chain (M8's direct production runs, with and without RunUMAP, which write the SCT assay too).
# The same runs without --export, in the same job, give the export's share. From the repository's
# top directory, after `cargo build --release -p scrust` and copying the binary to a file of its own:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m13_benchmark.sbatch <label> <object> [<object> ...]
# Outputs, data/clustering/<object>/rust/<label>/bench/:
#   t<threads>/           the run with --export (no stage dumps; summary.json's run section has
#                         each stage's time and peak RSS), time_v_t<threads>.txt and stderr beside it
#   n<threads>/           the same without --export (only summary.json and time_v are kept)
#   threads.txt           every file of t1 and t<many> byte-identical, the export included
#                         (summary.json: its result section and its manifest's entries), and the
#                         same files as the validation run's (VALIDATION_LABEL, when given: its run/
#                         minus stages/)
# The job reserves 32 CPUs and 96 GiB of the node; --exclusive is not used, so other jobs may share it.
#
#SBATCH --job-name=scrust-m13-bench
#SBATCH -A informatics
#SBATCH --partition=test
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 2 ] || { echo "usage: m13_benchmark.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
export OPENBLAS_CORETYPE=auto
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": "\([A-Za-z0-9]*\)".*/\1/p')
[ "$core" = SapphireRapids ] || { echo "the production profile needs SapphireRapids, this node gives $core" >&2; exit 1; }
many=${SLURM_CPUS_PER_TASK:-32}
status=0

result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }
entries_section() { python3 -I -c 'import json, sys; print(json.dumps(json.load(open(sys.argv[1]))["files"]["entries"], indent=1, sort_keys=True))' "$1"; }

# same files (summary.json: the result section, and with entries=1 the manifest's entries), <a>
# <b> <what> [entries]; stages/ left out (README.txt too against a run with stage dumps: it lists them)
same_files() {
    local a=$1 b=$2 differ=0 n=0
    while IFS= read -r f; do
        n=$((n + 1))
        [ "$f" = summary.json ] && continue
        [ "${4:-0}" = 0 ] && [ "$f" = README.txt ] && continue
        cmp -s "$a/$f" "$b/$f" || { echo "    $f differs ($3)" >&2; differ=$((differ + 1)); }
    done < <(cd "$a" && find . -type f -not -path './stages/*' | sed 's|^\./||' | sort)
    cmp -s <(result_section "$a/summary.json") <(result_section "$b/summary.json") || { echo "    result section differs ($3)" >&2; differ=$((differ + 1)); }
    if [ "${4:-0}" = 1 ]; then
        cmp -s <(entries_section "$a/summary.json") <(entries_section "$b/summary.json") || { echo "    the manifest's entries differ ($3)" >&2; differ=$((differ + 1)); }
    fi
    echo "$3: $n files, $differ differing"
    [ "$differ" = 0 ]
}

for o in "$@"; do
    d=data/clustering/$o
    out=$d/rust/$label/bench
    if [ -e "$out" ]; then echo "$o: $out exists, skipped"; continue; fi
    work=$out.partial.${SLURM_JOB_ID:-local}
    rm -rf "$work"
    mkdir -p "$work"
    for t in 1 "$many"; do
        echo "== $o, $t threads, --export"
        /usr/bin/time -v -o "$work/time_v_t$t.txt" "$bin" cluster --counts "$d/input" --preset workflow --export --out "$work/t$t" --threads "$t" 2> "$work/stderr_t$t.txt"
        tail -2 "$work/stderr_t$t.txt"
        grep -E "Elapsed \(wall|Maximum resident" "$work/time_v_t$t.txt"
        echo "== $o, $t threads, without --export"
        /usr/bin/time -v -o "$work/time_v_n$t.txt" "$bin" cluster --counts "$d/input" --preset workflow --out "$work/n$t" --threads "$t" 2> "$work/stderr_n$t.txt"
        grep -E "Elapsed \(wall|Maximum resident" "$work/time_v_n$t.txt"
    done
    same_files "$work/t1" "$work/t$many" "1 and $many threads" 1 > "$work/threads.txt" 2>&1 || status=1
    # keep only the record of the runs without --export (their metadata.tsv lacks nCount_SCT and
    # nFeature_SCT, and their result the corrected counts' summary: they are M12's files)
    find "$work/n1" "$work/n$many" -type f ! -name summary.json -delete
    if [ -n "${VALIDATION_LABEL:-}" ] && [ -d "$d/rust/$VALIDATION_LABEL/production/run" ]; then
        same_files "$work/t$many" "$d/rust/$VALIDATION_LABEL/production/run" "the validation run ($VALIDATION_LABEL)" >> "$work/threads.txt" 2>&1 || status=1
    fi
    cat "$work/threads.txt"
    mv "$work" "$out"
done
exit $status
