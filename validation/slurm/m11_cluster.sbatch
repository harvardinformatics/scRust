#!/bin/bash
# Milestone M11 (docs/PLAN.md C5: K1-K8, Threads; C6): the clustering chain from the counts,
# `scrust cluster --counts --preset workflow`, against every R run of M8 in the same OpenBLAS
# profile, with kernel-gated grading (A5b mechanism 3). From the repository's top directory, after
# `cargo build --release -p scrust` and copying the binary to a file of its own:
#   SCRUST_BIN=<copy> sbatch [--partition=...] [--cpus-per-task=N] [--mem=...] \
#       validation/slurm/m11_cluster.sbatch <label> portable|production <object> [<object> ...]
#
# Profiles (A5b): portable runs with OPENBLAS_CORETYPE=Haswell on any node; production lets
# OpenBLAS detect the kernel and must run on a Sapphire Rapids node (partition sapphire or test,
# or --constraint=sapphirerapids), where it detects SapphireRapids, the kernel of every workflow
# job (M8). The script checks the kernel the binary reports before it runs anything.
#
# For each object (data/clustering/<object>/), one Rust run, with the stage dumps at the level of
# M8's instrumented runs (validation/clustering_datasets.tsv: full on tier T, standard or light
# on tier P; --keep-stages on tier T), then `scrust validate cluster`:
#   - against each R run of the profile (runs/*-<profile>-*: instrumented, direct, without RunUMAP,
#     the v5 path; Sapphire Rapids and, portable, Cascade Lake), with the instrumented optimizer
#     (optimizer/), the saved object (saved/; its silhouette widths graded under its own kernel,
#     SapphireRapids) and, for tier P, the workflow job's logged numbers
#     (validation/clustering_panel.tsv): every one must pass;
#   - against the other profile's instrumented run on Sapphire Rapids: a deliberately mismatched
#     pair of kernels, graded in tolerance mode against the R-vs-R floor. Reported (vs_other_*),
#     not counted in the exit status;
#   - THREADS_CHECK=1 (the default on tier T): a second run at one thread, every file (stage dumps
#     included) byte-identical and summary.json's result section identical (threads.txt).
# Output, data/clustering/<object>/rust/<label>/<profile>/:
#   run/                  the Rust run (scrust cluster's outputs, stages/)
#   vs_<R run>.json, .md  each comparison; vs_other_<R run>.* the mismatched pairs
#   fingerprint.json, time_v.txt, stderr.txt, threads.txt
# An existing output directory is never overwritten: the object is skipped. Work in progress goes
# to <dir>.partial.<job id>, renamed when complete. Exit status 1 if a same-profile comparison
# fails or a thread check differs.
#
#SBATCH --job-name=scrust-m11-cl
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-08:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 3 ] || { echo "usage: m11_cluster.sbatch <label> portable|production <object> [<object> ...]" >&2; exit 2; }
label=$1
profile=$2
shift 2
case $profile in
    portable) export OPENBLAS_CORETYPE=Haswell; want=Haswell; other=production ;;
    production) export OPENBLAS_CORETYPE=auto; want=SapphireRapids; other=portable ;;
    *) echo "profile: portable or production, not $profile" >&2; exit 2 ;;
esac
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": "\([A-Za-z0-9]*\)".*/\1/p')
[ "$core" = "$want" ] || { echo "profile $profile needs OpenBLAS kernel $want, this node gives $core" >&2; exit 1; }
threads=${SLURM_CPUS_PER_TASK:-1}
job=${SLURM_JOB_ID:-local}
status=0

# summary.json's result section (to_file_text's layout: two-space indentation, "run" follows)
result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }

for o in "$@"; do
    d=data/clustering/$o
    [ -d "$d/input" ] || { echo "$o: no $d/input" >&2; status=1; continue; }
    row=$(awk -F'\t' -v o="$o" '$1 == o' validation/clustering_datasets.tsv)
    tier=$(echo "$row" | cut -f2)
    level=$(echo "$row" | cut -f4)
    out=$d/rust/$label/$profile
    if [ -e "$out" ]; then echo "$o: $out exists, skipped"; continue; fi
    work=$out.partial.$job
    rm -rf "$work"
    mkdir -p "$work"
    "$bin" fingerprint > "$work/fingerprint.json"
    keep=()
    [ "$tier" = T ] && keep=(--keep-stages)
    echo "== $o (tier $tier, dumps $level, $profile, $threads threads)"
    if ! /usr/bin/time -v -o "$work/time_v.txt" "$bin" cluster --counts "$d/input" --preset workflow --out "$work/run" \
            --dump-stages "$level" "${keep[@]}" --threads "$threads" 2> "$work/stderr.txt"; then
        echo "$o: scrust cluster failed" >&2
        tail -5 "$work/stderr.txt" >&2
        status=1
        continue
    fi
    tail -2 "$work/stderr.txt"
    panel=()
    [ "$tier" = P ] && panel=(--panel validation/clustering_panel.tsv --object "$o")
    for ref in "$d"/runs/*-"$profile"-*; do
        r=$(basename "$ref")
        if "$bin" validate cluster --reference "$ref" --candidate "$work/run" --optimizer "$d/optimizer" --saved "$d/saved" "${panel[@]}" \
                --json "$work/vs_$r.json" --markdown "$work/vs_$r.md" 2> "$work/vs_$r.log"; then
            echo "  vs $r: $(tail -1 "$work/vs_$r.log")"
        else
            echo "  vs $r: FAIL" >&2
            grep -E " differ | missing " "$work/vs_$r.log" | head -20 >&2 || true
            status=1
        fi
    done
    # the deliberately mismatched pair: the other profile's instrumented run on Sapphire Rapids
    ref=$d/runs/instrumented-$other-sapphirerapids
    if [ -d "$ref" ]; then
        r=$(basename "$ref")
        if "$bin" validate cluster --reference "$ref" --candidate "$work/run" --optimizer "$d/optimizer" \
                --json "$work/vs_other_$r.json" --markdown "$work/vs_other_$r.md" 2> "$work/vs_other_$r.log"; then
            echo "  vs $r (other kernel): $(tail -1 "$work/vs_other_$r.log")"
        else
            echo "  vs $r (other kernel): FAIL (reported, not counted)"
        fi
    fi
    if [ "${THREADS_CHECK:-$([ "$tier" = T ] && echo 1 || echo 0)}" = 1 ]; then
        "$bin" cluster --counts "$d/input" --preset workflow --out "$work/run_t1" --dump-stages "$level" "${keep[@]}" --threads 1 2> "$work/stderr_t1.txt"
        differ=0
        while IFS= read -r f; do
            [ "$f" = summary.json ] && continue
            cmp -s "$work/run/$f" "$work/run_t1/$f" || { echo "    $f differs at 1 thread" >&2; differ=$((differ + 1)); }
        done < <(cd "$work/run" && find . -type f | sed 's|^\./||' | sort)
        [ "$(cd "$work/run" && find . -type f | wc -l)" = "$(cd "$work/run_t1" && find . -type f | wc -l)" ] || differ=$((differ + 1))
        cmp -s <(result_section "$work/run/summary.json") <(result_section "$work/run_t1/summary.json") || { echo "    result section differs" >&2; differ=$((differ + 1)); }
        n=$(cd "$work/run" && find . -type f | wc -l)
        echo "threads $threads and 1: $n files, $differ differing (summary.json: its result section)" | tee "$work/threads.txt"
        [ "$differ" = 0 ] || status=1
        rm -rf "$work/run_t1"
    fi
    mv "$work" "$out"
done
exit $status
