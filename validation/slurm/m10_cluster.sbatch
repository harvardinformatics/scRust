#!/bin/bash
# Milestone M10 (docs/PLAN.md C5 K7, C6): `scrust cluster --snn` on R's own SNN graphs (teacher
# forcing, plan A5b), each run compared with its R run by `scrust validate cluster`. From the
# repository's top directory, after `cargo build --release -p scrust`:
#   sbatch [--cpus-per-task=N] validation/slurm/m10_cluster.sbatch <object> [<object> ...]
# For each object:
#   - the workflow's settings (--seed 12345) on every M8 instrumented run's own graph
#     (data/clustering/<object>/runs/instrumented-*/graphs/snn.mtx.gz), into
#     data/clustering/<object>/rust/m10/<run>/, against the run, the instrumented optimizer
#     (data/clustering/<object>/optimizer/) and, for tier P, the workflow job's logged numbers
#     (validation/clustering_panel.tsv);
#   - every FindClusters reference setting (cluster_references.sbatch settings,
#     data/clustering/<object>/cluster_refs/<setting>/) on the saved graph, with the setting's
#     arguments, into data/clustering/<object>/rust/m10/refs/<setting>/, against that run and its
#     optimizer/. Where R stopped with an error, Rust must stop with R's message
#     (refs/<setting>/error_check.txt). The `resolutions` setting runs again at one thread, and
#     every file must be byte-identical (summary.json: its result section).
# Each output directory holds Rust's files (with --keep-stages: every start's partition),
# time_v.txt and validate.{json,md}. SCRUST_BIN names the binary (default target/release/scrust;
# a copy keeps a long run safe from rebuilds). SKIP_REFS=1 runs only the instrumented runs.
# An existing output directory is never overwritten: its run is skipped. Work in progress goes to
# <dir>.partial.<job id>, renamed when complete. Exit status 1 if any comparison differs.
#
# No BLAS is involved, so any CPU and any thread count are fine; the results depend on neither.
#
#SBATCH --job-name=scrust-m10-cl
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary $bin sha256 $(sha256sum "$bin" | cut -d' ' -f1)"
threads=${SLURM_CPUS_PER_TASK:-1}
job=${SLURM_JOB_ID:-local}
status=0

# FindClusters' arguments of each reference setting (validation/R/cluster_references.R), as flags
setting_flags() {  # <setting> <tier>
    case $1 in
        workflow) echo "--seed 12345" ;;
        algorithm2) echo "--seed 12345 --algorithm 2" ;;
        algorithm3) echo "--seed 12345 --algorithm 3" ;;
        modularity2) echo "--seed 12345 --modularity 2" ;;
        modularity2_resolutions)
            if [ "$2" = T ]; then echo "--seed 12345 --modularity 2 --resolution 0.05,0.1,0.2"
            else echo "--seed 12345 --modularity 2 --resolution 0.002,0.005,0.01"; fi ;;
        resolutions) echo "--seed 12345 --resolution 0.2,0.4,0.8,1.2,2.0" ;;
        nogroup) echo "--seed 12345 --no-group-singletons" ;;
        nogroup_resolutions) echo "--seed 12345 --no-group-singletons --resolution 0.8,2.0" ;;
        starts1) echo "--seed 12345 --starts 1" ;;
        starts20) echo "--seed 12345 --starts 20" ;;
        iter1) echo "--seed 12345 --iterations 1" ;;
        iter50) echo "--seed 12345 --iterations 50" ;;
        seed0) echo "--seed 0" ;;
        seed1) echo "--seed 1" ;;
        *) return 1 ;;
    esac
}

# one Rust run against one R run: <graph dir> <R dir> <out> <optimizer dir> <panel object or -> <flags...>
one() {
    local graph=$1 ref=$2 out=$3 opt=$4 panel_obj=$5
    shift 5
    if [ -e "$out" ]; then
        echo "skip $out: exists"
        return 0
    fi
    local tmp=$out.partial.$job
    rm -rf "$tmp"
    mkdir -p "$(dirname "$out")"
    local err
    err=$(sed -n 's/^    "error": "\(.*\)",\{0,1\}$/\1/p' "$ref/summary.json" | head -n 1)
    if [ -n "$err" ]; then
        # R stopped: Rust must stop too, with R's message
        mkdir -p "$tmp"
        local rc=0
        "$bin" cluster --snn "$graph" --out "$tmp/out" --stats none --threads "$threads" "$@" 2> "$tmp/stderr.txt" || rc=$?
        if [ "$rc" -ne 0 ] && grep -qF -- "$err" "$tmp/stderr.txt"; then
            printf 'identical: R and Rust stop with "%s" (exit %s)\n' "$err" "$rc" > "$tmp/error_check.txt"
        else
            printf 'DIFFERS: R stopped with "%s"; Rust exit %s: %s\n' "$err" "$rc" "$(head -n 1 "$tmp/stderr.txt")" > "$tmp/error_check.txt"
            status=1
        fi
        cat "$tmp/error_check.txt"
        mv "$tmp" "$out"
        return 0
    fi
    /usr/bin/time -v -o "$tmp.time_v.txt" \
        "$bin" cluster --snn "$graph" --out "$tmp" --keep-stages --stats none --threads "$threads" "$@" 2> "$tmp.stderr.txt" \
        || { echo "FAILED: $out"; cat "$tmp.stderr.txt"; status=1; return 0; }
    mv "$tmp.time_v.txt" "$tmp/time_v.txt"
    mv "$tmp.stderr.txt" "$tmp/stderr.txt"
    local extra=()
    [ -d "$opt" ] && extra+=(--optimizer "$opt")
    [ "$panel_obj" != - ] && extra+=(--panel validation/clustering_panel.tsv --object "$panel_obj")
    if ! "$bin" validate cluster --reference "$ref" --candidate "$tmp" "${extra[@]}" \
        --json "$tmp/validate.json" --markdown "$tmp/validate.md" 2> "$tmp/validate.log"; then
        echo "DIFFERS: $out"
        grep -v ' identical \| n/a ' "$tmp/validate.log" | head -n 20
        status=1
    fi
    mv "$tmp" "$out"
}

for obj in "$@"; do
    base=data/clustering/$obj
    tier=$(awk -F'\t' -v o="$obj" '$1 == o { print $2 }' validation/clustering_datasets.tsv)
    [ -n "$tier" ] || { echo "unknown object $obj" >&2; status=1; continue; }
    panel=-
    [ "$tier" = P ] && panel=$obj
    for r in "$base"/runs/instrumented-*/; do
        r=${r%/}
        run=$(basename "$r")
        case $run in *.partial.*) continue ;; esac
        one "$r/graphs" "$r" "$base/rust/m10/$run" "$base/optimizer" "$panel" --seed 12345
    done
    if [ -d "$base/cluster_refs" ] && [ "${SKIP_REFS:-0}" != 1 ]; then
        for d in "$base"/cluster_refs/*/; do
            d=${d%/}
            s=$(basename "$d")
            [ -f "$d/summary.json" ] || continue
            flags=$(setting_flags "$s" "$tier") || { echo "unknown setting $s" >&2; status=1; continue; }
            # shellcheck disable=SC2086
            one "$base/saved/graphs" "$d" "$base/rust/m10/refs/$s" "$d/optimizer" - $flags
            if [ "$s" = resolutions ] && [ ! -e "$base/rust/m10/refs/resolutions.t1" ] && [ -f "$base/rust/m10/refs/resolutions/summary.json" ]; then
                t1=$base/rust/m10/refs/resolutions.t1
                # shellcheck disable=SC2086
                "$bin" cluster --snn "$base/saved/graphs" --out "$t1" --keep-stages --stats none --threads 1 $flags 2> /dev/null
                same=1
                for f in "$base"/rust/m10/refs/resolutions/*.tsv "$base"/rust/m10/refs/resolutions/*.txt "$base"/rust/m10/refs/resolutions/stages/*; do
                    g=$t1/${f#"$base"/rust/m10/refs/resolutions/}
                    case $f in */time_v.txt|*/stderr.txt|*/validate.*) continue ;; esac
                    cmp -s "$f" "$g" || { echo "threads: $f differs at 1 thread"; same=0; }
                done
                a=$(sed -n '/^  "result": {/,/^  },$/p' "$base/rust/m10/refs/resolutions/summary.json" | sha256sum)
                b=$(sed -n '/^  "result": {/,/^  },$/p' "$t1/summary.json" | sha256sum)
                [ "$a" = "$b" ] || { echo "threads: summary.json result differs at 1 thread"; same=0; }
                if [ "$same" = 1 ]; then
                    echo "identical at $threads threads and at 1 thread" > "$t1/threads_check.txt"
                else
                    echo "DIFFERS between $threads threads and 1 thread" > "$t1/threads_check.txt"
                    status=1
                fi
                cat "$t1/threads_check.txt"
            fi
        done
    fi
done
echo "status $status"
exit $status
