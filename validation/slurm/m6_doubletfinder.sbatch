#!/bin/bash
# Milestone M6 (docs/PLAN.md B5 and B6): the whole `scrust doubletfinder`, from the counts, the
# cluster labels and the seed, against every DoubletFinder reference run of the given objects
# and their saved workflow output. From the repository's top directory, after
# `cargo build --release -p scrust` and copying the binary to a file of its own:
#   SCRUST_BIN=<copy> sbatch [--cpus-per-task=N] [--mem=...] validation/slurm/m6_doubletfinder.sbatch <label> <object> [<object> ...]
#
# For each object, one Rust run per OpenBLAS profile its reference runs use (portable:
# OPENBLAS_CORETYPE=Haswell; production: the variable unset, so that OpenBLAS detects
# SapphireRapids; this job must run on a `sapphire` node), with stage dumps at the level of the
# profile's instrumented R run (light where there is none), then:
#   - against each instrumented R run of the profile (Sapphire Rapids and, for the portable
#     profile, Cascade Lake: the kernel, not the CPU, sets the bits): the final files, the
#     result section, every stage array block S1-S7, every RNG checkpoint, every scalar;
#   - against each direct R run of the profile: the final files and the result section;
#   - against the saved workflow output: the doublet set, the kept cells' pANN, the columns.
# Output, data/doubletfinder/<object>/rust/<label>/<profile>/:
#   run/                      the Rust run: cells.tsv, sweep_stats.tsv, bcmvn.tsv, summary.json,
#                             stages/
#   vs_<R run>.json, .md      each comparison (`scrust validate doubletfinder`)
#   vs_saved.json, .md        the saved workflow output
#   fingerprint.json          scrust fingerprint under the profile
#   time_v.txt, stderr.txt    /usr/bin/time -v and the run's messages
# An existing output directory is never overwritten: its profile is skipped. Work in progress
# goes to <dir>.partial.<job id>, renamed when complete. Exit status 1 if anything differs.
#
# SCRUST_BIN (environment) names the scrust binary to run (default target/release/scrust): a
# pinned copy, so that later builds do not change it. OpenBLAS runs one thread; the job's CPUs
# are scrust's threads, which change nothing but the time.
#
#SBATCH --job-name=scrust-m6-df
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
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
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -ge 2 ] || { echo "usage: m6_doubletfinder.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
threads=${SLURM_CPUS_PER_TASK:-1}
status=0

for obj in "$@"; do
    base=data/doubletfinder/$obj
    [ -d "$base/runs" ] || { echo "no runs for $obj" >&2; status=1; continue; }
    runs=$(cd "$base/runs" && ls -d */ 2>/dev/null | tr -d / | grep -v '\.partial\.' || true)
    for profile in portable production; do
        mine=$(printf '%s\n' $runs | grep -- "-$profile-" || true)
        [ -n "$mine" ] || continue
        out=$base/rust/$label/$profile
        if [ -e "$out" ]; then
            echo "skip $out: exists"
            continue
        fi
        case $profile in
            portable) core=Haswell; envset=(env OPENBLAS_CORETYPE=Haswell) ;;
            production) core=SapphireRapids; envset=(env -u OPENBLAS_CORETYPE) ;;
        esac
        # the dump level of the profile's instrumented R run (Sapphire Rapids first)
        level=light
        for r in instrumented-$profile-sapphirerapids $(printf '%s\n' $mine | grep '^instrumented-' || true); do
            j=$base/runs/$r/stages/stages.json
            if [ -f "$j" ]; then
                level=$(sed -n 's/.*"level": *"\([a-z]*\)".*/\1/p' "$j" | head -n 1)
                break
            fi
        done
        tmp=$out.partial.${SLURM_JOB_ID:-local}
        mkdir -p "$tmp"
        "${envset[@]}" "$bin" fingerprint > "$tmp/fingerprint.json"
        got=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
        if [ "$got" != "$core" ]; then
            echo "$obj $profile: OpenBLAS reports core '$got', expected $core; skipped" >&2
            status=1
            rm -rf "$tmp"
            continue
        fi
        echo "== $obj $profile ($core), dumps $level, $threads threads"
        if ! "${envset[@]}" /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" doubletfinder \
                --counts "$base/input" --clusters "$base/input/seurat_clusters.tsv" \
                --out "$tmp/run" --dump-stages "$level" --threads "$threads" 2> "$tmp/stderr.txt"; then
            cat "$tmp/stderr.txt" >&2
            echo "$obj $profile: scrust doubletfinder failed" >&2
            status=1
            continue
        fi
        cat "$tmp/stderr.txt" >&2
        for r in $mine; do
            echo "-- $obj $profile against $r"
            if ! "$bin" validate doubletfinder --reference "$base/runs/$r" --candidate "$tmp/run" \
                    --json "$tmp/vs_$r.json" --markdown "$tmp/vs_$r.md" --threads "$threads"; then
                status=1
            fi
        done
        echo "-- $obj $profile against the saved workflow output"
        if ! "$bin" validate doubletfinder --saved "$base/saved" --candidate "$tmp/run" \
                --json "$tmp/vs_saved.json" --markdown "$tmp/vs_saved.md"; then
            status=1
        fi
        mv "$tmp" "$out"
    done
done
exit $status
