#!/bin/bash
# Milestone M15b: scrust cluster --counts with SCTransform's vars.to.regress, one array task per
# row of validation/m15b_cases.tsv, in one profile, each run compared by scrust validate cluster
# with the case's R reference runs (validation/slurm/m15b_references.sbatch): the production R run
# (exact in the production profile; kernel-gated, so graded against [clustering.floor.kernel], in
# the portable one) and, when there is one, the portable R run (exact in the portable profile).
# The Rust options are the R run's, but where R's pattern finds no feature: there scrust's run
# gives R's all-zero percent.mt through a list of features without counts (the case's
# rust_mt_genes), and the task also checks that scrust refuses --vars-to-regress percent.mt by the
# pattern, naming --mt-genes. --keep-stages where R kept scale.data (levels full and standard), so
# it is compared value by value too.
#   SCRUST_BIN=<copy> sbatch --array=1-134 validation/slurm/m15b_cluster.sbatch <label> production|portable
# Output, never overwritten: data/clustering/m15b/rust/<label>/<object>/<case>/<profile>/: run/,
# validate_<R profile>.json and .md, refusal.txt, time_v.txt, stderr.txt, fingerprint.json.
#
#SBATCH --job-name=scrust-m15b-cl
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
root=$PWD
[ $# -eq 2 ] || { echo "usage: m15b_cluster.sbatch <label> production|portable" >&2; exit 2; }
label=$1 profile=$2
case $profile in
    portable) ;;
    production) export OPENBLAS_CORETYPE=auto ;;
    *) echo "profile: portable or production" >&2; exit 2 ;;
esac
source ./env.sh
bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' kernel=$core binary=$bin $(sha256sum "$bin" | cut -c1-16)"
if [ "$profile" = production ] && [ "$core" != SapphireRapids ]; then echo "production needs SapphireRapids, got '$core'" >&2; exit 3; fi
if [ "$profile" = portable ] && [ "$core" != Haswell ]; then echo "portable needs Haswell, got '$core'" >&2; exit 3; fi
threads=${SLURM_CPUS_PER_TASK:-1}

row=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n' validation/m15b_cases.tsv)
[ -n "$row" ] || { echo "no case $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
field() { cut -f"$1" <<< "$row"; }
obj=$(field 1); input=$(field 3); level=$(field 5); cs=$(field 6); vars=$(field 7); pm=$(field 8)
mt=$(field 9); md=$(field 10); rust_mt=$(field 11)
out=data/clustering/m15b/rust/$label/$obj/$cs/$profile
[ -e "$out" ] && { echo "skip $out: exists"; exit 0; }
tmp=$out.partial.${SLURM_JOB_ID:-$$}
rm -rf "$tmp"
mkdir -p "$tmp"
"$bin" fingerprint > "$tmp/fingerprint.json"
args=(cluster --counts "$input" --preset workflow --out "$tmp/run" --threads "$threads" --vars-to-regress "$vars")
[ "$pm" = yes ] && args+=(--percent-mt)
[ -z "$rust_mt" ] || args+=(--mt-genes "$rust_mt")
[ -z "$md" ] || args+=(--metadata "$md")
[ "$level" = light ] || args+=(--keep-stages)
echo "== $obj $cs ($profile, $core, $threads threads): scrust ${args[*]}"
/usr/bin/time -v -o "$tmp/time_v.txt" "$bin" "${args[@]}" 2> "$tmp/stderr.txt" || { cat "$tmp/stderr.txt" >&2; exit 1; }
tail -n 3 "$tmp/stderr.txt"
# the refusal: percent.mt by the pattern where it finds no feature
if [ "$cs" = pm ] && [ -z "$mt" ] && [ -n "$rust_mt" ]; then
    if "$bin" cluster --counts "$input" --preset workflow --out "$tmp/refused" --vars-to-regress percent.mt > /dev/null 2> "$tmp/refusal.txt"; then
        echo "refusal: NOT REFUSED" | tee -a "$tmp/refusal.txt" >&2
        exit 1
    fi
    grep -q -- "--mt-genes" "$tmp/refusal.txt" || { echo "refusal: no remedy named" >&2; exit 1; }
    rm -rf "$tmp/refused"
    echo "refusal: $(cat "$tmp/refusal.txt" | head -n 1)"
fi
status=0
for rp in production portable; do
    ref=$(ls -d data/clustering/m15b/runs/$obj/$cs/direct-$rp-* 2>/dev/null | grep -v partial | head -n 1 || true)
    [ -n "$ref" ] || continue
    if "$bin" validate cluster --reference "$ref" --candidate "$tmp/run" --json "$tmp/validate_$rp.json" --markdown "$tmp/validate_$rp.md" > "$tmp/validate_$rp.txt" 2>&1; then
        echo "against R ($rp): $(tail -n 1 "$tmp/validate_$rp.txt")"
    else
        echo "against R ($rp): $(tail -n 1 "$tmp/validate_$rp.txt")" >&2
        status=1
    fi
done
mv "$tmp" "$out"
exit $status
