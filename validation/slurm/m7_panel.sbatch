#!/bin/bash
# Milestone M7 (docs/PLAN.md B4-B6): the full panel. `scrust doubletfinder` on each of the 150
# production DoubletFinder objects of validation/doubletfinder_panel.tsv, from its export's
# counts and cluster labels (validation/slurm/m7_export.sbatch), compared with the workflow's
# saved output by `scrust validate doubletfinder --saved`: the doublet set, the kept cells' pANN,
# pK and nExp (the column names), and the singlets' counts (singlets/ against the saved object's
# RNA counts). R is not run.
#
# One array task per object (task i: the table's data row i), in the production profile: on the
# CPU type the object's workflow job ran on (the table's workflow_cpu, from sacct and the node's
# SLURM features; all 150 ran on Sapphire Rapids), with OPENBLAS_CORETYPE unset so that OpenBLAS
# detects the CPU, as it did in the workflow. The task stops before the run if the node lacks
# that CPU feature or OpenBLAS reports another core. No stage dumps. From the repository's top
# directory, with a pinned copy of the binary (cp target/release/scrust
# data/doubletfinder/logs/bin/scrust-<label>), and the exports done:
#   SCRUST_BIN=<copy> sbatch --array=1-150 [--output=<dir>/%x-%A_%a.sbatch.log] \
#       validation/slurm/m7_panel.sbatch <label>
# Output, never overwritten (an existing one is skipped), data/doubletfinder/<object>/rust/<label>/production/:
#   run/                     cells.tsv, sweep_stats.tsv, bcmvn.tsv, summary.json, singlets/
#   vs_saved.json, .md       the comparison with the saved workflow output
#   fingerprint.json         scrust fingerprint under the profile (node, CPU, OpenBLAS core)
#   node_features.txt        the node and its SLURM features
#   time_v.txt, stderr.txt   /usr/bin/time -v and the run's messages
# Work in progress goes to <dir>.partial.<job id>, renamed when complete. Exit status 1 unless
# every check of the comparison is made and identical. The job's CPUs are scrust's threads,
# which change nothing but the time; OpenBLAS runs one thread.
#
#SBATCH --job-name=scrust-m7-df
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
# production profile: env.sh unsets OPENBLAS_CORETYPE for "auto", so OpenBLAS detects the CPU
export OPENBLAS_CORETYPE=auto
source ./env.sh

# Environment fingerprint, for the log (A5b)
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -eq 1 ] || { echo "usage: m7_panel.sbatch <label>" >&2; exit 2; }
label=$1
threads=${SLURM_CPUS_PER_TASK:-1}

panel=validation/doubletfinder_panel.tsv
i=${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch --array=1-150 validation/slurm/m7_panel.sbatch <label>}
# the object and its workflow job's CPU type, by header name
read -r obj cpu < <(awk -F'\t' -v i="$i" 'NR == 1 { for (c = 1; c <= NF; c++) h[$c] = c; next }
    NR == i + 1 { print $h["object"], $h["workflow_cpu"] }' "$panel") || true
[ -n "${obj:-}" ] || { echo "$panel has no data row $i" >&2; exit 2; }
case $cpu in
    sapphirerapids) core=SapphireRapids ;;   # what OpenBLAS 0.3.34 detects there (M2)
    *) echo "$obj: no OpenBLAS core known for CPU type '$cpu'; run it where that CPU is" >&2; exit 2 ;;
esac
node=${SLURMD_NODENAME:-$(hostname -s)}
features=$(scontrol show node "$node" | sed -n 's/.*AvailableFeatures=\([^ ]*\).*/\1/p')
case ",$features," in
    *",$cpu,"*) ;;
    *) echo "$obj: node $node ($features) is not $cpu, the CPU type of the workflow's job" >&2; exit 1 ;;
esac

base=data/doubletfinder/$obj
for f in input/export.json saved/saved.json saved/singlets/export.json; do
    [ -f "$base/$f" ] || { echo "$obj: no $base/$f; export first (m7_export.sbatch)" >&2; exit 1; }
done
out=$base/rust/$label/production
if [ -e "$out" ]; then
    echo "skip $out: exists"
    exit 0
fi
tmp=$out.partial.${SLURM_JOB_ID:-local}
rm -rf "$tmp"   # this job's own incomplete output, if it was requeued
mkdir -p "$tmp"
echo "$node $features" > "$tmp/node_features.txt"
"$bin" fingerprint > "$tmp/fingerprint.json"
got=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
if [ "$got" != "$core" ]; then
    echo "$obj: OpenBLAS reports core '$got' on $node, expected $core" >&2
    rm -rf "$tmp"
    exit 1
fi

echo "== task $i: $obj on $node ($cpu, OpenBLAS $got), $threads threads"
if ! /usr/bin/time -v -o "$tmp/time_v.txt" "$bin" doubletfinder \
        --counts "$base/input" --clusters "$base/input/seurat_clusters.tsv" \
        --out "$tmp/run" --threads "$threads" 2> "$tmp/stderr.txt"; then
    cat "$tmp/stderr.txt" >&2
    echo "$obj: scrust doubletfinder failed; partial output left in $tmp" >&2
    exit 1
fi
cat "$tmp/stderr.txt" >&2
sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$tmp/time_v.txt"

status=0
echo "-- $obj against the saved workflow output"
"$bin" validate doubletfinder --saved "$base/saved" --candidate "$tmp/run" \
    --json "$tmp/vs_saved.json" --markdown "$tmp/vs_saved.md" || status=1
# Stricter than the exit status, which fails only on a difference: all nine checks made (the
# column names, pK, nExp, the doublet set, the kept cells' pANN, and the singlets' barcodes,
# features and matrix) and every one identical.
made=$(grep -c '"status":' "$tmp/vs_saved.json" || true)
same=$(grep -c '"status": "identical"' "$tmp/vs_saved.json" || true)
echo "checks identical: $same of $made"
[ "$made" -eq 9 ] && [ "$same" -eq 9 ] || status=1
mv "$tmp" "$out"
exit $status
