#!/bin/bash
# Milestone M15 (docs/PLAN.md C3h, C5's K12 and the singlets' re-clustering): `scrust doubletfinder`
# on each of the 150 production DoubletFinder objects of validation/doubletfinder_panel.tsv, from
# M7's export of its input, in the production profile on the CPU type of the object's workflow
# job, as M7 ran (m7_panel.sbatch). Five runs per object, each in a working directory of its own
# with the same relative paths, so the paths the runs record are the same:
#   given/       scrust doubletfinder --counts <input> --clusters <input labels> --out run
#                (M7's command: check (c), Step 1 unchanged), then the two-command route's second
#                step: scrust cluster --counts run/singlets --preset workflow --out run/singlets_clusters
#   hook/        scrust doubletfinder --counts <input> --out run: clustering on the fly (K12)
#   inprocess/   scrust doubletfinder --counts <input> --clusters <input labels> --out run
#                --recluster: the whole rule in one process
#   alone/       scrust cluster --counts <input> --preset workflow --stats none --out run/clusters:
#                the command whose output hook/run/clusters/ is, pca/ and neighbors/ aside (the
#                hook's --cluster-stats is none by default; the workflow preset's [stats] is not)
# The checks, each into checks/<name>.json (and printed):
#   c_step1          given/run against M7's run (label m7-b8c8a8c): cells.tsv, the sweep tables,
#                    singlets/, byte for byte, and summary.json's result
#   k12_labels       hook/run/clusters/seurat_clusters.tsv and its levels against the input's
#                    exported seurat_clusters, byte for byte
#   k12_files        hook/run against M7's run, as c_step1, clusters/ and result.clustering aside
#   k12_saved        hook/run against the saved workflow output (scrust validate doubletfinder
#                    --saved: M7's nine checks)
#   k12_logged       the hook's clustering against the numbers the input's workflow job logged
#                    (validation/clustering_panel.tsv)
#   k12_alone        hook/run/clusters/ against alone/run/clusters/: every file the nested output
#                    has, byte for byte, and its summary.json as JSON, the timings, the command line
#                    and the manifest's entries of pca/ and neighbors/ (not written there) aside
#   a_files          inprocess/run against given/run: DoubletFinder's files, and singlets_clusters/ against
#                    the two-command route's, byte for byte, summary.json's timings and command
#                    line aside (check (a) without the export; with it on T and P, m15_recluster.sbatch)
#   b_saved          inprocess/run/singlets_clusters against the saved object's re-clustering (the export of
#                    m15_export.sbatch): labels, levels, silhouette widths, purities, the RNA
#                    columns, and the numbers its workflow job logged
#   manifests        scrust validate manifest on hook/run and inprocess/run (nested outputs)
# Memory and time: /usr/bin/time -v of every run (time_*.txt).
#
# One array task per object (task i: the table's data row i). From the repository's top directory,
# with a pinned copy of the binary, M7's runs, and M15's exports (m15_export.sbatch) done:
#   SCRUST_BIN=<copy> sbatch --array=1-150%50 [--output=<dir>/%x-%A_%a.sbatch.log] \
#       validation/slurm/m15_panel.sbatch <label>
# Output, never overwritten (an existing one is skipped):
# data/doubletfinder/<object>/rust/<label>/production/: given/, hook/, inprocess/, alone/, checks/,
# time_*.txt, stderr_*.txt, fingerprint.json, node_features.txt, row.tsv (this object's line of
# the summary table). Work in progress goes to <dir>.partial.<job id>, renamed when complete. Exit
# status 1 unless every check is made and identical.
#
#SBATCH --job-name=scrust-m15-df
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
root=$PWD
# production profile: env.sh unsets OPENBLAS_CORETYPE for "auto", so OpenBLAS detects the CPU
export OPENBLAS_CORETYPE=auto
source ./env.sh
export PYTHONNOUSERSITE=1
py="$root/.envs/loaders-py/bin/python -I"

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
[ $# -eq 1 ] || { echo "usage: m15_panel.sbatch <label>" >&2; exit 2; }
label=$1
threads=${SLURM_CPUS_PER_TASK:-1}
m7=m7-b8c8a8c

panel=validation/doubletfinder_panel.tsv
i=${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch --array=1-150 validation/slurm/m15_panel.sbatch <label>}
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

base=$root/data/doubletfinder/$obj
inp=$base/input
m7run=$base/rust/$m7/production/run
saved_cl=$root/data/clustering/doubletfinder_$obj/saved
for f in "$inp/export.json" "$inp/seurat_clusters.tsv" "$base/saved/saved.json" "$m7run/summary.json" "$saved_cl/metadata.tsv"; do
    [ -f "$f" ] || { echo "$obj: no $f (M7's export and run, M15's export)" >&2; exit 1; }
done
out=$base/rust/$label/production
if [ -e "$out" ]; then
    echo "skip $out: exists"
    exit 0
fi
tmp=$out.partial.${SLURM_JOB_ID:-local}
rm -rf "$tmp"
mkdir -p "$tmp"/{given,hook,inprocess,alone,checks}
echo "$node $features" > "$tmp/node_features.txt"
"$bin" fingerprint > "$tmp/fingerprint.json"
got=$(sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' "$tmp/fingerprint.json" | head -n 1)
if [ "$got" != "$core" ]; then
    echo "$obj: OpenBLAS reports core '$got' on $node, expected $core" >&2
    rm -rf "$tmp"
    exit 1
fi

echo "== task $i: $obj on $node ($cpu, OpenBLAS $got), $threads threads"
# run <name> <working dir> <scrust arguments...>: timed, stderr kept
run() {
    local name=$1 dir=$2
    shift 2
    if ! (cd "$dir" && /usr/bin/time -v -o "$tmp/time_$name.txt" "$bin" "$@" 2> "$tmp/stderr_$name.txt"); then
        cat "$tmp/stderr_$name.txt" >&2
        echo "$obj: run $name failed; partial output left in $tmp" >&2
        exit 1
    fi
    tail -n 1 "$tmp/stderr_$name.txt"
    sed -n 's/^\s*\(Elapsed (wall clock) time\|Maximum resident set size\).*: /&/p' "$tmp/time_$name.txt"
}
run given "$tmp/given" doubletfinder --counts "$inp" --clusters "$inp/seurat_clusters.tsv" --out run --threads "$threads"
run recluster2 "$tmp/given" cluster --counts run/singlets --preset workflow --out run/singlets_clusters --threads "$threads"
run hook "$tmp/hook" doubletfinder --counts "$inp" --out run --threads "$threads"
run inprocess "$tmp/inprocess" doubletfinder --counts "$inp" --clusters "$inp/seurat_clusters.tsv" --out run --recluster --threads "$threads"
run alone "$tmp/alone" cluster --counts "$inp" --preset workflow --stats none --out run/clusters --threads "$threads"

status=0
cmp_runs() { $py "$root/validation/python/compare_runs.py" "$@"; }
check() {
    local name=$1
    shift
    echo "-- $name"
    if "$@"; then echo "$name: identical"; else echo "$name: DIFFERS" >&2; status=1; fi
}
check c_step1 cmp_runs "$m7run" "$tmp/given/run" --result-only --only-right singlets_clusters --json "$tmp/checks/c_step1.json"
k12_labels() {
    cmp "$tmp/hook/run/clusters/seurat_clusters.tsv" "$inp/seurat_clusters.tsv" &&
        cmp "$tmp/hook/run/clusters/seurat_clusters_levels.txt" "$inp/seurat_clusters_levels.txt"
}
check k12_labels k12_labels
check k12_files cmp_runs "$m7run" "$tmp/hook/run" --result-only --only-right clusters --ignore summary.json:result.clustering --json "$tmp/checks/k12_files.json"
k12_saved() {
    "$bin" validate doubletfinder --saved "$base/saved" --candidate "$tmp/hook/run" \
        --json "$tmp/checks/k12_saved.json" --markdown "$tmp/checks/k12_saved.md" || return 1
    local made same
    made=$(grep -c '"status":' "$tmp/checks/k12_saved.json" || true)
    same=$(grep -c '"status": "identical"' "$tmp/checks/k12_saved.json" || true)
    echo "checks identical: $same of $made"
    [ "$made" -eq 9 ] && [ "$same" -eq 9 ]
}
check k12_saved k12_saved
check k12_logged $py "$root/validation/python/compare_saved_clustering.py" "$tmp/hook/run/clusters" - \
    --panel "$root/validation/clustering_panel.tsv" --object "$obj" --json "$tmp/checks/k12_logged.json"
check k12_alone cmp_runs "$tmp/alone/run/clusters" "$tmp/hook/run/clusters" --only-left pca --only-left neighbors \
    --json "$tmp/checks/k12_alone.json"
check a_files cmp_runs "$tmp/given/run" "$tmp/inprocess/run" --ignore summary.json:run --ignore summary.json:files \
    --ignore summary.json:result.recluster --json "$tmp/checks/a_files.json"
check b_saved $py "$root/validation/python/compare_saved_clustering.py" "$tmp/inprocess/run/singlets_clusters" "$saved_cl" \
    --panel "$root/validation/clustering_panel.tsv" --object "doubletfinder_$obj" --json "$tmp/checks/b_saved.json"
manifests() { "$bin" validate manifest --dir "$tmp/hook/run" && "$bin" validate manifest --dir "$tmp/inprocess/run"; }
check manifests manifests

# this object's line of the summary table: the cells, the checks, peak RSS (KiB) and wall time (s)
# of each run, and M7's
peak() { sed -n 's/^\s*Maximum resident set size (kbytes): //p' "$1"; }
wall() { sed -n 's/^\s*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$1" | awk -F: '{ s = 0; for (k = 1; k <= NF; k++) s = s * 60 + $k; print s }'; }
cells=$(zcat "$inp/barcodes.tsv.gz" | wc -l)
singlets=$(zcat "$tmp/given/run/singlets/barcodes.tsv.gz" | wc -l)
{
    printf 'object\tcells\tsinglets\tstatus\tm7_peak_kib\tm7_wall_s\tgiven_peak_kib\tgiven_wall_s\thook_peak_kib\thook_wall_s\trecluster2_peak_kib\trecluster2_wall_s\tinprocess_peak_kib\tinprocess_wall_s\talone_peak_kib\talone_wall_s\tnode\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$obj" "$cells" "$singlets" "$([ $status -eq 0 ] && echo pass || echo FAIL)" \
        "$(peak "$base/rust/$m7/production/time_v.txt")" "$(wall "$base/rust/$m7/production/time_v.txt")" \
        "$(peak "$tmp/time_given.txt")" "$(wall "$tmp/time_given.txt")" "$(peak "$tmp/time_hook.txt")" "$(wall "$tmp/time_hook.txt")" \
        "$(peak "$tmp/time_recluster2.txt")" "$(wall "$tmp/time_recluster2.txt")" "$(peak "$tmp/time_inprocess.txt")" "$(wall "$tmp/time_inprocess.txt")" \
        "$(peak "$tmp/time_alone.txt")" "$(wall "$tmp/time_alone.txt")" "$node"
} > "$tmp/row.tsv"
cat "$tmp/row.tsv"
mv "$tmp" "$out"
exit $status
