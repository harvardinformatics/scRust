#!/bin/bash
# Milestone M15: the singlets' re-clustering in process (scrust doubletfinder --recluster) against
# the two-command route, with and without the export, on the five tier-T and five tier-P
# DoubletFinder objects, in one profile (portable: OPENBLAS_CORETYPE=Haswell, any node;
# production: OpenBLAS's own kernel on Sapphire Rapids, as the saved objects were made). Check (a):
# the in-process route's files equal the two-command route's (scrust doubletfinder, then
# scrust cluster --counts <out>/singlets --preset workflow [--export]), byte for byte, apart from
# the timings and the command line in summary.json; each route runs in a working directory of its
# own with the same relative paths, so the paths they record are the same. Check (b), against the
# saved DoubletFinder object's re-clustering (M8's exports on T, M15's on P): its labels and
# levels, silhouette widths, purities, the RNA columns and, with the export, nCount_SCT and
# nFeature_SCT; on P also the numbers its workflow job logged. And clustering on the fly with
# both scores (--cluster-stats silhouette,purity): its clusters/ against scrust cluster --counts
# <input> --stats silhouette,purity's output (every file it has, byte for byte, and its
# summary.json, pca/ and neighbors/ aside; Rust against Rust, so held to (a)'s standard), and
# against the saved input object's (M8's exports of the inputs on T and P), and on P the numbers
# the input's job logged. The saved objects were
# made under the SapphireRapids kernel: in the portable profile, (b) and the scores are reported,
# their differences graded by the floor elsewhere (M11, M12), not here.
#
#   SCRUST_BIN=<copy> sbatch --array=1-10 validation/slurm/m15_recluster.sbatch <label> portable|production
#
# Tasks 1-5 are tier T, 6-10 tier P (the objects below). Output, never overwritten:
# data/doubletfinder/<object>/rust/<label>/<profile>_recluster/: the routes' directories
# (inprocess/, inprocess_export/, two/, two_export/, hook/, alone/), checks/*.json, time_*.txt,
# stderr_*.txt, fingerprint.json. Exit status 1 unless (a) and hook_alone hold, and in the
# production profile every check.
#
#SBATCH --job-name=scrust-m15-rc
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-03:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
root=$PWD
[ $# -eq 2 ] || { echo "usage: m15_recluster.sbatch <label> portable|production" >&2; exit 2; }
label=$1 profile=$2
case $profile in
    portable) ;;                              # env.sh's default kernel, Haswell
    production) export OPENBLAS_CORETYPE=auto ;;
    *) echo "profile: portable or production, not '$profile'" >&2; exit 2 ;;
esac
source ./env.sh
export PYTHONNOUSERSITE=1
py="$root/.envs/loaders-py/bin/python -I"

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS" \
     "$(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
threads=${SLURM_CPUS_PER_TASK:-1}

OBJECTS=(soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test
         soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2)
i=${SLURM_ARRAY_TASK_ID:?submit as an array job: --array=1-10}
obj=${OBJECTS[$((i - 1))]:?no object $i}
tier=T
[ "$i" -gt 5 ] && tier=P

core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
if [ "$profile" = production ] && [ "$core" != SapphireRapids ]; then
    echo "$obj: the production profile needs the SapphireRapids kernel; this node gives '$core'" >&2
    exit 1
fi
base=$root/data/doubletfinder/$obj
inp=$base/input
# the saved objects' exports: the DoubletFinder output's re-clustering, and the input's clustering
saved_out=$root/data/clustering/doubletfinder_$obj/saved
saved_in=$root/data/clustering/$obj/saved
for f in "$inp/seurat_clusters.tsv" "$saved_out/metadata.tsv" "$saved_in/metadata.tsv"; do
    [ -f "$f" ] || { echo "$obj: no $f" >&2; exit 1; }
done
out=$base/rust/$label/${profile}_recluster
[ -e "$out" ] && { echo "skip $out: exists"; exit 0; }
tmp=$out.partial.${SLURM_JOB_ID:-local}
rm -rf "$tmp"
mkdir -p "$tmp"/{inprocess,inprocess_export,two,two_export/run,hook,alone,checks}
"$bin" fingerprint > "$tmp/fingerprint.json"
echo "== task $i: $obj (tier $tier) on $(hostname -s), $profile profile (OpenBLAS $core), $threads threads"

run() {
    local name=$1 dir=$2
    shift 2
    if ! (cd "$dir" && /usr/bin/time -v -o "$tmp/time_$name.txt" "$bin" "$@" 2> "$tmp/stderr_$name.txt"); then
        cat "$tmp/stderr_$name.txt" >&2
        echo "$obj: run $name failed; partial output left in $tmp" >&2
        exit 1
    fi
    tail -n 1 "$tmp/stderr_$name.txt"
}
given=(--counts "$inp" --clusters "$inp/seurat_clusters.tsv" --out run --threads "$threads")
run inprocess "$tmp/inprocess" doubletfinder "${given[@]}" --recluster
run inprocess_export "$tmp/inprocess_export" doubletfinder "${given[@]}" --recluster --recluster-export
run two_df "$tmp/two" doubletfinder "${given[@]}"
run two_cluster "$tmp/two" cluster --counts run/singlets --preset workflow --out run/singlets_clusters --threads "$threads"
# the export's second command reads the same singlets/ (a link), at the same relative path
ln -s ../../two/run/singlets "$tmp/two_export/run/singlets"
run two_cluster_export "$tmp/two_export" cluster --counts run/singlets --preset workflow --export --out run/singlets_clusters --threads "$threads"
run hook "$tmp/hook" doubletfinder --counts "$inp" --out run --cluster-stats silhouette,purity --threads "$threads"
run alone "$tmp/alone" cluster --counts "$inp" --preset workflow --stats silhouette,purity --out run/clusters --threads "$threads"

status=0 bstatus=0
cmp_runs() { $py "$root/validation/python/compare_runs.py" "$@"; }
saved_cmp() { $py "$root/validation/python/compare_saved_clustering.py" "$@"; }
check() {
    local name=$1 kind=$2
    shift 2
    echo "-- $name"
    if "$@"; then
        echo "$name: identical"
    else
        echo "$name: DIFFERS" >&2
        if [ "$kind" = a ]; then status=1; else bstatus=1; fi
    fi
}
check a_files a cmp_runs "$tmp/two/run" "$tmp/inprocess/run" --ignore summary.json:run --ignore summary.json:files \
    --ignore summary.json:result.recluster --json "$tmp/checks/a_files.json"
check a_export a cmp_runs "$tmp/two_export/run/singlets_clusters" "$tmp/inprocess_export/run/singlets_clusters" --json "$tmp/checks/a_export.json"
check hook_alone a cmp_runs "$tmp/alone/run/clusters" "$tmp/hook/run/clusters" --only-left pca --only-left neighbors --json "$tmp/checks/hook_alone.json"
manifests() { "$bin" validate manifest --dir "$tmp/inprocess/run" && "$bin" validate manifest --dir "$tmp/inprocess_export/run" && "$bin" validate manifest --dir "$tmp/hook/run"; }
check manifests a manifests
logged=()
[ "$tier" = P ] && logged=(--panel "$root/validation/clustering_panel.tsv" --object "doubletfinder_$obj")
check b_saved b saved_cmp "$tmp/inprocess_export/run/singlets_clusters" "$saved_out" "${logged[@]}" --json "$tmp/checks/b_saved.json"
hook_logged=()
[ "$tier" = P ] && hook_logged=(--panel "$root/validation/clustering_panel.tsv" --object "$obj")
check hook_saved b saved_cmp "$tmp/hook/run/clusters" "$saved_in" "${hook_logged[@]}" --json "$tmp/checks/hook_saved.json"
labels() { cmp "$tmp/hook/run/clusters/seurat_clusters.tsv" "$inp/seurat_clusters.tsv" && cmp "$tmp/hook/run/clusters/seurat_clusters_levels.txt" "$inp/seurat_clusters_levels.txt"; }
check hook_labels b labels

mv "$tmp" "$out"
echo "status: (a) $([ $status -eq 0 ] && echo identical || echo DIFFERS); against the saved objects $([ $bstatus -eq 0 ] && echo identical || echo "differs ($profile profile)")"
[ "$profile" = production ] && [ $bstatus -ne 0 ] && status=1
exit $status
