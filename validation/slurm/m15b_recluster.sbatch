#!/bin/bash
# Milestone M15b: SCTransform's vars.to.regress in scrust doubletfinder's two clusterings, in
# process against the commands run one after the other, on the five tier-T and five tier-P
# DoubletFinder objects (M15's, validation/slurm/m15_recluster.sbatch), in one profile
# (portable: OPENBLAS_CORETYPE=Haswell; production: OpenBLAS's own kernel on Sapphire Rapids).
# The options go in the parameter files, as a user gives them:
#   - clustering on the fly (--cluster-params hook.toml): vars_to_regress = percent.mt and
#     nCount_RNA; percent.mt by the workflow's pattern on the mouse objects, which have
#     mitochondrial genes named mt-, else by a list of ids (data/clustering/m15b/inputs/<object>/:
#     tier T's stand-in list, mt_genes.txt, or for the two objects whose pattern finds nothing a
#     list of features without counts, zero_genes.txt: percent.mt all 0, the QR's rank-deficient
#     path);
#   - the re-clustering (--recluster-params re.toml, --recluster-metadata metadata.tsv):
#     vars_to_regress = percent.mt and a --metadata column (score), percent.mt by the list (the
#     mouse objects' mt_genes.txt).
# Checks, each byte for byte (summary.json as JSON without the timings and the command line,
# validation/python/compare_runs.py):
#   a_files, a_recluster: scrust doubletfinder --recluster against scrust doubletfinder, then
#     scrust cluster --counts <out>/singlets --params re.toml --metadata metadata.tsv;
#   hook_alone: clustering on the fly's clusters/ against scrust cluster --counts <input>
#     --params hook.toml (its pca/ and neighbors/ aside);
#   df_unchanged: DoubletFinder's own files (all but singlets_clusters/ and summary.json, whose
#     result is compared without result.recluster) against M15's final label's in-process run of
#     the same profile (m15-a5de5e3): the re-clustering's regression leaves DoubletFinder as it was;
#   manifests: scrust validate manifest on each output.
# Each route runs in a working directory of its own, with the parameter files at the same relative
# path (../inputs/), so the paths they record are the same.
#
#   SCRUST_BIN=<copy> sbatch --array=1-10 validation/slurm/m15b_recluster.sbatch <label> portable|production
#
# Output, never overwritten: data/doubletfinder/<object>/rust/<label>/<profile>_m15b/: inputs/,
# the routes' directories (inprocess/, two/, hook/, alone/), checks/*.json, time_*.txt,
# stderr_*.txt, fingerprint.json. Exit status 1 unless every check holds.
#
#SBATCH --job-name=scrust-m15b-rc
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
[ $# -eq 2 ] || { echo "usage: m15b_recluster.sbatch <label> portable|production" >&2; exit 2; }
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

core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
if [ "$profile" = production ] && [ "$core" != SapphireRapids ]; then
    echo "$obj: the production profile needs the SapphireRapids kernel; this node gives '$core'" >&2
    exit 1
fi
if [ "$profile" = portable ] && [ "$core" != Haswell ]; then
    echo "$obj: the portable profile needs the Haswell kernel; this node gives '$core'" >&2
    exit 1
fi
inp=$root/data/doubletfinder/$obj/input
given_inputs=$root/data/clustering/m15b/inputs/$obj
m15=$root/data/doubletfinder/$obj/rust/m15-a5de5e3/${profile}_recluster/inprocess/run
for f in "$inp/seurat_clusters.tsv" "$given_inputs/metadata.tsv" "$m15/cells.tsv"; do
    [ -f "$f" ] || { echo "$obj: no $f" >&2; exit 1; }
done
out=$root/data/doubletfinder/$obj/rust/$label/${profile}_m15b
[ -e "$out" ] && { echo "skip $out: exists"; exit 0; }
tmp=$out.partial.${SLURM_JOB_ID:-local}
rm -rf "$tmp"
mkdir -p "$tmp"/{inputs,inprocess,two,hook,alone,checks}
"$bin" fingerprint > "$tmp/fingerprint.json"
echo "== task $i: $obj on $(hostname -s), $profile profile (OpenBLAS $core), $threads threads"

# the parameter files, next to the lists they name by a relative path
cp "$given_inputs/metadata.tsv" "$tmp/inputs/metadata.tsv"
hook_mt= re_mt=
if [ -f "$given_inputs/mt_genes.txt" ]; then
    cp "$given_inputs/mt_genes.txt" "$tmp/inputs/mt_genes.txt"
    re_mt='mt_genes = "mt_genes.txt"'
    case $obj in
        *_test) hook_mt='mt_genes = "mt_genes.txt"' ;;   # tier T: no mitochondrial gene, the stand-in list
        *) ;;                                            # the mouse objects: the workflow's pattern
    esac
else
    cp "$given_inputs/zero_genes.txt" "$tmp/inputs/zero_genes.txt"
    hook_mt='mt_genes = "zero_genes.txt"' re_mt='mt_genes = "zero_genes.txt"'
fi
printf '[sctransform]\nvars_to_regress = ["percent.mt", "nCount_RNA"]\n%s\n' "$hook_mt" > "$tmp/inputs/hook.toml"
printf '[sctransform]\nvars_to_regress = ["percent.mt", "score"]\n%s\n' "$re_mt" > "$tmp/inputs/re.toml"
head -n 20 "$tmp"/inputs/*.toml

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
run inprocess "$tmp/inprocess" doubletfinder "${given[@]}" --recluster --recluster-params ../inputs/re.toml \
    --recluster-metadata ../inputs/metadata.tsv
run two_df "$tmp/two" doubletfinder "${given[@]}"
run two_cluster "$tmp/two" cluster --counts run/singlets --preset workflow --params ../inputs/re.toml \
    --metadata ../inputs/metadata.tsv --out run/singlets_clusters --threads "$threads"
run hook "$tmp/hook" doubletfinder --counts "$inp" --out run --cluster-params ../inputs/hook.toml \
    --cluster-stats silhouette,purity --threads "$threads"
run alone "$tmp/alone" cluster --counts "$inp" --preset workflow --params ../inputs/hook.toml \
    --stats silhouette,purity --out run/clusters --threads "$threads"

status=0
cmp_runs() { $py "$root/validation/python/compare_runs.py" "$@"; }
check() {
    local name=$1
    shift
    echo "-- $name"
    if "$@"; then
        echo "$name: identical"
    else
        echo "$name: DIFFERS" >&2
        status=1
    fi
}
check a_files cmp_runs "$tmp/two/run" "$tmp/inprocess/run" --ignore summary.json:run --ignore summary.json:files \
    --ignore summary.json:result.recluster --json "$tmp/checks/a_files.json"
check a_recluster cmp_runs "$tmp/two/run/singlets_clusters" "$tmp/inprocess/run/singlets_clusters" --json "$tmp/checks/a_recluster.json"
check hook_alone cmp_runs "$tmp/alone/run/clusters" "$tmp/hook/run/clusters" --only-left pca --only-left neighbors \
    --json "$tmp/checks/hook_alone.json"
df_unchanged() {
    local f n=0
    while IFS= read -r f; do
        cmp -s "$m15/$f" "$tmp/inprocess/run/$f" || { echo "differs: $f"; return 1; }
        n=$((n + 1))
    done < <(cd "$m15" && find . -type f -not -path './singlets_clusters/*' -not -name summary.json | sort)
    [ "$(cd "$tmp/inprocess/run" && find . -type f -not -path './singlets_clusters/*' -not -name summary.json | wc -l)" = "$n" ] \
        || { echo "the file lists differ"; return 1; }
    $py - "$m15/summary.json" "$tmp/inprocess/run/summary.json" "$n" > "$tmp/checks/df_unchanged.txt" <<'EOF'
import json, sys
a, b = (json.load(open(p))["result"] for p in sys.argv[1:3])
for r in (a, b):
    r.pop("recluster", None)
print(f"{sys.argv[3]} files byte for byte; summary.json's result without result.recluster", "identical" if a == b else "DIFFERS")
sys.exit(0 if a == b else 1)
EOF
    local s=$?
    cat "$tmp/checks/df_unchanged.txt"
    return $s
}
check df_unchanged df_unchanged
manifests() {
    "$bin" validate manifest --dir "$tmp/inprocess/run" && "$bin" validate manifest --dir "$tmp/hook/run" \
        && "$bin" validate manifest --dir "$tmp/two/run/singlets_clusters" && "$bin" validate manifest --dir "$tmp/alone/run/clusters"
}
check manifests manifests
# what the runs regressed, as summary.json records it
$py - "$tmp" <<'EOF' | tee "$tmp/checks/regression.txt"
import json, sys, os
t = sys.argv[1]
for name, p in (("hook", "hook/run/clusters/summary.json"), ("recluster", "inprocess/run/singlets_clusters/summary.json")):
    s = json.load(open(os.path.join(t, p)))["result"]["sctransform"]
    r, pm = s.get("regression", {}), s.get("percent_mt", {})
    print(name, "regressed:", r.get("vars_to_regress"), "rank", r.get("rank"), "of", r.get("n_columns"),
          "| percent.mt:", pm.get("definition"), pm.get("n_features"), "feature(s), cells at 0:", pm.get("cells_zero"))
EOF

mv "$tmp" "$out"
echo "status: $([ $status -eq 0 ] && echo "every check holds" || echo DIFFERS)"
exit $status
