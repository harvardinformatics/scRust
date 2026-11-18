#!/bin/bash
# Benchmark (scRust data/benchmarks): the Rust ports' runs for one object, one repeat.
#
#   sbatch --job-name=bench-rust-<tool> run_rust.sbatch cluster <object> <rep>
#     scrust cluster --counts data/clustering/<object>/input --preset workflow [--export],
#     at 1 and 32 threads: variants n1, x1, n32, x32 (n: without --export, x: with it), in
#     that order
#   sbatch --job-name=bench-rust-<tool> run_rust.sbatch doubletfinder <object> <rep>
#     scrust doubletfinder --counts data/doubletfinder/<object>/input --clusters .../seurat_clusters.tsv,
#     then scrust cluster --counts <its singlets/> --preset workflow --export (the clustering of
#     the singlets that the workflow's doubletfinder.R runs after DoubletFinder), at 1 and 32
#     threads: variants d1 + r1, d32 + r32
#   SCRUST_BIN=<binary> sbatch --job-name=bench-rust-<tool> run_rust.sbatch doubletfinder_m15 <object> <rep>
#     milestone M15, with one binary: d and r as above; h: scrust doubletfinder without --clusters
#     (clustering on the fly); p: the whole rule in one process, scrust doubletfinder --clusters ...
#     --recluster --recluster-export (the singlets re-clustered from memory, with the export as r
#     has it), at 1 and 32 threads: variants d1, r1, h1, p1, d32, r32, h32, p32, in that order
#
# Submit from the repository's top directory (the paths below are the repository's). The binary is
# $SCRUST_BIN, by default main's code as verified for M13 (data/benchmarks/bin/scrust-b6f6e4e,
# sha256 6b373bdc...), which made the first tables; in the production profile (OPENBLAS_CORETYPE
# unset: env.sh with OPENBLAS_CORETYPE=auto; the OpenBLAS kernel must be SapphireRapids;
# OPENBLAS_NUM_THREADS=1, scrust splits its own work).
# The job has 33 CPUs: each run gets 32 (as M13's benchmark: --threads 1 or 32 on 32 CPUs) and
# the sampler (procmon.py, /proc every 10 ms) the 33rd. Before each run, its inputs are evicted
# from this node's page cache (procmon --evict), so every run reads them from Lustre.
#
# Output: runs/rust/<tool>/<object>/rep<rep>/<variant>/ (.partial.<jobid> until the job ends):
#   procmon.*, stderr.txt, summary.json (the run's own record: run.phases), outputs.tsv (every
#   output file's size; the large ones are then deleted), same.txt (outputs against the
#   reference runs, byte for byte), env.txt (node, CPU, kernel, binary, OpenBLAS kernel)
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a[32608,24102,24105-24112,24201-24212,26105-26112,26201-26212]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=33
#SBATCH --mem=48G
#SBATCH --time=0-01:30:00
#SBATCH --output=data/benchmarks/logs/%x-%j.log

set -uo pipefail
source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

# The partition mixes Xeon Platinum 8480+ (holy8a24xxx, 26xxx: excluded above) and 8480CL nodes
# (holy8a28xxx, 30xxx, 32xxx, where the workflow's jobs ran); every run of this benchmark is on
# an 8480CL, so R and Rust share the CPU model
model=$(lscpu | sed -n 's/^Model name: *//p')
case $model in *8480CL*) ;; *) echo "$(hostname -s): CPU model '$model', not an 8480CL; not run" >&2; exit 3 ;; esac

[ $# -eq 3 ] || { echo "usage: run_rust.sbatch cluster|doubletfinder|doubletfinder_m15 <object> <rep>" >&2; exit 2; }
tool=$1 obj=$2 rep=$3
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
BIN=${SCRUST_BIN:-$B/bin/scrust-b6f6e4e}
PY="env PYTHONNOUSERSITE=1 /usr/bin/python3 -I"

final=$B/runs/rust/$tool/$obj/rep$rep
[ -e "$final" ] && { echo "$final exists" >&2; exit 0; }
run=$final.partial.${SLURM_JOB_ID:-local}
rm -rf "$run"
mkdir -p "$run"

# production profile (env.sh: OPENBLAS_CORETYPE=auto -> unset; OPENBLAS_NUM_THREADS=1)
export OPENBLAS_CORETYPE=auto
source "$S/env.sh" || { echo "env.sh failed" >&2; exit 1; }
core=$("$BIN" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
{
    echo "job=${SLURM_JOB_ID:-local} host=$(hostname -s) cpu=$(lscpu | sed -n 's/^Model name: *//p' | tr ' ' '_')"
    echo "kernel=$(uname -r) cpus=$(taskset -cp $$ | sed 's/.*: //') openblas_core=$core OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-unset} OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    echo "binary=$BIN sha256=$(sha256sum "$BIN" | cut -c1-16) $("$BIN" --version) libopenblas=$(ldd "$BIN" | awk '/libopenblas/ {print $3}')"
    echo "tool=$tool object=$obj rep=$rep start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$run/env.txt"
cat "$run/env.txt"
[ "$core" = SapphireRapids ] || { echo "the production profile needs SapphireRapids, this node gives '$core'" >&2; exit 1; }

status=0
# one measured run: <variant> <input dir to watch and evict> <out dir> -- <scrust args>
measure() {
    local v=$1 input=$2 out=$3
    shift 4
    local d=$run/$v
    mkdir -p "$d"
    cp "$run/env.txt" "$d/env.txt"
    echo "== $v: scrust $*"
    $PY "$S/validation/benchmark/procmon.py" --out "$d/procmon" --interval 0.01 --child-ncpus 32 \
        --watch "$input" --watch "$out" --evict "$input" -- "$BIN" "$@" > "$d/stdout.txt" 2> "$d/stderr.txt"
    local st=$?
    echo "exit=$st" >> "$d/env.txt"
    tail -n 2 "$d/stderr.txt"
    [ -f "$out/summary.json" ] && cp "$out/summary.json" "$d/summary.json"
    find "$out" -type f -printf '%s\t%P\n' | sort -k2 > "$d/outputs.tsv"
    gzip -f "$d/procmon.samples.tsv"
    $PY -c 'import json,sys; r=json.load(open(sys.argv[1])); print("  wall %.2f s, cpu %.2f s, maxrss %d KiB, samples %d, sampler cpu %.2f s" % (r["wall_s"], r["rusage"]["cpu_s"], r["rusage"]["maxrss_kib"], r["samples"], r["sampler_cpu_s"]))' "$d/procmon.json"
    [ "$st" = 0 ] || status=1
    return $st
}

# outputs against a reference run's, byte for byte: <variant> <out> <ref dir> <file> ...
same() {
    local v=$1 out=$2 ref=$3
    shift 3
    local f
    for f in "$@"; do
        if [ ! -e "$ref/$f" ]; then
            echo "$f: no reference" >> "$run/$v/same.txt"
        elif cmp -s "$out/$f" "$ref/$f"; then
            echo "$f: identical" >> "$run/$v/same.txt"
        elif [[ $f == sct/*.mtx.gz ]] && cmp -s <(zcat "$out/$f" | grep -v '^%') <(zcat "$ref/$f" | grep -v '^%'); then
            # the export's matrices name the command and the input's path in their comment lines
            echo "$f: identical data lines (comment lines name the command's paths)" >> "$run/$v/same.txt"
        else
            echo "$f: DIFFERENT" >> "$run/$v/same.txt"
            status=1
        fi
    done
    echo "  vs $ref: $(grep -c identical "$run/$v/same.txt") identical, $(grep -c DIFFERENT "$run/$v/same.txt") different"
    return 0
}

# keep summary.json (already copied) and small files; delete the rest of an output dir
prune() { rm -rf "$1"; }

scratch=$run/scratch
mkdir -p "$scratch"
case $tool in
cluster)
    input=$S/data/clustering/$obj/input
    ref=$S/data/clustering/$obj/rust/m13-9ce0adb/bench/t1
    for v in n1 x1 n32 x32; do
        t=${v:1}
        out=$scratch/$v
        extra=()
        [ "${v:0:1}" = x ] && extra=(--export)
        measure "$v" "$input" "$out" -- cluster --counts "$input" --preset workflow "${extra[@]}" --out "$out" --threads "$t"
        files=(seurat_clusters.tsv seurat_clusters_levels.txt clusters.tsv optimizer_ids.tsv nclusters.txt cluster_ids.txt pca/embeddings.tsv neighbors/knn_index.tsv neighbors/snn.mtx.gz)
        [ "${v:0:1}" = x ] && files+=(metadata.tsv sct/corrected_counts.mtx.gz sct/scale_data.npy counts/matrix.mtx.gz)
        same "$v" "$out" "$ref" "${files[@]}"
        prune "$out"
    done
    ;;
doubletfinder)
    input=$S/data/doubletfinder/$obj/input
    ref=$S/data/doubletfinder/$obj/runs/direct-production-sapphirerapids
    for t in 1 32; do
        out=$scratch/d$t
        measure "d$t" "$input" "$out" -- doubletfinder --counts "$input" --clusters "$input/seurat_clusters.tsv" --out "$out" --threads "$t"
        same "d$t" "$out" "$ref" cells.tsv sweep_stats.tsv bcmvn.tsv
        # the clustering of the singlets (doubletfinder.R after DoubletFinder), from singlets/
        out2=$scratch/r$t
        measure "r$t" "$out/singlets" "$out2" -- cluster --counts "$out/singlets" --preset workflow --export --out "$out2" --threads "$t"
        cp "$out2/nclusters.txt" "$run/r$t/nclusters.txt" 2>/dev/null
        cp "$out2/cluster_ids.txt" "$run/r$t/cluster_ids.txt" 2>/dev/null
        if [ "$t" = 32 ]; then
            # the 32-thread runs against the 1-thread runs (kept until here)
            same d32 "$out" "$scratch/d1" cells.tsv singlets/matrix.mtx.gz
            same r32 "$out2" "$scratch/r1" seurat_clusters.tsv metadata.tsv sct/corrected_counts.mtx.gz
            prune "$scratch/d1"; prune "$scratch/r1"; prune "$out"; prune "$out2"
        fi
    done
    ;;
doubletfinder_m15)
    input=$S/data/doubletfinder/$obj/input
    ref=$S/data/doubletfinder/$obj/runs/direct-production-sapphirerapids
    for t in 1 32; do
        out=$scratch/d$t
        measure "d$t" "$input" "$out" -- doubletfinder --counts "$input" --clusters "$input/seurat_clusters.tsv" --out "$out" --threads "$t"
        same "d$t" "$out" "$ref" cells.tsv sweep_stats.tsv bcmvn.tsv
        out2=$scratch/r$t
        measure "r$t" "$out/singlets" "$out2" -- cluster --counts "$out/singlets" --preset workflow --export --out "$out2" --threads "$t"
        cp "$out2/nclusters.txt" "$run/r$t/nclusters.txt" 2>/dev/null
        cp "$out2/cluster_ids.txt" "$run/r$t/cluster_ids.txt" 2>/dev/null
        # clustering on the fly: R's labels and levels, and so R's DoubletFinder files
        out3=$scratch/h$t
        measure "h$t" "$input" "$out3" -- doubletfinder --counts "$input" --out "$out3" --threads "$t"
        same "h$t" "$out3" "$ref" cells.tsv sweep_stats.tsv bcmvn.tsv
        same "h$t" "$out3/clusters" "$input" seurat_clusters.tsv seurat_clusters_levels.txt
        # the whole rule in one process: the two-command route's files
        out4=$scratch/p$t
        measure "p$t" "$input" "$out4" -- doubletfinder --counts "$input" --clusters "$input/seurat_clusters.tsv" --out "$out4" \
            --recluster --recluster-export --threads "$t"
        same "p$t" "$out4" "$out" cells.tsv sweep_stats.tsv bcmvn.tsv singlets/matrix.mtx.gz
        same "p$t" "$out4/recluster" "$out2" seurat_clusters.tsv metadata.tsv pca/embeddings.tsv neighbors/snn.mtx.gz \
            sct/corrected_counts.mtx.gz sct/scale_data.npy counts/matrix.mtx.gz
        prune "$out"; prune "$out2"; prune "$out3"; prune "$out4"
    done
    ;;
*) echo "unknown tool $tool" >&2; exit 2 ;;
esac
rm -rf "$scratch"
echo "end=$(date -u +%Y-%m-%dT%H:%M:%SZ) status=$status" >> "$run/env.txt"
if [ "$status" = 0 ]; then
    mv "$run" "$final"
else
    echo "FAILED or differing outputs: kept in $run" >&2
fi
exit $status
