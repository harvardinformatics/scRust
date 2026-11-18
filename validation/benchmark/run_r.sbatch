#!/bin/bash
# Benchmark (scRust data/benchmarks): one run of a preprocessing-workflow R script, as the
# workflow's job runs it, with phase markers (scripts/workflow/scripts/, made by
# make_workflow_copies.py: the workflow's scripts plus marker lines, nothing else changed).
#
#   sbatch --mem=<mem> --time=<time> --job-name=bench-r-<tool> run_r.sbatch <tool> <object> <rep>
#   (submit_r.sh takes <mem> and <time> from r_runs.tsv)
#
# As the workflow's job (profiles/slurm: cpus_per_task 1; Snakemake's shell() sets the six
# *_NUM_THREADS variables to the rule's threads, 1; OPENBLAS_CORETYPE unset; the production R
# env d3ca7f730a70e6f748bb65a12a5ed3f9_; SCRNASEQ_PREPROCESS_SEED=12345; cwd holding
# workflow/scripts/): R is pinned to one CPU. The job has a second CPU for the sampler
# (procmon.py: /proc of the R process every 50 ms), so sampling takes nothing from R's CPU.
# The inputs are evicted from this node's page cache first, so R reads them from Lustre, as a
# workflow job on a fresh node does.
#
# Output: runs/r/<tool>/<object>/rep<rep>/ (a .partial.<jobid> directory until it finishes):
#   markers.tsv          the R process's phase marks (bench_markers.R)
#   procmon.*            the sampler's record (samples gzipped)
#   env.txt, outputs.tsv node, CPU, OpenBLAS kernel; the outputs' sizes (the large ones are then
#                        deleted: only the small files are kept)
#   check.txt            nclusters / cluster ids against the workflow's saved run (panel log)
#   stderr.txt           R's messages (the workflow's log)
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a[32608,24102,24105-24112,24201-24212,26105-26112,26201-26212]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=data/benchmarks/logs/%x-%j.log

set -uo pipefail
source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

# The partition mixes Xeon Platinum 8480+ (holy8a24xxx, 26xxx: excluded above) and 8480CL nodes
# (holy8a28xxx, 30xxx, 32xxx, where the workflow's jobs ran); every run of this benchmark is on
# an 8480CL, so R and Rust share the CPU model
model=$(lscpu | sed -n 's/^Model name: *//p')
case $model in *8480CL*) ;; *) echo "$(hostname -s): CPU model '$model', not an 8480CL; not run" >&2; exit 3 ;; esac

[ $# -eq 3 ] || { echo "usage: run_r.sbatch <tool> <object> <rep>" >&2; exit 2; }
tool=$1 obj=$2 rep=$3
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
REnv=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
SCRUST=$B/bin/scrust-b6f6e4e
PY="env PYTHONNOUSERSITE=1 /usr/bin/python3 -I"

line=$(awk -F'\t' -v t="$tool" -v o="$obj" '$1 == t && $2 == o' "$S/validation/benchmark/r_runs.tsv")
[ -n "$line" ] || { echo "no run $tool $obj in r_runs.tsv" >&2; exit 2; }
IFS=$'\t' read -r _t _o rule script mem tlimit inputs args <<< "$line"

final=$B/runs/r/$tool/$obj/rep$rep
[ -e "$final" ] && { echo "$final exists" >&2; exit 0; }
run=$final.partial.${SLURM_JOB_ID:-local}
rm -rf "$run"
mkdir -p "$run/out" "$run/work"
ln -s "$S/validation/benchmark/workflow" "$run/work/workflow"
args=${args//\{out\}/$run/out}

# the OpenBLAS kernel these nodes give the workflow's libopenblas with OPENBLAS_CORETYPE unset
# (scRust's rust-dev env has the build-identical OpenBLAS 0.3.34; `scrust fingerprint` reports it)
unset OPENBLAS_CORETYPE
core=$(env -u OPENBLAS_CORETYPE OPENBLAS_NUM_THREADS=1 "$SCRUST" fingerprint 2>/dev/null | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
{
    echo "job=${SLURM_JOB_ID:-local} host=$(hostname -s) cpu=$(lscpu | sed -n 's/^Model name: *//p' | tr ' ' '_')"
    echo "kernel=$(uname -r) cpus=$(taskset -cp $$ | sed 's/.*: //') mem_request=${SLURM_MEM_PER_NODE:-?}M openblas_core=$core"
    echo "tool=$tool object=$obj rule=$rule script=$script rep=$rep"
    echo "R=$REnv/bin/Rscript $("$REnv/bin/R" --version | head -n 1)"
    echo "copies: $(sed -n 1p "$S/validation/benchmark/workflow/SOURCES.tsv")"
    echo "start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$run/env.txt"
cat "$run/env.txt"
[ "$core" = SapphireRapids ] || echo "WARNING: OpenBLAS core '$core', not SapphireRapids" | tee -a "$run/env.txt"

watch=()
evict=()
IFS=';' read -r -a ins <<< "$inputs"
for p in "${ins[@]}"; do watch+=(--watch "$p"); evict+=(--evict "$p"); done

cd "$run/work"
# shellcheck disable=SC2086
$PY "$S/validation/benchmark/procmon.py" --out "$run/procmon" --interval 0.05 --child-ncpus 1 \
    "${watch[@]}" --watch "$run/out" --watch "$run/work" "${evict[@]}" -- \
    env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB -u OPENBLAS_CORETYPE \
        PYTHONNOUSERSITE=1 PATH="$REnv/bin:$PATH" \
        OMP_NUM_THREADS=1 GOTO_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
        SCRNASEQ_PREPROCESS_SEED=12345 BENCH_MARKERS="$run/markers.tsv" \
        "$REnv/bin/Rscript" "workflow/scripts/$script" $args \
    > "$run/stderr.txt" 2>&1
status=$?
cd "$run"
echo "end=$(date -u +%Y-%m-%dT%H:%M:%SZ) exit=$status" >> "$run/env.txt"
tail -n 5 "$run/stderr.txt"
echo "exit=$status"

# the outputs: sizes (and Rplots.pdf, which SoupX's autoEstCont draws into the cwd)
find "$run/out" "$run/work" -type f -printf '%s\t%P\t%p\n' | sort -k3 > "$run/outputs.tsv"
cat "$run/outputs.tsv"

# check: the number of clusters and their ids against the workflow's run (clustering_panel.tsv
# log_final_clusters for the clustering rules; the saved DoubletFinder log for doubletfinder)
{
    ncl=$(cat "$run"/out/*_nclusters.txt 2>/dev/null)
    echo "nclusters=$ncl ids=$(paste -sd, "$run"/out/*_cluster_ids.txt 2>/dev/null)"
    if [ "$tool" = cluster ]; then
        want=$(awk -F'\t' -v o="$obj" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == o {print $c["log_final_clusters"]}' \
            $S/validation/clustering_panel.tsv)
        echo "workflow's final clusters (clustering_panel.tsv): $want"
        [ "$ncl" = "$want" ] && echo "nclusters: same" || echo "nclusters: DIFFERENT"
    fi
} > "$run/check.txt" 2>&1
cat "$run/check.txt"

# keep the small files; delete the large outputs (.rds, the emptyDrops matrix) after recording sizes
find "$run/out" -type f -size +1M -delete
find "$run/out" -depth -type d -empty -delete 2>/dev/null
gzip -f "$run/procmon.samples.tsv"
rm -f "$run/work/workflow"
cd "$B"
if [ "$status" = 0 ]; then
    mv "$run" "$final"
else
    echo "FAILED: kept in $run" >&2
fi
exit $status
