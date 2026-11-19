#!/bin/bash
# Benchmark (scRust data/benchmarks): calibration of one blind spot of the Rust write attribution.
# `scrust cluster` writes pca/ with pca::output::write_outputs, which formats each file into
# memory before opening it, so the sampler cannot see the formatting of the stage's first file
# (embeddings.tsv) as writing. `scrust pca` writes the same files with the same function and
# marks that call as its own "write" phase. This job runs `scrust pca` on each tier-P object's
# SCTransform output (M13's benchmark export, data/clustering/<object>/rust/m13-9ce0adb/bench/t1/sct,
# read only) three times, under the sampler, to measure the writer's CPU time per object.
#   sbatch pca_calib.sbatch
# Output: runs/calib/pca/<object>/rep<k>/ (procmon.*, summary.json, env.txt)
#SBATCH --job-name=bench-calib-pca
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a[32608,24102,24105-24112,24201-24212,26105-26112,26201-26212]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=33
#SBATCH --mem=32G
#SBATCH --time=0-01:00:00
#SBATCH --output=data/benchmarks/logs/%x-%j.log

set -uo pipefail
source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

# The partition mixes Xeon Platinum 8480+ (holy8a24xxx, 26xxx: excluded above) and 8480CL nodes
# (holy8a28xxx, 30xxx, 32xxx, where the workflow's jobs ran); every run of this benchmark is on
# an 8480CL, so R and Rust share the CPU model
model=$(lscpu | sed -n 's/^Model name: *//p')
case $model in *8480CL*) ;; *) echo "$(hostname -s): CPU model '$model', not an 8480CL; not run" >&2; exit 3 ;; esac
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
BIN=${SCRUST_BIN:-$B/bin/scrust-b6f6e4e}
PY="env PYTHONNOUSERSITE=1 /usr/bin/python3 -I"
export OPENBLAS_CORETYPE=auto
source "$S/env.sh" || exit 1
core=$("$BIN" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
echo "job=${SLURM_JOB_ID:-local} host=$(hostname -s) cpu=$(lscpu | sed -n 's/^Model name: *//p' | tr ' ' '_') openblas_core=$core"
for obj in soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain filtered_seurat_tenx_10xmousebrain soupx_emptydrops_4500871lane \
           soupx_emptydrops_sscovelli filtered_seurat_emptydrops_sscovelli cellbender_fromraw_ghirsutum2; do
    sct=$S/data/clustering/$obj/rust/m13-9ce0adb/bench/t1/sct
    for rep in 1 2 3; do
        d=$B/runs/calib/pca/$obj/rep$rep
        [ -e "$d/summary.json" ] && continue
        mkdir -p "$d"
        out=$d/out
        rm -rf "$out"
        echo "job=${SLURM_JOB_ID:-local} host=$(hostname -s) cpu=$(lscpu | sed -n 's/^Model name: *//p' | tr ' ' '_') openblas_core=$core" > "$d/env.txt"
        $PY "$S/validation/benchmark/procmon.py" --out "$d/procmon" --interval 0.01 --child-ncpus 32 --watch "$sct" --watch "$out" -- \
            "$BIN" pca --sct "$sct" --out "$out" --threads 32 2> "$d/stderr.txt"
        cp "$out/summary.json" "$d/summary.json"
        find "$out" -type f -printf '%s\t%P\n' | sort -k2 > "$d/outputs.tsv"
        cmp -s "$out/embeddings.tsv" "$S/data/clustering/$obj/rust/m13-9ce0adb/bench/t1/pca/embeddings.tsv" && echo "embeddings.tsv identical to the chain's" >> "$d/env.txt"
        gzip -f "$d/procmon.samples.tsv"
        rm -rf "$out"
        $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); p=s["run"]["phases"]; e=0; out=[]
for k,v in p.items(): out.append("%s %.3f s" % (k, v["elapsed_s"]-e)); e=v["elapsed_s"]
print(sys.argv[2], "; ".join(out))' "$d/summary.json" "$obj rep$rep"
    done
done
