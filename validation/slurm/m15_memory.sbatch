#!/bin/bash
# Milestone M15: the peak memory of clustering on the fly, with and without settling the allocator
# before DoubletFinder (cli::settle_memory: glibc's mmap threshold fixed at 128 KiB and the free
# pages returned; SCRUST_NO_SETTLE=1 skips it), against the same binary given the labels (M7's
# command), on the objects named, in the production profile on Sapphire Rapids, 8 threads as M7.
# Each mode runs <repeats> times; /usr/bin/time -v records the peak. Output:
# data/doubletfinder/<object>/rust/<label>/memory/<mode>_<k>/ (the run) and memory.tsv (object,
# mode, repeat, peak RSS KiB, wall s).
#   SCRUST_BIN=<copy> sbatch validation/slurm/m15_memory.sbatch <label> <repeats> <object> ...
#
#SBATCH --job-name=scrust-m15-mem
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
root=$PWD
export OPENBLAS_CORETYPE=auto
source ./env.sh
bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16)"
core=$("$bin" fingerprint | sed -n 's/.*"openblas_corename": *"\([^"]*\)".*/\1/p' | head -n 1)
[ "$core" = SapphireRapids ] || { echo "the production profile needs SapphireRapids, not '$core'" >&2; exit 1; }
[ $# -ge 3 ] || { echo "usage: m15_memory.sbatch <label> <repeats> <object> ..." >&2; exit 2; }
label=$1 reps=$2
shift 2
threads=${SLURM_CPUS_PER_TASK:-1}
for obj in "$@"; do
    inp=$root/data/doubletfinder/$obj/input
    out=$root/data/doubletfinder/$obj/rust/$label/memory
    mkdir -p "$out"
    [ -f "$out/memory.tsv" ] || printf 'object\tmode\trepeat\tpeak_rss_kib\twall_s\tnode\n' > "$out/memory.tsv"
    for k in $(seq 1 "$reps"); do
        for mode in given hook hook_nosettle; do
            d=$out/${mode}_$k
            rm -rf "$d"
            mkdir -p "$d"
            args=(doubletfinder --counts "$inp" --out "$d/run" --threads "$threads")
            [ "$mode" = given ] && args+=(--clusters "$inp/seurat_clusters.tsv")
            if [ "$mode" = hook_nosettle ]; then
                SCRUST_NO_SETTLE=1 /usr/bin/time -v -o "$d/time_v.txt" "$bin" "${args[@]}" 2> "$d/stderr.txt"
            else
                /usr/bin/time -v -o "$d/time_v.txt" "$bin" "${args[@]}" 2> "$d/stderr.txt"
            fi
            peak=$(sed -n 's/^\s*Maximum resident set size (kbytes): //p' "$d/time_v.txt")
            wall=$(sed -n 's/^\s*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$d/time_v.txt" | awk -F: '{ s = 0; for (i = 1; i <= NF; i++) s = s * 60 + $i; print s }')
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$obj" "$mode" "$k" "$peak" "$wall" "$(hostname -s)" | tee -a "$out/memory.tsv"
            rm -rf "$d/run/singlets"   # identical to the panel's; the rest is kept
        done
    done
done
