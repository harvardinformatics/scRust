#!/bin/bash
# Milestone M15: peak memory and CPU time at the in-process hand-offs of scrust doubletfinder, and
# what settling the allocator there does (scio::memory::settle: glibc's mmap threshold fixed at
# 128 KiB and the free pages returned), on the objects named, in the production profile on
# Sapphire Rapids, 8 threads as M7. The modes (MODES):
#   given              scrust doubletfinder --clusters: M7's command
#   hook               clustering on the fly, settled at the hand-off to DoubletFinder
#   hook_nosettle      the same, not settled: NOSETTLE_BIN, the commit without that call
#                      (validation/slurm/m15_trial_bins.sh)
#   recluster2         scrust cluster --counts on given's singlets/: the two-command route's second
#   inprocess          the whole rule in one process (--recluster)
#   inprocess_settled  the same, also settled between DoubletFinder and the re-clustering:
#                      SETTLE2_BIN, the commit with that call added (m15_trial_bins.sh)
# Default "given hook hook_nosettle"; for the re-clustering, "given recluster2 inprocess
# inprocess_settled". Each mode runs <repeats> times, interleaved; /usr/bin/time -v records the
# peak, wall, user and system time; each run's summary.json records the peak at the end of each
# phase, so where the peak was reached. Output: data/doubletfinder/<object>/rust/<label>/memory/
# <mode>_<k>/ (the run, but singlets/) and memory.tsv (object, mode, repeat, peak RSS KiB, wall s,
# node, user s, system s).
#   SCRUST_BIN=<copy> [NOSETTLE_BIN=<trial>] [SETTLE2_BIN=<trial>] [MODES=...] \
#       sbatch validation/slurm/m15_memory.sbatch <label> <repeats> <object> ...
# (Jobs 61960681 and 61965453 ran an earlier form, with scrust-m15-49e8be7, in which
# SCRUST_NO_SETTLE=1 skipped the settle: memory.tsv then had six columns.)
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
modes=${MODES:-given hook hook_nosettle}
# the binary of each mode: the trial builds for the trial modes
bin_of() {
    case $1 in
        hook_nosettle) realpath "${NOSETTLE_BIN:?hook_nosettle needs NOSETTLE_BIN (m15_trial_bins.sh)}" ;;
        inprocess_settled) realpath "${SETTLE2_BIN:?inprocess_settled needs SETTLE2_BIN (m15_trial_bins.sh)}" ;;
        *) echo "$bin" ;;
    esac
}
for mode in $modes; do
    b=$(bin_of "$mode")
    [ "$b" = "$bin" ] || echo "mode $mode: binary=$b sha256=$(sha256sum "$b" | cut -c1-16)"
done
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
    [ -f "$out/memory.tsv" ] || printf 'object\tmode\trepeat\tpeak_rss_kib\twall_s\tnode\tuser_s\tsystem_s\n' > "$out/memory.tsv"
    for k in $(seq 1 "$reps"); do
        for mode in $modes; do
            d=$out/${mode}_$k
            rm -rf "$d"
            mkdir -p "$d"
            args=(doubletfinder --counts "$inp" --out "$d/run" --threads "$threads")
            case $mode in
                given) args+=(--clusters "$inp/seurat_clusters.tsv") ;;
                recluster2) args=(cluster --counts "$out/given_$k/run/singlets" --preset workflow --out "$d/run" --threads "$threads") ;;
                inprocess | inprocess_settled) args+=(--clusters "$inp/seurat_clusters.tsv" --recluster) ;;
            esac
            /usr/bin/time -v -o "$d/time_v.txt" "$(bin_of "$mode")" "${args[@]}" 2> "$d/stderr.txt"
            peak=$(sed -n 's/^\s*Maximum resident set size (kbytes): //p' "$d/time_v.txt")
            wall=$(sed -n 's/^\s*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$d/time_v.txt" | awk -F: '{ s = 0; for (i = 1; i <= NF; i++) s = s * 60 + $i; print s }')
            user=$(sed -n 's/^\s*User time (seconds): //p' "$d/time_v.txt")
            sys=$(sed -n 's/^\s*System time (seconds): //p' "$d/time_v.txt")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$obj" "$mode" "$k" "$peak" "$wall" "$(hostname -s)" "$user" "$sys" | tee -a "$out/memory.tsv"
            # singlets/, identical to the panel's, once nothing reads it any more
            [ "$mode" = given ] && [[ " $modes " == *" recluster2 "* ]] || rm -rf "$d/run/singlets"
        done
        rm -rf "$out/given_$k/run/singlets"
    done
done
