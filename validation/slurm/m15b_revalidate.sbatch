#!/bin/bash
# Milestone M15b: validates scrust's validation runs again (validation/slurm/m15b_cluster.sbatch's
# run/ directories, label <runs label>) with another scrust binary's validator, without running
# them again: one array task per row of validation/m15b_cases.tsv, both profiles' runs, each against
# every R run of its case (validation/slurm/m15b_references.sbatch). For validator changes that
# change no output (the runs' binary and the validator's must give the same files: the task first
# checks that the two binaries' `scrust --version` agree and records both sha256).
#   SCRUST_BIN=<validator> sbatch --array=1-134%20 validation/slurm/m15b_revalidate.sbatch <runs label> <validator label>
# Output, never overwritten: data/clustering/m15b/revalidate/<validator label>/<object>/<case>/<profile>/
# validate_<R profile>.{json,md,txt}; validators.txt.
#
#SBATCH --job-name=scrust-m15b-reval
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%A_%a.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -eq 2 ] || { echo "usage: m15b_revalidate.sbatch <runs label> <validator label>" >&2; exit 2; }
runs=$1 vlabel=$2
source ./env.sh
bin=$(realpath "${SCRUST_BIN:-target/release/scrust}")
[ -x "$bin" ] || { echo "no $bin" >&2; exit 1; }
row=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n' validation/m15b_cases.tsv)
[ -n "$row" ] || { echo "no case $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
obj=$(cut -f1 <<< "$row") cs=$(cut -f6 <<< "$row")
status=0
for profile in production portable; do
    cand=data/clustering/m15b/rust/$runs/$obj/$cs/$profile/run
    [ -d "$cand" ] || { echo "no $cand" >&2; status=1; continue; }
    out=data/clustering/m15b/revalidate/$vlabel/$obj/$cs/$profile
    [ -e "$out" ] && { echo "skip $out: exists"; continue; }
    tmp=$out.partial.${SLURM_JOB_ID:-$$}
    rm -rf "$tmp"
    mkdir -p "$tmp"
    {
        echo "validator=$bin sha256=$(sha256sum "$bin" | cut -d' ' -f1) $("$bin" --version)"
        echo "runs=$cand $(sed -n 's/.*"scrust_version": *"\([^"]*\)".*/scrust \1/p' "$(dirname "$cand")/fingerprint.json" | head -n 1)"
    } > "$tmp/validators.txt"
    for rp in production portable; do
        ref=$(ls -d data/clustering/m15b/runs/$obj/$cs/direct-$rp-* 2>/dev/null | grep -v partial | head -n 1 || true)
        [ -n "$ref" ] || continue
        if "$bin" validate cluster --reference "$ref" --candidate "$cand" --json "$tmp/validate_$rp.json" \
                --markdown "$tmp/validate_$rp.md" > "$tmp/validate_$rp.txt" 2>&1; then
            echo "$obj $cs $profile against R ($rp): $(tail -n 1 "$tmp/validate_$rp.txt")"
        else
            echo "$obj $cs $profile against R ($rp): $(tail -n 1 "$tmp/validate_$rp.txt")" >&2
            status=1
        fi
    done
    mv "$tmp" "$out"
done
exit $status
