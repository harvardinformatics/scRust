#!/bin/bash
# Milestone M15b: without the new options, every output is M15's. The M15b binary's runs of
# validation/slurm/m15_recluster.sbatch (the five tier-T and five tier-P DoubletFinder objects, both
# profiles: DoubletFinder given the labels; the re-clustering in process with and without the
# export; the two-command route with and without it; clustering on the fly; scrust cluster
# --counts) against M15's final label's runs of the same wrapper (m15-a5de5e3), route by route, with
# validation/python/compare_runs.py: every file byte for byte, summary.json as JSON without the
# timings, the command line, and what records the node (run.fingerprint, run.slurm_job_id), and
# (--nested-summaries) without the digest and size its manifest records for a nested output's
# summary.json, which holds that run's timings (the nested summary.json itself is compared as JSON).
#   sbatch validation/slurm/m15b_unchanged.sbatch <new label> [<old label>]
# Output: data/clustering/m15b/unchanged/<new label>/<object>_<profile>_<route>.json and
# summary.tsv.
#
#SBATCH --job-name=scrust-m15b-same
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -uo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
new=${1:?usage: m15b_unchanged.sbatch <new label> [<old label>]}
old=${2:-m15-a5de5e3}
py="env PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I"
out=data/clustering/m15b/unchanged/$new
[ -e "$out" ] && { echo "$out exists; give another label or move it" >&2; exit 2; }
mkdir -p "$out"
printf 'object\tprofile\troute\tstatus\n' > "$out/summary.tsv"
fails=0
for obj in soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test \
           soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2; do
    for profile in production portable; do
        a=data/doubletfinder/$obj/rust/$old/${profile}_recluster
        b=data/doubletfinder/$obj/rust/$new/${profile}_recluster
        for route in inprocess inprocess_export two two_export hook alone; do
            if [ ! -d "$a/$route" ] || [ ! -d "$b/$route" ]; then
                st=missing
            elif $py validation/python/compare_runs.py "$a/$route" "$b/$route" --ignore '*:run.fingerprint' \
                    --ignore '*:run.slurm_job_id' --nested-summaries --json "$out/${obj}_${profile}_${route}.json" > "$out/${obj}_${profile}_${route}.txt" 2>&1; then
                st=identical
            else
                st=DIFFERS
            fi
            [ "$st" = identical ] || fails=$((fails + 1))
            printf '%s\t%s\t%s\t%s\n' "$obj" "$profile" "$route" "$st" | tee -a "$out/summary.tsv"
        done
    done
done
echo "routes not identical: $fails"
[ "$fails" = 0 ]
