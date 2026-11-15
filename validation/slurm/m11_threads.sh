#!/bin/bash
# Milestone M11 (docs/PLAN.md C5, Threads): two validation runs of the same objects and profile at
# different thread counts (m11_cluster.sbatch with two labels, e.g. 8 threads and 1), every file
# compared byte for byte (the stage dumps included; summary.json: its result section). From the
# repository's top directory:
#   bash validation/slurm/m11_threads.sh <label a> <label b> <profile> <object> [<object> ...]
# Prints one line per object, and exits 1 if any file differs.
set -euo pipefail
a=$1 b=$2 profile=$3
shift 3
result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }
status=0
for o in "$@"; do
    da=data/clustering/$o/rust/$a/$profile/run
    db=data/clustering/$o/rust/$b/$profile/run
    differ=0
    n=0
    while IFS= read -r f; do
        n=$((n + 1))
        [ "$f" = summary.json ] && continue
        cmp -s "$da/$f" "$db/$f" || { echo "  $o: $f differs" >&2; differ=$((differ + 1)); }
    done < <(cd "$da" && find . -type f | sed 's|^\./||' | sort)
    [ "$n" = "$(cd "$db" && find . -type f | wc -l)" ] || { echo "  $o: different file lists" >&2; differ=$((differ + 1)); }
    cmp -s <(result_section "$da/summary.json") <(result_section "$db/summary.json") || { echo "  $o: result section differs" >&2; differ=$((differ + 1)); }
    ta=$(grep -m1 '"threads"' "$da/summary.json" | tr -dc 0-9)
    tb=$(grep -m1 '"threads"' "$db/summary.json" | tr -dc 0-9)
    echo "$o $profile: threads $ta and $tb: $n files, $differ differing (summary.json: its result section)"
    [ "$differ" = 0 ] || status=1
done
exit $status
