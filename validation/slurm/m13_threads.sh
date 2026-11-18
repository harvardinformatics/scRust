#!/bin/bash
# Milestone M13 (docs/PLAN.md C5, Threads): M11's m11_threads.sh for runs with the export: two
# validation runs of the same objects and profile at different thread counts (m13_cluster.sbatch
# with two labels, e.g. 8 threads and 1), every file compared byte for byte (the stage dumps and
# the export included: counts/, sct/ with the corrected counts and their log1p, README.txt;
# summary.json: its result section and its manifest's entries, every file's size and SHA-256). From
# the repository's top directory:
#   bash validation/slurm/m13_threads.sh <label a> <label b> <profile> <object> [<object> ...]
# Prints one line per object, and exits 1 if any file differs.
set -euo pipefail
a=$1 b=$2 profile=$3
shift 3
result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }
entries_section() { python3 -I -c 'import json, sys; print(json.dumps(json.load(open(sys.argv[1]))["files"]["entries"], indent=1, sort_keys=True))' "$1"; }
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
    cmp -s <(entries_section "$da/summary.json") <(entries_section "$db/summary.json") || { echo "  $o: the manifest's entries differ" >&2; differ=$((differ + 1)); }
    ta=$(grep -m1 '"threads"' "$da/summary.json" | tr -dc 0-9)
    tb=$(grep -m1 '"threads"' "$db/summary.json" | tr -dc 0-9)
    echo "$o $profile: threads $ta and $tb: $n files, $differ differing (summary.json: its result section and its manifest's entries)"
    [ "$differ" = 0 ] || status=1
done
exit $status
