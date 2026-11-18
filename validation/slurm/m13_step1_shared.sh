#!/bin/bash
# Milestone M13: Step 1 unchanged, the objects M6 and M7 share. m13_doubletfinder.sbatch runs each
# object and profile once per label: in its m6 set with M6's options (stage dumps), in its m7 set
# with M7's (none), into the same directory, so an object of both sets (five tier-P objects) is run
# with M6's options only and its m7 run is skipped. This compares those runs with M7's runs too:
# every file M7's run has (cells.tsv, sweep_stats.tsv, bcmvn.tsv, singlets/) byte for byte, and
# summary.json's result section (the dump level is in the run section). The extra files of the
# M6-option run (stages/) are not M7's to compare. Reads plain files only. From the repository's
# top directory:
#   bash validation/slurm/m13_step1_shared.sh <label> > data/doubletfinder/compare/m13_step1_shared_<label>.tsv
# One row per object; exit status 1 if any file differs or is missing.
set -euo pipefail
label=${1:?label}
old_label=m7-b8c8a8c
status=0
result_section() { awk '/^  "result": /{on=1} /^  "run": /{on=0} on' "$1"; }
printf 'object\tprofile\tcompared_with\tnew_dump_stages\tfiles\tidentical\tdiffering\tresult_section\n'
for old in data/doubletfinder/*/rust/$old_label/production/run; do
    obj=$(echo "$old" | cut -d/ -f3)
    new=data/doubletfinder/$obj/rust/$label/production/run
    [ -d "$new" ] || continue
    level=$(sed -n 's/.*"dump_stages": *"\([a-z]*\)".*/\1/p' "$new/summary.json" | head -n 1)
    [ -n "$level" ] && [ "$level" != none ] || continue  # run with M7's options: m13_doubletfinder.sbatch compared it
    n=0
    d=0
    while IFS= read -r f; do
        [ "$f" = summary.json ] && continue
        n=$((n + 1))
        if [ -f "$new/$f" ] && cmp -s "$new/$f" "$old/$f"; then :; else d=$((d + 1)); echo "$obj: $f differs or is missing" >&2; fi
    done < <(cd "$old" && find . -type f | sed 's|^\./||' | sort)
    rs=identical
    cmp -s <(result_section "$new/summary.json") <(result_section "$old/summary.json") || { rs=differs; status=1; }
    [ "$d" = 0 ] || status=1
    printf '%s\tproduction\t%s\t%s\t%s\t%s\t%s\t%s\n' "$obj" "$old_label" "$level" "$n" "$((n - d))" "$d" "$rs"
done
exit $status
