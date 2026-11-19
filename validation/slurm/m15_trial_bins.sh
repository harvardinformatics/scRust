#!/bin/bash
# Milestone M15: two trial builds of a commit, for measuring the allocator's settling at in-process
# hand-offs (scio::memory::settle; docs/validation/clustering.md, "Memory and time"). Never
# shipped, never committed: each is the commit's source with one call removed or added.
#   <bin>-nosettle   without the settle at clustering on the fly's hand-off to DoubletFinder
#   <bin>-settle2    with a settle also between DoubletFinder and the singlets' re-clustering
# From the repository's top directory, after source env.sh:
#   bash validation/slurm/m15_trial_bins.sh <commit>
# writes data/doubletfinder/logs/bin/scrust-m15-<commit>-nosettle and -settle2 (the commit's
# source is unpacked and built under data/doubletfinder/logs/bin/trial-<commit>/, removed after).
# validation/slurm/m15_memory.sbatch runs them as NOSETTLE_BIN and SETTLE2_BIN.
set -euo pipefail
commit=$(git rev-parse --short "${1:?commit}")
root=$PWD
work=$root/data/doubletfinder/logs/bin/trial-$commit
rm -rf "$work"
mkdir -p "$work"
for variant in nosettle settle2; do
    src=$work/$variant
    mkdir -p "$src"
    git archive "$commit" | tar -x -C "$src"
    PYTHONNOUSERSITE=1 python3 -I - "$src/crates/scrust/src/cmd_doubletfinder.rs" "$variant" <<'EOF'
import sys
path, variant = sys.argv[1], sys.argv[2]
s = open(path).read()
hook = "        scio::memory::settle();\n        phases.borrow_mut().mark(\"clustering: settle\");\n"
recluster = "        // DoubletFinder's own stages are gone; only its result's tables remain\n        drop(output);\n"
if s.count(hook) != 1 or s.count(recluster) != 1:
    sys.exit(f"{path}: the call sites are not where this script expects them")
if variant == "nosettle":
    s = s.replace(hook, "        phases.borrow_mut().mark(\"clustering: settle\");\n")
else:
    s = s.replace(recluster, recluster + "        scio::memory::settle();\n        phases.mark(\"recluster: settle\");\n")
open(path, "w").write(s)
EOF
    (cd "$src" && CARGO_TARGET_DIR=$work/target cargo build --release -p scrust --quiet)
    cp "$work/target/release/scrust" "$root/data/doubletfinder/logs/bin/scrust-m15-$commit-$variant"
    echo "scrust-m15-$commit-$variant: $(sha256sum "$root/data/doubletfinder/logs/bin/scrust-m15-$commit-$variant" | cut -c1-16)"
done
rm -rf "$work"
