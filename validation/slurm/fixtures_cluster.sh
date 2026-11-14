#!/bin/bash
# Copies M10's committed edge-case fixture (docs/PLAN.md C7) from data/clustering/edge/, made by
# validation/R/cluster_references.R --mode edge, into tests/fixtures/clustering/edge/:
#   cases.tsv           each case and setting: size, FindClusters' arguments, R's error (if R
#                       stopped), what the case exercises
#   <case>/snn.mtx.gz, <case>/barcodes.tsv   the graph, byte for byte
#   expected.tsv        the bytes and sha256 of every file R wrote for each setting (the labels,
#                       levels, the workflow's files, the K7 stage arrays, the RNG checkpoints, and
#                       the instrumented optimizer's starts, partitions and ids), except the
#                       summary.json files, whose run sections hold timings
#   printed.tsv         what FindClusters printed (but the elapsed time), its messages, R's error,
#                       and GroupSingletons' singletons, per setting and call
# crates/cluster/tests/edge_cases.rs recomputes every file from Rust's result and compares the
# digests; R's files themselves stay in data/clustering/edge/ (validation/manifest.toml). Plain
# files only, so a login node is enough. From the repository's top directory:
#   bash validation/slurm/fixtures_cluster.sh
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/../.."
src=data/clustering/edge
dst=tests/fixtures/clustering/edge
[ -f "$src/cases.tsv" ] || { echo "no $src/cases.tsv: run cluster_references.sbatch edge first" >&2; exit 1; }
rm -rf "$dst"
mkdir -p "$dst"
cp --preserve=mode "$src/cases.tsv" "$dst/"
printf 'case\tsetting\tfile\tbytes\tsha256\n' > "$dst/expected.tsv"
tail -n +2 "$src/cases.tsv" | cut -f1,2 | while read -r case setting; do
    if [ ! -d "$dst/$case" ]; then
        mkdir -p "$dst/$case"
        cp --preserve=mode "$src/$case/input/snn.mtx.gz" "$src/$case/input/barcodes.tsv" "$dst/$case/"
    fi
    (cd "$src/$case/$setting" && find . -type f ! -name summary.json ! -name optimizer.json ! -name stdout.txt ! -name stages.json -printf '%P\n' | LC_ALL=C sort) |
        while read -r f; do
            p=$src/$case/$setting/$f
            printf '%s\t%s\t%s\t%s\t%s\n' "$case" "$setting" "$f" "$(stat -c %s "$p")" "$(sha256sum "$p" | cut -d' ' -f1)"
        done >> "$dst/expected.tsv"
done
python3 -I - "$src" "$dst/printed.tsv" <<'EOF'
import csv, json, sys
src, out = sys.argv[1], sys.argv[2]
rows = [("case", "setting", "kind", "call", "text")]
with open(f"{src}/cases.tsv") as f:
    cases = list(csv.DictReader(f, delimiter="\t"))
for c in cases:
    d = f"{src}/{c['case']}/{c['setting']}"
    s = json.load(open(f"{d}/summary.json"))
    run = s["run"]
    for line in run["printed"]:
        if not line.startswith("Elapsed time"):
            rows.append((c["case"], c["setting"], "printed", "", line))
    for m in run["messages"]:
        rows.append((c["case"], c["setting"], "message", "", m.strip()))
    if run.get("error"):
        rows.append((c["case"], c["setting"], "error", "", run["error"]))
    for k, call in enumerate(s["result"]["calls"], start=1):
        rows.append((c["case"], c["setting"], "singletons", str(k), ",".join(call["singletons"])))
with open(out, "w") as f:
    for r in rows:
        if any("\t" in x or "\n" in x for x in r):
            sys.exit(f"tab or newline in {r}")
        f.write("\t".join(r) + "\n")
EOF
du -sb "$dst"
