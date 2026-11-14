#!/bin/bash
# Copies M9's committed edge-case fixture (docs/PLAN.md C7) from data/neighbors/edge/, made by
# validation/R/neighbors_references.R --mode edge, into tests/fixtures/neighbors/edge/: the cases
# (cases.tsv: size, parameters, the input's sha256, what each case exercises), every distinct input
# once, byte for byte, named by its sha256 (inputs/<sha256>.npy: several cases share one input),
# and the bytes and sha256 of every file R wrote (expected.tsv). crates/neighbors/tests/edge_cases.rs compares Rust's files with those digests;
# R's files themselves stay in data/neighbors/edge/ (validation/manifest.toml). Plain files only,
# so a login node is enough. From the repository's top directory:
#   bash validation/slurm/fixtures_neighbors.sh
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/../.."
src=data/neighbors/edge
dst=tests/fixtures/neighbors/edge
[ -f "$src/expected.tsv" ] || { echo "no $src/expected.tsv: run neighbors_references.sbatch edge first" >&2; exit 1; }
mkdir -p "$dst/inputs"
cp --preserve=mode "$src/cases.tsv" "$src/expected.tsv" "$dst/"
col=$(head -1 "$src/cases.tsv" | tr '\t' '\n' | grep -n '^embeddings_sha256$' | cut -d: -f1)
tail -n +2 "$src/cases.tsv" | cut -f1,"$col" | while read -r case sha; do
    [ "$(sha256sum "$src/$case/input/embeddings.npy" | cut -d' ' -f1)" = "$sha" ] || { echo "$case: input differs from cases.tsv" >&2; exit 1; }
    [ -e "$dst/inputs/$sha.npy" ] || cp --preserve=mode "$src/$case/input/embeddings.npy" "$dst/inputs/$sha.npy"
done
du -sb "$dst"
