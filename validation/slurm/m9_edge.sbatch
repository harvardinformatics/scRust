#!/bin/bash
# Milestone M9 (docs/PLAN.md C7): `scrust neighbors` on every constructed edge case of
# validation/R/neighbors_references.R --mode edge (data/neighbors/edge/), with each case's
# parameters and all three metrics, compared with R's references by `scrust validate neighbors`.
# cargo test checks the same cases through the library (crates/neighbors/tests/edge_cases.rs);
# this runs the command line on them and keeps the build's counts for the report. Seconds of
# work: `bash` on a login node is enough, or sbatch. From the repository's top directory:
#   bash validation/slurm/m9_edge.sbatch
# Outputs: data/neighbors/edge-rust/<case>/<metric>/ (scrust neighbors' files, summary.json,
# validate.json, validate.md; the index file is deleted once found identical, its sha256 kept),
# and data/neighbors/edge-rust/status.tsv. A case where FindNeighbors itself stops in R (one cell)
# must stop scrust neighbors too; status.tsv records both.
#
#SBATCH --job-name=scrust-m9-edge
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=0-00:30:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then
    source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh
fi
cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source ./env.sh

bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary $bin sha256 $(sha256sum "$bin" | cut -d' ' -f1) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
src=data/neighbors/edge
out=data/neighbors/edge-rust
[ -f "$src/cases.tsv" ] || { echo "no $src/cases.tsv" >&2; exit 1; }
[ -e "$out" ] && { echo "$out exists" >&2; exit 1; }
tmp=$out.partial.${SLURM_JOB_ID:-local}
rm -rf "$tmp"
mkdir -p "$tmp"
printf 'case\tmetric\trust\tr_graphs\tvalidate\n' > "$tmp/status.tsv"
status=0
# columns of cases.tsv: case n f k_param n_trees prune embeddings_sha256 what
while IFS=$'\t' read -r case n f k trees prune _sha _what; do
    for m in euclidean cosine manhattan; do
        o=$tmp/$case/$m
        mkdir -p "$(dirname "$o")"
        r_graphs=yes
        [ -f "$src/$case/$m/graphs/snn.mtx.gz" ] || r_graphs=no
        if "$bin" neighbors --pca "$src/$case/input" --out "$o" --dims "1:$f" --k "$k" --trees "$trees" \
                --metric "$m" --prune "$prune" --annoy yes --threads 4 2> "$tmp/$case.$m.stderr"; then
            rust=ok
            v=pass
            "$bin" validate neighbors --reference "$src/$case/$m" --candidate "$o" \
                --json "$o/validate.json" --markdown "$o/validate.md" 2>> "$tmp/$case.$m.stderr" || v=FAIL
            sha256sum "$o/annoy.ann" | sed "s|$o/||" > "$o/annoy.ann.sha256"
            [ "$v" = pass ] && rm "$o/annoy.ann"
        else
            rust=refused
            v=none
        fi
        # a refusal must match R's own failure to build the graphs, and only then
        if { [ "$rust" = refused ] && [ "$r_graphs" = yes ]; } || { [ "$rust" = ok ] && [ "$r_graphs" = no ]; } || [ "$v" = FAIL ]; then
            status=1
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$case" "$m" "$rust" "$r_graphs" "$v" >> "$tmp/status.tsv"
    done
done < <(tail -n +2 "$src/cases.tsv")
mv "$tmp" "$out"
column -t "$out/status.tsv"
echo "status $status"
exit $status
