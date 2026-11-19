#!/bin/bash
# Milestone M15b: the committed fixtures of the vars.to.regress tests (crates/scrust/tests/m15b.rs):
# for two tier-T objects (soupx_tenx_test, a v5 object; filtered_seurat_emptydrops_test, the v3
# path) and the four cases of validation/m15b_cases.tsv, the portable profile's R run
# (data/clustering/m15b/runs/<object>/<case>/direct-portable-*, OpenBLAS's Haswell kernel):
# metadata.tsv (gzipped), the levels and summary.json; and each case's inputs
# (data/clustering/m15b/inputs/<object>/: mt_genes.txt, zero_genes.txt, metadata.tsv). The counts
# are the committed fixtures' (identical to data/clustering/<object>/input). Copies only; no R.
# From the repository's top directory:
#   bash validation/slurm/fixtures_m15b.sh
set -euo pipefail
cd "${SCRUST_ROOT:-$PWD}"
for obj in soupx_tenx_test filtered_seurat_emptydrops_test; do
    dst=tests/fixtures/clustering/m15b/$obj
    rm -rf "$dst"
    mkdir -p "$dst/inputs"
    for f in mt_genes.txt zero_genes.txt metadata.tsv; do
        cp "data/clustering/m15b/inputs/$obj/$f" "$dst/inputs/$f"
    done
    for cs in pm ids ids_ncount meta; do
        src=$(ls -d data/clustering/m15b/runs/$obj/$cs/direct-portable-* | grep -v partial | head -n 1)
        [ -d "$src" ] || { echo "no portable R run for $obj $cs" >&2; exit 1; }
        grep -q '"openblas_corename": "Haswell"' "$src/summary.json" || { echo "$src: not the Haswell kernel" >&2; exit 1; }
        mkdir -p "$dst/$cs"
        gzip -n -9 -c "$src/metadata.tsv" > "$dst/$cs/metadata.tsv.gz"
        cp "$src/seurat_clusters_levels.txt" "$src/SCT_snn_res.0.8_levels.txt" "$src/summary.json" "$dst/$cs/"
        echo "$dst/$cs <- $src"
    done
done
du -sh tests/fixtures/clustering/m15b
