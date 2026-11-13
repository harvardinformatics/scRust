#!/bin/bash
# Copies M8's committed tier-T fixtures (docs/PLAN.md A1 and C4) from data/clustering/ into
# tests/fixtures/clustering/ (see its README.md for what each file is). Reads plain files only, so a
# login node is enough. From the repository's top directory:
#   bash validation/slurm/fixtures_clustering.sh
# Files are copied byte for byte, or gzipped (gzip -n -9, so the bytes depend on the content only)
# where that saves space; validation/R/manifest_clustering.R checks each against its source.
#
# The selection keeps the committed fixtures near A1's ~20 MB: the K0 vectors; for the seven
# objects cargo test covers (the two base objects and the five DoubletFinder inputs, whose counts
# are M2's tests/fixtures/doubletfinder/<object>/input/), the expected labels and scores of the
# portable profile (reference/) and of the saved object (saved/), and the instrumented
# optimizer's per-start results; for the two base objects also their counts and the stage
# references M9 and M10 test from R's own embeddings (the embeddings, the Annoy index, every
# cell's neighbours, both graphs, the optimizer's and GroupSingletons' ids, silhouette's and
# purity's intermediates). Everything else stays in data/ and the manifest.
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/../.."
src_root=data/clustering
dst_root=tests/fixtures/clustering
primary=runs/instrumented-portable-sapphirerapids

copy() {  # source, destination
    mkdir -p "$(dirname "$2")"
    cp --preserve=mode "$1" "$2"
}
gz() {  # source, destination (.gz)
    mkdir -p "$(dirname "$2")"
    gzip -n -9 -c "$1" > "$2"
}

rm -rf "$dst_root/k0"
mkdir -p "$dst_root/k0"
for f in "$src_root"/k0/*.tsv; do copy "$f" "$dst_root/k0/$(basename "$f")"; done

all="filtered_seurat_tenx_test filtered_seurat_emptydrops_test soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test"
base_objects="filtered_seurat_tenx_test filtered_seurat_emptydrops_test"
for obj in $all; do
    s=$src_root/$obj
    d=$dst_root/$obj
    [ -f "$s/$primary/summary.json" ] || { echo "no portable instrumented run of $obj" >&2; exit 2; }
    rm -rf "$d"
    copy "$s/input/assay.txt" "$d/input/assay.txt"
    gz "$s/saved/metadata.tsv" "$d/saved/metadata.tsv.gz"
    for f in seurat_clusters_levels.txt SCT_snn_res.0.8_levels.txt; do copy "$s/saved/$f" "$d/saved/$f"; done
    gz "$s/$primary/metadata.tsv" "$d/reference/metadata.tsv.gz"
    for f in seurat_clusters_levels.txt SCT_snn_res.0.8_levels.txt idents_levels.txt cluster_ids.txt nclusters.txt; do
        copy "$s/$primary/$f" "$d/reference/$f"
    done
    copy "$s/$primary/stages/stages.json" "$d/reference/stages/stages.json"
    for f in starts.tsv ids.txt; do copy "$s/optimizer/$f" "$d/optimizer/$f"; done
    gz "$s/optimizer/partitions.npy" "$d/optimizer/partitions.npy.gz"
done
for obj in $base_objects; do
    s=$src_root/$obj
    d=$dst_root/$obj
    for f in matrix.mtx.gz barcodes.tsv.gz features.tsv.gz; do copy "$s/input/$f" "$d/input/$f"; done
    copy "$s/$primary/pca/embeddings.npy" "$d/reference/pca/embeddings.npy"
    copy "$s/$primary/graphs/nn.mtx.gz" "$d/reference/graphs/nn.mtx.gz"
    copy "$s/$primary/graphs/snn.mtx.gz" "$d/reference/graphs/snn.mtx.gz"
    gz "$s/$primary/stages/annoy.ann" "$d/reference/stages/annoy.ann.gz"
    for f in K5_nn_index K5_nn_dist K7_optimizer_ids K7_final_ids K8_silhouette_width K8_silhouette_other \
             K8_centroids K8_cluster_var K9_dist50 K9_neighbor_count K9_neighbors_ptr K9_neighbors_idx \
             K9_weights K9_totals K9_purity K9_maximum; do
        gz "$s/$primary/stages/$f.npy" "$d/reference/stages/$f.npy.gz"
    done
done
du -sb "$dst_root"
