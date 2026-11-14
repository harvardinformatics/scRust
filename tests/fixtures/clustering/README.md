# tests/fixtures/clustering

The committed part of milestone M8's reference data for Step 2, the clustering crate (docs/PLAN.md C4 and C5; results in [`docs/validation/clustering_m8.md`](../../../docs/validation/clustering_m8.md)). Every file is a copy of a file under `data/clustering/`, byte for byte or gzipped (`gzip -n -9`; the name ends in `.gz` where the source's does not), made by [`validation/slurm/fixtures_clustering.sh`](../../../validation/slurm/fixtures_clustering.sh) and listed in `validation/manifest.toml` with its source. The rest of M8's data (all 27 tier-T objects, the 7 tier-P objects, every run in both profiles) stays in `data/clustering/`, listed in the manifest.

## k0/: the generators (stage K0)

Made from the original code, compiled with the production R env's compiler and Makeconf flags (`validation/cpp/build.sh`):

- `javarandom_nextint.tsv`, `javarandom_permutation.tsv`: Seurat 5.5.1's `JavaRandom` (`ModularityOptimizer.cpp`): `nextInt(n)` from fresh generators (seeds as the `int` FindClusters passes, including negative ones; n from 1 to 2^31 − 1, powers of two included), and successive `generateRandomPermutation` draws from one generator, as the optimizer's passes make them.
- `kiss64_kiss.tsv`, `kiss64_index.tsv`, `kiss64_flip.tsv`, `kiss64_mixed.tsv`: RcppAnnoy 0.0.23's `Kiss64Random` (`kissrandom.h`): `kiss()`, `index(n)`, `flip()`, and one stream through mixed calls as an index build makes them. Seed 1234567890987654321 is the default seed every Seurat build uses.
- `groupsingletons_sample.tsv`: R 4.5.3's `set.seed(1); sample(x, 1)` for n tied names, the draw GroupSingletons makes, with `.Random.seed[2]` afterwards (the draw count).

Lines starting with `#` name the columns; lists are comma-separated.

## \<object\>/: seven tier-T objects

The two base objects (`filtered_seurat_tenx_test`, `filtered_seurat_emptydrops_test`) and the five DoubletFinder inputs (`soupx_tenx_test`, `soupx_emptydrops_test`, `decontx_raw_tenx_test`, `decontx_raw_emptydrops_test`, `cellbender_fromraw_test`), whose counts are M2's fixtures: `tests/fixtures/doubletfinder/<object>/input/`. For each:

- `input/assay.txt`: `v5`, or `v3` for the emptyDrops object (SCTransform then dispatches to `SCTransform.Assay` directly). The base objects also have their counts here (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`: the RNA counts layer as saved, entries exactly as stored).
- `saved/`: the saved object's `metadata.tsv.gz` (barcode, nCount_RNA, nFeature_RNA, nCount_SCT, nFeature_SCT, SCT_snn_res.0.8, seurat_clusters, silhouette_width, neighborhood_purity; doubles with 17 significant digits) and the levels of `seurat_clusters` and `SCT_snn_res.0.8`. Made in the production profile (OpenBLAS's SapphireRapids kernel): its silhouette widths differ from the portable profile's in the last bits.
- `reference/`: the portable profile's run (`OPENBLAS_CORETYPE=Haswell`, data/clustering/<object>/runs/instrumented-portable-sapphirerapids/, whose final files equal the direct run's): `metadata.tsv.gz`, the levels, `idents_levels.txt`, `cluster_ids.txt` and `nclusters.txt` (as the workflow writes them), and `stages/stages.json` (every scalar: the Annoy and optimizer settings, the purity radius, correct()'s median log_umi, ...).
- `optimizer/`: the instrumented modularity optimizer on the saved SNN graph: `starts.tsv` (each start's iterations and modularity, `%.17g`), `partitions.npy.gz` (each start's own cluster ids, starts × cells) and `ids.txt` (the final partition, numbered by size, before GroupSingletons).

For the two base objects, `reference/` also holds the stage references M9 and M10 check from R's own embeddings:

- `pca/embeddings.npy`: the PCA embeddings (cells × 50); Annoy's items are the first 30 columns rounded to float;
- `stages/annoy.ann.gz`: Annoy's index as `save()` writes it (K4); the first n nodes are the items (K3);
- `stages/K5_nn_index.npy.gz`, `stages/K5_nn_dist.npy.gz`: every cell's 20 neighbours (1-based) and float distances;
- `graphs/nn.mtx.gz`, `graphs/snn.mtx.gz`: `SCT_nn` and `SCT_snn` as stored (K6);
- `stages/K7_optimizer_ids.npy.gz`, `stages/K7_final_ids.npy.gz`: the optimizer's result and GroupSingletons';
- `stages/K8_*.npy.gz`: approxSilhouette's width, nearest other cluster, centroids and within-cluster variances;
- `stages/K9_*.npy.gz`: neighborPurity's per-cell 50th distance, neighbour counts and lists (CSR: `neighbors_ptr`, 1-based `neighbors_idx`), weights, totals, purity and maximum.

Index conventions are in [`validation/R/lib/clustering_hooks.R`](../../../validation/R/lib/clustering_hooks.R). `.npy` files are NumPy format 1.0, little-endian, Fortran order for matrices; read the `.gz` ones through `scio::gz`.

## edge/: constructed edge cases for the cluster crate (milestone M10)

Made by `validation/R/cluster_references.R --mode edge` (Seurat 5.5.1's FindClusters and the instrumented optimizer, in the pinned env) and copied by [`validation/slurm/fixtures_cluster.sh`](../../../validation/slurm/fixtures_cluster.sh); R's files themselves stay in `data/clustering/edge/` (listed in the manifest). `crates/cluster/tests/edge_cases.rs` reads:

- `cases.tsv`: each case and setting (53): cells, edges, FindClusters' arguments (resolutions as `%.17g`), R's error where R stopped, what the case exercises;
- `<case>/snn.mtx.gz`, `<case>/barcodes.tsv`: the graph (stored order, values `%.17g`) and its cells;
- `expected.tsv`: the bytes and sha256 of every file R wrote for each setting: `metadata.tsv`, the levels files, `cluster_ids.txt`, `nclusters.txt`, `stages/K7_optimizer_ids.npy`, `stages/K7_final_ids.npy`, `stages/rng.npy` and `rng.tsv` (`.Random.seed` at each GroupSingletons call's entry and exit, after `set.seed(20261114)`), and the instrumented optimizer's `optimizer/<res.r>/starts.tsv`, `partitions.npy` and `ids.txt`. The test recomputes each of them from Rust's result in R's format;
- `printed.tsv`: what FindClusters printed (but the elapsed time), its messages, R's error, and each call's singletons in GroupSingletons' order.

