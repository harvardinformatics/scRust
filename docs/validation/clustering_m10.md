# FindClusters on R's own SNN graphs (milestone M10)

M10 ports Seurat 5.5.1's FindClusters to Rust, as the new `cluster` crate (docs/PLAN.md C1 step 5, C2, C3, C3c, C3d, C6): the modularity optimizer that `RunModularityClusteringCpp` runs (Seurat's C++ translation of ModularityOptimizer 1.3.0, with its Java-style generator), the Louvain algorithm, Louvain with multilevel refinement and smart local moving, both modularity functions, the best of the random starts and the numbering of its clusters, GroupSingletons, and the columns, levels and `seurat_clusters` that FindClusters.Seurat makes. It runs on R's own SNN graphs (*teacher forcing*, plan A5b): each R run's graph feeds the Rust stage, and everything after it is Rust's. None of it uses BLAS, so every check here is exact on any machine.

## Outcome

| Done when (plan C6's M10 row, and the brief) | Result |
|---|---|
| K0 (Java) exact | **Yes.** `JavaRandom` reproduces every vector Seurat's own `ModularityOptimizer.cpp` printed (M8): `nextInt(n)` for 43 values of n from 1 to 2^31 − 1 from 6 seeds (powers of two, negative seeds and the rejection loop's wrap-around region included), 100 draws each, and 15 streams of successive permutations. R's `set.seed(1); sample(x, 1)`, GroupSingletons' draw, matches for all 213 tie counts from 1 to 100,000, with `.Random.seed[2]` afterwards. |
| K7 exact on T and P from R's own SNN graphs, the objects with singletons included | **Yes, in all 98 instrumented runs** of the 27 tier-T and 7 tier-P objects, both profiles (Sapphire Rapids and Cascade Lake, and the v5 runs): the optimizer's ids before GroupSingletons, the final ids, the labels and levels, every line FindClusters printed but the elapsed time, GroupSingletons' singletons in order and its message, every random start's iterations, modularity (bitwise) and partition against the instrumented optimizer, the best start, `seurat_clusters` with its levels, `cluster_ids.txt` and `nclusters.txt`; for the four tier-P objects with singletons (1, 18, 22 and 40), R's `.Random.seed` after GroupSingletons; for tier P, the workflow jobs' logged numbers. |
| Other settings, against new R references | **Yes, 465 of 465.** FindClusters with algorithms 2 and 3, modularity function 2, five resolutions in one call (0.2, 0.4, 0.8, 1.2, 2.0), `group.singletons = FALSE`, starts 1 and 20, iterations 1 and 50 and seeds 0 and 1, on every tier-T and tier-P object's saved graph, each with the instrumented optimizer's per-start results: identical in every check in 461; in the other 4 (modularity function 2, where every cell ends a singleton), R stops with an error, and so does Rust, with R's message. |
| Constructed edge cases (C7) | **Yes.** 20 cases, 53 settings against R: isolated cells (ties of 2 and of 10), singletons whose earlier merges decide the later ones, networks of 63, 64, 65, 128, 256 and 1,024 cells, one and two cells, disconnected components, equal gains and equal-sized clusters, every cell a singleton, negative and extreme seeds. cargo test recomputes 717 of R's files from Rust's result and compares their sha256; all identical, and where R stops (3), Rust stops with R's message. |
| Threads | **Yes.** Resolutions run in parallel. Five resolutions at 1 and 32 threads on both 79,409-cell objects, and at 1 and 8 threads on all 34 objects: every file byte-identical, and `summary.json`'s result section; cargo test compares 1, 2, 3, 4 and 7 threads. |
| Time and memory on the 79,409-cell objects | **Measured** (one Xeon Platinum 8480CL node): R's FindClusters call 19.8–20.8 s at the workflow's settings, its process peak 1.10–1.13 GiB; Rust 14.1–14.8 s for the same clustering, the whole command 15.2–15.8 s, peak 153–157 MiB. Five resolutions: Rust 66–70 s at one thread, 15.0–15.3 s at 32 (R: 138–147 s, on another node). |
| `cargo test --workspace` passes, debug and release, no warnings | **Yes**: see [Tests](#tests). |

## What was built

| Where | M10 contents |
|---|---|
| `crates/cluster` (new) | `java_random.rs`: `JavaRandom` (`nextInt`, `generateRandomPermutation`). `network.rs`: `Network`, `matrixToNetwork` (the SNN graph's strictly lower triangle, walked as Eigen walks it), `createReducedNetwork`, `createSubnetworks`. `clustering.rs`: `Clustering`, `orderClustersByNNodes` (stable), `mergeClusters`. `vos.rs`: `calcQualityFunction`, local moving, Louvain, Louvain with multilevel refinement, smart local moving. `optimizer.rs`: `RunModularityClusteringCpp` (the argument checks with Seurat's messages, the rescaled resolution, the starts, the best start, what it prints). `singletons.rs`: GroupSingletons. `find_clusters.rs`: FindClusters over several resolutions, the columns and their levels, `seurat_clusters`; `find_clusters()`, the library entry point for M11's chain and M15's clustering on the fly. `params.rs`, `input.rs` (the graph from `scrust neighbors`' output or an R run's), `output.rs` (the files and `summary.json`'s result), `validate.rs` (K7). |
| `crates/scrust` | `scrust cluster --snn` (C3's clustering flags, `--graph-name`, `--keep-stages`, `--threads`), `scrust validate cluster`; flags without a value in the argument parser (`--no-group-singletons`, `--keep-stages`). |
| `validation/` | `R/cluster_references.R` (other settings, edge cases, R's benchmark), `R/report_clustering_m10.R`, `R/manifest_cluster.R`; `slurm/cluster_references.sbatch`, `m10_cluster.sbatch`, `m10_benchmark.sbatch`, `fixtures_cluster.sh`. |
| `tests/fixtures/clustering/edge/` (new, 0.74 MB) | the edge cases' graphs, settings, and R's file digests and printed lines. |

The workspace still uses only the Rust standard library. `rcompat` needed no addition: GroupSingletons uses its `set.seed`, `sample` and long-double `mean`, and the `res.` names its `as.character`, all exact since M1.

**The port.** Every function is translated statement by statement from Seurat's C++ and R code:

- **The network** is the strictly lower triangle of the SNN graph in stored order (column by column, rows ascending), each entry appended to both nodes' lists, so every list is ascending. Node weights (function 1) and the total edge weight are `std::accumulate` sums, left to right from 0.0.
- **Every sum in Seurat's order.** Cluster weights change by `+=` and `−=` in the order nodes are visited; a reduced network's node and edge weights are summed node by node, then list by list; the quality function adds each node's internal edge weights, then subtracts each cluster's `K·K·γ` in cluster order, then divides by `2W + self links`. Two starts that end at the same partition numbered differently therefore get modularities that differ in the last bits, and 1,364 of the 7,924 starts compared here did so (below): exact sums were necessary, not a refinement.
- **Every tie as Seurat breaks it**: a strictly larger gain wins and an equal one only for a lower cluster index; with no positive gain a node takes the last unused cluster; the first start with the maximum wins; `stable_sort` keeps equal-sized clusters in their order; GroupSingletons draws with R's `sample()` after `set.seed(1)`.
- **`nextInt`** takes Java's multiply-and-shift branch for powers of two and otherwise rejects with 32-bit wrap-around, as the compiled `Seurat.so` does. The seed is FindClusters' `int`, sign-extended to the generator's `uint64_t`.
- **GroupSingletons** computes each candidate's connectivity with R's `mean` (long double, with its second pass) over the singleton's SNN row at the candidate's current cells, zeros included; `SNN[i.cells, j.cells]` drops to a plain vector there, so `mean` is R's path (checked in R: `subsnn_is_object = false`).
- **What R prints** is reproduced line for line: the banner, nodes, edges, the algorithm, the maximum modularity as `%.4f` (or `Modularity:` for one start), communities, GroupSingletons' message.

## How the comparison works

For each R run, `scrust cluster --snn <graph> --keep-stages` reads the run's own SNN graph (`graphs/snn.mtx.gz`, entries in stored order, values parsed exactly from `%.17g`) and writes `seurat_clusters.tsv` and its levels, `clusters.tsv` and each column's levels, `cluster_ids.txt`, `nclusters.txt`, `metadata.tsv`, `optimizer_ids.tsv`, every start's partition (`stages/`) and `summary.json`. `scrust validate cluster --reference <R run> --candidate <out> --optimizer <instrumented optimizer>` then compares, resolution by resolution ([`validation/slurm/m10_cluster.sbatch`](../../validation/slurm/m10_cluster.sbatch); reports in `data/clustering/<object>/rust/m10/<run>/validate.{json,md}` and `rust/m10/refs/<setting>/`):

- **the settings** R recorded against the candidate's;
- **the optimizer's ids** (`stages/K7_optimizer_ids.npy`), **the final ids** (`K7_final_ids.npy`), **the labels** (`metadata.tsv`) and **the levels** (`<column>_levels.txt`, byte for byte);
- **what R printed** (`summary.json`, `run.printed`: every line but the elapsed time) and **GroupSingletons' messages**; the singletons in GroupSingletons' order (`stages.json`);
- **every random start** against the instrumented optimizer on the same graph (`starts.tsv`: iterations and modularity as a double, bitwise; the best start; `ids.txt`; `partitions.npy`, each start's partition);
- **R's `.Random.seed` after GroupSingletons** (`stages/rng.npy`): where Rust drew, R's state after the call must be Rust's; where it did not, R's state must not have changed;
- **`seurat_clusters`**, its levels, `idents_levels.txt`, `cluster_ids.txt`, `nclusters.txt`; and, for tier P, the numbers the workflow job logged (`validation/clustering_panel.tsv`).

**The R runs.**

- **The workflow's settings:** M8's 98 instrumented runs, whose hooks dumped the optimizer's ids, the final ids, the singletons and the RNG checkpoints, with M8's instrumented optimizer runs on the saved graphs (`data/clustering/<object>/optimizer/`; every run's graph equals the saved one, M8).
- **Other settings:** [`validation/R/cluster_references.R`](../../validation/R/cluster_references.R) runs FindClusters on each object's saved graph held as `SCT_snn` in a minimal Seurat object, with passive `trace()` hooks on `RunModularityClustering` and `GroupSingletons`, and the instrumented optimizer with the same arguments for every resolution, in the layout of an instrumented run (`data/clustering/<object>/cluster_refs/<setting>/`). The saved graphs hold their values as `%.17g` text, which R's own parser does not always read back exactly; each distinct value is corrected to the double whose text it is, and the matrix is checked to write back the same text (no object needed a correction: the 18 SNN values s/(40 − s) all parse exactly). The workflow setting reproduces every saved object's labels and levels and M8's printed numbers and messages (checked in each run's `summary.json`), so the rebuilt graph is the one FindClusters clustered.
- **Edge cases:** the same script with `--mode edge` builds the constructed graphs.

The real-data comparisons ran with the build of commit 5ddbd97, pinned as `data/clustering/bin/scrust-m10-5ddbd97` (sha256 bca7e5e4…). The benchmark ran with commit e0a40a9 (`scrust-m10-e0a40a9`, sha256 51fd8b2e…), whose clustering is the same: the only change after it is the validator's RNG check (below; `git diff e0a40a9 5ddbd97 -- crates/cluster/src` touches `validate.rs` only).

## Results

### Every comparison, by tier and kind

| tier | kind | comparisons | identical | objects | resolutions | starts | partition values | optimizer ids | labels | singletons | RNG states | R stopped | differing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T | instrumented | 82 | 82 | 27 | 82 | 820 | 179,790 | 17,979 | 35,958 | 0 | 0 | 0 | 0 |
| T | reference | 378 | 378 | 27 | 563 | 5,657 | 1,238,964 | 123,306 | 205,284 | 5,797 | 90 | 4 | 0 |
| P | instrumented | 16 | 16 | 7 | 16 | 160 | 6,255,350 | 625,535 | 1,251,070 | 184 | 9 | 0 | 0 |
| P | reference | 87 | 87 | 7 | 128 | 1,287 | 47,231,139 | 4,695,963 | 7,980,841 | 1,403 | 64 | 0 | 0 |

*identical*: every check identical, or (R stopped) Rust stopped with R's message. *starts*: random starts compared with the instrumented optimizer (iterations and modularity); *partition values*: their partitions, cell by cell; *labels*: the `res.` columns and `seurat_clusters`, cell by cell; *singletons*: singletons GroupSingletons handled (all resolutions); *RNG states*: R's `.Random.seed` compared after a GroupSingletons call that drew. One row per comparison, with the settings, clusters, maximum modularity, best start and timings: [`clustering_m10.tsv`](clustering_m10.tsv).

### Tier P, the workflow's settings

| object | run | cells | Q | communities | singletons | clusters | best start | distinct start modularities | K7 | log | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | portable-cascadelake | 3,108 | 0.8777 | 15 | 0 | 15 | 10 | 10 | = | = | 0.2 |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | 3,108 | 0.8777 | 15 | 0 | 15 | 10 | 10 | = | = | 0.2 |
| soupx_tenx_SRR18256555 | production-sapphirerapids | 3,108 | 0.8777 | 15 | 0 | 15 | 10 | 10 | = | = | 0.2 |
| soupx_tenx_10xmousebrain | portable-sapphirerapids | 11,831 | 0.8920 | 21 | 0 | 21 | 2 | 10 | = | = | 1.0 |
| soupx_tenx_10xmousebrain | production-sapphirerapids | 11,831 | 0.8920 | 21 | 0 | 21 | 2 | 10 | = | = | 1.0 |
| filtered_seurat_tenx_10xmousebrain | portable-sapphirerapids | 11,831 | 0.8898 | 21 | 0 | 21 | 2 | 10 | = | = | 1.0 |
| filtered_seurat_tenx_10xmousebrain | production-sapphirerapids | 11,831 | 0.8898 | 21 | 0 | 21 | 2 | 10 | = | = | 1.0 |
| soupx_emptydrops_4500871lane | portable-sapphirerapids | 26,598 | 0.9015 | 32 | 1 | 31 | 1 | 10 | = | = | 4.3 |
| soupx_emptydrops_4500871lane | production-sapphirerapids | 26,598 | 0.9015 | 32 | 1 | 31 | 1 | 10 | = | = | 4.5 |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | 79,409 | 0.8818 | 54 | 18 | 36 | 8 | 10 | = | = | 18.9 |
| soupx_emptydrops_sscovelli | production-sapphirerapids | 79,409 | 0.8818 | 54 | 18 | 36 | 8 | 10 | = | = | 19.2 |
| filtered_seurat_emptydrops_sscovelli | portable-sapphirerapids | 79,409 | 0.8722 | 59 | 22 | 37 | 9 | 10 | = | = | 14.7 |
| filtered_seurat_emptydrops_sscovelli | production-sapphirerapids | 79,409 | 0.8722 | 59 | 22 | 37 | 9 | 10 | = | = | 14.7 |
| filtered_seurat_emptydrops_sscovelli | production-sapphirerapids-v5 | 79,409 | 0.8722 | 59 | 22 | 37 | 9 | 10 | = | = | 17.6 |
| cellbender_fromraw_ghirsutum2 | portable-sapphirerapids | 59,323 | 0.8567 | 62 | 40 | 22 | 4 | 10 | = | = | 10.6 |
| cellbender_fromraw_ghirsutum2 | production-sapphirerapids | 59,323 | 0.8567 | 62 | 40 | 22 | 4 | 10 | = | = | 10.3 |

*Q*: the maximum modularity as printed (`%.4f`); *communities*: the optimizer's, as printed; *clusters*: after GroupSingletons; *K7*: every check identical; *log*: the workflow job's nodes, edges, maximum modularity, communities, singletons and final clusters (`validation/clustering_panel.tsv`) identical; *Rust s*: building the network and the resolution, one thread, on the `shared` partition's Cascade Lake nodes.

Production and portable runs of an object have the same SNN graph (M8), so their results are the same; each is checked against its own R run.

### Other settings

Every reference setting, summed over the objects (27 tier-T and 7 tier-P; modularity function 2 on tier T and, at lower resolutions, the three tier-P objects up to 11,831 cells):

| setting | FindClusters' arguments | comparisons | identical | resolutions | singletons | RNG states | R stopped | starts | Rust s | R s |
|---|---|---|---|---|---|---|---|---|---|---|
| workflow | `random.seed = 12345` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 53.4 | 106.0 |
| algorithm2 | `algorithm = 2` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 52.2 | 91.4 |
| algorithm3 | `algorithm = 3` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 198.2 | 273.3 |
| modularity2 | `modularity.fxn = 2` (tier T) | 27 | 27 | 23 | 4,940 | 23 | 4 | 230 | 0.3 | 8.5 |
| modularity2_resolutions | `modularity.fxn = 2`, 3 resolutions (T: 0.05, 0.1, 0.2; P: 0.002, 0.005, 0.01) | 30 | 30 | 90 | 175 | 49 | 0 | 900 | 2.6 | 15.4 |
| resolutions | `resolution = c(0.2, 0.4, 0.8, 1.2, 2.0)` | 34 | 34 | 170 | 762 | 46 | 0 | 1,700 | 48.7 | 435.2 |
| nogroup | `group.singletons = FALSE` | 34 | 34 | 34 | 81 | 0 | 0 | 340 | 53.1 | 60.8 |
| nogroup_resolutions | `group.singletons = FALSE`, `resolution = c(0.8, 2.0)` | 34 | 34 | 68 | 513 | 0 | 0 | 680 | 52.8 | 118.9 |
| starts1 | `n.start = 1` | 34 | 34 | 34 | 81 | 4 | 0 | 34 | 6.0 | 26.2 |
| starts20 | `n.start = 20` | 34 | 34 | 34 | 81 | 4 | 0 | 680 | 99.0 | 125.0 |
| iter1 | `n.iter = 1` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 38.0 | 58.7 |
| iter50 | `n.iter = 50` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 52.0 | 72.6 |
| seed0 | `random.seed = 0` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 46.4 | 70.8 |
| seed1 | `random.seed = 1` | 34 | 34 | 34 | 81 | 4 | 0 | 340 | 46.8 | 73.2 |

Every other argument is the workflow's (`random.seed = 12345` except in seed0 and seed1). *Rust s* and *R s* sum the clustering itself over the objects (Rust: the network and every resolution at 8 threads, on Cascade Lake nodes; R: the FindClusters call, on other `shared` nodes), not a benchmark. Tier P, each cell giving the clusters of each resolution and the singletons (all resolutions):

| setting | soupx_tenx_SRR18256555 | soupx_tenx_10xmousebrain | filtered_seurat_tenx_10xmousebrain | soupx_emptydrops_4500871lane | soupx_emptydrops_sscovelli | filtered_seurat_emptydrops_sscovelli | cellbender_fromraw_ghirsutum2 |
|---|---|---|---|---|---|---|---|
| workflow | = (15; 0) | = (21; 0) | = (21; 0) | = (31; 1) | = (36; 18) | = (37; 22) | = (22; 40) |
| algorithm2 | = (15; 0) | = (23; 0) | = (23; 0) | = (32; 1) | = (38; 18) | = (35; 22) | = (24; 40) |
| algorithm3 | = (15; 0) | = (24; 0) | = (24; 0) | = (32; 1) | = (38; 18) | = (37; 22) | = (25; 40) |
| modularity2_resolutions | = (15,19,28; 3) | = (37,73,115; 9) | = (40,75,113; 14) |  |  |  |  |
| resolutions | = (12,14,15,17,21; 0) | = (12,17,21,27,40; 0) | = (12,18,21,26,38; 0) | = (20,23,31,36,46; 5) | = (19,26,36,44,56; 90) | = (18,25,37,45,55; 110) | = (10,15,22,30,38; 200) |
| nogroup | = (15; 0) | = (21; 0) | = (21; 0) | = (32; 1) | = (37; 18) | = (38; 22) | = (23; 40) |
| nogroup_resolutions | = (15,21; 0) | = (21,40; 0) | = (21,38; 0) | = (32,47; 2) | = (37,57; 36) | = (38,56; 44) | = (23,39; 80) |
| starts1 | = (16; 0) | = (22; 0) | = (22; 0) | = (31; 1) | = (36; 18) | = (35; 22) | = (21; 40) |
| starts20 | = (15; 0) | = (21; 0) | = (21; 0) | = (31; 1) | = (36; 18) | = (37; 22) | = (24; 40) |
| iter1 | = (15; 0) | = (21; 0) | = (22; 0) | = (31; 1) | = (36; 18) | = (36; 22) | = (22; 40) |
| iter50 | = (15; 0) | = (21; 0) | = (21; 0) | = (31; 1) | = (36; 18) | = (37; 22) | = (22; 40) |
| seed0 | = (15; 0) | = (24; 0) | = (23; 0) | = (31; 1) | = (38; 18) | = (36; 22) | = (23; 40) |
| seed1 | = (15; 0) | = (24; 0) | = (22; 0) | = (31; 1) | = (37; 18) | = (38; 22) | = (21; 40) |

With `group.singletons = FALSE` the singletons are a cluster of their own, `singleton`: hence one cluster more than the workflow's for the objects with singletons.

**Where R stops.** With modularity function 2 at the grid's resolution 0.8, nearly every cell of a test object is a singleton (the resolution is not scaled by the total edge weight), and on 4 of the 27 (doubletfinder_cellbender_fromraw_test and three post-hoc objects) every cell is. GroupSingletons then has no cluster to join them to, `max()` of nothing warns, and `sample()` of no candidate stops R with "invalid first argument". `scrust cluster` stops there too, with R's message; on the other 23, the 4,940 singletons merge exactly as in R.

### The near-ties

1,364 of the 7,924 random starts compared, in 220 comparisons of 25 objects, ended within 1e-12 of their best start's modularity without equalling it. Each such start reached the best start's partition, numbered differently, so the quality function's sums ran over the clusters in another order (M8's finding). Every one of them is reproduced bitwise, and with them every best start: the first start at the maximum, by strict `>`. A sum taken in another order than Seurat's would have chosen another start in some of these, and with it other cluster numbers.

### Singletons

- **On real data at the workflow's settings,** every singleton is a cell without SNN edges: all 81 of the four tier-P objects (1, 18, 22, 40). Each ties at connectivity 0 with every cluster (31, 36, 37 and 22 candidates), so R's `sample()` places it, and every draw starts from `set.seed(1)`. Rust's state after the draw equals R's in all 9 runs with singletons.
- **Over all comparisons,** 7,384 singletons: 5,453 merged by their connectivity (13 of them tied between two or more clusters), 1,337 cells without edges tied with every cluster (up to 56), and 594 left as `singleton` (`group.singletons = FALSE`). Tier T's high resolutions and modularity function 2 produce singletons with edges whose earlier merges change the connectivities of the later ones: the merges reproduce R's in every case.
- **R's `.Random.seed`** after GroupSingletons was compared in 163 calls that drew, including draws whose rejection sampling took a second uniform (`.Random.seed[2]` = 2 after the call), and found unchanged in every call that did not draw.

### `nextInt`'s draws

| tier | comparisons | draws | power-of-two draws | comparisons with power-of-two draws |
|---|---|---|---|---|
| T | 456 | 4,797,452 | 48,230 | 377 |
| P | 103 | 210,836,800 | 106,350 | 51 |

No object has a power-of-two number of cells, so every power-of-two draw comes from a reduced network or subnetwork of 2^k nodes, the case C2 anticipated. Real data takes Java's other branch 154,580 times, and each time as R does.

## Constructed edge cases

`validation/R/cluster_references.R --mode edge` (job 60501674, under a minute) builds 20 graphs, deterministic or from a fixed seed, and runs FindClusters and the instrumented optimizer on each with one to four settings, into `data/clustering/edge/`. The committed fixture holds the graphs, the settings, and the size and sha256 of each of the 717 files R wrote, with what FindClusters printed (`tests/fixtures/clustering/edge/`, its [README](../../tests/fixtures/clustering/README.md)). cargo test (`crates/cluster/tests/edge_cases.rs`) runs every setting, recomputes each of R's files from Rust's result in R's format (the labels, levels and workflow files; `K7_optimizer_ids.npy`, `K7_final_ids.npy`, `rng.npy`, `rng.tsv`; the optimizer's `starts.tsv`, `partitions.npy` and `ids.txt`) and compares the digests, then compares the printed lines, messages and singletons, and R's errors. All identical:

| case | cells | edges | settings | what it exercises |
|---|---|---|---|---|
| isolated_two_clusters | 30 | 132 | workflow, nogroup, seed1 | two weighted cliques of 12 interleaved, 6 cells without edges at the end: ties of 2 |
| isolated_ten_clusters | 65 | 150 | workflow, nogroup | ten cliques of 6 and 5 cells without edges: ties of 10, which `sample()` resolves to the 9th |
| isolated_twelve_clusters | 64 | 120 | workflow | twelve cliques of 5 and 4 cells without edges; 64 nodes, a power of two |
| singleton_order | 6 | 5 | modularity2, modularity2_nogroup | pairs A = {1, 2}, B = {3, 4}; singleton 5 (to A 0.5, to 6 0.75) joins A, and then singleton 6 (to B 0.4) joins A through 5 |
| singleton_order_reversed | 6 | 5 | modularity2 | the same with cells 5 and 6 swapped |
| one_cell | 1 | 0 | workflow | R stops: "Matrix contained no network data.  Check format." |
| two_cells | 2 | 1 | workflow | two connected cells |
| two_cells_apart | 2 | 0 | workflow | R stops: no network |
| components | 20 | 42 | workflow, algorithm2, algorithm3 | cliques of 7, 5 and 5 (equal sizes keep their order), a cell without edges, a pair |
| ring8 | 8 | 8 | workflow, algorithm3, resolution 2 | a ring of 8 equal edges: every gain ties |
| lattice6 | 36 | 60 | workflow, starts1 | a 6 × 6 lattice of equal edges |
| all_singletons | 10 | 9 | modularity2, modularity2_nogroup | a path of weak edges: no node joins another; R stops ("invalid first argument") unless `group.singletons = FALSE` |
| snn63, snn64, snn65, snn128, snn256, snn1024 | 63–1,024 | 1,295–47,740 | workflow, three resolutions (0.3, 1, 3), algorithm2, algorithm3 | SNN graphs from FindNeighbors (k = 20) on 4-component Gaussian mixtures: powers of two and their neighbours |
| snn300_high_resolution | 300 | 7,336 | resolution 5, nogroup at 2, resolutions in the order 1.2, 0.4, 2.0, modularity2 at 0.05, 0.2, 0.5 | many singletons with edges; levels sorted as text with `singleton`; `seurat_clusters` from the last resolution given |
| seeds | 200 | 4,898 | seeds −1, 2147483647, −12345 | the `int` seed sign-extended to `uint64_t` |

The test also checks that the cases reach what they are for: 5 or more singletons tied with more than 8 clusters, the second singleton of singleton_order following the first, and more than 100 power-of-two draws on reduced networks where the input is not a power of two. The instrumented optimizer, given negative seeds as text, wraps them as R's conversion does (`std::stoull`), and equals R's optimizer in all 51 settings where R's optimizer finished.

## Benchmark

Both 79,409-cell objects, on their saved graphs, on one node (Xeon Platinum 8480CL, 32 CPUs, job 60510816, [`m10_benchmark.sbatch`](../../validation/slurm/m10_benchmark.sbatch)):

- **R**: the workflow's FindClusters call alone (`cluster_references.R --mode bench`).
- **Rust**: `scrust cluster --snn`, the whole command. Its *clustering s* is building the network and every resolution (the optimizer and GroupSingletons).

Wall time and peak RSS are from `/usr/bin/time -v`.

| object | side | resolutions | threads | clustering s | whole process | peak RSS MiB | note |
|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | R | 1 | 1 | 19.8 | 0:36.52 | 1,122 | VmHWM unchanged by the call (1,122 MiB: R, Seurat, the graph) |
| soupx_emptydrops_sscovelli | Rust | 1 | 1 | 14.1 | 0:15.19 | 157 | reading the graph and writing every output included |
| soupx_emptydrops_sscovelli | Rust | 5 | 1 | 66.4 | 1:07.48 | 156 | |
| soupx_emptydrops_sscovelli | Rust | 5 | 32 | 15.0 | 0:16.05 | 435 | |
| filtered_seurat_emptydrops_sscovelli | R | 1 | 1 | 20.8 | 0:35.84 | 1,153 | VmHWM 1,117 → 1,153 MiB across the call |
| filtered_seurat_emptydrops_sscovelli | Rust | 1 | 1 | 14.8 | 0:15.80 | 153 | reading the graph and writing every output included |
| filtered_seurat_emptydrops_sscovelli | Rust | 5 | 1 | 70.3 | 1:11.29 | 153 | |
| filtered_seurat_emptydrops_sscovelli | Rust | 5 | 32 | 15.3 | 0:16.39 | 416 | |

- **Threads.** Every output file, at 1 and 32 threads, is byte-identical, and so is `summary.json`'s result section (`threads.txt`).
- **Time.** One resolution is sequential, as in Seurat: the ten starts share one generator. Rust takes 14.1–14.8 s where R's call takes 19.8–20.8 s on the same node (31.2–36.0 s in M8's chain, on other nodes). The local moving visits nodes in random order, so the time goes to memory latency; removing the bounds checks of the inner loop saved about 8% in a trial, and was not kept. Several resolutions run in parallel: five take 15.0–15.3 s at 32 threads, the slowest resolution's time. R runs them in turn: 138.2 and 146.9 s in the `resolutions` references (another node).
- **Memory.** Rust holds the graph (3.2 million entries), the network (3.1 million list entries) and, per resolution in progress, the working arrays and reduced networks: 153–157 MiB at one thread, 416–435 MiB with five resolutions at once. R's process holds 1.1 GiB before the call (R, Seurat, the graph), and the call itself adds 0–36 MiB.

## Tests

`cargo test --workspace` passes in both profiles: 222 tests in 37 suites, none failing or ignored, in debug and in release, with no compiler warnings. 33 of them are new in M10:

- `crates/cluster`, unit tests (25): `JavaRandom` (Java's values, sign extension, permutations, draw counts), the network (lists, weights, the lower triangle in Eigen's order, reduced networks and subnetworks), `orderClustersByNNodes` and `mergeClusters`, the three algorithms on a small graph, isolated nodes, a single node drawing nothing, the argument checks with Seurat's messages, what R prints, smart local moving running every iteration, GroupSingletons (order of appearance, earlier merges counting, ties, `group.singletons = FALSE`, every cluster a singleton), the levels, the graph checks, the column names, parsing;
- `tests/k0_vectors.rs` (3): the K0 vectors;
- `tests/tier_t.rs` (3): K7 on both base objects' committed references (M8), every start bitwise with its partition, and five resolutions at 1, 3 and 7 threads;
- `tests/edge_cases.rs` (1): the 53 edge-case settings against R's 717 file digests, printed lines and errors;
- `crates/scrust/tests/cli.rs` (1): `scrust cluster` and `scrust validate cluster` end to end on a base object, five resolutions at 1 and 4 threads (identical files), `scrust neighbors`' own output as the input (the same labels as from R's graph), and the refusals (Leiden, another starting point, duplicate resolutions, Seurat's argument checks).

## Findings for the port

- **The near-ties are everywhere.** 17% of the starts compared ended within 1e-12 of their best start without equalling it: same partition, another numbering, another order of the quality function's sums. Only exact sums in Seurat's order reproduce which start wins.
- **Real data's singletons are cells without SNN edges.** All 81 in the tier-P objects have connectivity 0 with every cluster, so their clusters come from R's draw, not from the graph. The draw is exact, so this costs nothing; but the singletons' labels are arbitrary in a way the workflow's users may not expect.
- **`group.singletons = FALSE` changes `seurat_clusters`' level order.** The `res.` column keeps numbers first, sorted by value, then `singleton`; but `seurat_clusters` passes its levels through `as.numeric()`, which warns on `singleton`, and Seurat then sorts them as text: `0, 1, 10, 11, …, 2, 20, …, singleton`. GroupSingletons also returns before its draw and its message.
- **Modularity function 2 does not scale the resolution.** At the grid's 0.8 nearly every cell is a singleton, and where every cell is, R stops ("invalid first argument" from `sample()`). Resolutions that give clusters depend on the object's size: 0.05–0.2 for the 300-cell test objects, 0.002–0.01 for 3,000–12,000 cells, below 1e-3 for the larger objects (where R's GroupSingletons, at singletons × clusters × cells, would take hours).
- **One draw can leave `.Random.seed` where it was.** When several resolutions draw in GroupSingletons, each after `set.seed(1)`, a call that consumes as many uniforms as the previous one ends in the state it started from. Whether R drew must be read from the singletons, not from the RNG states; the validator first did the latter and flagged 28 checks in 21 comparisons as differing (fixed in 5ddbd97, before the final runs).
- **R's parser does not round trip `%.17g` in general,** but the 18 SNN values do. The reference script corrects any value it misreads, and none needed it.
- **`SNN[i.cells, j.cells]` is a plain vector** in GroupSingletons (one singleton row, two or more cluster columns), so its connectivity is R's `mean()` with the second pass, not `sum()/(n·m)` (checked in R).
- **Seurat's parallel path drops `group.singletons`.** With a `future` plan of several workers, FindClusters calls GroupSingletons without it, so singletons are grouped even when the caller asks otherwise. The workflow runs the sequential path, which the port reproduces.

## Corrections to the plan

Corrected in [`docs/PLAN.md`](../PLAN.md) on this branch:

- **C1 step 5 and C3d, `seurat_clusters`' levels:** sorted as text when a level is `singleton` (`group.singletons = FALSE`), not numerically.
- **C2 and C3d, GroupSingletons:** no draw and no message with `group.singletons = FALSE`; R stops when every cluster is a singleton; every tier-P singleton is a cell without edges; Seurat's parallel path drops `group.singletons`.
- **C3d, column names:** two resolutions with the same `res.` name overwrite each other in R; the port refuses them.
- **C3, commands and outputs:** `scrust cluster --snn` exists, with `--graph-name`, `--keep-stages` and `--threads`; `optimizer_ids.tsv` and, with `--keep-stages`, every start's partition.
- **C3, components:** `rcompat` needed no addition for M10 (`set.seed`, `sample`, the long-double `mean` and the number formatting exist since M1).
- **C4, the parameter grid:** modularity function 2 at 0.8 is degenerate, and R stops on 4 of the 27 test objects.

Plan facts the data confirmed as written: the network of the lower triangle in stored order, with ascending lists; the resolution rescaled for function 1; one generator per resolution, seeded once; `nextInt`'s wrap-around and its power-of-two branch, reached on reduced networks; equal gains to the lowest cluster index; the last unused cluster when no gain is positive; smart local moving running all its iterations; the first start at the maximum by strict `>`; `stable_sort` by size; the singletons and candidates in order of first appearance; the long-double `mean` over the current cells, earlier merges included; `set.seed(1)` before every draw; the numeric levels; `seurat_clusters` from the last resolution; the `res.` names at 15 significant digits.

## Deviations and notes

- **Duplicate resolution names are refused.** FindClusters would overwrite the earlier column; the port refuses the call instead of reproducing a lost column.
- **Graphs the optimizer could misread are refused:** not square, entries out of stored order, non-finite values, or non-positive weights below the diagonal. Seurat's graphs never have them, and the optimizer's bookkeeping (`edgeWeightPerCluster[l] == 0` marks a cluster not yet seen) is only defined for positive weights; R would read past its buffers on some of them.
- **`--stats`** defaults to `silhouette,purity`, as the workflow computes them. This version clusters only: when the scores are asked for, it says on stderr and in `summary.json` (`run.stats`) that they are not computed yet (silhouette in M11, purity in M12) and writes the labels. `--pca`, `--purity-k`, `--purity-weighted`, `--counts`, `--sct`, `--preset`, `--params`, `--export` and `--metadata` are refused with a message naming M11–M13.
- **`--no-group-singletons` and `--keep-stages` take no value**: `scrust`'s argument parser now accepts such flags (M9 had to write `--annoy yes|no`).
- **The graph's name** (`--graph-name`, default `SCT_snn`) is FindClusters' `graph.name`: it only names the columns.
- **Not ported:** the optimizer's standalone `main()` and its edge-file reader (FindClusters' `edge.file.name`), `identifyComponents`, `removeCluster` and `removeSmallClusters`, which FindClusters never calls; Leiden (`algorithm = 4`, refused); `initial.membership` and `node.sizes` (Leiden's).
- **The R references are larger than M9's**: 465 settings (14 per tier-T object, 12 or 13 per tier-P object), because the brief asked for every tier-T object and as many tier-P objects as stayed cheap: all seven took 0.5 to 15.5 minutes of R each. They take 355 MB in `data/clustering/<object>/cluster_refs/`, the edge cases 6.4 MB.
- **The manifest** gains a "FindClusters" section from [`validation/R/manifest_cluster.R`](../../validation/R/manifest_cluster.R), before M9's "Neighbours" section: 468 run tables and 606 file entries. Its compact form lists one `SHA256SUMS` per reference directory (10,195 files in all; `sha256sum -c SHA256SUMS` checks them) instead of a table per file, which would have added about 3 MB. The committed fixture and the edge cases' graphs are listed file by file. Rust's validation outputs (`data/clustering/<object>/rust/m10/`) are not listed, as in M9.
- **Committed fixtures:** `tests/fixtures/` grows from 23.46 MB to 24.20 MB (the edge cases; 0.34 MB of it the 1,024-cell graph). Only digests of R's files are committed; the files are in `data/clustering/edge/`.
- **Licences:** NOTICE now lists the code ported from Seurat's `ModularityOptimizer.cpp` and `RModularityOptimizer.cpp` (MIT; Waltman and van Eck's ModularityOptimizer 1.3.0, translated by Seurat's authors) and from Seurat's R code (GroupSingletons, FindClusters). Each source file carries the port notice.

## Reproducing

From the repository's top directory, after `source env.sh`, `cargo build --release -p scrust` and `bash validation/cpp/build.sh`:

```bash
sbatch validation/slurm/cluster_references.sbatch settings soupx_tenx_test   # any objects; the tier picks the settings
sbatch validation/slurm/cluster_references.sbatch edge
bash validation/slurm/fixtures_cluster.sh
sbatch validation/slurm/m10_cluster.sbatch soupx_tenx_test soupx_emptydrops_sscovelli   # any objects
sbatch validation/slurm/m10_benchmark.sbatch soupx_emptydrops_sscovelli filtered_seurat_emptydrops_sscovelli
Rscript validation/R/report_clustering_m10.R docs/validation/clustering_m10.tsv   # through srun
Rscript validation/R/manifest_cluster.R                                           # through srun
```

One run by hand, on a committed fixture (a second, any machine):

```bash
f=tests/fixtures/clustering/filtered_seurat_tenx_test
target/release/scrust cluster --snn $f/reference/graphs --out /path/to/out --seed 12345 --keep-stages
target/release/scrust validate cluster --reference $f/reference --candidate /path/to/out --optimizer $f/optimizer
```

R runs with the hermetic flags of `validation/R/README.md`, through sbatch or srun. SLURM logs of the M10 jobs: `scrust-m10-*.sbatch.log` in the submit directory (gitignored; job ids above and in the manifest).
