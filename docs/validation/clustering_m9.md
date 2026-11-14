# The neighbour stage on R's own embeddings (milestone M9)

M9 ports Seurat 5.5.1's FindNeighbors path to Rust, as the new `neighbors` crate (docs/PLAN.md C1 step 4, C3a, C3b, C6): RcppAnnoy 0.0.23's `Kiss64Random`, the Annoy 1.17.3 index it builds for the Euclidean, Angular (`annoy.metric = "cosine"`) and Manhattan metrics, Annoy's file format, the search, the `SCT_nn` graph and `ComputeSNN` with pruning. It runs on R's own inputs (*teacher forcing*, plan A5b): each R run's PCA embeddings feed the Rust stage, and everything after them is Rust's. None of it uses BLAS, so every check here is exact on any machine.

## Outcome

| Done when (plan C6's M9 row, and the brief) | Result |
|---|---|
| K0 (Kiss64) exact | **Yes.** `Kiss64Random` reproduces every vector RcppAnnoy's own `kissrandom.h` printed (M8): `kiss()`, `index(n)` for 21 values of n from 1 to 2^32 + 1, `flip()`, and a stream of 1,449 mixed calls, from the default seed and three others. |
| K3–K6 exact on T and P from R's own embeddings | **Yes, in all 98 instrumented runs** of the 27 tier-T and 7 tier-P objects, both profiles (Sapphire Rapids and Cascade Lake, and the v5 runs): every float item (K3), every byte of the index file (K4; up to 440,701 nodes, 59.9 MB), every cell's 20 neighbour ids and float distances (K5), and both graphs' dimensions, stored entries in order and values bitwise (K6). 12.9 million neighbour ids, as many distances and 28.3 million SNN entries, none differing. |
| Angular and Manhattan | **Yes, 114 of 114.** New R references (FindNeighbors with `annoy.metric` cosine and manhattan, RcppAnnoy's saved index) on the embeddings of 57 instrumented runs: the 27 tier-T objects in both profiles, soupx_tenx_SRR18256555 in both, soupx_emptydrops_sscovelli (79,409 cells) in the portable profile. Every check identical, as for Euclidean. |
| Constructed edge cases | **Yes.** 25 cases × 3 metrics against R (duplicated cells, margins of exactly 0, failed splits kept imbalanced or randomised, n ≤ K, fewer cells than k, powers of two, overflow, zero vectors, other parameters): cargo test compares all 732 files R wrote by sha256, identical, and checks that every metric reaches every rare branch; the command line agrees on all 72 that R completes, and refuses the 3 one-cell cases where R's FindNeighbors itself fails. |
| Threads | **Yes.** At 1 and 32 threads on both 79,409-cell objects, every output file, `annoy.ann` included, is byte-identical; cargo test compares 1, 3, 4 and 7 threads on tier T. |
| Time and memory on the largest object | **Measured** (79,409 cells): R's FindNeighbors call 27.2–27.8 s (31.4–33.4 s inside M8's chain), its process peak 1.6 GiB; Rust 20.1–20.4 s at one thread and 4.2–4.4 s at 32 (build 3.7 s, sequential as in RcppAnnoy), the whole command 24.8 s and 6.8–7.3 s with reading and writing, peak RSS 262–452 MiB. |
| `cargo test --workspace` passes, debug and release, no warnings | **Yes**: see [Tests](#tests). |

No FFI fallback was needed: the port is byte-identical everywhere.

## What was built

| Where | M9 contents |
|---|---|
| `crates/neighbors` (new) | `kiss.rs`: `Kiss64Random`. `annoy.rs`: Annoy's node array in its exact memory layout (so the file is the array), `add_item`, `build` (`thread_build`, `_make_tree`, `create_split`, `two_means`, `side`, `_split_imbalance`, the root copies), `get_nns_by_vector` (`_get_all_nns`), `save`/`to_bytes`, `from_bytes` (Annoy's `load()`, roots found as Annoy finds them), and the build's counts of its rare branches. `stl.rs`: libstdc++'s `priority_queue` (`push_heap`, `pop_heap`) and `partial_sort` (`__heap_select`, `__sort_heap`), move for move. `graph.rs`: `SCT_nn` and `ComputeSNN`, parallel over columns. `params.rs`: FindNeighbors' arguments with Seurat's defaults, `dims` and `a/b` parsing. `lib.rs`: `find_neighbors`, the in-process entry point (for M11's chain and M15's clustering on the fly). `validate.rs`: the comparison with an R run. |
| `crates/scio` | `mtx`: plain Matrix Market files written gzipped (pieces formatted and compressed in parallel, the bytes independent of the thread count) and read in stored order; `npy::encode_f64`/`encode_i32`. |
| `crates/scrust` | `scrust neighbors` (C3's flags, plus `--annoy yes\|no` and `--threads`) and `scrust validate neighbors`. |
| `validation/` | `R/neighbors_references.R` (Angular and Manhattan references, edge cases, R's benchmark), `R/report_clustering_m9.R`, `R/manifest_neighbors.R`; `slurm/neighbors_references.sbatch`, `m9_neighbors.sbatch`, `m9_edge.sbatch`, `m9_benchmark.sbatch`, `fixtures_neighbors.sh`. |
| `tests/fixtures/neighbors/edge/` (new, 0.38 MB) | the edge cases' inputs, parameters and R's file digests. |

The workspace still uses only the Rust standard library.

**The port.** Every function is translated statement by statement from `annoylib.h`, `kissrandom.h` and Seurat's R and C++ code, with the arithmetic as compiled for R (`-O2 -march=nocona`, checked in `RcppAnnoy.so` with `objdump`):

- **float throughout**, with no FMA. GCC vectorizes only element-wise products, sums and differences (`mulps`, `subps`, `addps` across a vector's components); every reduction stays a sequential `addss` chain, so scalar Rust gives the same bits.
- **two places mix in double**, as written: `_split_imbalance` (the share in double, rounded to float, compared with the doubles 0.95 and 0.99) and the Angular distance `2.0 - 2.0 * pq / sqrt(ppqq)`, whose `sqrt` is the float one (`sqrtss`) and whose quotient and difference are double, rounded to float on return.
- **Node ids, leaves and the root copies** come from appending to the array in Annoy's order. Bytes Annoy never writes are zero, as `_reallocate_nodes` leaves them.
- **The search** uses libstdc++'s exact heap moves. For finite distances the order is total and any correct queue would do. The overflow edge case produces NaN distances, where only the same moves give R's order; the port matches it there too.

## How the comparison works

For each R run with FindNeighbors' stage dumps, `scrust neighbors --pca <run>/pca --dims 1:30 --metric <m> --annoy yes` reads the run's own embeddings and writes `knn_index.tsv`, `knn_distance.tsv`, `nn.mtx.gz`, `snn.mtx.gz`, `annoy.ann` and `summary.json`. `scrust validate neighbors --reference <R dir> --candidate <out>` then compares ([`validation/slurm/m9_neighbors.sbatch`](../../validation/slurm/m9_neighbors.sbatch); reports in `data/clustering/<object>/rust/m9/<run>/<metric>/validate.{json,md}`):

- **the parameters** R recorded (metric, dimensions, trees, k, search.k) against the candidate's;
- **K3**: R's `getItemsVector` for every cell against the first n nodes of the candidate's index, every float bitwise;
- **K4**: `annoy.ann` against R's saved index, byte for byte. On a difference the report names the first differing node and field; none occurred;
- **K5**: every cell's k ids and distances, bitwise;
- **K6**: `SCT_nn` and `SCT_snn` as R stored them: dimensions, every entry in stored order, values bitwise.

**The R runs.**

- **Euclidean:** M8's instrumented runs, whose hooks dumped the items, a second identical build's saved index, the neighbours and both graphs.
- **Angular and Manhattan:** [`validation/R/neighbors_references.R`](../../validation/R/neighbors_references.R) writes the same layout:
  - the index of Seurat's own `AnnoyBuildIndex`, saved;
  - the items;
  - `FindNeighbors(embeddings[, 1:30], annoy.metric = m)`'s neighbours (`return.neighbor = TRUE`) and graphs.

  It also checks that FindNeighbors' search equals Seurat's `AnnoySearch` on a second, unsaved build and on the saved index. Both held in every reference.
- **Edge cases:** the same script with `--mode edge` builds the constructed inputs.

The Rust index files are deleted once found identical (they duplicate R's), keeping their sha256.

The real-data comparisons and the benchmark ran with the build of commit db5e474, pinned as `data/clustering/bin/scrust-m9-db5e474` (sha256 f89b2993…). The edge cases ran through the command line with the build of commit 5b9b26f (`scrust-m9-5b9b26f`, sha256 3454691f…), which adds the validator's reading of NaN distances (below). Nothing after db5e474 changes a computation: the other changes are tests, scripts, documentation and two doc comments (`git diff db5e474 -- crates/neighbors/src`).

## Results

### Every comparison, by tier and metric

| tier | metric | comparisons | identical | objects | cells | items (K3) | index nodes (K4) | neighbour ids (K5) | distances (K5) | nn entries (K6) | snn entries (K6) | differing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T | euclidean | 82 | 82 | 27 | 17,979 | 539,370 | 99,999 | 359,580 | 359,580 | 359,580 | 1,891,171 | 0 |
| T | cosine | 54 | 54 | 27 | 11,808 | 354,240 | 66,476 | 236,160 | 236,160 | 236,160 | 673,732 | 0 |
| T | manhattan | 54 | 54 | 27 | 11,808 | 354,240 | 65,784 | 236,160 | 236,160 | 236,160 | 1,499,020 | 0 |
| P | euclidean | 16 | 16 | 7 | 625,535 | 18,766,050 | 3,467,858 | 12,510,700 | 12,510,700 | 12,510,700 | 26,422,715 | 0 |
| P | cosine | 3 | 3 | 2 | 85,625 | 2,568,750 | 562,583 | 1,712,500 | 1,712,500 | 1,712,500 | 4,615,961 | 0 |
| P | manhattan | 3 | 3 | 2 | 85,625 | 2,568,750 | 473,957 | 1,712,500 | 1,712,500 | 1,712,500 | 3,312,159 | 0 |

Euclidean tier T: 27 objects × the portable runs on Sapphire Rapids and Cascade Lake and the production run, plus filtered_seurat_emptydrops_test's v5 run. One row per comparison, with the build's counts and timings: [`clustering_m9.tsv`](clustering_m9.tsv).

### Tier P, run by run

Each row is identical to its R run in every check (K3–K6 `=`). The last columns count the build's rare branches; these are the branches C7 names, and real data reaches them:

| object | run | metric | cells | index nodes | node bytes | K3-K6 | snn entries | self missing | zero margins | retried | kept imbalanced | random sides | depth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | portable-sapphirerapids | cosine | 3,108 | 19,905 | 132 | = | 227,922 | 0 | 0 | 346 | 125 | 4 | 15 |
| soupx_tenx_SRR18256555 | production-sapphirerapids | cosine | 3,108 | 19,905 | 132 | = | 227,922 | 0 | 0 | 346 | 125 | 4 | 15 |
| soupx_tenx_SRR18256555 | portable-cascadelake | euclidean | 3,108 | 17,565 | 136 | = | 222,928 | 0 | 0 | 12 | 1 | 0 | 10 |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | euclidean | 3,108 | 17,565 | 136 | = | 222,928 | 0 | 0 | 12 | 1 | 0 | 10 |
| soupx_tenx_SRR18256555 | production-sapphirerapids | euclidean | 3,108 | 17,565 | 136 | = | 222,928 | 0 | 0 | 12 | 1 | 0 | 10 |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | manhattan | 3,108 | 17,464 | 136 | = | 227,906 | 0 | 0 | 3 | 0 | 0 | 9 |
| soupx_tenx_SRR18256555 | production-sapphirerapids | manhattan | 3,108 | 17,464 | 136 | = | 227,906 | 0 | 0 | 3 | 0 | 0 | 9 |
| soupx_tenx_10xmousebrain | portable-sapphirerapids | euclidean | 11,831 | 65,629 | 136 | = | 735,417 | 0 | 0 | 32 | 2 | 0 | 11 |
| soupx_tenx_10xmousebrain | production-sapphirerapids | euclidean | 11,831 | 65,629 | 136 | = | 735,417 | 0 | 0 | 32 | 2 | 0 | 11 |
| filtered_seurat_tenx_10xmousebrain | portable-sapphirerapids | euclidean | 11,831 | 65,696 | 136 | = | 736,123 | 0 | 1 | 30 | 4 | 0 | 12 |
| filtered_seurat_tenx_10xmousebrain | production-sapphirerapids | euclidean | 11,831 | 65,696 | 136 | = | 736,123 | 0 | 1 | 30 | 4 | 0 | 12 |
| soupx_emptydrops_4500871lane | portable-sapphirerapids | euclidean | 26,598 | 147,106 | 136 | = | 1,483,292 | 0 | 0 | 108 | 41 | 2 | 14 |
| soupx_emptydrops_4500871lane | production-sapphirerapids | euclidean | 26,598 | 147,106 | 136 | = | 1,483,292 | 0 | 0 | 108 | 41 | 2 | 14 |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | cosine | 79,409 | 522,773 | 132 | = | 4,160,117 | 1686 | 0 | 17029 | 18708 | 17254 | 35 |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | euclidean | 79,409 | 440,121 | 136 | = | 3,208,769 | 0 | 1 | 409 | 111 | 14 | 17 |
| soupx_emptydrops_sscovelli | production-sapphirerapids | euclidean | 79,409 | 440,121 | 136 | = | 3,208,769 | 0 | 1 | 409 | 111 | 14 | 17 |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | manhattan | 79,409 | 439,029 | 136 | = | 2,856,347 | 0 | 1 | 212 | 63 | 6 | 16 |
| filtered_seurat_emptydrops_sscovelli | portable-sapphirerapids | euclidean | 79,409 | 440,701 | 136 | = | 3,040,053 | 0 | 2 | 401 | 74 | 10 | 17 |
| filtered_seurat_emptydrops_sscovelli | production-sapphirerapids | euclidean | 79,409 | 440,701 | 136 | = | 3,040,053 | 0 | 2 | 401 | 74 | 10 | 17 |
| filtered_seurat_emptydrops_sscovelli | production-sapphirerapids-v5 | euclidean | 79,409 | 440,701 | 136 | = | 3,040,053 | 0 | 2 | 401 | 74 | 10 | 17 |
| cellbender_fromraw_ghirsutum2 | portable-sapphirerapids | euclidean | 59,323 | 327,978 | 136 | = | 2,153,285 | 0 | 221 | 173 | 20 | 0 | 16 |
| cellbender_fromraw_ghirsutum2 | production-sapphirerapids | euclidean | 59,323 | 327,978 | 136 | = | 2,153,285 | 0 | 221 | 173 | 20 | 0 | 16 |

Column key:

- *self missing*: cells not in their own neighbour list.
- *zero margins*: points whose margin was exactly 0 in a split, each costing a `flip()`.
- *retried*: splits accepted at the second or third attempt.
- *kept imbalanced*: three attempts at 0.95 or more, the last at 0.99 or less.
- *random sides*: the fallback.
- *depth*: the deepest split node.

Production and portable runs of an object differ in their embeddings' last bits and in some components' signs (M8). Each is checked against its own R run.

### The rare branches on real data

| tier | metric | comparisons | split nodes | zero margins | in runs | retried splits | kept imbalanced | random sides | cells missing from their own list |
|---|---|---|---|---|---|---|---|---|---|
| T | euclidean | 82 | 36,910 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | cosine | 54 | 24,634 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | manhattan | 54 | 24,288 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | euclidean | 16 | 1,420,481 | 452 | 9 | 2,743 | 581 | 62 | 0 |
| P | cosine | 3 | 243,002 | 0 | 0 | 17,721 | 18,958 | 17,262 | 1,686 |
| P | manhattan | 3 | 194,025 | 1 | 1 | 218 | 63 | 6 | 0 |

Tier T's shallow trees (depth 4–6) never leave the main path, and tier P reaches every branch. One miscounted draw would change every later tree (C7), and the index files are identical, so every draw was counted as R counts it.

- **Duplicated cells.** cellbender_fromraw_ghirsutum2's twin cells (M8) give 221 zero margins. Its 6 cells with a lower-numbered twin get the twin first and keep themselves second, as in R; 59,317 cells are first in their own list.
- **Angular on sscovelli.** 17,254 splits fall back to random sides. Their zero normals give both children the same priority, so 1,686 cells do not find themselves among the 1,000 candidates. R does exactly the same (K5 identical).

## Angular and Manhattan references

The R references are in `data/clustering/<object>/neighbors/<run>/<metric>/` (manifest section "Neighbours"). They come from the embeddings of the runs listed in the outcome: 56 runs with both metrics in one job (1.3 min, 0.6 GiB peak), and soupx_emptydrops_sscovelli in another (4 min, 2.9 GiB). Each confirms, in R, that FindNeighbors searches the same index a fresh `AnnoyBuildIndex` builds, and that the saved index searches identically.

- **Angular** builds a different forest from Euclidean (132-byte nodes, 4.2 million SNN entries on sscovelli against 3.2 million) and is checked as thoroughly: every node, neighbour and graph entry.
- **Seurat's cosine conversion is dead code.** Seurat means to convert Angular distances to cosine distances, but `res$dist <- 0.5 * (res$dist * res$dist)` adds a list element that is never read. The distances it returns are Annoy's Angular ones, and the port's are identical to them (plan C3a).

## Constructed edge cases

`validation/R/neighbors_references.R --mode edge` (job 60383021, 42 s) builds 25 inputs, deterministic or from a fixed seed, and runs FindNeighbors and RcppAnnoy on each with all three metrics, into `data/neighbors/edge/`. The committed fixture holds the inputs, the parameters and the size and sha256 of each of the 732 files R wrote (`tests/fixtures/neighbors/edge/`, its [README](../../tests/fixtures/neighbors/README.md)).

- **cargo test** (`crates/neighbors/tests/edge_cases.rs`) recomputes everything and compares the 732 digests. It also checks that each metric reaches each rare branch: per metric, at least one margin of exactly 0, a retried split, a split kept imbalanced, a random-side split and a cell missing from its own list.
- **The command line** ([`m9_edge.sbatch`](../../validation/slurm/m9_edge.sbatch)) runs `scrust neighbors` and `scrust validate neighbors` on every case, into `data/neighbors/edge-rust/`. The cells give the result and the build's counts:

| case | n | f | k | trees | prune | euclidean | cosine | manhattan |
|---|---|---|---|---|---|---|---|---|
| duplicates | 240 | 30 | 20 | 50 | 1/15 | = (0-margins 5,699, self missing 25) | = (0-margins 6,315, retried 4, self missing 25) | = (0-margins 5,363, self missing 25) |
| lattice | 225 | 2 | 20 | 50 | 1/15 | = | = (0-margins 7,509, retried 86, random 22, self missing 1) | = (0-margins 15) |
| two_sites | 150 | 3 | 20 | 50 | 1/15 | = (0-margins 32,637, retried 26, self missing 100) | = (0-margins 32,593, retried 24, self missing 100) | = (0-margins 32,462, retried 17, self missing 100) |
| ulp_outlier | 200 | 4 | 20 | 50 | 1/15 | = (0-margins 56,752, retried 10, self missing 179) | = (0-margins 56,797, retried 30, self missing 179) | = (0-margins 56,752, retried 10, self missing 179) |
| exp_spacing | 250 | 2 | 20 | 50 | 1/15 | = (0-margins 34,069, retried 267, kept 1075, random 4, self missing 27) | = (0-margins 25,679, retried 1007, kept 109, random 2465, self missing 191) | = (0-margins 79,695, retried 59, random 1) |
| exp_spacing_dup | 240 | 2 | 20 | 50 | 1/15 | = (0-margins 5,130, retried 224, kept 798, random 1) | = (0-margins 2,217, retried 470, kept 304, random 3173, self missing 154) | = (0-margins 4,886, retried 187, kept 628) |
| collinear | 150 | 6 | 20 | 50 | 1/15 | = | = (0-margins 34,759, retried 349, random 75, self missing 150) | = (0-margins 172) |
| zero_rows | 100 | 4 | 20 | 50 | 1/15 | = (0-margins 3,795, retried 1) | = (0-margins 5,127, retried 1) | = (0-margins 3,623, retried 1) |
| n1_f30 | 1 | 30 | 20 | 50 | 1/15 | refused, as R | refused, as R | refused, as R |
| n2_f30 | 2 | 30 | 20 | 50 | 1/15 | = | = | = |
| n3_f30 | 3 | 30 | 20 | 50 | 1/15 | = | = | = |
| n5_f30 | 5 | 30 | 20 | 50 | 1/15 | = | = | = |
| n32_f30 | 32 | 30 | 20 | 50 | 1/15 | = | = | = |
| n33_f30 | 33 | 30 | 20 | 50 | 1/15 | = | = | = |
| n4_f2 | 4 | 2 | 20 | 50 | 1/15 | = | = | = |
| n5_f2 | 5 | 2 | 20 | 50 | 1/15 | = | = | = |
| fewer_than_k | 12 | 30 | 20 | 50 | 1/15 | = | = | = |
| pow2_64 | 64 | 10 | 20 | 50 | 1/15 | = | = | = |
| pow2_256 | 256 | 10 | 20 | 50 | 1/15 | = | = | = |
| pow2_1024 | 1024 | 10 | 20 | 50 | 1/15 | = | = | = |
| trees10 | 300 | 30 | 20 | 10 | 1/15 | = | = | = |
| k5 | 300 | 30 | 5 | 50 | 1/15 | = | = | = |
| k50 | 300 | 30 | 50 | 50 | 1/15 | = | = | = |
| prune0 | 300 | 30 | 20 | 50 | 0 | = | = | = |
| prune02 | 300 | 30 | 20 | 50 | 0.2 | = | = | = |

What the cases exercise:

- **duplicates**: 40 points five times each, interleaved, and 40 more copies of one of them, so 45 identical cells, more than k. The 20 nearest are all at distance 0 and go to the lower ids, so 25 cells miss themselves, as C3i says.
- **lattice, two_sites, ulp_outlier**: exact midpoints, so margins of exactly 0. In ulp_outlier, one cell sits a float ulp from 199 copies of another, in a coordinate whose odd mantissa makes the midpoint round onto the outlier. Splits that cannot separate them are retried.
- **exp_spacing, exp_spacing_dup**: exponentially spaced points (2^-118.75 to 2^117.8). Splits peel off a few cells, so they are kept imbalanced or randomised. Squared distances underflow and overflow in float, and the Angular distances of overflowing points are NaN: they test the libstdc++ heap moves.
- **collinear**: rank-one cells. **zero_rows**: zero vectors, which Angular's two-means skips and whose distance is 2.
- **n ≤ K, the first split, fewer cells than k**:
  - n1–n32 with f = 30 (K = 32), and n4 with f = 2: every root is a leaf;
  - n33 and n5 with f = 2: the first split;
  - `fewer_than_k`: FindNeighbors lowers k to 11, with its warning.
- **powers of two**: no Annoy code path branches on them (unlike `JavaRandom.nextInt`, M10), so these only confirm that nothing does.
- **trees10, k5, k50, prune0, prune02**: non-default parameters.
- **n1_f30**: for a single cell, FindNeighbors lowers k to 0. Its search returns no neighbours, and `sparseMatrix` then fails (`'i' and 'j' must not contain NA`). `scrust neighbors` refuses that case too. The index and items R built for it are identical.

## Benchmark

Both 79,409-cell objects, on their production instrumented runs' embeddings, on one node (Xeon Platinum 8480CL, 32 CPUs, job 60383132, [`m9_benchmark.sbatch`](../../validation/slurm/m9_benchmark.sbatch)):

- **R**: the workflow's FindNeighbors call alone (`neighbors_references.R --mode bench`).
- **Rust**: `scrust neighbors` at 1 and 32 threads, the whole command. Its *FindNeighbors s* is build + search + graphs.

Wall time and peak RSS are from `/usr/bin/time -v`.

| object | side | threads | FindNeighbors s | build s | search s | graphs s | whole process | peak RSS MiB | note |
|---|---|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | R | 1 | 27.2 |  |  |  | 0:37.37 | 1,617 | VmHWM 453 → 1,617 MiB across the call; in M8's chain 31.4 s and 32.7 s (with and without RunUMAP) |
| soupx_emptydrops_sscovelli | Rust | 1 | 20.1 | 3.66 | 15.49 | 0.99 | 0:24.76 | 264 | reading the embeddings and writing every output included |
| soupx_emptydrops_sscovelli | Rust | 32 | 4.2 | 3.70 | 0.41 | 0.09 | 0:06.80 | 452 | reading the embeddings and writing every output included |
| filtered_seurat_emptydrops_sscovelli | R | 1 | 27.8 |  |  |  | 0:33.66 | 1,640 | VmHWM 455 → 1,640 MiB across the call; in M8's chain 33.4 s and 31.7 s (with and without RunUMAP) |
| filtered_seurat_emptydrops_sscovelli | Rust | 1 | 20.4 | 3.72 | 15.67 | 0.98 | 0:24.80 | 262 | reading the embeddings and writing every output included |
| filtered_seurat_emptydrops_sscovelli | Rust | 32 | 4.4 | 3.73 | 0.57 | 0.09 | 0:07.25 | 420 | reading the embeddings and writing every output included |

- **Threads.** All five outputs (`knn_index.tsv`, `knn_distance.tsv`, `nn.mtx.gz`, `snn.mtx.gz`, `annoy.ann`) are byte-identical at 1 and 32 threads, and identical with and without `--annoy yes`; the index's sha256 equals R's saved index's (`threads.txt`).
- **Time.** The build is sequential, 3.7 s, because one generator runs through all 50 trees, as in RcppAnnoy. The search is the one-thread cost (15.5 s, about 0.2 ms per cell): it is bound by memory latency, visiting nodes scattered over a 60 MB index. It runs 27–38 times faster at 32 threads. The graphs take 0.1–1 s. Writing the outputs adds 2.4–4.4 s, mostly compressing the 3.2 million SNN entries.
- **Memory.** Rust holds the embeddings, the 60 MB index, the neighbour lists and the graphs: 262–264 MiB at one thread. At 32 threads, 420–452 MiB, mostly SNN pieces compressed in parallel. R's process peaks at 1.6 GiB, its FindNeighbors call adding 1.2 GiB to the 0.45 GiB of R, Seurat and the embeddings.

## Tests

`cargo test --workspace` passes in both profiles: 189 tests in 32 suites, none failing or ignored, in debug and in release, with no compiler warnings. 29 of them are new in M9:

- `crates/neighbors`, unit tests (19): Kiss64, the libstdc++ heap moves (including NaN priorities) and `partial_sort`, node layouts, `_split_imbalance`'s rounding, the Angular distance's arithmetic, a small index and its reload, the SNN graph against A·Aᵀ by brute force (several prunes, 1 and 3 threads), parameter parsing and refusals, the item rounding, fewer cells than k, a single cell, and thread independence;
- `tests/k0_kiss64.rs` (4): the K0 vectors;
- `tests/tier_t.rs` (3): K3–K6 on both base objects' committed stage references (M8), at 1, 3 and 7 threads;
- `tests/edge_cases.rs` (1): the 25 edge cases × 3 metrics against R's 732 digests, with each metric reaching each rare branch;
- `crates/scio` (1): the Matrix Market writer and reader;
- `crates/scrust/tests/cli.rs` (1): `scrust neighbors` and `scrust validate neighbors` end to end on a base object at 1 and 4 threads (identical files), and the refusals (hamming, dims beyond the PCA).

## Findings for the port

- **Angular nodes are 12 + 4f bytes, not 16 + 4f.** Angular has no offset, and its squared norm shares four bytes with the first child (a union). An item therefore holds its norm where a split node holds a child. Leaves hold f + 2 ids with every metric.
- **The 95% and 99% thresholds are crossed in float.**
  - `_split_imbalance` rounds the larger share to float before comparing it with the doubles 0.95 and 0.99.
  - An exact 95% split (19 of 20) gives 0.94999999 and passes.
  - An exact 99% split (99 of 100) gives 0.99000001 and counts as above 0.99, so it falls back to random sides.
- **The Angular distance is computed half in double.** The float `sqrt` of `pp * qq` is converted to double for the quotient and the difference, and the result is rounded to float. Computing it all in float, or with a double `sqrt`, would change distances in their last bits.
- **Equal priorities pop by the larger node id.** The search also depends on libstdc++'s heap moves wherever the order is not total. The overflow edge case produces NaN distances, and its index, neighbours and graphs are identical to R's only because the port copies those moves.
- **FindNeighbors fails on a single cell** (k = 0). The workflow's guard stops at 50 cells or fewer anyway.
- **The saved index searches as the built one does.** Annoy's `load()` finds 2q − 1 roots in an index whose roots are leaves (n ≤ K), not q; such roots hold every item, so the results do not change. FindNeighbors searches the in-memory index anyway, and R's checks confirm the saved index gives the same neighbours (stages.json, every reference).
- **Real data reaches the rare branches** (tables above). The tier-P runs make 2,743 retried splits, 581 kept imbalanced and 62 random-side fallbacks with Euclidean, and many more with Angular. Tier T reaches none of them. Exactness on tier P therefore tests the draw counting C7 worries about.
- **The text output spells every NaN alike.** NaN distances arise only from overflow, in the exp_spacing case with cosine. `knn_distance.tsv` writes them as `NaN`, which loses the sign of x86's default NaN, so `scrust validate neighbors` counts NaN = NaN as identical and says so in its report. cargo test compares those bits, through the arrays, and they are identical.

## Corrections to the plan

Corrected in [`docs/PLAN.md`](../PLAN.md) on this branch:

- **C3a, nodes:** 16 + 4f bytes for Euclidean and Manhattan, 12 + 4f for Angular (the norm in a union with the first child); bytes Annoy never writes are zero.
- **C3a, the split:** the retry and fallback thresholds as `_split_imbalance` rounds them. An exact 95% split passes; an exact 99% split falls back to random sides.
- **C3a, metrics:** how Manhattan and Angular build. Angular's normalised centroids, its skipped zero vectors, its margin without an offset, and its distance's float `sqrt` with a double quotient and difference.
- **C1 step 4:** FindNeighbors fails on a single cell (k = 0).
- **C3i:** the port copies libstdc++'s heap moves, which matters for NaN distances.
- **C3, commands and outputs:** `scrust neighbors` exists, with `--annoy yes|no` and `--threads`. Its knn tables carry a header and the barcode.

Plan facts the data confirmed as written: the draw counts of C2; `Kiss64Random` from the default seed in one stream through all trees; K = f + 2; node ids (items, then each tree's nodes in creation order, smaller side first, then the root copies); `two_means`' weighted distances and update; the offset; the flip at a margin of exactly 0; the queue's tie order; `partial_sort`'s tie order and the missing self of duplicated cells; the candidates sorted, deduplicated and measured in float; the float conversion once per value; `SCT_nn` and `ComputeSNN` (s/(k + (k − s)) in double, `<` prune, the diagonal kept, rows sorted within each column); and the Angular distances Seurat returns.

## Deviations and notes

- **`--annoy yes|no`** takes a value, because `scrust`'s argument parser reads `--name value` pairs. The plan says only "on request".
- **The Angular and Manhattan references are deliberately small** (the brief's "small R references"): tier T in both profiles, the small P object in both, and the largest object in one profile. They take 0.48 GB in `data/`, 0.38 GB of it the largest object's; the edge cases 0.04 GB; M9's Rust outputs 0.8 GB (the benchmark's duplicates replaced by their checksums). The 98 Euclidean comparisons cover all 34 objects in both profiles, and the build code paths that depend on size are metric-independent.
- **No R run for the remaining P objects' other metrics**, and none for Hamming, which Seurat would run on integers; the port refuses `annoy.metric = "hamming"`.
- **Committed fixtures:** `tests/fixtures/` grows from 23.07 MB to 23.46 MB. Only digests are committed for the edge cases; R's files are in `data/neighbors/edge/`.
- **The manifest** gains a "Neighbours" section from [`validation/R/manifest_neighbors.R`](../../validation/R/manifest_neighbors.R): 2,886 files and 115 run tables. It sits before the clustering section, so M8's and M2's generators leave it alone. Rust's validation outputs are not listed, as in M3.
- **Licences:** NOTICE now lists the code ported from Annoy (Apache-2.0, Spotify AB; the licence text in `third_party/licenses/Apache-2.0.txt`, the changes stated), RcppAnnoy, Seurat and libstdc++. `annoy.rs` and `kiss.rs` carry the modification notice.

## Reproducing

From the repository's top directory, after `source env.sh` and `cargo build --release -p scrust`:

```bash
sbatch validation/slurm/neighbors_references.sbatch edge
sbatch validation/slurm/neighbors_references.sbatch real filtered_seurat_tenx_test/instrumented-portable-sapphirerapids   # and the other runs
bash validation/slurm/fixtures_neighbors.sh
sbatch validation/slurm/m9_neighbors.sbatch filtered_seurat_tenx_test soupx_emptydrops_sscovelli   # any objects
bash validation/slurm/m9_edge.sbatch
sbatch validation/slurm/m9_benchmark.sbatch soupx_emptydrops_sscovelli filtered_seurat_emptydrops_sscovelli
Rscript validation/R/report_clustering_m9.R docs/validation/clustering_m9.tsv
Rscript validation/R/manifest_neighbors.R
```

One run by hand, on a committed fixture (a second, any machine):

```bash
f=tests/fixtures/clustering/filtered_seurat_tenx_test/reference
target/release/scrust neighbors --pca $f/pca --out /path/to/out --dims 1:30 --annoy yes
target/release/scrust validate neighbors --reference $f --candidate /path/to/out
```

The report and manifest scripts read plain files only; R runs with the hermetic flags of `validation/R/README.md`. SLURM logs of the M9 jobs: `scrust-m9-*.sbatch.log` in the submit directory (gitignored; job ids above and in the manifest).
