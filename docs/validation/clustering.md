# Clustering in Rust (Step 2; milestones M11 and M12)

M11 runs the workflow's clustering chain end to end in Rust, from the counts to the labels and the silhouette widths (docs/PLAN.md C1, C3, C5, C6): SCTransform at seed 12345, the workflow's PCA guard, RunPCA with 50 PCs, FindNeighbors (Annoy, k = 20, dims 1:30), FindClusters (Louvain, resolution 0.8, seed 12345, singletons grouped) and bluster's `approxSilhouette`. It chains the `sctransform` and `pca` crates of Step 1 with M9's `neighbors` and M10's `cluster`, and adds the silhouette widths, `scrust cluster --counts` with its staged starts, parameters, preset and TOML file, and `scrust validate cluster` for the whole chain. It also adds kernel-gated grading (A5b mechanism 3) to both validators and a kernel check to the fixture tests.

M12 adds the workflow's second score, neighbourhood purity: bluster's `neighborPurity`, with its search done exhaustively with knncolle's arithmetic and a certificate of where KMKNN, the algorithm R uses, could differ. Its outcome comes first ([M12 outcome](#m12-outcome)), its details after M11's ([M12: neighbourhood purity](#m12-neighbourhood-purity)).

This is Step 2's main report, as [doubletfinder.md](doubletfinder.md) is Step 1's; M13 to M17 extend it. The milestones before it have their own: [clustering_m8.md](clustering_m8.md) (the R reference data and the R-vs-R floor), [clustering_m9.md](clustering_m9.md) (Annoy and the SNN graph), [clustering_m10.md](clustering_m10.md) (FindClusters).

## M12 outcome

| Done when (plan C6's M12 row, and the brief) | Result |
|---|---|
| K9 on T and P, teacher-forced | **Exact.** `neighborPurity` on R's own embeddings and labels, in the 98 comparisons with instrumented runs (all 34 objects, both profiles): 1,033 checks identical: each cell's 50th distance, the radius, the neighbour counts, the neighbour lists in R's order where R kept them, the weights, sums, totals, purities and maxima, and R's `neighborhood_purity` column. BLAS-free, so identical on the other kernel's embeddings as well (713 checks in the 68 mismatched pairs). |
| K9 end to end, from the counts, under the reference's kernel; the files against the direct runs | **Exact.** 222 comparisons with every R run of the profile, all in exact mode, all passing: 1,831 K9 checks identical (the stage dumps, the scalars, `neighborhood_purity` and `purity_maximum` in the files), none differing or missing; every direct run's purities identical, with and without RunUMAP ([M12 results](#m12-results)). |
| The saved purities | **Identical in the production profile**, bit for bit, on all 34 objects. In the portable profile they are graded under the saved objects' kernel: identical except in the 28 portable runs M8 found differing from their saved objects in one or two cells, the median pair, which stay within the floor. |
| Everything M11 checked stays identical | **Yes.** All 31,269 checks of M11's 290 reports keep their status and counts in M12's: the labels and levels, the silhouette widths, the files, K1 to K8. |
| The certificate, and the decision point | **Fired; decided.** It flags the median pair of every odd-sized object whenever that pair lies within the radius (17 objects; 23 of the 68 runs; 0 or 1 ulp below the radius squared), and nothing else; R's KMKNN kept every flagged pair. The user chose to keep the exhaustive result and report the flagged pairs, graded "flagged, reported" (passing) and labelled `median_pair` or `other` ([The decision point](#the-decision-point)). |
| Independent of the thread count | **Yes.** 8 and 1 threads on all 34 objects in the portable profile and on the 27 tier-T objects in the production profile, 1, 8 and 32 threads on the 7 tier-P objects in the production profile: every file byte-identical, the K9 dumps and the certificate included ([M12 threads](#m12-threads)). |
| Time measured on the largest objects; the benchmark | **C5's criterion met** ([M12 benchmark](#m12-benchmark)). The whole run with both scores takes 2.4 to 5.6 times less time than R's at one thread, 8.2 to 15.8 times less at 32 threads, and 8 to 39 times less memory: 3.15 GiB at one thread and 3.19 GiB at 32 on the 79,409-cell objects (bound: 10 GiB). Purity alone takes 80.7 and 83.5 s on those at one thread, against R's 180.3 and 192.6 s, and 2.3 and 2.8 s at 32 threads; it adds nothing to the peak memory. |
| Kernel-gated grading of K9 | **Yes.** The 68 mismatched-kernel pairs are in tolerance mode and pass: 276 purity checks within `[clustering.floor.kernel]`, none beyond, the BLAS-free ones identical. |
| Documentation | This report; `[clustering.K9]`, the K9 entries of `[clustering.floor.kernel]` and M12's measurements in `validation/tolerances.toml`; M12's section of `validation/manifest.toml`; bluster's `neighborPurity`, BiocNeighbors and knncolle in `NOTICE`; the plan's corrections ([M12 corrections](#m12-corrections-to-the-plan)). |
| `cargo test --workspace` passes, debug and release, no warnings | **Yes**: 267 tests, 14 of them new ([M12 tests](#m12-tests)). Under the SapphireRapids kernel, the 14 tests that compare a BLAS-dependent result with a portable fixture stop at the kernel check, and the other 253 pass. |

## M11 outcome

| Done when (plan C6's M11 row, and the brief) | Result |
|---|---|
| K1–K8 met on T and P in both profiles | **Yes.** 222 comparisons of the Rust chain from the counts with every R run of its profile (the 27 tier-T and 7 tier-P objects; instrumented runs, direct runs with and without RunUMAP, the v5 runs; R on Sapphire Rapids and on Cascade Lake), all in exact mode, all passing: 20,383 checks identical, none differing or missing. K1 (SCTransform at seed 12345, through both assay paths) 4,925 checks identical; K2 (RunPCA with 50 PCs, `work = 66` up to 79,409 cells) 2,926; K3–K6 836; K7 3,696; K8 (the silhouette, `metadata.tsv`, `pca/`, the `result` numbers) 5,852. |
| The direct runs' files and the saved labels identical | **Yes.** Every direct run's files, with and without RunUMAP: `metadata.tsv`'s columns, the levels, `cluster_ids.txt`, `nclusters.txt`, `pca/`, both graphs, the `result` numbers. The labels and both columns' levels are identical to the saved objects' in all 68 runs (34 objects, both profiles), and so are the silhouette widths, bit for bit, in the production profile; in the portable profile they are graded under the saved objects' kernel and stay within the floor (118 checks, at most 3.3e-13 apart). |
| Independent of the thread count | **Yes.** 8 and 1 threads on all 34 objects in the portable profile and on the 27 tier-T objects in the production profile; 1, 8 and 32 threads on the 7 tier-P objects in the production profile: every file byte-identical, the stage dumps included, and `summary.json`'s result section ([Threads](#threads)). cargo test compares 1, 2 and 5 threads. |
| A kernel mismatch puts both validators' reports in tolerance mode, shown on deliberately mismatched pairs | **Yes.** `scrust validate cluster`: 68 pairs (each Rust run against the other profile's instrumented R run), all in tolerance mode with the reason given, all passing: the BLAS-free checks identical, 1,028 BLAS-dependent checks within `[clustering.floor.kernel]`, none beyond it. `scrust validate doubletfinder`: 20 pairs of M6's runs, all in tolerance mode, all passing (674 BLAS-dependent checks within `[doubletfinder.floor.kernel]`). The floor is now measured per stage on all ten objects, the two large ones from instrumented production R runs made for it ([DoubletFinder's two large objects](#doubletfinders-two-large-objects-per-stage)); with M2's floor, 16 of the first 18 pairs passed. cargo test shows both validators switching, and that a run without a recorded kernel is never graded as exact. With matching kernels Step 1's results are unchanged: 99 reports of M6's and M7's runs re-graded, every verdict and every check's status as before, and the two new comparisons with the large objects' production runs identical stage by stage. |
| The fixture tests stop with a message naming both kernels when run under another kernel | **Yes.** Under the SapphireRapids kernel (job 61108077, this report's tree), `cargo test --workspace --release` gives 13 failures, each a kernel check naming both kernels, and the other 240 tests pass: the BLAS-free, teacher-forced and Rust-against-Rust tests run under any kernel ([Tests](#tests)). The check replaces the silent skips. |
| The benchmark in `docs/validation/clustering.md` | **Yes; C5's criterion met** ([Benchmark](#benchmark)). At one thread Rust takes 2.7 to 5.7 times less time than R's direct run without RunUMAP (2.0 to 2.6 times less than R's chain alone, without its purity and writing), at 32 threads 8.6 to 17.4 times less, and 8 to 39 times less memory: 3.15 GiB at one thread and 3.19 GiB at 32 on the 79,409-cell objects (bound: 10 GiB), against R's 57.5 and 62.1 GiB. |
| Silhouette, teacher-forced | **Exact.** `approxSilhouette` on R's own embeddings and labels: 490 checks in the 98 comparisons with instrumented runs, all identical. |
| `scrust cluster --counts`, the staged starts, the parameters, the preset and the TOML file | **Built.** The four starting points write byte-identical files (cargo test); every A7b parameter has a flag and a TOML key, with Seurat's defaults; the refused code paths stop with a message; purity is reported as not computed yet (M12) and `--export` is refused (M13). |
| Documentation | This report; `[clustering.K1]`, `[clustering.K2]`, `[clustering.K8]`, `[clustering.kernel_gating]`, `[clustering.kernel_gating_shown]`, `[clustering.threads]`, `[clustering.performance]` and the floor's new entries in `validation/tolerances.toml`; M11's section of `validation/manifest.toml`; bluster in `NOTICE`; the plan's corrections ([below](#corrections-to-the-plan)). |
| `cargo test --workspace` passes, debug and release, no warnings | **Yes**: 253 tests ([Tests](#tests)). |

## Step 2 so far

| Milestone | What it established | Report |
|---|---|---|
| M8 | the R reference data: 27 tier-T and 7 tier-P objects, direct and instrumented runs in the portable (Haswell) and production (SapphireRapids) OpenBLAS profiles, each production run reproducing its saved object; stage dumps; the R-vs-R floor; the marker references; the instrumented optimizer | [clustering_m8.md](clustering_m8.md) |
| M9 | `crates/neighbors`: Annoy (RcppAnnoy 0.0.23) and the SNN graph, byte-identical to R on R's own embeddings | [clustering_m9.md](clustering_m9.md) |
| M10 | `crates/cluster`: Seurat's modularity optimizer, GroupSingletons, FindClusters' columns and levels, exact on R's own SNN graphs | [clustering_m10.md](clustering_m10.md) |
| M11 | the chain from the counts, `approxSilhouette`, `scrust cluster --counts` and its staged starts, kernel-gated grading in both validators, the benchmark | this report |
| M12 | `neighborPurity`: the exhaustive search with knncolle's arithmetic and its certificate, the weights and sums, purity in `scrust cluster` and its validation; at the decision point the user chose to keep the exhaustive result and report the flagged pairs (the median pair of every odd-sized object) | this report |

## What was built

| Where | M11 contents |
|---|---|
| `crates/cluster` | `chain.rs`: the chain as one library entry point (C3h), from the counts (`chain::from_counts`), SCTransform's `scale.data` (`from_scale_data`), the embeddings (`from_embeddings`) or the SNN graph (`from_snn`), with the workflow's PCA guard; an `Observer` sees each stage's result as it is made, so nothing is kept that the next stage does not need. `chain_params.rs`: every parameter A7b lists, with Seurat's defaults, the `workflow` preset, the TOML file and the refusals. `chain_output.rs`: `scrust cluster`'s outputs and the stage dumps. `metadata.rs`: `nCount_RNA` and `nFeature_RNA` as `CreateSeuratObject` computes them, and `--metadata`'s columns. `silhouette.rs`: bluster 1.20.0's `approxSilhouette`. `validate.rs`: `scrust validate cluster` extended to K1-K8, R's final files, the saved object and teacher-forced silhouette, kernel-gated. |
| `crates/scio` | `toml.rs`: a TOML reader (standard library only), for `--params` and the floors in `validation/tolerances.toml`. `stages.rs`: the instrumented R runs' stage-dump layout, read and written (moved from `doubletfinder`, which now uses it, so the chain and DoubletFinder share one implementation). `grade.rs`: kernel-gated grading (the mode, the classes, the floor's bounds, sign-aligned array statistics) and the floor measurement. `sha256.rs`: a streaming SHA-256 (scale.data's digest). |
| `crates/sctransform`, `pca`, `neighbors` | `output.rs`: each command's file layout and `summary.json` sections, written by the command and by the chain, so both write the same bytes. `record.rs` (sctransform, pca): the S2 and S3 stage dumps, shared by DoubletFinder's runs and the chain. |
| `crates/doubletfinder` | `validate.rs`: kernel-gated grading. `record.rs`, `stages.rs`: on `scio::stages`. |
| `crates/rcompat` | `blas::fixture_kernel` and `blas::require_fixture_kernel`: the fixture tests' kernel check. |
| `crates/scrust` | `scrust cluster` (`--counts`, `--sct`, `--pca`, `--snn`; `--preset`, `--params`, every parameter's flag; `--dump-stages`, `--keep-stages`, `--metadata`, `--annoy`); `scrust validate cluster` (`--saved`, `--saved-kernel`); `scrust validate floor`; `scrust sctransform`, `pca` and `neighbors` carry the counts' per-cell columns forward. |
| `validation/slurm` | `m11_floor.sbatch`, `m11_cluster.sbatch`, `m11_doubletfinder.sbatch`, `m11_df_references.sbatch`, `m11_kernel_tests.sbatch`, `m11_benchmark.sbatch`, `m11_threads.sh`, `m11_report.sh`, `m11_manifest.sh`. |
| `tests/fixtures/kernels.tsv` (new, 3 KB) | the OpenBLAS kernel each BLAS-dependent fixture without its own run record was made under. |

The workspace still uses only the Rust standard library and the conda env's OpenBLAS.

### The chain

`scrust cluster --counts <10x dir> --preset workflow --out <dir>` runs, in one process:

1. **SCTransform** (`sctransform::start`, then `finish`): `set.seed(12345)` in each method along the call chain, vst's setup, fits, regularisation and residuals, the 3,000 variable features and `scale.data`, as in Step 1 with the seed changed. The counts are freed once the setup has kept the kept genes' rows; `nCount_RNA` and `nFeature_RNA` are taken before.
2. **The workflow's guard** (`require_min_cells_for_pca`, on with the preset): 50 cells or fewer stop the run with the workflow's message.
3. **RunPCA** (`pca::run_pca`): `set.seed(12345)`, irlba with `nv = 50` (`work = 66`, capped at min(m, n) − 2 for small objects), the embeddings, loadings and standard deviations. SCTransform's working data is freed first, `scale.data` after.
4. **FindNeighbors** (`neighbors::find_neighbors`, M9): dims 1:30, k = 20, 50 trees, prune 1/15.
5. **FindClusters** (`cluster::find_clusters`, M10): resolution 0.8, seed 12345, 10 starts, 10 iterations, GroupSingletons.
6. **Silhouette** (`approx_silhouette`) on the neighbour search's dims (1:30), as `add_silhouette_to_metadata` calls bluster.

RunUMAP is left out: it changes nothing downstream (C1; M8 checked it on all 34 objects). The random-number stream is R's (`rcompat::RRng`), reseeded by SCTransform and RunPCA as in R; Annoy and the optimizer have their own generators, and GroupSingletons reseeds with `set.seed(1)`.

**`approxSilhouette`** (bluster 1.20.0, as the workflow calls it on `emb[, 1:30]` and `seurat_clusters`): the clusters are the levels present, in level order; each centroid is R's `colMeans` and each within-cluster variance `sum(colMeans(sweep(x, 2, centroid)^2))`, both accumulated in 80-bit long double as R's `do_colsum` and `rsum` do (`rcompat::F80`); each cell's root-mean-square distance to cluster i is `sqrt(colSums((t(x) - centroid)^2) + var_i)`, again with `colSums` in long double; the closest other cluster wins by strict `<` in level order; the width is `(other - self) / pmax(other, self)`. With fewer than 2 cells or clusters the workflow leaves every width NA, and so does the port. No BLAS is involved.

### Parameters

Every parameter A7b lists has Seurat 5.5.1's default and can be set by a flag or in the TOML file given by `--params`. They apply in this order: the defaults, the preset (`--preset`, else the file's `preset`), the file's values, the flags. Numbers are read as R reads them (`as.numeric`, R's own decimal parser), and a fraction written as text (`"1/15"`) is divided in double.

| stage | parameter (Seurat) | default | flag | TOML (`[section] key`) |
|---|---|---|---|---|
| SCTransform | `seed.use` | 1448145 | `--sct-seed` | `[sctransform] seed` |
| | `ncells` (fitting cells) | 5000 | `--ncells` | `ncells` |
| | `variable.features.n` | 3000 | `--variable-features` | `variable_features` |
| | `clip.range` | ±sqrt(N/30) | `--clip-range lo,hi` | `clip_range = [lo, hi]` |
| RunPCA | `npcs` | 50 | `--npcs` | `[pca] npcs` |
| | `seed.use` | 42 | `--pca-seed` | `seed` |
| | `weight.by.var` | TRUE | `--weight-by-var` | `weight_by_var` |
| | the workflow's guard | off | `--pca-guard none\|true\|N` | `guard` |
| FindNeighbors | `dims` | 1:10 | `--dims` | `[neighbors] dims` |
| | `k.param` | 20 | `--k` | `k` |
| | `n.trees` | 50 | `--trees` | `trees` |
| | `annoy.metric` | euclidean | `--metric` | `metric` |
| | `prune.SNN` | 1/15 | `--prune` | `prune` |
| | AnnoySearch's `search.k` | −1 | `--search-k` | `search_k` |
| FindClusters | `resolution` (several allowed) | 0.8 | `--resolution` | `[cluster] resolution` |
| | `algorithm` | 1 | `--algorithm` | `algorithm` |
| | `n.start` | 10 | `--starts` | `starts` |
| | `n.iter` | 10 | `--iterations` | `iterations` |
| | `random.seed` | 0 | `--seed` | `seed` |
| | `modularity.fxn` | 1 | `--modularity` | `modularity` |
| | `group.singletons` | TRUE | `--no-group-singletons` | `group_singletons` |
| | `graph.name` | SCT_snn | `--graph-name` | `graph_name` |
| scores | which (`add_silhouette_to_metadata`) | silhouette, purity | `--stats` | `[stats] compute` |
| | purity's `k`, `weighted` | 50, TRUE | `--purity-k`, `--purity-weighted` | `purity_k`, `purity_weighted` |

**`--preset workflow`** sets the workflow's calls: seed 12345 in SCTransform, RunPCA and FindClusters, `dims = 1:30`, and the PCA guard. A file such as

```toml
preset = "workflow"
[neighbors]
k = 30
prune = "1/15"
[cluster]
resolution = [0.4, 0.8, 1.2]
```

runs the workflow's chain with k = 30 and three resolutions.

**Refused** with a message, as `SctParams::check` does: `vars.to.regress`, `vst.flavor` other than v2, `method` other than glmGamPoi_offset, `do.scale = TRUE`, `nn.method = "rann"` or `"hnsw"`, `l2.norm = TRUE`, `annoy.metric = "hamming"`, Leiden (`algorithm = 4`), `seed.use = NULL`, and any key the file has that is not a parameter. Values that keep Seurat's code path (another k, resolution, seed, metric, ...) are plumbed through, with their defaults and refusals tested; checking them against R is M16's parameter grid.

**Purity** (M12): computed by default, as in the workflow, with bluster's `k = 50` and weights ([M12](#m12-neighbourhood-purity); M11's binary reported it as not computed yet). **`--export`** (M13) is refused with a message naming M13.

### Outputs and routes

`scrust cluster` writes, as C3 lists them: `seurat_clusters.tsv` and `seurat_clusters_levels.txt`, `clusters.tsv` with each column's levels, `cluster_ids.txt`, `nclusters.txt`, `optimizer_ids.tsv` (M10); `metadata.tsv` with `barcode`, `nCount_RNA`, `nFeature_RNA`, each `res.` column, `seurat_clusters`, `silhouette_width` and `silhouette_other`, `neighborhood_purity` and `purity_maximum` (M12), then `--metadata`'s columns; `pca/` in `scrust pca`'s layout and `neighbors/` in `scrust neighbors`' (`annoy.ann` with `--annoy yes`); with `--keep-stages`, `sct/` in `scrust sctransform`'s layout and every start's partition; with `--dump-stages`, `stages/` in the instrumented R runs' layout (purity's neighbour lists with `--purity-lists`, M12). `summary.json`'s `result` holds FindClusters' result (M10), each stage's result, R's digests of `.Random.seed` after SCTransform and after RunPCA (`rng_after_sha256`) and of `scale.data` in R's layout (`scale_data_sha256`), and the silhouette's and the purity's summaries (the purity's with its certificate, M12); `run` holds the parameters, the inputs, the files, each stage's time and peak RSS, and the environment fingerprint.

**The staged starts.** `--sct` reads `scrust sctransform`'s output, `--pca` `scrust pca`'s (or an R run's `pca/`), `--snn` `scrust neighbors`' (or an R run's `graphs/`), with `--pca` for the scores. `scrust sctransform` now writes the counts' per-cell columns (`cell_metadata.tsv`) and `scrust pca` and `scrust neighbors` carry them forward, so every route writes the same `metadata.tsv`. A staged start copies the stage outputs it was given (`pca/`, `neighbors/`, and with `--keep-stages` `sct/`) byte for byte, with their `summary.json`'s result. cargo test runs the four routes on a tier-T object and finds every output file byte-identical, and FindClusters' and the silhouette's results identical ([Tests](#tests)).

**Memory.** Only what the next stage needs is kept: the counts until SCTransform's setup, the kept genes' rows until `scale.data` is made, `scale.data` until the PCA is done, then the embeddings and the graphs.

## How the comparison works

`scrust validate cluster --reference <R run> --candidate <Rust run> [--optimizer <dir>] [--saved <dir>] [--panel <tsv> --object <name>]` compares what both sides hold:

- **K1-K9 from the stage dumps**, when both have `stages/` (an instrumented R run; a Rust run with `--dump-stages`): every block of every array R dumped (S2a-S2d, S3, K3, K5, K7, K8, and since M12 K9, its neighbour lists cell by cell), every `.Random.seed` checkpoint both recorded, every scalar both recorded, and Annoy's index file (K4), byte for byte. For the light-level objects, which kept no copy of `scale.data`, its digest in R's layout. R's arrays of later milestones (K10: the corrected counts, M13) are listed as not compared; so are R's checkpoints after RunUMAP (its draws are not reproduced, and nothing later reads them: GroupSingletons reseeds) and the scalars of R's own structural checks.
- **K6 and the files**, against R's final files: both graphs (`graphs/` against `neighbors/`: dimensions, every stored entry in order, every value); `metadata.tsv`'s `nCount_RNA`, `nFeature_RNA`, `silhouette_width` and (M12) `neighborhood_purity`, and `silhouette_other` and `purity_maximum` against R's `K8_silhouette_other` and `K9_maximum`; `pca/` (the embeddings, the loadings and their features, the standard deviations, the total variance); the `result` section's numbers (cells, clusters and their sizes, `scale.data`'s shape and digest, the variable features, the optimizer's numbers, R's digests of `.Random.seed` after SCTransform and RunPCA). The columns of later milestones (`nCount_SCT`, `nFeature_SCT`: M13) are left out, and the report says so.
- **K7** (M10): the settings, the optimizer's ids, the final ids, the labels and levels, what FindClusters printed, the singletons and GroupSingletons' messages, R's `.Random.seed` after GroupSingletons, every start against the instrumented optimizer, `seurat_clusters` and the workflow's files, and for tier P the workflow job's logged numbers.
- **The saved object** (`--saved`, M8's export): the labels and both columns' levels, and the silhouette widths and (M12) purities, graded under the saved object's kernel (SapphireRapids: every workflow job ran on Sapphire Rapids, where OpenBLAS detects it; M8).
- **Silhouette, teacher-forced**: `approxSilhouette` on R's own embeddings (`pca/embeddings.npy`) and labels, against R's `K8_*` arrays and `silhouette_width` column. BLAS-free, so exact on any machine. Since M12 also `neighborPurity`, against R's `K9_*` arrays, its radius and `neighborhood_purity` column, with the certificate of that search.
- **The purity certificate** (M12): the candidate's flagged pairs, reported as "flagged, reported" ([M12](#m12-neighbourhood-purity)).

The report (JSON and markdown) gives every check's stage, class, criterion, status and counts, and the mode.

### Kernel-gated grading

OpenBLAS picks its kernel per CPU, and BLAS products differ in their last bits between kernels (A5b). So both validators now grade by the kernels the two runs recorded (`run.fingerprint.openblas_corename`):

- **Exact mode** when they match: every check must be identical.
- **Tolerance mode** when they differ, and whenever a run records none (a run without a recorded kernel is never graded as exact): each BLAS-free check stays exact; each BLAS-dependent check is graded against the R-vs-R floor, `[clustering.floor.kernel]` or `[doubletfinder.floor.kernel]` in `validation/tolerances.toml` (embedded in the binary when it is built). A check passes when it stays within every bound of its array's entry (the keys starting with `max_`, and `values_differing` and `cells_differing`, which are totals over the measured objects); an array the floor does not list stays exact, because the floor found it identical across kernels. PCA components may flip sign between kernels, so embeddings, loadings, Annoy's items and the silhouette's centroids are compared after aligning each column's sign, as the floor was measured. A digest (scale.data's, for the light objects) cannot be graded against bounds: in tolerance mode a differing digest is reported, not graded.
- **The report says** which mode, and why (for example "the reference recorded OpenBLAS kernel SapphireRapids, the candidate Haswell"), and gives each check's class and criterion.

**Classes.** A check is BLAS-dependent when any of its inputs passed through a BLAS product in the candidate's run:

| validator | candidate | BLAS-free | BLAS-dependent |
|---|---|---|---|
| cluster | from the counts | the settings, the barcodes, `nCount_RNA`, `nFeature_RNA`; S2a; SCTransform's random-number states and those at RunPCA's and irlba's entry (their draw counts are fixed); teacher-forced silhouette | everything from S2b on: S2b (glmGamPoi's Cox-Reid term calls `dgemv`, `nlminb` calls `ddot` and `dnrm2`), S2c, S2d and scale.data's digest, S3 and irlba's exit state, K3-K8, the labels, the graphs, the files, the `result` numbers, the saved object's checks |
| cluster | from `--sct` | as above up to SCTransform | RunPCA and everything after it |
| cluster | from R's embeddings or SNN graph | every check | none |
| doubletfinder | from the counts | S1 and S2a, and the random-number states of the sampling and of SCTransform's fits, in the contexts whose doublets come from a fresh `set.seed(12345)` (the first sweep iteration and the final run); S1 redrawn from R's checkpoints; S6 and S7 recomputed from R's S5 | S2b-S2d, S3, S4-S7, the other sweep iterations (their doublets are drawn after the previous irlba), irlba's random-number states, the final files, the `result` section, the saved output's checks |
| doubletfinder | teacher-forced (`--from-stages`) | every check | none |

The classification is recorded in `tolerances.toml` (`[clustering.kernel_gating]`, `[doubletfinder.kernel_gating]`) and in every report.

## Results

Label `m11-eedf376`: binary `data/clustering/bin/scrust-m11-eedf376` (sha256 `1bfbd2d1...`, built from commit eedf376), [`m11_cluster.sbatch`](../../validation/slurm/m11_cluster.sbatch), 8 threads. The portable profile ran with `OPENBLAS_CORETYPE=Haswell` on Cascade Lake nodes (Xeon Platinum 8268; jobs 60803803 and 60803842), and its references were made on Sapphire Rapids and Cascade Lake: the same kernel gives the same bits on both CPUs (M8). The production profile ran with the kernel detected, SapphireRapids, on Sapphire Rapids (Xeon Platinum 8480CL, partition `test`; jobs 60803843 and 60803844). Each run dumped its stages at the level of M8's instrumented runs (full on tier T, standard or light on tier P) and was compared with every R run of its profile, the instrumented optimizer, the saved object and, on tier P, the workflow job's logged numbers. One row per comparison: [`clustering_m11.tsv`](clustering_m11.tsv) ([`m11_report.sh`](../../validation/slurm/m11_report.sh)).

### Every comparison

| tier | profile | comparison | comparisons | passed | mode | checks identical | within floor | differing | missing |
|---|---|---|---|---|---|---|---|---|---|
| T | portable | the profile's R runs | 108 | 108 | exact | 10,470 | 108 | 0 | 0 |
| T | production | the profile's R runs | 82 | 82 | exact | 6,707 | 0 | 0 | 0 |
| P | portable | the profile's R runs | 10 | 10 | exact | 1,316 | 10 | 0 | 0 |
| P | production | the profile's R runs | 22 | 22 | exact | 1,890 | 0 | 0 | 0 |
| T | portable | the other kernel's instrumented run | 27 | 27 | tolerance | 3,462 | 396 | 0 | 0 |
| T | production | the other kernel's instrumented run | 27 | 27 | tolerance | 3,462 | 396 | 0 | 0 |
| P | portable | the other kernel's instrumented run | 7 | 7 | tolerance | 901 | 118 | 0 | 0 |
| P | production | the other kernel's instrumented run | 7 | 7 | tolerance | 901 | 118 | 0 | 0 |

*R runs*: per tier-T object, portable: direct and instrumented on Sapphire Rapids and Cascade Lake (4); production: direct, direct without RunUMAP and instrumented (3), and for the emptyDrops object the v5 run (4). Per tier-P object, portable: the instrumented run (1), and for the small object also the direct run and both runs on Cascade Lake (4); production: as tier T. *Checks* exclude those that do not apply (a direct run has no stage dumps; later milestones' columns).

**222 comparisons with the profile's own R runs, all exact, 20,383 checks identical, 1.11 billion values compared.** The only checks that are not identical are the portable runs' silhouette widths against the saved objects (118): the saved objects were made under SapphireRapids, so these checks are graded under that kernel, in tolerance mode, and every one stays within the floor: at most 2.4e-14 apart on tier T and 3.3e-13 on tier P, as M8 found between R's own portable runs and the saved objects. In the production profile every check against the saved object is identical, the silhouette widths bit for bit, on all 34 objects.

By check, over the 222 comparisons: K1 (SCTransform) 4,925 identical; K2 (RunPCA) 2,926; K3-K6 (Annoy's items and index, the neighbours, both graphs) 836; K7 (FindClusters, the workflow's files and logged numbers) 3,696; K8 and the files (silhouette, metadata, `pca/`, the `result` numbers; the validation binary files the stage dumps' Annoy, neighbour and optimizer scalars here too) 5,852; the saved objects 1,214 identical and 118 within the floor; silhouette teacher-forced on R's embeddings 490 (every comparison with an instrumented run); the settings and the barcodes 444.

### Tier T

| object | cells | clusters | portable: R runs | portable: checks identical / within floor | production: R runs | production: checks identical | threads (8 and 1; portable; production) | mismatched kernels (both directions) |
|---|---|---|---|---|---|---|---|---|
| posthocfilt_threshold_doubletfinder_decontx_raw_tenx_test | 80 | 2 | 4 | 386 / 4 | 3 | 242 | identical; identical | 2 of 2 pass, 28 checks within floor |
| posthocfilt_threshold_scdblfinder_decontx_raw_tenx_test | 87 | 2 | 4 | 386 / 4 | 3 | 242 | identical; identical | 2 of 2 pass, 28 checks within floor |
| posthocfilt_threshold_doubletfinder_decontx_raw_emptydrops_test | 92 | 3 | 4 | 386 / 4 | 3 | 242 | identical; identical | 2 of 2 pass, 28 checks within floor |
| posthocfilt_threshold_scdblfinder_decontx_raw_emptydrops_test | 112 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 36 checks within floor |
| posthocfilt_threshold_doubletfinder_cellbender_fromraw_test | 135 | 2 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 32 checks within floor |
| posthocfilt_threshold_doubletfinder_soupx_emptydrops_test | 137 | 2 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| posthocfilt_threshold_doubletfinder_soupx_tenx_test | 142 | 2 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 32 checks within floor |
| posthocfilt_threshold_scdblfinder_soupx_tenx_test | 154 | 2 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 32 checks within floor |
| posthocfilt_threshold_scdblfinder_soupx_emptydrops_test | 161 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 32 checks within floor |
| posthocfilt_threshold_scdblfinder_cellbender_fromraw_test | 163 | 2 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| doubletfinder_decontx_raw_emptydrops_test | 237 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 30 checks within floor |
| doubletfinder_soupx_emptydrops_test | 237 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 30 checks within floor |
| doubletfinder_cellbender_fromraw_test | 247 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 30 checks within floor |
| scdblfinder_decontx_raw_emptydrops_test | 262 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 32 checks within floor |
| scdblfinder_soupx_emptydrops_test | 264 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| decontx_raw_emptydrops_test | 267 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| filtered_seurat_emptydrops_test | 267 | 4 | 4 | 388 / 4 | 4 | 392 | identical; identical | 2 of 2 pass, 28 checks within floor |
| soupx_emptydrops_test | 267 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| doubletfinder_decontx_raw_tenx_test | 273 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 30 checks within floor |
| doubletfinder_soupx_tenx_test | 276 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| cellbender_fromraw_test | 278 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| scdblfinder_cellbender_fromraw_test | 278 | 4 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| scdblfinder_soupx_tenx_test | 292 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| scdblfinder_decontx_raw_tenx_test | 296 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| decontx_raw_tenx_test | 300 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| filtered_seurat_tenx_test | 300 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |
| soupx_tenx_test | 300 | 3 | 4 | 388 / 4 | 3 | 243 | identical; identical | 2 of 2 pass, 28 checks within floor |

*Checks* sum over the R runs. *Threads*: the run again at one thread, every file byte-identical (62 files: the outputs, `sct/`, every start's partition and the stage dumps), and `summary.json`'s result section identical.

### Tier P

| object | cells | clusters | portable: R runs, checks | production: R runs, checks | portable run (8 threads, Cascade Lake) | production run (8 threads, Sapphire Rapids) | mismatched kernels: portable, production (checks within floor) |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 15 | 4: 398 identical, 4 within floor | 3: 248 identical | 4.9 s, 0.18 GiB | 3.7 s, 0.18 GiB | pass (17), pass (17) |
| soupx_tenx_10xmousebrain | 11,831 | 21 | 1: 153 identical, 1 within floor | 3: 248 identical | 19.9 s, 1.44 GiB | 16.1 s, 1.44 GiB | pass (16), pass (16) |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 21 | 1: 153 identical, 1 within floor | 3: 248 identical | 19.6 s, 1.44 GiB | 16.0 s, 1.44 GiB | pass (17), pass (17) |
| soupx_emptydrops_4500871lane | 26,598 | 31 | 1: 153 identical, 1 within floor | 3: 248 identical | 37.3 s, 1.48 GiB | 28.3 s, 1.48 GiB | pass (18), pass (18) |
| soupx_emptydrops_sscovelli | 79,409 | 36 | 1: 153 identical, 1 within floor | 3: 248 identical | 105.4 s, 3.16 GiB | 72.1 s, 3.16 GiB | pass (17), pass (17) |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 37 | 1: 153 identical, 1 within floor | 4: 402 identical | 104.3 s, 3.16 GiB | 73.7 s, 3.16 GiB | pass (17), pass (17) |
| cellbender_fromraw_ghirsutum2 | 59,323 | 22 | 1: 153 identical, 1 within floor | 3: 248 identical | 107.4 s, 2.66 GiB | 69.4 s, 2.66 GiB | pass (16), pass (16) |

The runs' time and peak RSS include their stage dumps (standard: scale.data; light: none of the residual-sized matrices). Every tier-P object exercises a code path tier T cannot: SCTransform's 2,000-gene sampling (all seven) and 5,000-cell sampling (all but the small object), irlba's working dimension of 66 with more restarts, Annoy's deeper trees, and GroupSingletons on 1, 18, 22 and 40 cells without SNN edges (soupx_emptydrops_4500871lane, soupx_emptydrops_sscovelli, filtered_seurat_emptydrops_sscovelli, cellbender_fromraw_ghirsutum2), whose `sample()` draws after `set.seed(1)` match R's `.Random.seed`. The workflow jobs' logged numbers (nodes, edges, maximum modularity, communities, singletons, final clusters) are identical for all seven. **K1** on the three light-level objects compares R's digest of scale.data, equal in both profiles. **The v3 assay path**: the emptyDrops objects' instrumented runs (v3) and their v5 runs both compare identical, as M8 found them identical to each other.

### Kernel-gated grading, shown on mismatched kernels

Each Rust run was also compared with the other profile's instrumented R run: a Haswell Rust run against R's SapphireRapids run, and a SapphireRapids Rust run against R's Haswell run, 68 pairs. Every report says "tolerance mode: the reference recorded OpenBLAS kernel SapphireRapids, the candidate Haswell" (or the reverse), and every pair passes: the BLAS-free checks are identical (S2a, SCTransform's random-number states, the settings, the counts' columns, teacher-forced silhouette), and of the BLAS-dependent ones 1,028 differ within `[clustering.floor.kernel]` and the rest are identical (the labels and levels, the neighbour ids, both graphs, the optimizer's ids and every start's modularity, the printed numbers, GroupSingletons' states). Within the floor, per pair: the S2d residual statistics (tier P), the PCA's embeddings, loadings, standard deviations and irlba's d, irlba's iterations (5 objects), Annoy's items and index file, a few neighbour distances, the silhouette widths, centroids and variances, and their files (`pca/`, `metadata.tsv`'s `silhouette_width`). On the light-level objects scale.data's digests differ between kernels; the reports list the digest as not gradable.

**DoubletFinder.** [`m11_doubletfinder.sbatch`](../../validation/slurm/m11_doubletfinder.sbatch) re-graded Step 1's own runs with the kernel-gated validator, and compared each new report with the old one check by check. It ran twice. The first run used M2's floor (job 60804010, binary `scrust-m11-eedf376`, `data/doubletfinder/compare/m11_verdicts_m11-eedf376.tsv`). The second used the floor measured on all ten objects, once the two large objects' instrumented production R runs existed (job 61104652, binary `scrust-m11-63e7b3b`, `m11_verdicts_m11-63e7b3b.tsv`). The second run:

| runs | comparisons | mode | verdicts | checks | checks whose status changed |
|---|---|---|---|---|---|
| M6 (m6-bc1b48e), every R run of the run's profile, 10 objects | 49 | exact | 49 pass: 47 as before, and the 2 new comparisons with the large objects' production instrumented runs | 15,995 | 0 |
| M6, the saved outputs | 20 | (not gated) | 20 pass, as before | 140 | 0 |
| M7 (m7-b8c8a8c), the saved outputs, 30 objects (6 samples × 5 arms, sscovelli and amellifera included) | 30 | (not gated) | 30 pass, as before | 270 | 0 |
| M6, the other kernel's instrumented run | 20 | tolerance | 20 pass | 11,896 | (new) |

With the same kernel, every verdict and every check's status is the one Step 1's validator gave. The kernel check itself is identical as before; in tolerance mode it is listed, not failed. The new saved-output reports of M6's runs have one check more than M6's own, "singlets: missing": M6b added that check after M6's label, whose runs wrote no `singlets/`. A missing check does not fail a DoubletFinder report, then or now. In tolerance mode, 674 BLAS-dependent checks differ within `[doubletfinder.floor.kernel]`: S2d's and S3's arrays, and S4's digest in cellbender_fromraw_ghirsutum2's final run.

With M2's floor, the first run gave the same verdicts except for two pairs: Rust's production runs of the two large objects against R's portable runs differed beyond that floor in three checks.

- soupx_emptydrops_sscovelli, `sweep_pN0.1`: the S3 embeddings by 3.665e-11 and the loadings by 8.370e-14, against bounds of 3.6e-11 and 8.1e-14.
- cellbender_fromraw_ghirsutum2, the final run: S4's digest for 3 cells, against "identical".

These are R's own differences between the kernels on these objects, which M2 had not measured per stage ([below](#doubletfinders-two-large-objects-per-stage)). With the floor measured on them, the three checks are within it and nothing else changes: every other check of the first run's 85 M6 reports has the same status in the second. The two new pairs, Rust's portable runs of the large objects against the new production R runs, pass too.

### Threads

Every output is byte-identical at one thread and at many:

| tier | profile | runs | threads | files compared per run | differing |
|---|---|---|---|---|---|
| T | portable | 27 | 8 and 1 (m11_cluster.sbatch's second run) | 62: the outputs, `sct/`, every start's partition, the stage dumps | 0 |
| T | production | 27 | 8 and 1 | 62 | 0 |
| P | portable | 7 | 8 (label m11-eedf376) and 1 (m11-eedf376-t1, job 60805276), [`m11_threads.sh`](../../validation/slurm/m11_threads.sh) | 55-56: the outputs and the stage dumps | 0 |
| P | production | 7 | 1 and 32 (the benchmark), and 8 (the validation run) | 17: the outputs | 0 |

`summary.json` is compared by its `result` section (its `run` section records the run). The tier-P portable runs at one thread were validated as the 8-thread runs were, and their 17 reports are identical to the 8-thread runs' check by check, counts and details included. cargo test compares 1, 2 and 5 threads on a tier-T object, and the four routes at 2 threads with the one-process route.

Where the work is split, each result is placed by index, never by finishing order: SCTransform's fitting bins and residuals (by gene), irlba's two large products (only where OpenBLAS's kernels keep R's bits, M5), Annoy's search (by cell; the build is sequential, as in RcppAnnoy), the SNN graph (by cell), and FindClusters' resolutions. The silhouette is sequential (under 4 s on 79,409 cells).


## The R-vs-R floor, measured again

Tolerance mode needs a bound for every BLAS-dependent array it may meet. M8's floor table listed the arrays it summarised; M11 measured the floor again with the statistics the grader itself uses (`scrust validate floor`, [`m11_floor.sbatch`](../../validation/slurm/m11_floor.sbatch), job 60798602): the instrumented R runs of all 34 clustering objects and of DoubletFinder's eight objects with runs in both profiles, Haswell against SapphireRapids on Sapphire Rapids, block by block (`data/clustering/floor_kernel_m11/`, `data/doubletfinder/compare/floor_kernel_m11/`). DoubletFinder's was then measured on all ten objects, the two large ones included ([below](#doubletfinders-two-large-objects-per-stage)).

- **It reproduces both earlier floors.** DoubletFinder's, to every printed digit of M2's table (S2d's residual statistics 5.68e-14 and 1.36e-16, scale.data 2.13e-14 and 1.18e-15, the embeddings 3.55e-11 and 2.41e-13, the loadings 8.03e-14 and 2.97e-13, irlba's d 7.50e-12 and 3.80e-15, the standard deviations 6.75e-14 and 3.84e-15); the clustering entries M8 listed likewise.
- **It adds what M8 did not tabulate,** now in `[clustering.floor.kernel]`: irlba's d (6.94e-12, 3.56e-14 of the largest value), the loadings (1.2e-7 absolute, on component 50 of the smallest tier-T objects; 1.6e-12 on tier P), the silhouette's centroids (6.68e-11, after aligning signs) and within-cluster variances (1.75e-9, 1.67e-13 scaled), Annoy's index file (its bytes differ in all 34 objects; the neighbour ids and both graphs are identical in all 34), and irlba's iterations and products (below). DoubletFinder's floor gains irlba's d, from M2's own measurement.
- **Two of M8's rounded values were below the data,** by less than their last digit: S2d's residual statistics' scaled bound (3.2e-17, measured 3.204e-17) and the neighbour distances' (1.9e-6, measured 1.907e-6). They are now 3.21e-17 and 1.91e-6; every entry is the maximum over the objects, rounded up.
- **Identical across kernels in every object:** everything else, including the total variance and every scalar of SCTransform, Annoy and the optimizer, the neighbour ids, both graphs and the labels.

**irlba's restarts depend on the kernel.** In 5 of the 34 clustering objects (small post-hoc tier-T objects, 80 to 163 cells), irlba runs one or two restarts more or fewer under one kernel than the other (`irlba_iter` differs by up to 2, `irlba_mprod` by up to 12). Its draws do not change (no restart draws a new vector, so `.Random.seed` is identical), and its results agree to the floor's bounds. M8's floor compared stage arrays and random-number states, not these scalars; the floor now lists them.

### DoubletFinder's two large objects, per stage

M2 measured DoubletFinder's floor per stage on eight objects. The two large ones had only a portable instrumented run and a production direct run, so their floor was measured end to end only: their final files are identical across kernels. These are soupx_emptydrops_sscovelli (79,409 cells) and cellbender_fromraw_ghirsutum2 (59,323). M11's first re-grading in tolerance mode (job 60804010) flagged both against a floor never measured at their size. M11 then made the two missing R runs, checked them as references, and measured the floor again on all ten objects.

**The runs.** [`reference_doubletfinder.sbatch`](../../validation/slurm/reference_doubletfinder.sbatch) `production instrumented` made them at the light dump level of the objects' portable instrumented runs, one object per job, the two in parallel on one `test` node (Xeon Platinum 8480CL, holy8a26501). `OPENBLAS_CORETYPE` was unset, and the job's `scrust fingerprint` reported the SapphireRapids kernel before R started.

| object | job | wall | peak RSS | the portable instrumented run (M2) |
|---|---|---|---|---|
| soupx_emptydrops_sscovelli | 61097286 | 35 min 2 s | 127.5 GiB | 36 min 59 s, 127.5 GiB |
| cellbender_fromraw_ghirsutum2 | 61097287 | 35 min 25 s | 153.5 GiB | 38 min 48 s, 153.1 GiB |

**Checked as references** ([`m11_df_references.sbatch`](../../validation/slurm/m11_df_references.sbatch), job 61102374, binary `scrust-m11-eedf376`):

- **Against the direct production run.** Each run was checked against its object's direct production run, as M2 checked its instrumented runs: `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv` are byte-identical, and so is the `result` section of `summary.json`, the digest of R's final `.Random.seed` included. The run's own hook checks are clean (`data/doubletfinder/compare/<object>/instr_vs_direct_production-sapphirerapids.tsv`).
- **M6's production Rust runs against them,** stage by stage, in exact mode: every check is identical. It is Step 1's first stage-by-stage comparison of these two objects in the production profile (M6 compared them with the direct runs and the saved outputs). In the format of [doubletfinder.md](doubletfinder.md#instrumented-runs)'s table:

| object | R run | dumps | S1 | S2a | S2b | S2c | S2d | S3 | S4 | S5 | S6 | S7 | RNG checkpoints | files | result | values compared | differing | n/a | missing | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | instrumented-production-sapphirerapids | light | = (35) | = (126) | = (42) | = (63) | = (49) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 11,481,122 | 0 | 28 | 0 | yes |
| cellbender_fromraw_ghirsutum2 | instrumented-production-sapphirerapids | light | = (35) | = (126) | = (42) | = (63) | = (49) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 13,658,206 | 0 | 28 | 0 | yes |

**The floor on ten objects.** `m11_floor.sbatch doubletfinder floor_kernel_m11_ten_objects` (job 61102448, `data/doubletfinder/compare/floor_kernel_m11_ten_objects/`) compares each object's two kernels' R runs, block by block, with the statistics the grader uses:

| array | eight objects (M2's table, M11's first measurement) | soupx_emptydrops_sscovelli | cellbender_fromraw_ghirsutum2 | the bound before | the bound now |
|---|---|---|---|---|---|
| S2d residual statistics: max abs, scaled | 5.684e-14, 1.364e-16 | 1.066e-14, 3.267e-17 | 3.553e-15, 1.139e-16 | 5.7e-14, 1.4e-16 | unchanged |
| S3 embeddings: max abs, scaled | 3.548e-11, 2.407e-13 | **3.665e-11, 3.133e-13** | 1.531e-11, 1.116e-13 | 3.6e-11, 2.5e-13 | **3.7e-11, 3.2e-13** |
| S3 loadings: max abs, scaled | 8.034e-14, 2.970e-13 | **8.370e-14**, 2.674e-13 | 7.752e-14, 1.496e-13 | 8.1e-14, 3.0e-13 | **8.4e-14**, 3.0e-13 |
| S3 standard deviations: max rel | 7.211e-15 | 6.232e-15 | 4.573e-15 | 7.3e-15 | unchanged |
| irlba's d: max abs, scaled | 7.503e-12, 3.799e-15 | 5.912e-12, 2.057e-15 | 2.956e-12, 2.732e-15 | 7.6e-12, 3.8e-15 | unchanged |
| S4 digest: cells (values) differing | 0 | 0 | **3 (6)**, the final run | identical (not listed) | **`cells_differing = 3, values_differing = 6`** |

Everything else the dumps hold is identical across kernels on both large objects, as on the eight: S1, S2a, S2b, S2c, every RNG checkpoint, every scalar (irlba's iterations and products included), S5, S6 and S7. So are the final files, as M2 found end to end.

- **The PCA moves slightly more on soupx_emptydrops_sscovelli** than on any of the eight, in one sweep iteration (`sweep_pN0.1`). These are exactly the differences the first re-grading reported between Rust's production run and R's portable run: Rust's production run is R's production run bit for bit (above), so under mismatched kernels a Rust-against-R pair shows R's own kernel difference.
- **R itself orders 3 cells' neighbours differently** in cellbender_fromraw_ghirsutum2's final run (k = 2,373). For each of the three cells, both kernels give the same k + 1 = 2,374 nearest positions, and two or three of them swap ranks (ranks 13 to 15, 81 and 82, 84 and 85). Their distances are equal to within 1.1e-15, 4 to 15 units in the last place, and two pairs are exact ties under the portable kernel, ordered by index. This was recomputed from each run's own final embeddings with `fields::rdist`'s arithmetic and R's `order()`, and the recomputed digests equal R's under both kernels. pANN and the calls depend only on which cells are among the first k, so S5 and S7 stay identical, and so do the final files. A digest has no magnitude to bound, so the floor's entry counts cells and values, as `[clustering.floor.kernel]` does for the neighbour distances and purities; the other nine objects keep S4 identical. Plan B5's S4 criterion holds as well: 99.995% of the cells are identical (all but 3 of 59,323), and each difference is a near-tie.

`[doubletfinder.floor.kernel]` takes the new maxima, rounded up as before, and the S4 entry. A unit test (`floor_bounds_the_large_objects`) grades the measured maxima against the floor as built.

## Benchmark

[`m11_benchmark.sbatch`](../../validation/slurm/m11_benchmark.sbatch) (job 60804011) ran `scrust cluster --counts <object>/input --preset workflow`, without stage dumps, on each tier-P object at 1 thread and at 32, in the production profile on one Sapphire Rapids node (holy8a24301, Xeon Platinum 8480CL, partition `test`, 32 CPUs and 96 GiB reserved), with the validation runs' binary. R's chain is M8's direct production runs on the same input directory (`reference_clustering.R`: R 4.5.3 with OpenBLAS at one thread, on the same CPU type), without RunUMAP, and with it for context. Time and peak RSS are the whole process's (`/usr/bin/time -v`); the stages' times come from each run's `summary.json` (`run.phases`, on both sides). One row per run: [`clustering_m11_bench.tsv`](clustering_m11_bench.tsv) (`m11_report.sh <label> bench`).

**Whole runs.** *R's chain*: from the timer's start (after R has loaded its packages) to the end of the silhouette, reading the input included: R's like for like with Rust's run, which reads the counts and writes its outputs but computes no purity yet (M12's benchmark adds it: [M12 benchmark](#m12-benchmark)). *Purity* (M12) and *writing* (R's run's files in the export's layout, the SCT assay's corrected counts and scale.data among them; the workflow's export is M13's) are R's shares that Rust does not have yet. *R ÷ Rust*: R's whole run, then R's chain, against Rust at 1 thread; then R's whole run against Rust at 32 threads; memory, against Rust at 1 thread.

| object | cells | R without RunUMAP: time, peak RSS | R's chain | R's purity, writing | R with RunUMAP: time | Rust, 1 thread: time, peak RSS | Rust, 32 threads: time, peak RSS | R ÷ Rust, time | R ÷ Rust, memory |
|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 47.1 s, 2.50 GiB | 30.3 s | 0.2 s, 8.1 s | 71.4 s (9.5 s RunUMAP) | **11.7 s, 0.15 GiB** | **3.0 s, 0.18 GiB** | 4.0×, 2.6×; 15.8× | 16× |
| soupx_tenx_10xmousebrain | 11,831 | 215.0 s, 12.02 GiB | 91.9 s | 1.0 s, 113.6 s | 295.0 s (15.7 s RunUMAP) | **38.0 s, 1.44 GiB** | **14.7 s, 1.44 GiB** | 5.7×, 2.4×; 14.6× | 8× |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 211.5 s, 12.08 GiB | 87.3 s | 0.9 s, 115.3 s | 226.4 s (12.7 s RunUMAP) | **37.2 s, 1.44 GiB** | **14.9 s, 1.44 GiB** | 5.7×, 2.3×; 14.2× | 8× |
| soupx_emptydrops_4500871lane | 26,598 | 204.8 s, 23.36 GiB | 153.5 s | 9.7 s, 33.5 s | 240.8 s (32.1 s RunUMAP) | **75.2 s, 1.48 GiB** | **23.9 s, 1.48 GiB** | 2.7×, 2.0×; 8.6× | 16× |
| soupx_emptydrops_sscovelli | 79,409 | 808.4 s, 62.13 GiB | 474.6 s | 192.6 s, 131.9 s | 883.2 s (107.7 s RunUMAP) | **205.2 s, 3.15 GiB** | **59.8 s, 3.19 GiB** | 3.9×, 2.3×; 13.5× | 20× |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 778.0 s, 57.52 GiB | 456.8 s | 180.3 s, 130.7 s | 1,010.8 s (117.1 s RunUMAP) | **203.0 s, 3.15 GiB** | **61.0 s, 3.19 GiB** | 3.8×, 2.3×; 12.8× | 18× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 905.9 s, 104.13 GiB | 604.8 s | 74.4 s, 202.8 s | 1,144.7 s (82.9 s RunUMAP) | **255.5 s, 2.65 GiB** | **52.1 s, 2.69 GiB** | 3.5×, 2.4×; 17.4× | 39× |

**By stage**, in seconds (R: `run.phases`; Rust: `run.phases`, each stage with the writing of its own files and stage records, `write` being the labels and `metadata.tsv`):

| object | run | load | SCTransform | RunPCA | FindNeighbors | FindClusters | silhouette | purity | write | whole run | peak RSS GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | R, no RunUMAP | 4.0 | 20.7 | 4.3 | 0.85 | 0.22 | 0.26 | 0.18 | 8.1 | 47.1 | 2.50 |
| soupx_tenx_SRR18256555 | Rust, 1 thread | 0.36 | 6.9 | 3.6 | 0.59 | 0.17 | 0.04 | 0.18† | 0.01 | 11.7 | 0.15 |
| soupx_tenx_SRR18256555 | Rust, 32 threads | 0.36 | 1.4 | 0.70 | 0.34 | 0.17 | 0.04 | 0.018† | 0.01 | 3.0 | 0.18 |
| soupx_tenx_10xmousebrain | R, no RunUMAP | 23.3 | 58.8 | 6.0 | 2.5 | 1.1 | 0.22 | 1.0 | 113.6 | 215.0 | 12.02 |
| soupx_tenx_10xmousebrain | Rust, 1 thread | 5.2 | 21.2 | 7.7 | 2.6 | 1.0 | 0.21 | 1.88† | 0.0 | 38.0 | 1.44 |
| soupx_tenx_10xmousebrain | Rust, 32 threads | 5.2 | 5.6 | 1.7 | 1.0 | 1.0 | 0.21 | 0.10† | 0.0 | 14.7 | 1.44 |
| filtered_seurat_tenx_10xmousebrain | R, no RunUMAP | 21.0 | 56.2 | 6.1 | 2.6 | 1.2 | 0.23 | 0.94 | 115.3 | 211.5 | 12.08 |
| filtered_seurat_tenx_10xmousebrain | Rust, 1 thread | 5.2 | 21.4 | 6.9 | 2.4 | 1.1 | 0.21 | 1.87† | 0.0 | 37.2 | 1.44 |
| filtered_seurat_tenx_10xmousebrain | Rust, 32 threads | 5.2 | 5.6 | 1.7 | 1.1 | 1.1 | 0.21 | 0.11† | 0.0 | 14.9 | 1.44 |
| soupx_emptydrops_4500871lane | R, no RunUMAP | 24.2 | 102.1 | 15.5 | 7.3 | 4.2 | 0.32 | 9.7 | 33.5 | 204.8 | 23.36 |
| soupx_emptydrops_4500871lane | Rust, 1 thread | 5.5 | 42.6 | 16.0 | 6.4 | 4.0 | 0.67 | 9.23† | 0.00 | 75.2 | 1.48 |
| soupx_emptydrops_4500871lane | Rust, 32 threads | 5.4 | 8.3 | 3.1 | 2.5 | 3.9 | 0.67 | 0.37† | 0.0 | 23.9 | 1.48 |
| soupx_emptydrops_sscovelli | R, no RunUMAP | 42.0 | 294.7 | 71.7 | 32.7 | 32.5 | 1.1 | 192.6 | 131.9 | 808.4 | 62.13 |
| soupx_emptydrops_sscovelli | Rust, 1 thread | 8.4 | 109.6 | 49.2 | 20.5 | 15.0 | 2.4 | 83.5† | 0.08 | 205.2 | 3.15 |
| soupx_emptydrops_sscovelli | Rust, 32 threads | 8.2 | 19.5 | 8.4 | 6.2 | 15.1 | 2.3 | 2.27† | 0.07 | 59.8 | 3.19 |
| filtered_seurat_emptydrops_sscovelli | R, no RunUMAP | 36.9 | 284.6 | 68.8 | 31.7 | 33.7 | 1.1 | 180.3 | 130.7 | 778.0 | 57.52 |
| filtered_seurat_emptydrops_sscovelli | Rust, 1 thread | 8.2 | 108.7 | 47.8 | 20.2 | 15.6 | 2.4 | 80.7† | 0.06 | 203.0 | 3.15 |
| filtered_seurat_emptydrops_sscovelli | Rust, 32 threads | 8.2 | 19.8 | 8.8 | 6.3 | 15.4 | 2.4 | 2.81† | 0.06 | 61.0 | 3.19 |
| cellbender_fromraw_ghirsutum2 | R, no RunUMAP | 37.1 | 477.1 | 49.3 | 20.9 | 19.3 | 1.0 | 74.4 | 202.8 | 905.9 | 104.13 |
| cellbender_fromraw_ghirsutum2 | Rust, 1 thread | 9.3 | 181.3 | 39.8 | 14.3 | 9.6 | 1.1 | 43.6† | 0.0 | 255.5 | 2.65 |
| cellbender_fromraw_ghirsutum2 | Rust, 32 threads | 9.2 | 20.8 | 7.2 | 4.5 | 9.2 | 1.1 | 1.57† | 0.00 | 52.1 | 2.69 |

† Purity: M12's benchmark (job 61652420, [M12 benchmark](#m12-benchmark)), on another node of the same type; M11's binary did not compute it, and the other columns are M11's runs.

- **C5's criterion is met.** At one thread, as R runs, Rust takes less time than R's whole run on every object (2.7 to 5.7 times less) and than R's chain alone (2.0 to 2.6 times less); at 32 threads, 8.6 to 17.4 times less than R's whole run. It peaks at 0.15 to 3.19 GiB against R's 2.50 to 104.13 GiB (8 to 39 times less); on the two 79,409-cell objects at 3.15 GiB at one thread and 3.19 GiB at 32, under the 10 GiB bound.
- **Where the time goes.** SCTransform is 53 to 71% of Rust's one-thread run and 2.4 to 3.0 times faster than R's. RunPCA at one thread is about R's speed (0.8 to 1.5 times): both spend it in irlba's products, in the same OpenBLAS kernel. At 32 threads SCTransform runs 3.8 to 8.7 times faster than at one, RunPCA 4.1 to 5.9 times, FindNeighbors 1.7 to 3.3 times (Annoy's build is sequential, as RcppAnnoy's); FindClusters not at all with one resolution (the starts run in turn, M10), and the silhouette is sequential. So on the 79,409-cell objects at 32 threads FindClusters is a quarter of the run.
- **Where the memory goes.** R's peak is set in SCTransform, which holds the residuals of every kept gene, dense, on five objects, and in RunPCA, 1.9 to 2.9 GiB above SCTransform's, on the two 79,409-cell objects (M8). Rust's is set in SCTransform on every object: it holds the kept genes' counts, two bins' products (500 genes each, 0.6 GiB at 79,409 cells) and the 3,000 rows that become scale.data (1.8 GiB), and keeps no other residuals. Reading the counts peaks at 1.83 GiB on the 79,409-cell objects.
- **The same files.** The benchmark runs' 17 output files are byte-identical at 1 and 32 threads, and identical to the validation runs' (8 threads, with stage dumps), on all seven objects ([Threads](#threads)).
- **RunUMAP**, for context: R's runs with it spent 9.5 to 117.1 s in RunUMAP. M8's R runs shared their nodes (partition `test`), so differences of a minute between R runs are noise (M8).
- **The run record's offset.** The validation binary timed the chain's stages from the making of its output writer but placed them after the end of reading, so in `run.phases` SCTransform's time includes the interval between the two (computing `nCount_RNA` and `nFeature_RNA`: at most 0.12 s in these runs, from each record's own `sctransform_s`) and `write`'s lacks it, which is why `write` can read 0.0 above. Fixed after these runs (`scrust cluster` now puts the chain's times on its own clock, and cargo test checks the phases are in order); the whole-run times come from `/usr/bin/time` and are not affected.

## Tests

`cargo test --workspace` passes in debug and in release: 253 tests in 38 suites, none failing or ignored, with no compiler warnings, under env.sh's kernel (Haswell). 31 of them are new in M11:

- `crates/rcompat` (1): `fixture_kernel` reads a fixture's kernel from its run record and from `tests/fixtures/kernels.tsv`, and the check passes under the fixture's kernel and stops under another;
- `crates/scio` (11): the TOML reader (values and tables, strings and inline tables, arrays of tables, the errors TOML requires, special floats, and the repository's own `tolerances.toml`), the stage-dump layout written and read back, kernel-gated grading (the modes, sign-aligned statistics, grading against a floor's bounds), the streaming SHA-256;
- `crates/doubletfinder`, unit test (1): `floor_bounds_the_large_objects`. The floor as built holds the two large objects' R-vs-R maxima: S3 on soupx_emptydrops_sscovelli, and S4's 3 cells on cellbender_fromraw_ghirsutum2, where one cell more fails it. Under matching kernels, and for an array the floor does not list, the grading stays exact;
- `crates/cluster`, unit tests (11): `approxSilhouette` (a case by hand, NA with one cluster, absent levels and ties to the first cluster, a long-double sum that double precision would lose), the parameters (Seurat's defaults and the preset, a TOML file under flags, every refusal), the PCA guard, `CreateSeuratObject`'s columns and `--metadata`'s, R's `.Random.seed` digest;
- `crates/cluster/tests/m11_chain.rs` (4): the chain from the counts on all seven committed tier-T objects against M8's portable runs (every label, level and workflow file, every silhouette width bitwise, the counts' columns, 26 scalars of SCTransform, irlba, Annoy and the optimizer; for the two base objects also the embeddings, the K5, K7 and K8 arrays, Annoy's index file and both graphs); the same files at 1, 2 and 5 threads; silhouette teacher-forced on R's embeddings and labels; the stages an observer sees;
- `crates/scrust/tests/cli.rs` (3, and one extended): the four routes (`--counts`, `--sct`, `--pca`, `--snn` with `--pca`, after `scrust sctransform`, `pca` and `neighbors`) give byte-identical files and the same clustering and silhouette in `summary.json`; a TOML file under the preset and a flag after it, the values recorded, the refusals and the PCA guard; `scrust validate cluster` against a reference that records no kernel runs in tolerance mode and says so; and `scrust validate doubletfinder` on a run recorded under another kernel runs in tolerance mode, with the reason and each check's class.

**The fixture tests' kernel check.** Every test that compares a BLAS-dependent result with a portable fixture now calls `rcompat::blas::require_fixture_kernel` first. It reads the kernel in use from `openblas_get_corename()` (not from `OPENBLAS_CORETYPE`) and the fixture's from its run record (`reference/summary.json`'s fingerprint) or from `tests/fixtures/kernels.tsv`, and stops with both kernels named when they differ. This replaces the silent skips (`skipped: needs OPENBLAS_CORETYPE=Haswell`), which looked only at the variable. The inventory:

| test | compares | BLAS | kernel check |
|---|---|---|---|
| `rcompat` `nlminb_against_r` | nlminb on R's vectors | PORT's `ddot`, `dnrm2` | added |
| `pca` `irlba_against_r` | irlba on R's vectors | `dgemm`, `dgemv`, `dgesdd` | replaces the skip |
| `doubletfinder` `m4_tier_t`: `s2b_fits_are_exact` | S2b teacher-forced | Cox-Reid `dgemv`, nlminb | added |
| `m4_tier_t`: `s2d_residuals_are_exact`, `s2_chained_run_is_exact` | S2d, S2 chained | the residuals' `dgemm` | replaces the skip |
| `m5_tier_t`: `s3_after_rust_sctransform_is_exact`, `s3_teacher_forced_on_data` | S3 | irlba | replaces the skip |
| `m6_tier_t`: both tests | the whole tool | all | replaces the skip |
| `scrust` `cli.rs`: `doubletfinder_runs_every_stage_from_the_counts`, `sctransform_matches_r_byte_for_byte`, `pca_matches_r_byte_for_byte`, `cluster_validate_without_a_recorded_kernel_is_not_exact` | against R's files | all | replaces the skip (the last is new) |
| `cluster` `m11_chain`: `the_chain_from_the_counts_is_r_s` | the chain | all | new |
| BLAS-free or Rust against Rust: `m4_tier_t` S2a and S2c, `m5_tier_t` `r_irlba_draws_only_its_start_vector`, `m2_fixtures`, `m3_tier_t`, `m4_vectors`, `kde`, `knn`, `sctransform` vectors, `neighbors`, `cluster` K0, K7 and edge cases, `m11_chain`'s other three tests, `cli.rs`'s routes, parameters and the other commands | | none | none needed |

**Under another kernel** (`validation/slurm/m11_kernel_tests.sbatch auto`, job 61108077, on this report's tree: a Sapphire Rapids node, `OPENBLAS_CORETYPE` unset, so OpenBLAS runs its SapphireRapids kernel), `cargo test --workspace --release` gives 240 passes and 13 failures, and every failure is a kernel check: "this result depends on OpenBLAS's kernel. The fixture (...) was made under kernel Haswell, but this process runs kernel SapphireRapids. Run the tests under the fixture's kernel: OPENBLAS_CORETYPE=Haswell (env.sh's default), on a CPU that supports it." The 13 are exactly the BLAS-dependent fixture tests above (DoubletFinder's M6 thread test among them: it also compares with R's files); every teacher-forced and BLAS-free test, and every test that compares Rust with Rust only (the chain's and M3's thread tests, the routes), passes, the kernel check's own test and the floor's included (`data/clustering/m11_kernel_tests/SapphireRapids/`). Earlier runs of the same job gave the same 13 failures: 60804256, before the kernel check's own test was added, and 60815653, before the floor's.


## Findings

- **The chain is exact end to end under the reference's kernel,** on all 34 objects in both profiles: SCTransform at seed 12345 through both assay paths (the emptyDrops objects' instrumented runs are v3, and their v5 runs compare the same), irlba with 50 PCs at production scale (`work = 66`, up to 79,409 cells), Annoy, the SNN graph, the optimizer, GroupSingletons (the 81 singletons of the four tier-P objects with singletons) and the silhouette widths, bit for bit. C7's two risks for M11, "irlba with 50 PCs at scale" and "SCTransform at seed 12345", did not materialise.
- **irlba's restarts depend on the kernel** (R-vs-R floor): one or two more or fewer in 5 of the 34 objects, with the same draws and results within the floor. M8 compared stage arrays and random-number states; these scalars were not in its floor, and are now.
- **DoubletFinder's floor now covers its two largest objects.** M2 measured the floor per stage on eight objects, and the two large ones only end to end. Against that floor, tolerance mode flagged Rust's production runs of soupx_emptydrops_sscovelli and cellbender_fromraw_ghirsutum2 against R's portable runs, and 16 of the 18 mismatched pairs passed.

  The instrumented production R runs of the two objects, made for this, show that R itself moves that much between the kernels:
  - sscovelli's S3 embeddings by 3.665e-11 and its loadings by 8.370e-14;
  - ghirsutum2's neighbour order for 3 cells of the final run: near-ties within the same neighbour set, with pANN and the calls identical.

  The floor measured on all ten objects takes them, and all 20 mismatched pairs pass. Under matching kernels, M6's production Rust runs are identical to the new R runs stage by stage. So under mismatched kernels each Rust-against-R pair differs exactly as R differs from itself.
- **Annoy's index file differs across kernels in every object** while the neighbour ids and both graphs are identical: a PCA component that changes sign flips one coordinate of every item and every split normal, which changes the stored bytes but no distance and no split (the Euclidean build is invariant under a coordinate's sign). The floor records the file as differing (`may_differ`); K5 and K6 carry the claim.
- **`nlminb` calls BLAS** (`ddot`, `dnrm2` in PORT), so S2b depends on BLAS twice over (A5b corrected). Across kernels S2b was identical in every measured object, so its floor is "identical".
- **The silhouette needs R's long double.** R's `colMeans`, `colSums` and `sum` accumulate in 80-bit long double; the port does so too (`rcompat::F80`), and its widths are bitwise R's in every comparison, teacher-forced or from the counts.
- **R's random-number checkpoints after RunUMAP cannot be reproduced and need not be.** RunUMAP draws from R's stream; the Annoy build, the search and the optimizer do not; GroupSingletons reseeds with `set.seed(1)`. So R's states between RunUMAP and GroupSingletons carry UMAP's draws and nothing reads them; the validator lists them as not compared, and compares R's state after GroupSingletons where it drew (M10's rule).
- **R's digest of scale.data is in R's layout:** the rows in the SCT assay's gene order (not `VariableFeatures()` order), the doubles column by column (M8's `raw_sha256`). Rust's `result.scale_data_sha256` is computed the same way, streaming, and equals R's in every same-kernel comparison, all 34 objects in both profiles; for the three light-level objects it is K1's only check of scale.data.

## Corrections to the plan

Corrected in [docs/PLAN.md](../PLAN.md) on this branch:

- **A5b mechanism 4: what in DoubletFinder depends on BLAS.** The plan said "only the Cox–Reid term in S2b, S2d and S3". S2b's `nlminb` calls BLAS too: R's PORT routines `dd7tpr` and `dv2nrm` are `ddot` and `dnrm2` (M4's manifest notes it). The classification of kernel-gated grading follows the corrected statement.
- **A5b mechanism 3:** the switch to tolerance mode and the fixture tests' kernel check are now described as built, not planned.
- **A1, C3 components:** `scio` also holds the TOML reader, the stage-dump layout (moved from `doubletfinder`) and kernel-gated grading.
- **C3, commands and parameters:** `scrust cluster`'s starting points write the same files because `scrust sctransform` writes the counts' per-cell columns and `scrust pca` and `scrust neighbors` carry them forward; every SCTransform, RunPCA and FindNeighbors parameter has a flag of its own (`--sct-seed`, `--pca-seed`, ...; `--seed` stays FindClusters'); `--dump-stages`; the preset turns on the workflow's PCA guard, which is off without it.
- **C7, irlba with 50 PCs at scale, and SCTransform at seed 12345:** their outcomes (K2 and K1 identical under the reference's kernel on every object, K1 through both assay paths; irlba's restarts depend on the kernel in 5 objects).
- **C6, M11's row: the floor DoubletFinder's tolerance mode uses.** The plan said "M2's". It is M2's floor measured per stage on all ten objects in M11: M2 had the two large objects end to end only.


## Deviations and notes

- **The validation binary's stage labels.** The validation and benchmark runs used `scrust-m11-eedf376`. After them, the validator's labels of a few stage-dump scalars were corrected: Annoy's (`annoy_*`), the neighbours' (`knn_*`) and the optimizer's (`optimizer_*`, `n_singletons`) scalars were filed under stage K8 in its reports, and are now under K4, K5 and K7. Only the labels changed; every check, its status and its counts are the same. The summary line printed to stderr also names the saved object's checks apart now.
- **Purity is not computed** (M12): asked for by default, as in the workflow, it is reported as not computed yet on stderr and in `summary.json`. The `neighborhood_purity` column is left out of the file comparisons.
- **`--dump-stages` needs `--counts`.** The stage dumps describe the whole chain; a staged start refuses them.
- **New files in Step 1's commands' outputs.** `scrust sctransform` writes `cell_metadata.tsv` (the counts' per-cell columns), and `scrust pca` and `scrust neighbors` copy it forward when their input has it. Every file they wrote before is unchanged, byte for byte (the CLI tests against R's files pass).
- **`scrust validate floor`** is new: it measures the R-vs-R floor between pairs of runs with the grader's own statistics.
- **The saved objects' kernel** is a parameter (`--saved-kernel`, default SapphireRapids): the exported saved objects record no kernel, and M8 established that every workflow job ran with SapphireRapids.
- **Digests in tolerance mode** are reported, not graded: a digest has no magnitude to hold against the floor's bounds.
- **`scrust validate doubletfinder --saved`** is not kernel-gated: its checks read only the calls, the pANN as written and the counts, which the floor found identical across kernels, so they are exact in either mode (each is classified BLAS-dependent in the report, criterion exact).
- **The DoubletFinder demonstration ran twice** (jobs 60804010 and 61104652). With M2's floor, 16 of 18 pairs passed. Once the two large objects' instrumented production R runs existed and the floor was measured on all ten objects, all 20 pairs passed; the two new pairs are the large objects' portable Rust runs against the new R runs. The cluster validator's demonstration passes on all 68 pairs.
- **Two new R runs.** The instrumented production R runs of DoubletFinder's two large objects (jobs 61097286 and 61097287) were made after this report's first version, at the user's request, to measure the floor on those objects. They are the only R runs M11 made. Like the M2 runs, they are listed in the manifest, in M11's section.
- **R's `colMeans` and `colSums`** are computed where the silhouette needs them, with `rcompat`'s exact long double (`F80`, tested bit for bit against the x87 since M1) and R's loop order, rather than as new `rcompat` functions: the silhouette is their only user so far. `rcompat` gains the fixture tests' kernel check, with its own test.
- **Not ported:** RunUMAP (A7b); Leiden and the other refused paths; `--export` and the corrected counts (M13).
- **Rust's validation outputs** (`data/clustering/<object>/rust/m11-eedf376/`, `m11-eedf376-t1/`) are not listed in the manifest, as in M9 and M10; `clustering_m11.tsv` has one row per comparison.
- **The run record's stage times.** In the validation binary's `run.phases`, SCTransform's time includes the computing of `nCount_RNA` and `nFeature_RNA` (at most 0.12 s) and `write`'s time lacks it: the chain's stages were timed from the making of the output writer but placed after the end of reading. Fixed after the runs; only the run record was affected ([Benchmark](#benchmark)).

## Reproducing

From the repository's top directory, after `source env.sh` and `cargo build --release -p scrust` (copy the binary to a file of its own, here `data/clustering/bin/scrust-m11-eedf376`):

```bash
export SCRUST_BIN=data/clustering/bin/scrust-m11-eedf376
sbatch validation/slurm/m11_floor.sbatch
sbatch --partition=shared validation/slurm/m11_cluster.sbatch m11-eedf376 portable <objects>
sbatch --partition=test validation/slurm/m11_cluster.sbatch m11-eedf376 production <objects>   # Sapphire Rapids
sbatch --partition=shared --cpus-per-task=1 --mem=48G validation/slurm/m11_cluster.sbatch m11-eedf376-t1 portable <tier-P objects>
bash validation/slurm/m11_threads.sh m11-eedf376 m11-eedf376-t1 portable <tier-P objects>
sbatch validation/slurm/m11_doubletfinder.sbatch m11-eedf376 <M7 objects>
sbatch --partition=test validation/slurm/m11_kernel_tests.sbatch auto
VALIDATION_LABEL=m11-eedf376 sbatch validation/slurm/m11_benchmark.sbatch m11-eedf376 <tier-P objects>
bash validation/slurm/m11_report.sh m11-eedf376 > docs/validation/clustering_m11.tsv
bash validation/slurm/m11_report.sh m11-eedf376 bench > docs/validation/clustering_m11_bench.tsv
```

DoubletFinder's two large objects (the R runs, their checks, the floor on all ten objects, and the re-grading with it; the binary `data/clustering/bin/scrust-m11-63e7b3b` is built from commit 63e7b3b, which holds that floor; the commit after it changes only comments in `tolerances.toml`, which the binary embeds):

```bash
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 sbatch --partition=test --mem=200G --time=0-03:00:00 validation/slurm/reference_doubletfinder.sbatch production instrumented soupx_emptydrops_sscovelli
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 sbatch --partition=test --mem=250G --time=0-03:00:00 validation/slurm/reference_doubletfinder.sbatch production instrumented cellbender_fromraw_ghirsutum2
SCRUST_BIN=data/clustering/bin/scrust-m11-eedf376 sbatch validation/slurm/m11_df_references.sbatch m11-references soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2
SCRUST_BIN=data/clustering/bin/scrust-m11-eedf376 sbatch validation/slurm/m11_floor.sbatch doubletfinder floor_kernel_m11_ten_objects
SCRUST_BIN=data/clustering/bin/scrust-m11-63e7b3b sbatch validation/slurm/m11_doubletfinder.sbatch m11-63e7b3b <M7 objects>
```

One run by hand, on a committed fixture (seconds, any machine with AVX2, under env.sh's Haswell kernel):

```bash
f=tests/fixtures/clustering/filtered_seurat_tenx_test
target/release/scrust cluster --counts $f/input --preset workflow --out /path/to/out --dump-stages full
target/release/scrust validate cluster --reference $f/reference --candidate /path/to/out --optimizer $f/optimizer
```

The committed reference records no kernel, so this comparison runs in tolerance mode and says so; every check is identical. Against `data/clustering/<object>/runs/instrumented-portable-sapphirerapids` (which records Haswell) it runs in exact mode.

## M12: neighbourhood purity

M12 computes the workflow's second score, bluster 1.20.0's `neighborPurity`, as `add_silhouette_to_metadata` calls it: on the PCA's dims 1:30 and `seurat_clusters`, with `k = 50`, `weighted = TRUE` and BiocNeighbors 2.4.0's default index, KMKNN (docs/PLAN.md C1 step 6, C3e). The search is exhaustive, with knncolle's arithmetic, and a certificate reports every pair that KMKNN's floating-point bounds could treat differently. The certificate fired, on the median pair of every odd-sized object, and at M12's decision point the user chose to keep the exhaustive result and report the flagged pairs ([The decision point](#the-decision-point)).

### What M12 built

| Where | M12 contents |
|---|---|
| `crates/knn/src/purity.rs` (new) | the search `neighborPurity` makes: each cell's k-th distance (`findDistance`), R's `median` as the radius, each cell's neighbours within it (`findNeighbors`), with knncolle's arithmetic, exhaustively, in parallel over tiles of cells; the certificate (the margin, the flagged pairs, each labelled by kind) |
| `crates/cluster/src/purity.rs` (new) | `neighborPurity`'s weights, per-cluster sums (`sum_neighbor_weights`), long-double totals (`rowSums`), purities and maxima (`max.col`); the workflow's guard (NA with fewer than 2 cells or clusters) |
| `crates/cluster/src/chain.rs`, `chain_output.rs` | purity in the chain's library entry point (C3h), on every route; `metadata.tsv`'s `neighborhood_purity` and `purity_maximum`; `summary.json`'s `result.purity` (the radius, the mean neighbour count, the purities' digest, the certificate with every flagged pair); the K9 stage dumps in R's layout, and the certificate's `K9_certificate_*` arrays, which R has not |
| `crates/cluster/src/validate.rs` | K9 from the stage dumps (the integer arrays the floor bounds graded as numbers; the neighbour lists cell by cell); `neighborhood_purity` and `purity_maximum` in the files; the saved object's purities, under its kernel; K9 teacher-forced on R's embeddings and labels; the certificate, whose flagged pairs are reported ("flagged, reported", a status that passes) |
| `crates/scrust` | `scrust cluster`: purity on every route, `--purity-lists auto\|yes\|no` (the neighbour lists in the stage dumps, as R's instrumented runs choose them), the radius and the certificate on stderr; `scrust validate cluster`: the verdict names reported certificates |
| `validation/tolerances.toml` | `[clustering.K9]` with the decided criterion; the K9 entries M11 measured but did not list (sums, totals, radius, mean neighbour count) and the neighbour lists measured cell by cell, in `[clustering.floor.kernel]`; K9 in `[clustering.kernel_gating]` |
| `validation/slurm` | `m12_certificate.sbatch` (the decision point's survey), `m12_cluster.sbatch`, `m12_benchmark.sbatch`, `m12_kernel_tests.sbatch`, `m12_lists_floor.sh`, `m12_report.sh`, `m12_manifest.sh` |
| `NOTICE`, `third_party/licenses/MIT-knncolle.txt` | bluster's `neighborPurity` and `sum_neighbor_weights`, BiocNeighbors' two searches, knncolle's arithmetic (MIT) |

The workspace still uses only the Rust standard library and the conda env's OpenBLAS; purity uses no BLAS.

### The search

`neighborPurity` makes two BiocNeighbors calls on one KMKNN index:

1. **`findDistance(idx, k = 50)`**: each cell's distance to its 50th nearest other cell. knncolle fills a queue with the 51 nearest cells, the cell itself among them, drops the cell itself (or, when more than 51 duplicates sit at distance 0, the farthest), and BiocNeighbors returns the last. Whatever the ties, that is the 51st smallest squared distance to all cells (the cell itself at 0), which is the 50th smallest to the other cells, and its square root. BiocNeighbors caps k at n − 1, with a warning, and returns 0 when no neighbour is left; the port does both.
2. **`median()`** of those distances is the radius, one for all cells: the middle value, or the long-double mean of the middle two when n is even (`rcompat`'s R `median`).
3. **`findNeighbors(idx, threshold = radius)`**: each cell's other cells whose squared distance is `<=` `radius * radius`, in double (knncolle's `denormalize`), listed by squared distance.

**Arithmetic.** knncolle's `EuclideanDistance::raw` is `out = 0; out += delta * delta` over the 30 dims in order, in double: the sum `fields::rdist` makes before its square root, which `knn` already computed for DoubletFinder. Searches compare these squared values. The port computes every pair this way, in tiles of 8 cells against blocks of 256 cells; each element sees the same IEEE subtract, multiply and add whatever the vector width, and a test checks the AVX-512, AVX2 and scalar paths against the pair-by-pair definition, bit for bit. Step 1 keeps each cell's 51 smallest squared distances; step 3 visits each cell's pairs within the radius (and within the certificate's bands) and counts them by cluster. Two passes over all pairs and no n × n matrix; tiles run in parallel and every result is placed by its cell, so the output does not depend on the thread count.

**Neighbour order.** BiocNeighbors lists a cell's neighbours by squared distance, and equal squared distances by KMKNN's internal index (the cells reordered by k-means cluster); the port orders those by the cell's index. Every list R kept (the 27 tier-T objects and soupx_tenx_SRR18256555) is identical to the port's, in order: no list holds two neighbours at the same squared distance.

### The weights and the sums

`neighborPurity`'s weights are `1/table(clusters)` (1 divided by the cluster's size, in double), or 1 with `weighted = FALSE`. `sum_neighbor_weights` sets each cell's own cluster's entry to its weight, then adds each neighbour's weight to its cluster's entry, in double. Every cell of a cluster has the same weight, so each entry is that weight added to 0 once per counted cell, in any order; the port takes it from a table of repeated sums per cluster. The totals are R's `rowSums`: one long-double accumulator per cell, the clusters added in level order from 0, then rounded to double (R 4.5.3's `do_colsum`). The purity is the cell's own cluster's entry divided by its total. The maximum is `max.col(ties.method = "first")` (`R_max_col`): the first cluster strictly larger than every earlier one, written as `purity_maximum`. With fewer than 2 cells or clusters every purity is NA, as the workflow's guard leaves them.

### The certificate

KMKNN (`knncolle_kmknn.hpp` as BiocNeighbors 2.4.0 builds it: `KNNCOLLE_KMKNN_USE_UPPER` undefined, so only the lower bound is used) keeps, for each cell `j` of a cluster with centre `c`, `r_j = sqrt(raw(x_j, c))`. For a query `t` it computes `D = sqrt(raw(t, c))` and skips the cells with `r_j < fl(D − S)` (`lower_bd`, by `std::lower_bound` on the sorted `r_j`, or the whole cluster when its largest `r_j` is below it), where `S` is the search's distance threshold: the radius in step 3; in step 1, the square root of the largest squared distance in its queue once the queue is full. It computes and compares every pair it visits exactly as the exhaustive search does, so it can only drop cells; and since `D ≤ d + r_j` (`d = |t − x_j|`), in exact arithmetic it never drops one that belongs. With rounding (`u = 2⁻⁵³`, `m` dims):

- a squared distance has a relative error of at most `g = γ(m + 2) = (m + 2)u / (1 − (m + 2)u)`: one rounding in the subtraction, two in the square, `m − 1` in the sum of non-negative terms; a distance, after `sqrt`, at most `e = g/2 + g²/2 + u + ug`;
- a skip therefore implies `S < d(1 + e) + κ·r_j`, with `κ = 2e + u + eu`;
- KMKNN's centres are means of cells (kmeans' `compute_centroids`, from which Hartigan–Wong's last centres always come) or cells (k-means++'s seeds and the edge cases), so each lies in the cells' bounding box, widened per dim by its rounding, `(γ(n + 2) + 2u)·max|x|`. Then `r_j ≤ B`, that box's diagonal.

**Step 3** (threshold `T`, the radius). A pair is in the result when its squared distance is `≤ fl(T·T)`. KMKNN can drop only a pair with `d > (T − κB)/(1 + e)`, so with a squared distance above `L = ((T − κB)/(1 + e))²(1 − g)`. Flagged: every pair within the radius above `L`. (A duplicate of the query, at 0, has `r_j = D` bit for bit and is never dropped.)

**Step 1** (the 50th distance). KMKNN's result can differ only if it drops one of the true 51 nearest while its threshold `θ`, the largest squared distance in its full queue, exceeds the true 51st smallest `v`: if `θ = v`, the queue already holds 51 values `≤ v` and the result is `v` whatever is dropped later. So `θ > v` is another cell's squared distance, and the drop needs `θ < H = ((sqrt(v/(1 − g))(1 + e) + κB)/(1 − u))²`. Flagged: every cell with another cell's squared distance in `(v, H)`. The 50th neighbour itself, at `v`, is not flagged: a tie at `v` cannot change the result.

Each constant is enlarged by 1% and each threshold moved 16 ulps outwards, so the flagged set can only be larger than the exact one. The margin grows with the data's extent `B`, not with the radius: on the reference objects `B` is 70 to 84 on tier T and 350 to 1,080 on tier P, and step 3's band is 220 to 12,760 ulps of the radius squared, not "a few units in the last place" (C3e, corrected).

**The median pair.** When the radius is a cell's own 50th distance, always with an odd number of cells (the median is then the middle value), that cell and its 50th neighbour lie on the radius by construction: their squared distance `v` gives `radius = sqrt(v)`, and the pair is within the radius exactly when `fl(radius · radius) ≥ v`, which the embeddings' last bits decide. The certificate flags such a pair, both ways, whenever it is within the radius, and labels it `median_pair`; every other flag is `other`. It cannot clear the pair without KMKNN's centres: KMKNN could drop it only if the far cell lay on the segment from the query to its k-means centre, with the triangle inequality tight to about `κB`, which on these objects means an angle of the order of 10⁻⁷ radians.

**Reported** in `summary.json` (`result.purity.certificate`: the counts by step and by kind, every flagged pair with its cells and barcodes, step, kind, squared distance, threshold, ulps from it and the band's end; the constants `g`, `e`, `κ`, `B`; the width of step 3's band and the gap to the nearest unflagged pair within the radius; the smallest step-1 gap in band widths), in the stage dumps (`K9_certificate_pairs`, `_step`, `_kind`, `_raw`, `_threshold`, and the scalars `purity_certificate_flagged*`), on stderr, and in `scrust validate cluster`'s reports, where a flagged certificate has the status "flagged, reported" and passes.

### The decision point

The certificate fired on the first odd-sized object. Before asking, the search ran on R's own embeddings and labels in both profiles for all 34 objects ([`m12_certificate.sbatch`](../../validation/slurm/m12_certificate.sbatch), binary `scrust-m12-survey`, jobs 61647015 and 61647016, `data/clustering/m12_certificate/m12-survey/`):

- 68 runs (34 objects × 2 profiles): the certificate flags the median pair of each of the 17 odd-sized objects whenever it lies within the radius (23 runs, at 0 or 1 ulp below the radius squared), and nothing else. On the even-sized objects the nearest pair within the radius lies at least 2.6·10⁸ ulps below the radius squared, against bands of 220 to 12,760 ulps; no cell has another squared distance within 985 band widths above its 50th.
- In every run the exhaustive search's K9 equals R's KMKNN's bit for bit (1,053 teacher-forced checks and 181 others: the 50th distances, the radius, the counts, the lists, the sums, the totals, the purities, the maxima), so R's KMKNN kept every flagged pair.

**Decided (the user, 2026-11-17): keep the exhaustive result.** The flagged pairs are reported (summary.json, the stage dumps, the reports and this report) and graded "flagged, reported", passing, not failing; each flag is labelled by kind (`median_pair` or `other`), so M17 can tell a structural flag from a new one; M17's comparison of the 840 saved objects' purities checks every flagged pair, so any R run that dropped one would show there. A KMKNN port stays optional ("Can follow"). The plan records the outcome in C3e, C5's K9 row, C6's M12 and M17 rows, "Can follow", C2 and C7 ([corrections](#m12-corrections-to-the-plan)).

### M12 results

Label `m12-4231635`: binary `data/clustering/bin/scrust-m12-4231635` (sha256 `f9d1f5bc...`, built from commit 4231635), [`m12_cluster.sbatch`](../../validation/slurm/m12_cluster.sbatch), 8 threads, M11's design: each Rust run from the counts dumps its stages at the level of M8's instrumented runs (full on tier T; standard on the small and medium tier-P objects, with purity's neighbour lists on the small one, `--purity-lists yes`, as M8's R run kept them; light on the three large ones), and is compared with every R run of its profile, the instrumented optimizer, the saved object and, on tier P, the workflow job's logged numbers; then with the other profile's instrumented run, a deliberately mismatched pair of kernels. The portable profile ran with `OPENBLAS_CORETYPE=Haswell` on Cascade Lake nodes (`shared`; jobs 61652407 and 61652409), the production profile with the kernel detected, SapphireRapids, on Sapphire Rapids (`test`; jobs 61652412 and 61652416; `sapphire` had 2,650 jobs waiting). One row per comparison: [`clustering_m12.tsv`](clustering_m12.tsv) ([`m12_report.sh`](../../validation/slurm/m12_report.sh)).

| tier | profile | comparison | comparisons | passed | mode | checks identical | within floor | differing | missing | certificates flagged, reported |
|---|---|---|---|---|---|---|---|---|---|---|
| T | portable | the profile's R runs | 108 | 108 | exact | 12,120 | 132 | 0 | 0 | 54 |
| T | production | the profile's R runs | 82 | 82 | exact | 7,645 | 0 | 0 | 0 | 36 |
| P | portable | the profile's R runs | 10 | 10 | exact | 1,530 | 14 | 0 | 0 | 8 |
| P | production | the profile's R runs | 22 | 22 | exact | 2,146 | 0 | 0 | 0 | 4 |
| T | portable | the other kernel's instrumented run | 27 | 27 | tolerance | 4,104 | 492 | 0 | 0 | 18 |
| T | production | the other kernel's instrumented run | 27 | 27 | tolerance | 4,104 | 492 | 0 | 0 | 18 |
| P | portable | the other kernel's instrumented run | 7 | 7 | tolerance | 1,038 | 160 | 0 | 0 | 5 |
| P | production | the other kernel's instrumented run | 7 | 7 | tolerance | 1,038 | 160 | 0 | 0 | 5 |

*Certificates flagged, reported*: per comparison, the candidate's certificate (`result.purity.certificate`) and, against an instrumented run, the certificate of the teacher-forced search on R's embeddings; each flags the median pair of an odd-sized object when it is within the radius.

**290 comparisons, all passing.** The 222 with the profile's own R runs are in exact mode: 23,441 checks identical, none differing or missing; the only checks within the floor (146) are against the saved objects in the portable profile, graded under the saved objects' kernel: the 118 silhouette widths M11 found, and 28 purities, those of the 28 portable R runs M8 found differing from their saved objects in one or two cells (the median pair, [Findings](#m12-findings)). In the production profile every check against the saved object is identical, the purities bit for bit, on all 34 objects.

**K9.** Over the 222 comparisons, 3,058 purity checks are identical: from the stage dumps (each cell's 50th distance, the neighbour counts, the neighbour lists cell by cell in R's order where R kept them, the weights, the sums, the totals, the purities, the maxima; the radius, the mean neighbour count and the other purity scalars), `metadata.tsv`'s `neighborhood_purity` and `purity_maximum` against the direct and instrumented runs, the saved objects' purities, and K9 teacher-forced on R's own embeddings and labels; 102 certificates are flagged and reported. **Teacher-forced**, BLAS-free: in the 98 comparisons with instrumented runs, 1,033 checks of every K9 array, the radius and R's `neighborhood_purity` column, all identical, and 32 certificates flagged and reported; in the 68 mismatched pairs, the same on the other kernel's embeddings (713 identical, 23 reported), as BLAS-free checks must be.

**Everything M11 checked is unchanged.** Every M11 report of these 290 comparisons was compared with M12's, check by check: all 31,269 checks M11 graded (its placeholders for later milestones aside) have the same status and the same counts in M12's report: the labels and levels, the silhouette widths, the files, K1 to K8.

**The certificate, object by object** (the Rust runs from the counts; R's own embeddings give the same, [The decision point](#the-decision-point)): on every even-sized object nothing is flagged, and the nearest pair within the radius lies at least 2.6·10⁸ ulps below the radius squared, against bands of 220 to 12,760 ulps; on every odd-sized object the median pair is flagged, both ways, whenever it is within the radius; no cell has another squared distance within 985 band widths above its 50th. Every flag is `median_pair`; none is `other`.

#### Tier T

Per object, both profiles' comparisons with their own R runs: the purity checks identical (K9 from the stage dumps, the files, the saved object, teacher-forced) and the certificates reported; the median pair flagged in each profile's Rust run (its cells, 1-based, and how many ulps below the radius squared its squared distance lies; "-" when it lies outside the radius or the object is even-sized); and the mismatched pairs' purity checks within the floor (both directions).

| object | cells | clusters | R runs (portable, production) | purity checks identical | certificates flagged, reported | median pair flagged (ulps below the radius squared): portable | production | mismatched kernels: purity checks within floor |
|---|---|---|---|---|---|---|---|---|
| posthocfilt_threshold_doubletfinder_decontx_raw_tenx_test | 80 (even) | 2 | 4, 3 | 99 | 0 | - | - | 4 |
| posthocfilt_threshold_scdblfinder_decontx_raw_tenx_test | 87 (odd) | 2 | 4, 3 | 91 | 4 | - | 66–80 (0) | 18 |
| posthocfilt_threshold_doubletfinder_decontx_raw_emptydrops_test | 92 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| posthocfilt_threshold_scdblfinder_decontx_raw_emptydrops_test | 112 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| posthocfilt_threshold_doubletfinder_cellbender_fromraw_test | 135 (odd) | 2 | 4, 3 | 89 | 10 | 32–93 (1) | 32–93 (1) | 4 |
| posthocfilt_threshold_doubletfinder_soupx_emptydrops_test | 137 (odd) | 2 | 4, 3 | 91 | 4 | - | 57–101 (1) | 18 |
| posthocfilt_threshold_doubletfinder_soupx_tenx_test | 142 (even) | 2 | 4, 3 | 99 | 0 | - | - | 4 |
| posthocfilt_threshold_scdblfinder_soupx_tenx_test | 154 (even) | 2 | 4, 3 | 99 | 0 | - | - | 4 |
| posthocfilt_threshold_scdblfinder_soupx_emptydrops_test | 161 (odd) | 3 | 4, 3 | 89 | 6 | 26–138 (0) | - | 18 |
| posthocfilt_threshold_scdblfinder_cellbender_fromraw_test | 163 (odd) | 2 | 4, 3 | 91 | 4 | - | 122–132 (1) | 18 |
| doubletfinder_decontx_raw_emptydrops_test | 237 (odd) | 4 | 4, 3 | 89 | 10 | 228–232 (1) | 228–232 (0) | 4 |
| doubletfinder_soupx_emptydrops_test | 237 (odd) | 3 | 4, 3 | 89 | 10 | 25–86 (1) | 25–86 (0) | 4 |
| doubletfinder_cellbender_fromraw_test | 247 (odd) | 3 | 4, 3 | 89 | 6 | 151–226 (0) | - | 18 |
| scdblfinder_decontx_raw_emptydrops_test | 262 (even) | 4 | 4, 3 | 99 | 0 | - | - | 4 |
| scdblfinder_soupx_emptydrops_test | 264 (even) | 4 | 4, 3 | 99 | 0 | - | - | 4 |
| decontx_raw_emptydrops_test | 267 (odd) | 4 | 4, 3 | 89 | 10 | 117–197 (0) | 117–197 (0) | 4 |
| filtered_seurat_emptydrops_test | 267 (odd) | 4 | 4, 4 | 118 | 6 | 106–136 (0) | - | 18 |
| soupx_emptydrops_test | 267 (odd) | 4 | 4, 3 | 89 | 10 | 141–212 (0) | 141–212 (0) | 4 |
| doubletfinder_decontx_raw_tenx_test | 273 (odd) | 3 | 4, 3 | 89 | 10 | 43–167 (0) | 43–167 (1) | 4 |
| doubletfinder_soupx_tenx_test | 276 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| cellbender_fromraw_test | 278 (even) | 4 | 4, 3 | 99 | 0 | - | - | 4 |
| scdblfinder_cellbender_fromraw_test | 278 (even) | 4 | 4, 3 | 99 | 0 | - | - | 4 |
| scdblfinder_soupx_tenx_test | 292 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| scdblfinder_decontx_raw_tenx_test | 296 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| decontx_raw_tenx_test | 300 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| filtered_seurat_tenx_test | 300 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |
| soupx_tenx_test | 300 (even) | 3 | 4, 3 | 99 | 0 | - | - | 4 |

#### Tier P

| object | cells | clusters | R runs (portable, production) | purity checks identical | certificates flagged, reported | median pair flagged (ulps below the radius squared): portable | production | mismatched kernels: purity checks within floor |
|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 (even) | 15 | 4, 3 | 99 | 0 | - | - | 4 |
| filtered_seurat_tenx_10xmousebrain | 11,831 (odd) | 21 | 1, 3 | 55 | 4 | - | 409–10819 (1) | 16 |
| soupx_tenx_10xmousebrain | 11,831 (odd) | 21 | 1, 3 | 58 | 2 | 1214–1698 (0) | - | 12 |
| soupx_emptydrops_4500871lane | 26,598 (even) | 31 | 1, 3 | 60 | 0 | - | - | 4 |
| cellbender_fromraw_ghirsutum2 | 59,323 (odd) | 22 | 1, 3 | 57 | 2 | 32068–47204 (0) | - | 16 |
| filtered_seurat_emptydrops_sscovelli | 79,409 (odd) | 37 | 1, 4 | 84 | 2 | 8542–38756 (0) | - | 16 |
| soupx_emptydrops_sscovelli | 79,409 (odd) | 36 | 1, 3 | 57 | 2 | 7745–60959 (1) | - | 16 |


### M12 threads

Every output is byte-identical at one thread and at many, purity's included:

| tier | profile | runs | threads | files compared per run | differing |
|---|---|---|---|---|---|
| T | portable | 27 | 8 and 1 (`m12_cluster.sbatch`'s second run) | 71, or 76 where a pair is flagged (the outputs, `sct/`, every start's partition, the stage dumps with K9 and the certificate's arrays) | 0 |
| T | production | 27 | 8 and 1 | 71 or 76 | 0 |
| P | portable | 7 | 8 (label m12-4231635) and 1 (m12-4231635-t1, job 61652418), [`m11_threads.sh`](../../validation/slurm/m11_threads.sh) | 63 to 68: the outputs and the stage dumps | 0 |
| P | production | 7 | 1 and 32 (the benchmark, job 61652420), and 8 (the validation run) | 17: the outputs | 0 |

`summary.json` is compared by its `result` section, where the purity's summary and certificate are. The tier-P portable runs at one thread were validated as the 8-thread runs were, and their 17 reports equal the 8-thread runs' check by check, details included. cargo test compares 1, 2, 5 and 16 threads ([M12 tests](#m12-tests)). The purity search splits its work by tiles of 8 cells and places every result by its cell; the sums are per cell.

### M12 benchmark

[`m12_benchmark.sbatch`](../../validation/slurm/m12_benchmark.sbatch) (job 61652420) ran M11's benchmark again with purity: `scrust cluster --counts <object>/input --preset workflow`, both scores, without stage dumps, on each tier-P object at 1 thread and at 32, in the production profile on one Sapphire Rapids node (holy8a24302, Xeon Platinum 8480CL, partition `test`, 32 CPUs and 96 GiB reserved), with the validation runs' binary, against M8's direct production R runs on the same inputs (R with OpenBLAS at one thread, the same CPU type), without RunUMAP and, for context, with it. Time and peak RSS are the whole process's (`/usr/bin/time -v`), each stage's time comes from `summary.json` (`run.phases`). One row per run: [`clustering_m12_bench.tsv`](clustering_m12_bench.tsv) (`m12_report.sh <label> bench`).

**Whole runs.** *R's chain, through purity*: from the timer's start (after R has loaded its packages) to the end of `neighborPurity`, reading the input included: R's like for like with Rust's run, which now computes both scores. *R's writing* is R's run's files in the export's layout, which Rust does not have yet (M13). *R ÷ Rust*: R's whole run, then R's chain, against Rust at 1 thread; then R's whole run against Rust at 32 threads; memory, against Rust at 1 thread.

| object | cells | R without RunUMAP: time, peak RSS | R's chain, through purity | R's purity, writing | R with RunUMAP: time | Rust, 1 thread: time, peak RSS | Rust, 32 threads: time, peak RSS | R ÷ Rust, time | R ÷ Rust, memory |
|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 47.1 s, 2.50 GiB | 30.5 s | 0.2 s, 8.1 s | 71.4 s (9.5 s RunUMAP) | **11.3 s, 0.15 GiB** | **3.2 s, 0.18 GiB** | 4.2×, 2.7×; 14.8× | 16× |
| soupx_tenx_10xmousebrain | 11,831 | 215.0 s, 12.02 GiB | 92.8 s | 1.0 s, 113.6 s | 295.0 s (15.7 s RunUMAP) | **38.6 s, 1.44 GiB** | **14.7 s, 1.44 GiB** | 5.6×, 2.4×; 14.6× | 8× |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 211.5 s, 12.08 GiB | 88.2 s | 0.9 s, 115.3 s | 226.4 s (12.7 s RunUMAP) | **39.6 s, 1.44 GiB** | **14.7 s, 1.44 GiB** | 5.3×, 2.2×; 14.4× | 8× |
| soupx_emptydrops_4500871lane | 26,598 | 204.8 s, 23.36 GiB | 163.2 s | 9.7 s, 33.5 s | 240.8 s (32.1 s RunUMAP) | **86.2 s, 1.48 GiB** | **25.1 s, 1.48 GiB** | 2.4×, 1.9×; 8.2× | 16× |
| soupx_emptydrops_sscovelli | 79,409 | 808.4 s, 62.13 GiB | 667.2 s | 192.6 s, 131.9 s | 883.2 s (107.7 s RunUMAP) | **314.4 s, 3.15 GiB** | **66.9 s, 3.19 GiB** | 2.6×, 2.1×; 12.1× | 20× |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 778.0 s, 57.52 GiB | 637.1 s | 180.3 s, 130.7 s | 1010.8 s (117.1 s RunUMAP) | **297.9 s, 3.15 GiB** | **67.3 s, 3.19 GiB** | 2.6×, 2.1×; 11.6× | 18× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 905.9 s, 104.13 GiB | 679.2 s | 74.4 s, 202.8 s | 1144.7 s (82.9 s RunUMAP) | **305.9 s, 2.65 GiB** | **57.4 s, 2.69 GiB** | 3.0×, 2.2×; 15.8× | 39× |

**Purity alone**, R's `neighborPurity` (KMKNN, one thread; M8's direct runs) against Rust's purity stage (both searches, the certificate, the sums):

| object | cells | mean neighbours | R's purity (KMKNN, 1 thread) | Rust, 1 thread | Rust, 32 threads | R ÷ Rust, 1 thread | R ÷ Rust, 32 threads |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 340.4 | 0.18 s | 0.18 s | 0.018 s | 1.0× | 10× |
| soupx_tenx_10xmousebrain | 11,831 | 139.5 | 0.97 s | 1.88 s | 0.102 s | 0.5× | 10× |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 139.6 | 0.94 s | 1.87 s | 0.110 s | 0.5× | 9× |
| soupx_emptydrops_4500871lane | 26,598 | 573.4 | 9.74 s | 9.23 s | 0.368 s | 1.1× | 26× |
| soupx_emptydrops_sscovelli | 79,409 | 448.5 | 192.58 s | 83.46 s | 2.273 s | 2.3× | 85× |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 468.1 | 180.34 s | 80.66 s | 2.807 s | 2.2× | 64× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 386.8 | 74.44 s | 43.64 s | 1.566 s | 1.7× | 48× |

- **C5's criterion is met.** At one thread Rust takes less time than R's whole run on every object (2.4 to 5.6 times less) and than R's chain through purity (1.9 to 2.7 times less); at 32 threads, 8.2 to 15.8 times less than R's whole run. It peaks at 0.15 to 3.19 GiB against R's 2.50 to 104.13 GiB (8 to 39 times less); on the two 79,409-cell objects at 3.15 GiB at one thread and 3.19 GiB at 32, under the 10 GiB bound.
- **Purity's time.** The exhaustive search makes both of KMKNN's searches over every pair, so at one thread it takes about as long as R's KMKNN on the 3,108- and 26,598-cell objects, twice as long on the two 11,831-cell ones (1.9 s against 1.0 s), and 1.7 to 2.3 times less on the 59,323- and 79,409-cell ones (44 to 83 s against 74 to 193 s), where R's KMKNN spends the longest. It runs in parallel, which R's call does not (`BPPARAM = NULL`): at 32 threads it takes 0.02 to 2.8 s, 9 to 85 times less than R's. Purity is 14 to 27% of Rust's one-thread run on the 59,323- and 79,409-cell objects, 3 to 4% at 32 threads.
- **Purity's memory.** It never raises the process's peak: in every run of the benchmark the peak RSS recorded after the purity stage equals the one recorded before it (`run.phases`). It holds the 30 columns it searches (19 MB on 79,409 cells), the per-cell sums (cells × clusters), and per thread a tile of 8 cells' distances to 256 cells.
- **The same files.** The benchmark runs' 17 output files are byte-identical at 1 and 32 threads, and identical to the validation runs' (8 threads, with stage dumps), on all seven objects.
- **The rest of the chain** is M11's code: its stages take 0.86 to 1.24 times their time in M11's benchmark (another node of the same type, shared with other jobs). M11's by-stage table now holds these runs' purity times ([Benchmark](#benchmark)).

### M12 tests

`cargo test --workspace` passes in debug and in release: 267 tests in 39 suites, none failing or ignored, with no compiler warnings, under env.sh's kernel (Haswell). 14 are new in M12:

- `crates/knn` (7): `every_code_path_gives_knncolle_s_bits` (the tiled kernel's AVX-512, AVX2 and scalar paths against the pair-by-pair definition, bit for bit, and its sums equal to `fields::rdist`'s before the square root); `kth_distance_and_radius_by_brute_force` (the k-th distance, the radius and the neighbour lists against a full sort, k capped at n − 1, k = 0, sorted and unsorted lists); `duplicates_ties_and_the_self_exclusion`; `independent_of_the_thread_count` (1, 2, 5 and 16 threads, the certificate included); `the_median_cell_s_kth_neighbour_sits_on_the_radius` (flagged both ways, as `median_pair`, exactly when within the radius); `a_pair_just_above_the_kth_distance_is_flagged` (a step-1 flag, `other`); `the_margin_bounds_the_rounding` (the raw distance's error bound against a long-double computation);
- `crates/cluster`, unit tests (4): a case by hand, R's repeated sums and `max.col`'s ties to the first cluster, unweighted sums, NA with one cluster and absent levels skipped, R's long-double `rowSums` (`1 + 2⁻⁵³ + 2⁻⁵³`);
- `crates/cluster/tests/m12_purity.rs` (3): `neighborPurity` teacher-forced on the two base objects' R embeddings and labels (every K9 array bit for bit, the neighbour lists in BiocNeighbors' order, the radius and the mean neighbour count, the `neighborhood_purity` column, at 1 and 3 threads, and the certificate: the median pair on the odd-sized object, nothing on the even one); purity from the counts on the seven committed tier-T objects against M8's portable runs (the column, `purity_maximum`, the certificate's flags and their kinds, under the fixtures' kernel); other settings (k = 5, 20 and 400, unweighted, dims 1:10) against the definitions written as plain loops.

Extended: `m11_chain`'s thread test now covers purity's columns and stage dumps (1, 2 and 5 threads, every file byte-identical), and its stage-order test the `purity` stage; `cli.rs`'s routes test finds `neighborhood_purity` and `purity_maximum` byte-identical by all four routes, and the purity's summary identical; its validation test on an odd-sized object (267 cells) finds the certificate reported ("2 pair(s) flagged (2 the median pair, 0 other)") and every other check identical; and `--snn` without `--pca` says that neither score is computed.

**Under another kernel** ([`m12_kernel_tests.sbatch`](../../validation/slurm/m12_kernel_tests.sbatch) `auto`, job 61652585, commit 4231635, a Sapphire Rapids node, OpenBLAS's SapphireRapids kernel): 253 passes and 14 failures, each a kernel check naming both kernels: M11's 13 and M12's purity from the counts. The teacher-forced purity tests and every BLAS-free test pass under any kernel (`data/clustering/m12_kernel_tests/SapphireRapids/`).

### M12 findings

- **The median pair is purity's whole dependence on the kernel beyond the last bits.** M8 found 22 neighbour counts in 11 objects, and 19 purities in 10, differing between the kernels. Each is a median pair: on an odd-sized object the median cell's 50th neighbour lies exactly on the radius, and whether `fl(radius · radius)` reaches its squared distance is decided by the embeddings' last bits. Of the 17 odd-sized objects, 11 have the pair inside the radius under one kernel and outside under the other, and they are exactly the 11 objects of M8's 22 counts (M11's floor measurement, `K9_neighbor_count`); the other 6 have it inside under both. The neighbour lists, measured cell by cell, differ across kernels in exactly those cells (12 cells in the 6 tier-T objects whose runs kept the lists, each list one cell longer or shorter; `floor_kernel_m12_lists.tsv`).
- **R's KMKNN keeps every pair the certificate flags.** In every comparison, with R's own embeddings and from the counts, the exhaustive search's K9 equals R's bit for bit, the flagged pairs included: 23 flagged runs among R's 68 instrumented runs, and every direct run.
- **The certificate's margin is set by the data's extent, not by the radius:** `B`, the bound on a cell's distance to a k-means centre, is 3 to 6 times the radius on tier T and 31 to 210 times on tier P, so step 3's band is hundreds to thousands of ulps of the radius squared. Even so, on even-sized objects the nearest pair within the radius lies at least 2.6·10⁸ ulps below it, and no cell's 50th distance has a neighbour within 985 band widths above it: away from the median pair, purity is far from any rounding boundary.
- **No two neighbours of a cell lie at the same squared distance** in any list R kept, so BiocNeighbors' order (equal distances by KMKNN's internal index) and the port's (by cell index) coincide everywhere.
- **The exhaustive search is faster than R's KMKNN** at one thread, and runs in parallel ([M12 benchmark](#m12-benchmark)); its peak memory is below the chain's.

### M12 corrections to the plan

Corrected in [docs/PLAN.md](../PLAN.md) on this branch:

- **C3e, "Where KMKNN could differ":** KMKNN computes every pair it visits as the exhaustive search does, so it can only drop cells; the margin within which its bounds can drop one is not "a few units in the last place" but grows with the data's extent (`B` times about `2(dims + 3)·2⁻⁵³`; 220 to 12,760 ulps of the radius squared on the reference objects), with its derivation here and in `crates/knn/src/purity.rs`; step 1's condition concerns another cell's distance just above the 50th, not the 50th neighbour itself; the median pair; M12's decision point, settled.
- **C3e, steps and memory:** BiocNeighbors caps k at n − 1; the search runs in tiles of 8 cells against blocks of 256 and keeps each cell's 51 smallest squared distances, rather than a buffer of all of a cell's distances.
- **C2:** KMKNN's generator needs no port; the certificate does flag pairs at a rounding boundary (the median pair), which M12's decision settles.
- **C3, commands:** `scrust cluster --purity-lists auto|yes|no`.
- **C3i:** the certificate's flagged pairs are reported.
- **C5, K9:** the checks (the neighbour lists, the weights, sums, totals and maxima; the same from the counts and in the files) and the criterion: flagged pairs reported, not failed, each labelled by kind (M12's decision), in place of "no flagged pair".
- **C6, M12's row:** the saved purities are identical in the production profile, the portable runs' differences graded by the floor, as M8 measured and as A5b, K8 and C3i say (the row said "identical on T and P"); the decision point's outcome. **M17's row:** its comparison of the saved purities covers the flagged pairs. **"Can follow":** the KMKNN port no longer waits on the certificate never firing.
- **C7, purity at a rounding boundary:** the outcome.

### M12 deviations and notes

- **Two binaries.** The decision point's survey used `scrust-m12-survey`, built from M12's tree before its first commit; its flags had no kinds yet and it graded a flagged certificate as differing. The validation and benchmark runs used `scrust-m12-4231635`, with the kinds and the decided grading. Their certificates agree object by object and profile by profile.
- **No new R runs:** M12 compares with M8's runs, and measures the neighbour lists' floor from them (`m12_lists_floor.sh`).
- **`--purity-lists`** is new: `auto` keeps purity's neighbour lists in the stage dumps at level full, as M8's instrumented runs did; `yes` at any level (M8's run of the small tier-P object kept them at level standard).
- **Arrays and scalars R has not.** The Rust run's dumps add `K9_certificate_pairs`, `_step`, `_kind`, `_raw` and `_threshold` (when a pair is flagged) and the scalars `purity_search`, `purity_k_used` and `purity_certificate_flagged*`; the validator compares only what R dumped. R's `purity_index_class` (`KmknnIndex`) is R's own scalar, listed as not compared.
- **A new status.** `scrust validate cluster`'s reports gain "flagged, reported", which passes, and a column for it per stage; the verdict line names the certificates reported.
- **R's `.Random.seed` at the chain's end** (after `neighborPurity`) stays uncompared: `neighborPurity` draws nothing from R's stream (Rcpp's `rng = false`; KMKNN has its own generator), so that state is R's after GroupSingletons, compared under K7 where it drew, or after RunUMAP, which is not ported.
- **The column order of `metadata.tsv`** keeps M11's `silhouette_width`, `silhouette_other`, then `neighborhood_purity` and `purity_maximum`; C3 lists the workflow's columns first.
- **The teacher-forced purity in the validator** uses every CPU the job has; its results do not depend on it.
- **Rust's validation outputs** (`data/clustering/<object>/rust/m12-4231635/`, `m12-4231635-t1/`) and the survey's runs are not listed in the manifest, as in M9 to M11; `clustering_m12.tsv` and the survey's tables summarise them.

### Reproducing M12

From the repository's top directory, after `source env.sh` and `cargo build --release -p scrust` (copy the binary to a file of its own, here `data/clustering/bin/scrust-m12-4231635`):

```bash
export SCRUST_BIN=data/clustering/bin/scrust-m12-4231635
sbatch validation/slurm/m12_certificate.sbatch <label> <objects>
sbatch --partition=shared --cpus-per-task=8 --mem=24G validation/slurm/m12_cluster.sbatch m12-4231635 portable <tier-T objects>
sbatch --partition=shared --cpus-per-task=8 --mem=64G validation/slurm/m12_cluster.sbatch m12-4231635 portable <tier-P objects>
sbatch --partition=test --cpus-per-task=8 --mem=24G validation/slurm/m12_cluster.sbatch m12-4231635 production <tier-T objects>
sbatch --partition=test --cpus-per-task=8 --mem=64G validation/slurm/m12_cluster.sbatch m12-4231635 production <tier-P objects>
sbatch --partition=shared --cpus-per-task=1 --mem=48G validation/slurm/m12_cluster.sbatch m12-4231635-t1 portable <tier-P objects>
bash validation/slurm/m11_threads.sh m12-4231635 m12-4231635-t1 portable <tier-P objects>
VALIDATION_LABEL=m12-4231635 sbatch validation/slurm/m12_benchmark.sbatch m12-4231635 <tier-P objects>
sbatch --partition=test validation/slurm/m12_kernel_tests.sbatch auto
bash validation/slurm/m12_lists_floor.sh > data/clustering/floor_kernel_m12_lists.tsv
bash validation/slurm/m12_report.sh m12-4231635 > docs/validation/clustering_m12.tsv
bash validation/slurm/m12_report.sh m12-4231635 bench > docs/validation/clustering_m12_bench.tsv
```

One run by hand, on a committed fixture (seconds, under env.sh's Haswell kernel): as for M11 ([Reproducing](#reproducing)); `summary.json`'s `result.purity` holds the radius and the certificate, and on `filtered_seurat_emptydrops_test` (267 cells) the certificate flags its median pair, cells 106 and 136, both ways.
