# Clustering in Rust (Step 2; milestones M11, M12, M13 and M14)

M11 runs the workflow's clustering chain end to end in Rust, from the counts to the labels and the silhouette widths (docs/PLAN.md C1, C3, C5, C6): SCTransform at seed 12345, the workflow's PCA guard, RunPCA with 50 PCs, FindNeighbors (Annoy, k = 20, dims 1:30), FindClusters (Louvain, resolution 0.8, seed 12345, singletons grouped) and bluster's `approxSilhouette`. It chains the `sctransform` and `pca` crates of Step 1 with M9's `neighbors` and M10's `cluster`, and adds the silhouette widths, `scrust cluster --counts` with its staged starts, parameters, preset and TOML file, and `scrust validate cluster` for the whole chain. It also adds kernel-gated grading (A5b mechanism 3) to both validators and a kernel check to the fixture tests.

M12 adds the workflow's second score, neighbourhood purity: bluster's `neighborPurity`, with its search done exhaustively with knncolle's arithmetic and a certificate of where KMKNN, the algorithm R uses, could differ. Its outcome comes after M13's ([M12 outcome](#m12-outcome)), its details after M11's ([M12: neighbourhood purity](#m12-neighbourhood-purity)).

M13 adds the rest of the SCT assay, the corrected counts of sctransform's `correct()`, and the export: the input counts, SCTransform's normalised data and the per-cell metadata as plain files that no 10x reader takes for raw counts, with a manifest of every file in every command's `summary.json`. Its outcome comes after M14's ([M13 outcome](#m13-outcome)), its details after M12's ([M13: the corrected counts and the export](#m13-the-corrected-counts-and-the-export)).

M14 adds the loaders: an R package (Seurat, SingleCellExperiment) and a Python package (AnnData) that read any scrust output through its manifest and put each matrix only into a slot meant for its scale. Each loader's object is compared with the saved workflow object (K11), and the workflow's own marker steps run on the Seurat loader's object (K13). Its outcome comes first ([M14 outcome](#m14-outcome)), its details last ([M14: the loaders](#m14-the-loaders)).

This is Step 2's main report, as [doubletfinder.md](doubletfinder.md) is Step 1's; M13 to M17 extend it. The milestones before it have their own: [clustering_m8.md](clustering_m8.md) (the R reference data and the R-vs-R floor), [clustering_m9.md](clustering_m9.md) (Annoy and the SNN graph), [clustering_m10.md](clustering_m10.md) (FindClusters).

## M14 outcome

| Done when (plan C6's M14 row, and the brief) | Result |
|---|---|
| K11, the R loaders (Seurat, SingleCellExperiment), against the saved `.rds` on T and P | **Production profile: identical; portable profile: within the floor.** 136 comparisons (34 objects, 2 profiles, 2 loaders), all passing. In the production profile every component of the Seurat object (1,834) and of the SingleCellExperiment (1,462) is identical to the saved `.rds`: the RNA counts, the SCT layers, the SCTModel (`cells_step1` also on its own) and the SCT assay (also as whole objects), the chain's metadata columns and levels, the identities, the PCA and both graphs. In the portable profile 249 and 201 components differ within `[clustering.floor.kernel]` (the PCA, silhouette widths, a few purities, and on tier P scale.data and the model's residual statistics: exactly M8's measured floor), none beyond it, and every other component is identical ([K11 results](#k11-results)). |
| K11, the Python loader (AnnData), against M8's exports of the saved objects on T and P | **The same.** 68 comparisons, all passing: in the production profile 1,619 components identical to M8's exports of the saved objects; in the portable profile 198 within the floor, none beyond; the scale.data of the three light-level objects, of which M8 kept only a digest, reported as differing from it in the portable profile (graded by the R loaders against the `.rds`). |
| `cells_step1`, drawn again (the coordinating session's decision) | **Identical on every object, and checked against the export.** The R loaders take the seed and `ncells` from the manifest and the cells in `counts/` order, and leave the caller's random-number state as it was. `SCTModel cells_step1`, compared on its own, is identical in all 136 R comparisons: every cell on the 28 objects with at most 5,000 cells, the redraw of 5,000 on the six above them (24 comparisons: both profiles, both R loaders). The self-check passes on every load of those 12 exports, in K11, K13 and M17: vst's step-1 gene draw, made again from the redrawn cells, gives the list the export holds. It costs 2 to 3 s on the largest objects. M17's comparison includes the row. The tests find the caller's state untouched under three generator settings, the draw equal to scrust's own record of it, and tampered exports refused ([notes](#m14-deviations-and-notes)). |
| The exports load exactly what their files hold | **Yes, both profiles.** Every value the loaders read was checked against the files: 3,627,555,688 values in R (formatted again with `%.17g` against the files' text, `.npy` as bytes) and 3,627,049,492 in Python (parsed again with Python's `float()`), none differing. |
| `Read10X` and `read_10x_mtx` fail on `sct/` and on the top level of every export loaded | **Yes, on all 68 exports.** Seurat 5.5.1's `Read10X` fails on every one of the 340 directories but the 68 `counts/`, which it reads as the input; scanpy 1.13.0's `read_10x_mtx`, with every prefix (52,088 attempts), reads only `counts/` (directly, and from the top directory with `prefix="counts/"`), as the input. |
| A file offered to the wrong slot is refused, every combination of role, scale and slot | **Yes.** The R loaders' 18 slots and the Python loader's 9 each take one role, one scale and one format; their tests offer every other one of schema 1 (8,496 and 4,248 wrong combinations) and every one is refused. Tampered manifests (the corrected counts declared raw counts; only their scale changed; the residuals declared raw counts; the counts and log1p layers swapped) are refused by all three builders, by the slot check or by the comment lines' "NOT raw counts" ([M14 tests](#m14-tests)). |
| What is not rebuilt, and why | Listed in every comparison as not rebuilt, since the export does not hold them: the `umap` reduction (not part of the chain), the `commands` (scrust records its run in `summary.json`), other workflow steps' metadata columns (`percent.mt`; scDblFinder's and DoubletFinder's columns on their outputs; the emptyDrops objects' `as.Seurat()` columns), and the RNA assay's v3 class on the emptyDrops path (the counts are identical); in the AnnData, the SCTModel's `cells_step1` above 5,000 cells. The export adds bluster's `silhouette_other` and `purity_maximum` ([what the loaders do not rebuild](#what-the-loaders-do-not-rebuild)). |
| K13 on T and P | **Byte-identical, 68 of 68.** The workflow's marker steps on the Seurat loader's object reproduce the saved `<prefix>_markergenes.csv` of all 34 objects from both profiles' exports, the clusters combined in the saved table's order: 82 to 327 rows on tier T, 22,943 to 81,448 on tier P, up to 5.7 minutes and 10.1 GiB ([K13](#k13-the-markers)). |
| The one-line commands' `.rds` and `.h5ad` reload into identical objects | **Yes, all 264:** every `.rds` from `save_rds()` (Seurat and SingleCellExperiment) and `.h5ad` from `write_h5ad()` (`python -m scrustio h5ad`) rereads into an object identical to the loader's, for the 68 exports and the 20 DoubletFinder outputs (up to 2.1 GB of `.rds` and 3.8 GB of `.h5ad`). |
| Other commands' outputs | **Loaded as the saved objects hold them.** The ten `scrust doubletfinder` outputs (five tier-T and five tier-P inputs, both profiles) load into objects whose singlets' raw counts and DoubletFinder columns, named as DoubletFinder names them, equal the saved DoubletFinder objects' (60 comparisons, R and Python). `scrust sctransform --corrected-counts` loads as the SCT assay alone, `add_scrust()` adds calls and labels to an existing object (the packages' tests). |
| The R loader ready for M17 | **Ready.** `loader_compare.R --loaders seurat --files no --round-trip no` compares one export with its saved `.rds` and writes `k11.tsv` and `compare.json`; it passes on both assay paths, the emptyDrops objects, the objects with singletons and the 99,452-feature object, and takes at most 3.1 minutes (the 99,452-feature object; 4.0 in the first run, on another node of the same type) and 25.8 GiB (the 79,409-cell objects); 6 s and 0.4 GiB on tier T ([the R loader for M17](#the-r-loader-for-m17)). |
| Packaging and envs | The R package installs into the repository's `.rlib/` (gitignored), never into the R env, and checks Seurat 5.5.1 and SeuratObject 5.4.0 (warning with others); the Python package is put on the path by its tests, not installed; `envs/loaders-py.yml` gains pytest 9.2.0, re-solved, the pin recorded and `.envs/loaders-py` recreated from it. |
| Tests | `cargo test --workspace` passes in debug and in release: 281 tests in 41 suites each, none failing or ignored, no compiler warnings (M14 changes no Rust code but two README sentences). The R loader's tests: 72, all passing (13 of them on `cells_step1`); the Python loader's: 12, all passing (in the production R env and in `.envs/loaders-py` rebuilt from the new pin) ([M14 tests](#m14-tests)). |
| Documentation | This report; each package's README; `[clustering.K11]` and `[clustering.K13]` in `validation/tolerances.toml`; M14's section of `validation/manifest.toml`; the loaders in `NOTICE`; `docs/output_manifest.md`'s "For loaders"; the plan's corrections and the M14b milestone ([M14 corrections](#m14-corrections-to-the-plan)). |

## M13 outcome

| Done when (plan C6's M13 row, and the brief) | Result |
|---|---|
| K10 on T and P in both profiles | **Exact.** 222 comparisons of `scrust cluster --counts --export` with every R run of its profile (the 27 tier-T and 7 tier-P objects; instrumented runs, direct runs with and without RunUMAP, the v5 runs; R on Sapphire Rapids and on Cascade Lake), all in exact mode, all passing: 4,200 K10 checks identical, none differing or missing (the corrected counts and their log1p on their data lines, or by R's digest where R kept no copy; R's digest of the counts; the genes and cells; `nCount_SCT` and `nFeature_SCT`; the K10 stage dumps: `correct()`'s genes, median log_umi, residuals beyond the clip range, entries, settings and `.Random.seed` at its entry) ([M13 results](#m13-results)). |
| K10 teacher-forced, under the reference's kernel | **Exact.** `correct()` on R's own model parameters, cell attributes and genes and the counts, in the 98 comparisons with instrumented runs: 964 checks identical (every corrected count and log1p value against R's files, or R's digests; `nCount_SCT`, `nFeature_SCT`; the median log_umi, the residuals beyond the clip range, the entries). |
| The files against the direct runs; the saved objects (`nCount_SCT`, `nFeature_SCT`) | **Identical.** Every direct run's corrected counts on their data lines (and log1p data, where R kept them), with and without RunUMAP; against the saved objects, 2,632 K10 checks identical, none differing, in both profiles on all 34 objects: `nCount_SCT`, `nFeature_SCT`, the corrected counts on their data lines, the genes and cells, R's digests of both layers (graded under the saved objects' kernel: M8 found the corrected counts identical across kernels). |
| The mismatched-kernel pairs pass in tolerance mode; K10's floor recorded | **Yes.** The 68 pairs (each Rust run against the other profile's instrumented R run) are in tolerance mode and pass, K10 identical in every one, teacher-forced too (1,890 and 666 checks): `[clustering.floor.kernel]` records K10 from M8's measurement, identical, with zero bounds. |
| Everything M11 and M12 checked keeps its status | **Yes.** Of the 37,015 checks of M12's 290 reports, 36,337 keep their status and counts in M13's (the labels and levels, the silhouette widths, the purities and certificates, the files, K1 to K9) and the other 678, M12's placeholders for M13 (`nCount_SCT`, `nFeature_SCT`, the K10 arrays), are now compared and identical; none changed. M13's validator re-grading M12's own runs gives the same (290 passing, 36,337 unchanged) ([M13 results](#m13-results)). Without `--export`, M13's `scrust cluster` writes M12's files byte for byte: all 34 objects' validation runs made again, 2,388 files, none differing ([M13 deviations and notes](#m13-deviations-and-notes)). |
| Every file the manifest lists present and described | **Yes.** Every command writes `summary.json`'s `files` section (schema `scrust-files` 1, [`docs/output_manifest.md`](../output_manifest.md)) and checks it before writing it: every file listed once with its role, scale, format, shape, identifiers, size and SHA-256, every listed file present, the naming rules. `scrust validate manifest` passes on all 75 validation runs' exports (the 68 at 8 threads and the 7 at one thread: 6,090 files, 34.6 GB, every checksum recomputed) and on the 45 DoubletFinder runs of Step 1's check; cargo test checks the manifests of `scrust sctransform`, `pca`, `neighbors`, `cluster` and `doubletfinder`. |
| `Read10X` and `read_10x_mtx` fail on the export as designed | **Yes** ([The guard tests](#the-guard-tests)). On three exports (two tier-T objects and the small tier-P object) and on `scrust sctransform --corrected-counts`' outputs, Seurat 5.5.1's `Read10X` fails on every directory but `counts/`, which it reads identical to the input (18 directories, all as required); scanpy 1.13.0's `read_10x_mtx`, with no prefix and with every beginning of every file's path as `prefix=` (2,889 attempts), reads only `counts/` (the directory, and the top directory with `prefix="counts/"`), identical to the input. |
| The export's memory and time measured; C5's target | **Met** ([M13 benchmark](#m13-benchmark)). With `--export`, Rust at one thread takes 2.0 to 3.3 times less time than R's whole run (which makes the same SCT assay and writes its files), at 32 threads 8.6 to 15.3 times less, and peaks at 0.16 to 3.41 GiB against R's 2.50 to 104.13 GiB: 3.34 and 3.35 GiB at one thread and 3.41 GiB at 32 on the 79,409-cell objects (bound: 10 GiB). The export's share: 3.3 to 87.2 s at one thread (24 to 30% of the run; 107% on the two 11,831-cell objects, whose corrected counts are the densest) and 0.5 to 15.8 s at 32 (16 to 36%), nearly all `correct()` and writing `sct/`; copying `counts/` takes 0.2 to 0.7 s and the manifest's checksums 0.1 to 2.0 s. It raises the peak by at most 0.41 GiB at one thread and 0.71 GiB at 32. |
| Independent of the thread count, the export included | **Yes** ([M13 threads](#m13-threads)). 8 and 1 threads on all 34 objects in the portable profile and on the 27 tier-T objects in the production profile, 1, 8 and 32 threads on the 7 tier-P objects in the production profile: every file byte-identical, the export (`counts/`, both SCT layers, `metadata.tsv`, the README) and the K10 dumps included, and `summary.json`'s result section and manifest entries; the 17 one-thread tier-P reports equal the 8-thread runs' check by check. |
| Step 1 unchanged | **Yes** ([Step 1 unchanged](#step-1-unchanged)). `scrust doubletfinder` with M13's binary gives M6b's files byte for byte on M6's ten objects in both profiles (cells.tsv, the sweep tables, `singlets/`, the stage dumps; `summary.json`'s result section), and M7's on all 30 objects M11 re-graded (cells.tsv, the sweep tables, `singlets/`); each of the 45 runs' manifests passes. Step 1's 119 earlier reports re-graded with M13's validator: 15,995 + 140 + 270 checks with matching kernels and 11,896 in tolerance mode, none changing status, as in M11's re-grading. |
| The guard tests' env | `envs/loaders-py.yml` (scanpy 1.13.0, anndata 0.13.4, Python 3.13.18), its solve in `envs/loaders-py.linux-64.pin.txt`, created in `.envs/loaders-py` and used with `PYTHONNOUSERSITE=1` and `python -I`; M14 extends it. |
| Documentation | This report; [`docs/output_manifest.md`](../output_manifest.md), the manifest's schema; `[clustering.K10]`, K10's entries in `[clustering.floor.kernel]` and M13's measurements in `validation/tolerances.toml`; M13's section of `validation/manifest.toml`; sctransform's `correct()` and Seurat's SCT assay in `NOTICE`; the plan's corrections ([M13 corrections](#m13-corrections-to-the-plan)). |
| `cargo test --workspace` passes, debug and release, no warnings | **Yes**: 281 tests, 14 of them new ([M13 tests](#m13-tests)). Under the SapphireRapids kernel, the 17 tests that compare a BLAS-dependent result with a portable fixture stop at the kernel check, and the other 264 pass. |

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
| M13 | sctransform's `correct()` (the corrected counts, their log1p, `nCount_SCT`, `nFeature_SCT`), `scrust sctransform --corrected-counts` and `scrust cluster --export`, the guards against reading normalised data as raw counts, the output manifest of every command (`docs/output_manifest.md`) | this report |
| M14 | the loaders, `loaders/R/scrustio` (Seurat, SingleCellExperiment) and `loaders/python/scrustio` (AnnData): the manifest's rules, the slots and their guard, the SCT assay and its SCTModel as SCTransform makes them, the one-line `.rds` and `.h5ad` commands; K11 and K13 on T and P in both profiles; the R loader's comparison for M17 | this report |

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

## M13: the corrected counts and the export

M13 adds SCTransform's corrected counts, the last part of the SCT assay the workflow's objects hold, and the export: the input counts, SCTransform's normalised data and the per-cell metadata, as plain files that cannot be taken for raw counts (docs/PLAN.md C3f). It also gives every scrust command's `summary.json` a manifest of the files it wrote, in one versioned schema that M14's loaders read ([`docs/output_manifest.md`](../output_manifest.md)).

### What M13 built

| Where | M13 contents |
|---|---|
| `crates/sctransform/src/correct.rs` (new), `residuals.rs`, `vst.rs` | sctransform 0.4.3's `correct()`, made in S2d's pass (`residuals_corrected`, `finish_with`): the corrected counts kept sparse per gene and transposed into compressed columns; their `log1p`, `nCount_SCT`, `nFeature_SCT`, R's digests of the values; the SCT layers written as Matrix Market with their comment lines |
| `crates/sctransform/src/output.rs` | the SCT files' comment lines, the manifest's entries for `scrust sctransform`'s layout, the corrected counts' summary in `result` |
| `crates/scio/src/manifest.rs` (new) | the output manifest: its schema (`scrust-files` 1), the entries and provenance, the digests (parallel over files), the check of a directory against its manifest, the naming rules, the README's file list |
| `crates/scio/src/copy.rs` (new), `sha256.rs`, `mtx.rs`, `npy.rs`, `tenx.rs` | copying files with their checksums; SHA-256 of files, on the SHA extensions where the CPU has them; Matrix Market written from computed values, and compared on data lines; `.npy` shapes from the header |
| `crates/cluster/src/chain.rs`, `chain_output.rs`, `validate.rs` | `correct` in the chain's options; `sct/` with the corrected counts, `nCount_SCT` and `nFeature_SCT` in `metadata.tsv`, the K10 stage dumps in R's layout; the manifest's entries of every file the run wrote; K10 in `scrust validate cluster` (stage dumps, files, digests, the saved object, teacher-forced) and the manifest check |
| `crates/pca`, `crates/neighbors` (`output.rs`) | the manifest's entries of their layouts |
| `crates/scrust` | `scrust sctransform --corrected-counts`; `scrust cluster --export` (from `--counts` only); the manifest in `scrust sctransform`, `pca`, `neighbors`, `cluster` and `doubletfinder`; `scrust validate cluster --counts`; `scrust validate manifest` |
| `envs/loaders-py.yml`, `.linux-64.pin.txt` (new) | the hermetic Python env of the guard test (scanpy 1.13.0), which M14 extends for the Python loader |
| `validation/R/guard_read10x.R`, `validation/python/guard_read10x.py` (new) | the guard tests: Seurat's `Read10X`, scanpy's `read_10x_mtx` with every prefix |
| `validation/slurm` | `m13_cluster.sbatch`, `m13_benchmark.sbatch`, `m13_doubletfinder.sbatch`, `m13_guard.sbatch`, `m13_kernel_tests.sbatch`, `m13_regrade.sbatch`, `m13_noexport.sbatch`, `m13_threads.sh`, `m13_unchanged.sh`, `m13_step1_shared.sh`, `m13_report.sh`, `m13_manifest.sh`, `fixtures_m13.sh` |
| `docs/output_manifest.md` (new) | the manifest's schema, for M14's loaders |

The workspace still uses only the Rust standard library and the conda env's OpenBLAS (the SHA extensions through `std::arch`).

### The corrected counts

Seurat's SCT assay holds, as its `counts` layer, the corrected counts of sctransform 0.4.3's `correct()`, which `vst()` itself calls on its own residuals when Seurat passes `return_corrected_umi = TRUE` (`do.correct.umi`, the default), before it clips them (M8). The port follows `correct()` operation by operation:

1. **The residuals** are vst's own: the Pearson residuals of every model gene (detected in at least 5 cells) and every cell, `(y − μ)/√max(μ + μ²/θ, min_variance)`, as S2d computes them, read before vst clips them at ±√N.
2. **The regressors** are each cell's `[1, log_umi]` with `log_umi` replaced by R's `median()` of it over all cells (`apply(..., 2, function(x) rep(median(x), length(x)))`).
3. **For each bin of 500 genes**, in the counts' order (vst's bins): μ = `exp(tcrossprod(model_pars_fit[genes_bin, -1], regressor_data))`, the OpenBLAS `dgemm` R calls, with R's arguments and shapes; the variance `μ + μ²/θ` (μ where θ = ∞; `μ²` is `μ * μ`, as R computes `x^2`), not floored; the corrected value `μ + residual * sqrt(variance)`.
4. **`round(x, 0)`**, which R 4.5.3's `fround` turns into `nearbyint()`: half to even. Negatives become 0, and `make.sparse()` drops the zeros: positive integers, stored sparse.
5. **The rest of the assay:** `data` is their `log1p` (glibc's, through `rcompat::libm`); adding the assay computes `nCount_SCT` (`colSums`, exact for integers) and `nFeature_SCT` (the entries above 0, every stored entry).

**In S2d's pass, bin by bin.** R holds vst's residual matrix and `correct()`'s output, genes × cells, both dense: 104.1 GiB on cellbender_fromraw_ghirsutum2 (M8). The port makes the corrected counts in the same pass as the residuals, from each residual before it is clipped. For each bin the thread that makes the products first runs R's `dgemm` at the median depth into the free product buffer, keeps it as runs of equal values per gene, and then makes the bin's own product in the same buffer: no second buffer. (Every cell has the same regressors, so in exact arithmetic each gene's row holds one value; the kernel's blocks of columns may differ in their last bits, and the runs keep R's bits wherever they change. Each run's `exp`, variance and square root are computed once.) Each gene's corrected counts are kept as (cell, count) pairs of 32-bit integers, and transposed into compressed columns, Seurat's `dgCMatrix` layout, once every bin is done. No dense genes × cells matrix is made, and the result is placed by gene, never by the order the threads finish, so it does not depend on the thread count. A value that is not a count the port stores (not finite, or above 2³² − 1) stops the run with the gene and the cell named; R would keep it, and its own writer would refuse it.

**Only when asked for.** Clustering does not need them, so they are made only by `scrust cluster --export` and `scrust sctransform --corrected-counts`; every other run, DoubletFinder's included, runs S2d as before.

### The export

`scrust cluster --counts <dir> --export` writes, besides `scrust cluster`'s own files:

```
<out>/
  README.txt                       what each file holds and how to load it with standard readers
  summary.json                     result, run and files: the manifest (docs/output_manifest.md)
  metadata.tsv                     + nCount_SCT and nFeature_SCT, after nFeature_RNA (Seurat's order)
  counts/                          the input's matrix, barcodes and features, copied byte for byte with their
                                   checksums: the export's only 10x directory, and its only raw counts
  sct/                             scrust sctransform's layout (cell_attr.tsv, gene_attr.tsv, model_pars.tsv,
                                   model_pars_fit.tsv, variable_features.tsv, scale_data.npy), and
    corrected_counts.mtx.gz        Seurat's SCT counts layer (integer Matrix Market)
    log1p_corrected_counts.mtx.gz  Seurat's SCT data layer (real Matrix Market)
    sct_features.tsv.gz            the model genes, one name per line, in the counts' order (both layers' rows)
    sct_cells.tsv.gz               the cells, one barcode per line (their columns)
  pca/, neighbors/, the labels' files, optimizer_ids.tsv   as without --export
```

The two layers are written as R's reference runs wrote theirs (`mtx_write_plain`): entries in stored (column-major) order, values as `%.17g` writes them, the field `integer` or `real`. `scrust sctransform --corrected-counts` writes the SCT part alone: its usual files, the four above and `README.txt`.

**Which routes export.** Only a run from the counts. The export copies the input counts, and the corrected counts need vst's residuals of every model gene before clipping, which only SCTransform's own pass has: `scrust sctransform` keeps the variable features' clipped residuals (`scale_data.npy`), and the later stages keep nothing of SCTransform. So `scrust cluster --sct`, `--pca` and `--snn` refuse `--export`, with that reason; `scrust sctransform --corrected-counts` writes the SCT part alone, and a sweep over k or the resolution reuses one SCTransform as before.

### The guards

1. **No 10x names for normalised data.** Only `counts/` holds files named as a 10x reader looks for (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, `genes.tsv`, gzipped or not), and it holds only those; no other file's name ends in `matrix.mtx(.gz)`. The manifest's check (below) enforces both rules on every command's output before `summary.json` is written.
2. **Each normalised matrix says what it is.** Right after the Matrix Market banner:

   ```
   %%MatrixMarket matrix coordinate integer general
   % scrust: SCTransform's CORRECTED counts (Seurat's SCT assay, counts layer). NOT raw counts: do not use them as raw UMI counts.
   % They are sctransform 0.4.3's correct() as Seurat 5.5.1's SCTransform() makes them: every cell's counts corrected to the median sequencing depth (log10 UMIs 2.737986600594553) and rounded.
   % Rows: the model genes, sct_features.tsv.gz; columns: the cells, sct_cells.tsv.gz. The raw counts are elsewhere (an export's counts/).
   % Made by scrust 0.1.0: scrust cluster --counts data/clustering/soupx_tenx_test/input --preset workflow --export
   % Input sha256: 4b8be3ff...  data/clustering/soupx_tenx_test/input/matrix.mtx.gz
   % Input sha256: 53ca5d7a...  data/clustering/soupx_tenx_test/input/barcodes.tsv.gz
   % Input sha256: 9c62c608...  data/clustering/soupx_tenx_test/input/features.tsv.gz
   % Described in summary.json (files) and README.txt.
   ```

   The command leaves out the options that change no value of the file (`--out`, `--threads`, and for `scrust cluster` `--keep-stages`, `--dump-stages`, `--purity-lists`, `--annoy`, `--metadata`), so that the files are the same bytes at any thread count and in any output directory; the manifest keeps the command line as given.
3. **The manifest gives every file a role and a scale**, and the loaders (M14) put a matrix only into a slot meant for its scale.

### The manifest

Every command that writes `summary.json` (`scrust sctransform`, `pca`, `neighbors`, `cluster`, `doubletfinder`) adds its `files` section, the manifest: schema `scrust-files`, version 1, documented in [`docs/output_manifest.md`](../output_manifest.md) for M14's loaders.

- **Per file:** its path, role, scale, format (and whether it is gzipped and has a header), shape (and stored entries), where its row and column identifiers are (a file and its column), a description, its size and its SHA-256. `summary.json` itself is listed without size or checksum.
- **Provenance:** scrust's version, the command line as given, the command as the comment lines give it, each input (the option, the path, every file read with its size and SHA-256), the R functions and packages reproduced, the parameters and the seeds.
- **Completeness:** every file of the output directory, in any subdirectory, is listed once; every listed file exists with its size and checksum; every identifier file is listed; the roles, scales and formats are known; the naming rules hold. Each command checks this before it writes `summary.json`, and fails if anything is off; `scrust validate manifest --dir <dir>` checks any output again, its checksums included.

The checksums of a large export cost seconds: the files are digested in parallel, largest first, and on the SHA extensions where the CPU has them (Sapphire Rapids; not Cascade Lake).

### M13 results

Label `m13-9ce0adb`: binary `data/clustering/bin/scrust-m13-9ce0adb` (sha256 `6eb5c58c...`, built from commit 9ce0adb), [`m13_cluster.sbatch`](../../validation/slurm/m13_cluster.sbatch), 8 threads, M12's design with `--export`: each Rust run from the counts writes the export and dumps its stages at the level of M8's instrumented runs (full on tier T, standard or light on tier P), its manifest is checked (`scrust validate manifest`), and it is compared with every R run of its profile, the instrumented optimizer, the saved object and, on tier P, the workflow job's logged numbers; then with the other profile's instrumented run, a deliberately mismatched pair of kernels. The portable profile ran with `OPENBLAS_CORETYPE=Haswell` on Cascade Lake nodes (`shared`; jobs 61798840 and 61798841, and 61798842 at one thread), the production profile with the kernel detected, SapphireRapids, on Sapphire Rapids (`test`; jobs 61798838 and 61798839). Each validator's teacher-forced K10 runs under the job's kernel. One row per comparison: [`clustering_m13.tsv`](clustering_m13.tsv); one row per object and profile for K10: [`clustering_m13_k10.tsv`](clustering_m13_k10.tsv) ([`m13_report.sh`](../../validation/slurm/m13_report.sh)).

| tier | profile | comparison | comparisons | passed | mode | checks identical | within floor | differing | missing | K10 identical | K10 teacher-forced identical | certificates flagged, reported |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T | portable | the profile's R runs | 108 | 108 | exact | 16,224 | 132 | 0 | 0 | 2,160 | 540 | 54 |
| T | production | the profile's R runs | 82 | 82 | exact | 10,423 | 0 | 0 | 0 | 1,432 | 280 | 36 |
| P | portable | the profile's R runs | 10 | 10 | exact | 1,960 | 14 | 0 | 0 | 238 | 72 | 8 |
| P | production | the profile's R runs | 22 | 22 | exact | 2,852 | 0 | 0 | 0 | 370 | 72 | 4 |
| T | portable | the other kernel's instrumented run | 27 | 27 | tolerance | 5,157 | 492 | 0 | 0 | 756 | 270 | 18 |
| T | production | the other kernel's instrumented run | 27 | 27 | tolerance | 5,157 | 492 | 0 | 0 | 756 | 270 | 18 |
| P | portable | the other kernel's instrumented run | 7 | 7 | tolerance | 1,297 | 160 | 0 | 0 | 189 | 63 | 5 |
| P | production | the other kernel's instrumented run | 7 | 7 | tolerance | 1,297 | 160 | 0 | 0 | 189 | 63 | 5 |

*K10 identical*: the checks of stage K10, from the counts (the stage dumps, the files, the columns, R's digests); the saved objects' are in the K10 table below. *K10 teacher-forced*: `correct()` on R's own model, against instrumented runs only. *Certificates*: M12's (`result.purity.certificate`), unchanged.

**290 comparisons, all passing.** The 222 with the profile's own R runs are in exact mode: 31,459 checks identical, none differing or missing. The only checks within the floor are M12's 146, against the saved objects in the portable profile, graded under the saved objects' kernel: the 118 silhouette widths M11 found and the 28 purities M12 found. The 68 mismatched pairs are in tolerance mode and pass, with M12's 1,304 checks within `[clustering.floor.kernel]`, none beyond it: K10 adds no check within the floor.

**K10.** Over the 222 comparisons, 4,200 K10 checks are identical, none differing or missing: from the stage dumps (`correct()`'s genes, its median log_umi, the residuals beyond the clip range, the entries, its settings, `.Random.seed` at its entry), the files (`sct/corrected_counts.mtx.gz` on its data lines against every R run; `sct/log1p_corrected_counts.mtx.gz` on its data lines where R kept it, tier T, else R's digest of the values; R's digest of the counts' values; `sct_features.tsv.gz`, `sct_cells.tsv.gz`; the layers' summary) and `metadata.tsv`'s `nCount_SCT` and `nFeature_SCT`, against the instrumented runs, the direct runs with and without RunUMAP and the v5 runs. **Against the saved objects**, 2,632 K10 checks identical, none differing, in both profiles on all 34 objects: `nCount_SCT`, `nFeature_SCT`, the corrected counts on their data lines, the genes and cells, R's digests of both SCT layers and the layers' summary (graded under the saved objects' kernel, as M8 found the corrected counts identical across kernels). **Teacher-forced**, in the 98 comparisons with instrumented runs: 964 checks identical (every corrected count and log1p value against R's files or digests, `nCount_SCT`, `nFeature_SCT`, the median, the residuals beyond the clip range, the entries), under the job's kernel, R's. **On the other kernel**, in the 68 mismatched pairs: 1,890 K10 checks identical, and 666 teacher-forced, graded in tolerance mode (the validator's kernel is not R's) and identical: C7's claim that the corrected counts do not depend on the kernel holds on all 34 objects in both directions.

**Everything M12 checked keeps its status.** Every M12 report of these 290 comparisons was compared with M13's, check by check ([`m13_unchanged.sh`](../../validation/slurm/m13_unchanged.sh), `data/clustering/compare/m13_unchanged_m13-9ce0adb.tsv`): of the 37,015 checks M12 graded, 36,337 have the same status and counts in M13's report (the labels and levels, the silhouette widths, the purities and certificates, the files, K1 to K9), and the other 678 are the placeholders M12 left for M13 (`nCount_SCT`, `nFeature_SCT` and the K10 arrays, "not compared"), now compared and identical; none changed, none went missing. **M13's validator grades M12's runs as M12's did.** [`m13_regrade.sbatch`](../../validation/slurm/m13_regrade.sbatch) (job 61799048) graded M12's 290 runs again with M13's binary (`data/clustering/compare/m13_regrade_m13-9ce0adb.tsv`): all 290 pass, the same 36,337 checks unchanged, none changed or missing. M12's runs wrote no corrected counts, so their K10 checks from the counts are n/a (K10 is compared when the run made it: `--export`); K10 teacher-forced, which needs only R's model and the counts, is identical there too (1,630 checks, as in M13's own runs).

#### K10, object by object

The Rust runs from the counts (8 threads); the size of the corrected counts, the median and the residuals beyond the clip range are the same in both profiles, and R's. *R runs*: the comparisons with the profile's R runs and their K10 checks; *teacher-forced*: the K10 checks of `correct()` on R's model, in the comparisons with the profile's instrumented runs (two where the profile has two: on Sapphire Rapids and Cascade Lake in the portable profile, on tier T and soupx_tenx_SRR18256555; with and without the v5 path in the production profile, on filtered_seurat_emptydrops_test and filtered_seurat_emptydrops_sscovelli); *saved object*: the K10 checks against it; *mismatched kernels*: each profile's run against the other's instrumented run, both directions together.

| object | tier | cells | model genes | corrected entries | median log_umi | residuals beyond the clip range | portable: K10 identical (R runs; teacher-forced) | production: K10 identical (R runs; teacher-forced) | saved object: K10 identical (portable, production) | mismatched kernels: K10 identical (teacher-forced), both directions |
|---|---|---|---|---|---|---|---|---|---|---|
| posthocfilt_threshold_doubletfinder_decontx_raw_tenx_test | T | 80 | 694 | 21,613 | 2.8211580513463046 | 91 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_scdblfinder_decontx_raw_tenx_test | T | 87 | 702 | 23,375 | 2.8109042806687006 | 88 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_doubletfinder_decontx_raw_emptydrops_test | T | 92 | 718 | 25,027 | 2.8463366714932414 | 87 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_scdblfinder_decontx_raw_emptydrops_test | T | 112 | 738 | 30,280 | 2.8331452389810172 | 95 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_doubletfinder_cellbender_fromraw_test | T | 135 | 756 | 36,181 | 2.8943160626844384 | 88 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_doubletfinder_soupx_emptydrops_test | T | 137 | 766 | 36,607 | 2.909020854211156 | 83 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_doubletfinder_soupx_tenx_test | T | 142 | 768 | 36,789 | 2.8890140914225233 | 87 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_scdblfinder_soupx_tenx_test | T | 154 | 776 | 40,437 | 2.8978948373479394 | 95 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_scdblfinder_soupx_emptydrops_test | T | 161 | 787 | 43,211 | 2.9057958803678687 | 96 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| posthocfilt_threshold_scdblfinder_cellbender_fromraw_test | T | 163 | 793 | 44,270 | 2.904715545278681 | 94 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| doubletfinder_decontx_raw_emptydrops_test | T | 237 | 786 | 47,488 | 2.6263403673750423 | 83 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| doubletfinder_soupx_emptydrops_test | T | 237 | 799 | 52,325 | 2.767897616018091 | 92 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| doubletfinder_cellbender_fromraw_test | T | 247 | 796 | 53,357 | 2.747411807886423 | 77 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| scdblfinder_decontx_raw_emptydrops_test | T | 262 | 791 | 53,634 | 2.641463923211238 | 88 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| scdblfinder_soupx_emptydrops_test | T | 264 | 818 | 59,847 | 2.79098790652668 | 90 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| decontx_raw_emptydrops_test | T | 267 | 794 | 54,219 | 2.6314437690131722 | 87 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| filtered_seurat_emptydrops_test | T | 267 | 821 | 61,367 | 2.8000293592441343 | 89 | 4 runs: 80; 20 | 4 runs: 80; 20 | 48, 48 | 56 (20) |
| soupx_emptydrops_test | T | 267 | 818 | 60,659 | 2.795184589682424 | 89 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| doubletfinder_decontx_raw_tenx_test | T | 273 | 776 | 51,301 | 2.5877109650189114 | 70 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| doubletfinder_soupx_tenx_test | T | 276 | 802 | 56,694 | 2.7205709510981366 | 79 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| cellbender_fromraw_test | T | 278 | 821 | 62,493 | 2.7792354816302822 | 82 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| scdblfinder_cellbender_fromraw_test | T | 278 | 821 | 62,493 | 2.7792354816302822 | 82 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| scdblfinder_soupx_tenx_test | T | 292 | 812 | 61,303 | 2.737986600594553 | 88 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| scdblfinder_decontx_raw_tenx_test | T | 296 | 781 | 55,821 | 2.5860240170901276 | 75 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| decontx_raw_tenx_test | T | 300 | 782 | 56,279 | 2.5814911852927165 | 70 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| filtered_seurat_tenx_test | T | 300 | 821 | 65,546 | 2.760799148289472 | 90 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| soupx_tenx_test | T | 300 | 815 | 62,883 | 2.737986600594553 | 87 | 4 runs: 80; 20 | 3 runs: 52; 10 | 48, 36 | 56 (20) |
| soupx_tenx_SRR18256555 | P | 3,108 | 15,256 | 2,229,910 | 3.175076696935747 | 644 | 4 runs: 76; 18 | 3 runs: 49; 9 | 44, 33 | 54 (18) |
| filtered_seurat_tenx_10xmousebrain | P | 11,831 | 20,948 | 41,522,923 | 3.987487506133093 | 3,068 | 1 run: 27; 9 | 3 runs: 49; 9 | 11, 33 | 54 (18) |
| soupx_tenx_10xmousebrain | P | 11,831 | 20,926 | 41,030,052 | 3.9824521513849898 | 3,120 | 1 run: 27; 9 | 3 runs: 49; 9 | 11, 33 | 54 (18) |
| soupx_emptydrops_4500871lane | P | 26,598 | 22,831 | 7,373,505 | 2.5728716022004803 | 229 | 1 run: 27; 9 | 3 runs: 49; 9 | 11, 33 | 54 (18) |
| cellbender_fromraw_ghirsutum2 | P | 59,323 | 56,106 | 56,742,959 | 3.0916669575956846 | 50 | 1 run: 27; 9 | 3 runs: 49; 9 | 11, 33 | 54 (18) |
| filtered_seurat_emptydrops_sscovelli | P | 79,409 | 22,025 | 28,919,769 | 2.718501688867274 | 1,067 | 1 run: 27; 9 | 4 runs: 76; 18 | 11, 44 | 54 (18) |
| soupx_emptydrops_sscovelli | P | 79,409 | 22,014 | 27,565,772 | 2.6972293427597176 | 5,455 | 1 run: 27; 9 | 3 runs: 49; 9 | 11, 33 | 54 (18) |

*Residuals beyond the clip range*: residuals vst clips at ±√N that `correct()` uses unclipped; the port reads them before clipping, as R does. The corrected counts hold 21,613 to 56,742,959 entries (206.7 million over the 34 objects).

### The guard tests

[`m13_guard.sbatch`](../../validation/slurm/m13_guard.sbatch) (job 61799098, a Cascade Lake node of `shared`) tried the readers people use on three exports of the validation runs (production profile: filtered_seurat_tenx_test, filtered_seurat_emptydrops_test, soupx_tenx_SRR18256555) and, for each object, on a `scrust sctransform --seed 12345 --corrected-counts` output made in the job:

- **Seurat 5.5.1's `Read10X()`**, in the production R env ([`guard_read10x.R`](../../validation/R/guard_read10x.R); no `.rds` read), on every directory of each output (stage dumps aside): it must fail everywhere but the export's `counts/`, and read `counts/` identical to the input read the same way (dimensions, names, every entry and value).
- **scanpy 1.13.0's `read_10x_mtx()`** (anndata 0.13.4, Python 3.13.18), in `.envs/loaders-py` with `PYTHONNOUSERSITE=1` and `python -I` ([`guard_read10x.py`](../../validation/python/guard_read10x.py)), on the same directories, without `prefix=` and with every beginning of every file's path below the directory as `prefix=`. `read_10x_mtx(path, prefix=p)` reads `<p>matrix.mtx.gz` (or the legacy names), so a prefix can name a file only if it begins that file's path: the attempts are exhaustive.

| export | `Read10X`: directories as required | read | `read_10x_mtx`: attempts | read | as required |
|---|---|---|---|---|---|
| filtered_seurat_tenx_test | 6 of 6 | `counts/` only, identical to the input | 963 | 2: `counts/`, and the top directory with `prefix="counts/"`, both identical to the input | 963 |
| filtered_seurat_emptydrops_test | 6 of 6 | `counts/` only, identical to the input | 963 | 2, as above | 963 |
| soupx_tenx_SRR18256555 | 6 of 6 | `counts/` only, identical to the input | 963 | 2, as above | 963 |

The six directories are the export's top directory, `counts/`, `sct/`, `pca/` and `neighbors/`, and the top directory of `scrust sctransform --corrected-counts`' output. Every other attempt fails: `Read10X` with "Barcode file missing", `read_10x_mtx` with "Did not find file .../matrix.mtx.gz". Per directory and prefix: `data/clustering/<object>/rust/m13-9ce0adb/guard/read10x_seurat.tsv` and `read10x_scanpy.tsv`. cargo test checks the same rules on every committed object's export without either reader: the naming rules, and both readers' file lookup with every prefix ([M13 tests](#m13-tests)).

### Step 1 unchanged

M13 changed code `scrust doubletfinder` runs: SCTransform's residual pass, which now also serves `correct()`, and the manifest every command writes. So [`m13_doubletfinder.sbatch`](../../validation/slurm/m13_doubletfinder.sbatch) ran it again with M13's binary, with the options of the run it is compared with, and compared every file with that run's byte for byte (`summary.json` by its result section; the run section records the run, and M13 adds the files section):

| runs | objects | profile | compared with | files per run | files differing | result sections differing | manifests passing | jobs |
|---|---|---|---|---|---|---|---|---|
| M6's options (stage dumps at the R run's level) | 10 (5 tier T, 5 tier P) | portable and production | M6b's runs (`m6b-13592a6`): cells.tsv, sweep_stats.tsv, bcmvn.tsv, `singlets/`, `stages/` | 42 to 60 | 0 | 0 | 20 of 20 | 61798844 |
| M7's options (no dumps) | 25 of M7's 30 (6 samples × 5 arms) | production | M7's runs (`m7-b8c8a8c`): cells.tsv, sweep_stats.tsv, bcmvn.tsv, `singlets/` | 6 | 0 | 0 | 25 of 25 | 61798931 |
| M6's options, compared with M7's runs as well | the other 5 of M7's 30, in M6's set too | production | M7's runs: every file they have | 6 | 0 | 0 | (above) | [`m13_step1_shared.sh`](../../validation/slurm/m13_step1_shared.sh) |

One row per run: `data/doubletfinder/compare/m13_step1_m13-9ce0adb.tsv`, and `m13_step1_shared_m13-9ce0adb.tsv` for the five objects of both sets (soupx_tenx_SRR18256555, soupx_tenx_10xmousebrain, soupx_emptydrops_4500871lane, soupx_emptydrops_sscovelli, cellbender_fromraw_ghirsutum2): the script runs an object and profile once per label, so these ran with M6's options, and their files and result sections are compared with M7's runs too, every file M7's run has (the stage dumps are the extra files of M6's options). The production runs ran on Sapphire Rapids (`test`), with OpenBLAS's SapphireRapids kernel checked before each.

**Earlier reports re-grade unchanged.** [`m11_doubletfinder.sbatch`](../../validation/slurm/m11_doubletfinder.sbatch) with M13's binary (job 61798855) re-graded the same 119 reports as M11, Step 1's own runs against every R run, the saved outputs and the other kernel's instrumented runs, each check against M6's and M7's reports:

| runs | comparisons | mode | verdicts | checks | checks whose status changed |
|---|---|---|---|---|---|
| M6 (`m6-bc1b48e`), every R run of the run's profile, 10 objects | 49 | exact | 49 pass | 15,995 | 0 |
| M6, the saved outputs | 20 | (not gated) | 20 pass | 140 | 0 |
| M7 (`m7-b8c8a8c`), the saved outputs, 30 objects | 30 | (not gated) | 30 pass | 270 | 0 |
| M6, the other kernel's instrumented run | 20 | tolerance | 20 pass, 674 checks within `[doubletfinder.floor.kernel]` | 11,896 | (no earlier report) |

These are M11's numbers exactly (`data/doubletfinder/compare/m11_verdicts_m13-9ce0adb.tsv`).

### M13 threads

Every output is byte-identical at one thread and at many, the export included:

| tier | profile | runs | threads | files compared per run | differing |
|---|---|---|---|---|---|
| T | portable | 27 | 8 and 1 (`m13_cluster.sbatch`'s second run) | 80, or 85 where a pair is flagged (the export: `counts/`, `sct/` with both layers, `metadata.tsv`; every start's partition; the stage dumps with K10 and the certificate's arrays) | 0 |
| T | production | 27 | 8 and 1 | 80 or 85 | 0 |
| P | portable | 7 | 8 (label m13-9ce0adb) and 1 (m13-9ce0adb-t1, job 61798842), [`m13_threads.sh`](../../validation/slurm/m13_threads.sh) | 78 to 83: the export and the stage dumps | 0 |
| P | production | 7 | 1 and 32 (the benchmark, job 61798843), and 8 (the validation run) | 31: the export (no stage dumps); the 8-thread run's files less its dumps | 0 |

`summary.json` is compared by its `result` section and its manifest's entries: every file's path, role, scale, shape, size and SHA-256 (the manifest's provenance holds the command line, whose `--threads` and `--out` differ). The tier-P portable runs at one thread were validated as the 8-thread runs were, and their 17 reports equal the 8-thread runs' check by check, details included. cargo test compares 1 and 3 threads, with and without the stage dumps ([M13 tests](#m13-tests)). The corrected counts are placed by gene and cell, never in the order the threads finish, and the files are written from them in a fixed order; each matrix's gzip stream is cut into pieces of whole columns (about 131,072 entries each, fixed by the matrix alone), each compressed on its own and written in order.

### M13 benchmark

[`m13_benchmark.sbatch`](../../validation/slurm/m13_benchmark.sbatch) (job 61798843) ran `scrust cluster --counts <object>/input --preset workflow --export`, and the same without `--export`, without stage dumps, on each tier-P object at 1 thread and at 32, in the production profile on one Sapphire Rapids node (holy8a24301, Xeon Platinum 8480CL, partition `test`, 32 CPUs and 96 GiB reserved), with the validation runs' binary, against M8's direct production R runs on the same inputs (R with OpenBLAS at one thread, the same CPU type). Time and peak RSS are the whole process's (`/usr/bin/time -v`), each stage's time comes from `summary.json` (`run.phases`). One row per run: [`clustering_m13_bench.tsv`](clustering_m13_bench.tsv) (`m13_report.sh <label> bench`).

**Whole runs with the export, against R's.** R's direct run computes the whole SCT assay, the corrected counts included (`do.correct.umi`), and writes its files in the export's layout (*R's writing*, part of its time): R's like for like with `--export`. *R ÷ Rust*: R's run without RunUMAP against Rust's with `--export`.

| object | cells | R without RunUMAP: time, peak RSS | R's writing | Rust with `--export`, 1 thread: time, peak RSS | Rust with `--export`, 32 threads: time, peak RSS | R ÷ Rust, time: 1 thread; 32 threads | R ÷ Rust, memory: 1 thread; 32 threads |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 47.1 s, 2.50 GiB | 8.1 s | **14.4 s, 0.16 GiB** | **3.1 s, 0.28 GiB** | 3.3×; 15.3× | 15×; 9× |
| soupx_tenx_10xmousebrain | 11,831 | 215.0 s, 12.02 GiB | 113.6 s | **77.6 s, 1.44 GiB** | **16.9 s, 1.75 GiB** | 2.8×; 12.7× | 8×; 7× |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 211.5 s, 12.08 GiB | 115.3 s | **78.1 s, 1.44 GiB** | **17.1 s, 1.76 GiB** | 2.7×; 12.3× | 8×; 7× |
| soupx_emptydrops_4500871lane | 26,598 | 204.8 s, 23.36 GiB | 33.5 s | **102.0 s, 1.48 GiB** | **23.7 s, 1.56 GiB** | 2.0×; 8.6× | 16×; 15× |
| soupx_emptydrops_sscovelli | 79,409 | 808.4 s, 62.13 GiB | 131.9 s | **335.3 s, 3.34 GiB** | **60.9 s, 3.41 GiB** | 2.4×; 13.3× | 19×; 18× |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 778.0 s, 57.52 GiB | 130.7 s | **329.6 s, 3.35 GiB** | **61.3 s, 3.41 GiB** | 2.4×; 12.7× | 17×; 17× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 905.9 s, 104.13 GiB | 202.8 s | **373.9 s, 3.06 GiB** | **59.8 s, 3.40 GiB** | 2.4×; 15.1× | 34×; 31× |

**The export's share**, each run with `--export` against the same run without it, on the same node. *SCTransform's stage* gains `correct()`, made in S2d's pass, and writing `sct/` (without `--export` or `--keep-stages` the chain writes no `sct/`): the model's tables, `scale_data.npy` and both layers; *`counts/`*: copying the input with its checksums; *the manifest*: the SHA-256 of the export's files (every command digests every file it writes; the export adds `counts/` and `sct/`).

| object | the export: Rust's output, size | 1 thread: without `--export`, with it (difference) | the difference: SCTransform's stage (`correct()`, `sct/`), `counts/`, the manifest | 32 threads: without, with (difference) | the difference: SCTransform's stage, `counts/`, the manifest | peak RSS without, with: 1 thread; 32 threads |
|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 0.11 GB | 11.1 s, 14.4 s (+3.3 s, 30%) | 3.0 s, 0.2 s, 0.1 s | 2.6 s, 3.1 s (+0.5 s, 20%) | 0.4 s, 0.1 s, 0.1 s | 0.15, 0.16 GiB; 0.18, 0.28 GiB |
| soupx_tenx_10xmousebrain | 0.72 GB | 37.4 s, 77.6 s (+40.2 s, 107%) | 39.4 s, 0.4 s, 0.6 s | 12.6 s, 16.9 s (+4.3 s, 34%) | 4.0 s, 0.4 s, 0.2 s | 1.44, 1.44 GiB; 1.44, 1.75 GiB |
| filtered_seurat_tenx_10xmousebrain | 0.72 GB | 37.6 s, 78.1 s (+40.5 s, 108%) | 39.9 s, 0.4 s, 0.6 s | 12.6 s, 17.1 s (+4.5 s, 36%) | 4.1 s, 0.3 s, 0.2 s | 1.44, 1.44 GiB; 1.44, 1.76 GiB |
| soupx_emptydrops_4500871lane | 0.89 GB | 81.2 s, 102.0 s (+20.9 s, 26%) | 19.5 s, 0.4 s, 0.7 s | 20.4 s, 23.7 s (+3.3 s, 16%) | 3.0 s, 0.3 s, 0.5 s | 1.48, 1.48 GiB; 1.48, 1.56 GiB |
| soupx_emptydrops_sscovelli | 2.47 GB | 266.6 s, 335.3 s (+68.7 s, 26%) | 61.3 s, 0.6 s, 2.0 s | 51.6 s, 60.9 s (+9.3 s, 18%) | 8.5 s, 0.5 s, 1.6 s | 3.15, 3.34 GiB; 3.19, 3.41 GiB |
| filtered_seurat_emptydrops_sscovelli | 2.48 GB | 264.9 s, 329.6 s (+64.6 s, 24%) | 56.8 s, 0.6 s, 2.0 s | 51.0 s, 61.3 s (+10.3 s, 20%) | 8.4 s, 0.5 s, 1.6 s | 3.15, 3.35 GiB; 3.19, 3.41 GiB |
| cellbender_fromraw_ghirsutum2 | 2.18 GB | 286.7 s, 373.9 s (+87.2 s, 30%) | 80.1 s, 0.7 s, 1.8 s | 44.0 s, 59.8 s (+15.8 s, 36%) | 14.0 s, 0.6 s, 1.2 s | 2.65, 3.06 GiB; 2.69, 3.40 GiB |

- **C5's target is met.** With the export, Rust at one thread takes 2.0 to 3.3 times less time than R's whole run, at 32 threads 8.6 to 15.3 times less, and peaks at 0.16 to 3.41 GiB against R's 2.50 to 104.13 GiB (7 to 34 times less); on the two 79,409-cell objects at 3.34 and 3.35 GiB at one thread and 3.41 GiB at 32, under the 10 GiB bound.
- **The export's cost is SCTransform's stage:** `correct()` in S2d's pass and writing `sct/` (its `scale_data.npy`, 0.07 to 1.9 GB, is written only with `--export` or `--keep-stages`): 3.0 to 80.1 s at one thread, 0.4 to 14.0 s at 32. Copying `counts/` with its checksums takes 0.2 to 0.7 s, the manifest's digests 0.1 to 2.0 s more. At one thread the export adds 24 to 30% to the run, and 107% on the two 11,831-cell objects, whose corrected counts are the densest: 41 million entries, 17% of the model genes × cells (1.2 to 4.7% on the others), written as 270 MB of gzipped text. At 32 threads it adds 16 to 36%. The rest of each difference (at most 7 s at one thread, on the large objects) is the other stages' variation between runs, mostly purity's. At one thread the whole export, the corrected counts included, takes less time than R's writing of its files alone (3.3 to 87.2 s against 8.1 to 202.8 s).
- **Memory.** The export raises the peak by at most 0.41 GiB at one thread and 0.71 GiB at 32, on cellbender_fromraw_ghirsutum2 (56.7 million corrected entries): the corrected counts are held as 32-bit (cell, count) pairs per gene and transposed once into compressed columns, and each writer holds a batch of compressed pieces per thread; `scale_data.npy` is written from the rows already in memory. At one thread on the 11,831- and 26,598-cell objects the peak does not move: another stage sets it.
- **The same files.** On all seven objects the benchmark's 31 output files are byte-identical at 1 and 32 threads, and identical to the validation runs' (8 threads, with stage dumps) ([M13 threads](#m13-threads)).
- **Without `--export`** the run is M12's chain, byte for byte ([M13 deviations and notes](#m13-deviations-and-notes)); it took 0.85 to 0.98 times M12's benchmark time at one thread and 0.76 to 0.86 times at 32 (another node of the same type).

### M13 tests

`cargo test --workspace` passes in debug and in release: 281 tests in 41 suites, none failing or ignored, with no compiler warnings, under env.sh's kernel (Haswell). 14 are new in M13:

- `crates/scio` (7): `copy::copies_byte_for_byte_through_a_link` (a 9 MiB file and a symbolic link copied byte for byte with their checksums, an existing target refused); `manifest::the_naming_rules` (a normalised matrix under a 10x name, beside one, or ending in `matrix.mtx.gz`, and a 10x name on another file, each caught); `manifest::round_trip_and_checks` (the manifest written and read back, the digests at 1 and 3 threads, an unlisted, a missing and a changed file each caught); `manifest::the_canonical_command`; `mtx::computed_values_and_data_lines` (the writer from computed values gives the same data lines as from a matrix, the comparison on data lines); `sha256::the_sha_extensions_give_the_scalar_digest` (every length from 0 to 300 bytes and up to 1 MiB, the SHA-extension path against the scalar one where the CPU has it); `sha256::a_file_s_digest`;
- `crates/sctransform` (2): `correct::transposition_and_column_statistics` (the transposition, `nCount_SCT`, `nFeature_SCT`, R's digests, `log1p`'s table); `residuals::corrected_counts_match_the_plain_computation_for_any_number_of_threads` (the fused pass against `correct()` written out literally, the full residual matrix and the full product at the median depth element by element, for an odd and an even number of cells, θ = ∞ genes, a one-gene bin, at 1, 2, 3 and 7 threads and three gene groupings; the residuals unchanged);
- `crates/cluster/tests/m13_correct.rs` (1): `correct()` teacher-forced on R's model, cell attributes and genes for the two base objects, against R's corrected counts and log1p data (the SHA-256 of their data lines, formatted independently of the writer), R's digests, `nCount_SCT`, `nFeature_SCT`, the median log_umi, the residuals beyond the clip range and the entries, at 1 and 3 threads;
- `crates/scrust/tests/export.rs` (4): `the_export_reproduces_r_s_corrected_counts` (`scrust cluster --export` on the seven committed objects: K10 against M8's portable runs, `nCount_SCT` and `nFeature_SCT` against R's `metadata.tsv` and in Seurat's column order, `counts/` the input byte for byte with the manifest's checksums equal to the inputs', the manifest valid, the README, the naming rules and the 10x readers' own file lookup: Seurat's `Read10X` and scanpy's `read_10x_mtx` with every prefix of every file's path, only `counts/` readable); `the_export_does_not_depend_on_the_thread_count` (1 and 3 threads, and with the stage dumps: every file byte-identical but the README with the dumps, which lists them; the manifest's entries and the command identical; K10 in the stage dumps, `correct_entry` between `fits_exit` and `vst_exit` as R records it); `only_a_run_from_the_counts_exports`; `sctransform_writes_the_sct_part` (`--corrected-counts` at seed 12345 against R's K10, the files shared with a run without it byte-identical, `--from-stages` refused).

Extended: `cli.rs`' DoubletFinder, `scrust sctransform`, `pca` and `neighbors` tests check each output's manifest (`scrust validate manifest`, and the library's own check with the checksums); the routes test records the new `manifest` phase; the cluster refusals test refuses `--export` from `--snn`. The release suite's first run under M13 found two of M12's tests out of date in that way (the phases; R's graphs copied by `--snn` without `knn_index.tsv`, whose cells the graphs' entries now name from `metadata.tsv`), both fixed before the validation binary was built.

**Under another kernel** ([`m13_kernel_tests.sbatch`](../../validation/slurm/m13_kernel_tests.sbatch) `auto`, job 61798854, commit 9ce0adb, a Sapphire Rapids node, OpenBLAS's SapphireRapids kernel): 264 passes and 17 failures, each a kernel check naming both kernels: M12's 14 and M13's three that compare a BLAS-dependent result with a portable fixture (`scrust cluster --export` and `scrust sctransform --corrected-counts` against R's K10, and `correct()` teacher-forced: its μ is a `dgemm`). The thread test, the refusal and every BLAS-free test pass under any kernel (`data/clustering/m13_kernel_tests/SapphireRapids/`).

### M13 corrections to the plan

Corrected in [docs/PLAN.md](../PLAN.md) on this branch:

- **A2, the loaders' Python env:** M13 created `envs/loaders-py.yml` (scanpy 1.13.0, for the guard test), with its solve in `envs/loaders-py.linux-64.pin.txt` as `rust-dev`'s; M14 extends it.
- **A3:** `summary.json` lists every output file with its role and scale, in the schema of `docs/output_manifest.md`.
- **C3, components (`scio`):** the Matrix Market writers compute values as they write them (the `log1p` layer needs no second matrix); the manifest is `scio::manifest`, with the naming rules; `counts/` copies the input's three 10x files with their checksums, not a whole directory (M8's input directories also hold `assay.txt` and `export.json`).
- **C3f, the corrected counts:** computed only when asked for (`--export`, `--corrected-counts`), in S2d's pass, bin by bin, the product at the median depth in the free product buffer and kept as its runs of equal values; no dense genes × cells matrix and no third product buffer.
- **C3f, the guards:** the comment lines' command leaves out the options that change no value of the file, so that the files are the same at any thread count and in any output directory; the scale of the `data` layer is `log1p_corrected_counts` (the plan wrote `log1p`), and the vocabulary is the schema document's.
- **C3f, the manifest:** versioned (`scrust-files` 1), documented in `docs/output_manifest.md`, checked by each command before it writes it and by `scrust validate manifest`.
- **C3f, which routes export (new):** only a run from the counts, and why; `scrust sctransform --corrected-counts` writes the SCT part alone.
- **C3g:** the loaders read the manifest's schema, `docs/output_manifest.md`.
- **C5, K10:** the checks (the K10 stage values, R's digests where R kept no copy, the saved objects) and the criterion of K10 teacher-forced: graded under the validator's own kernel against R's.
- **C7, corrected counts near .5, and the export's size:** their outcomes.

### M13 deviations and notes

- **The validation binary.** Every M13 run used `scrust-m13-9ce0adb`, built from commit 9ce0adb. The commits after it add this report, the manifest's section, the plan's corrections and, in `validation/tolerances.toml`, which the binary embeds, M13's measurements (`[clustering.K10]`'s `settled`, `[clustering.kernel_gating_shown]`'s `m13`, `measured_m13`, `compared_m13`); the floor's entries the grading reads are the binary's. The Python guard script's version line changed after the binary (it reads scanpy's version through `importlib.metadata`, scanpy 1.13 deprecating `__version__`); it is not part of the binary.
- **No new R runs.** M13 compares with M8's runs; the guard test runs Seurat's `Read10X` on the exports, in the production R env, through sbatch, and reads no `.rds`.
- **Values R would keep and the port does not.** R's `correct()` keeps any value its arithmetic gives; a corrected value that is not finite or exceeds 2³² − 1 stops the port with the gene and cell named (R's own writer refuses non-finite values too). None occurs in the 34 objects.
- **`sct_features.tsv.gz` holds the genes' names**, as R's reference runs write it and as Seurat names the SCT assay's rows; the ids are in `gene_attr.tsv`, row for row.
- **The README lists the directory's files,** so it is the export's only file that differs between runs with and without stage dumps; the comparisons of such runs leave it out, and every other file is the same bytes.
- **Without `--export`, `scrust cluster` writes M12's files.** `metadata.tsv` gains `nCount_SCT` and `nFeature_SCT` after `nFeature_RNA` only with `--export`. [`m13_noexport.sbatch`](../../validation/slurm/m13_noexport.sbatch) (job 61802055) made M12's validation run of each of the 34 objects again with M13's binary and M12's options, in the portable profile: every file byte-identical to M12's run (2,388 files: 70 or 75 per tier-T object, 62 to 67 per tier-P object, the stage dumps included), `summary.json`'s result section identical, each manifest passing (`data/clustering/compare/m13_noexport_m13-9ce0adb.tsv`).
- **`scrust sctransform --corrected-counts`** refuses `--from-stages` (DoubletFinder's validation contexts, whose corrected counts no reference holds).
- **K10 teacher-forced in the validator** reads the counts from `--counts`, else from the R run's `run.input`, else from the candidate's `counts/`; it needs the R run's `stages/vst/` (instrumented runs).
- **The hermetic env** was created with conda's package cache on lab storage and `CONDA_REGISTER_ENVS=false`, so nothing was written under the home directory (checked after the creation).
- **Rust's validation outputs** (`data/clustering/<object>/rust/m13-9ce0adb/`, `m13-9ce0adb-t1/`, `m13-9ce0adb-noexport/`, `data/doubletfinder/<object>/rust/m13-9ce0adb/`) are not listed in the manifest, as in M9 to M12; `clustering_m13.tsv`, `clustering_m13_k10.tsv` and `clustering_m13_bench.tsv` summarise them.

### Reproducing M13

From the repository's top directory, after `source env.sh` and `cargo build --release -p scrust` (copy the binary to a file of its own, here `data/clustering/bin/scrust-m13-9ce0adb`), and, once, the guard test's Python env from its pinned solve (with conda's package cache off the home directory, as `envs/loaders-py.yml` says):

```bash
module load python
```

```bash
CONDA_REGISTER_ENVS=false conda create -p .envs/loaders-py --file envs/loaders-py.linux-64.pin.txt
```

```bash
export SCRUST_BIN=data/clustering/bin/scrust-m13-9ce0adb
sbatch --partition=shared --cpus-per-task=8 --mem=24G validation/slurm/m13_cluster.sbatch m13-9ce0adb portable <tier-T objects>
sbatch --partition=shared --cpus-per-task=8 --mem=64G validation/slurm/m13_cluster.sbatch m13-9ce0adb portable <tier-P objects>
sbatch --partition=test --cpus-per-task=8 --mem=24G validation/slurm/m13_cluster.sbatch m13-9ce0adb production <tier-T objects>
sbatch --partition=test --cpus-per-task=8 --mem=64G validation/slurm/m13_cluster.sbatch m13-9ce0adb production <tier-P objects>
sbatch --partition=shared --cpus-per-task=1 --mem=48G validation/slurm/m13_cluster.sbatch m13-9ce0adb-t1 portable <tier-P objects>
bash validation/slurm/m13_threads.sh m13-9ce0adb m13-9ce0adb-t1 portable <tier-P objects>
VALIDATION_LABEL=m13-9ce0adb sbatch validation/slurm/m13_benchmark.sbatch m13-9ce0adb <tier-P objects>
sbatch validation/slurm/m13_guard.sbatch m13-9ce0adb production filtered_seurat_tenx_test filtered_seurat_emptydrops_test soupx_tenx_SRR18256555
sbatch --partition=test validation/slurm/m13_kernel_tests.sbatch auto
sbatch --partition=test validation/slurm/m13_doubletfinder.sbatch m13-9ce0adb m6 <M6's ten objects>
sbatch --partition=test validation/slurm/m13_doubletfinder.sbatch m13-9ce0adb m7 <M7's 30 objects>
sbatch --partition=shared validation/slurm/m11_doubletfinder.sbatch m13-9ce0adb <M11's objects>
sbatch validation/slurm/m13_regrade.sbatch m13-9ce0adb <the 34 objects>
sbatch validation/slurm/m13_noexport.sbatch m13-9ce0adb <the 34 objects>
bash validation/slurm/m13_step1_shared.sh m13-9ce0adb > data/doubletfinder/compare/m13_step1_shared_m13-9ce0adb.tsv
bash validation/slurm/m13_unchanged.sh m12-4231635 m13-9ce0adb > data/clustering/compare/m13_unchanged_m13-9ce0adb.tsv
bash validation/slurm/m13_unchanged.sh m12-4231635 m13-9ce0adb regrade > data/clustering/compare/m13_regrade_m13-9ce0adb.tsv
bash validation/slurm/m13_report.sh m13-9ce0adb > docs/validation/clustering_m13.tsv
bash validation/slurm/m13_report.sh m13-9ce0adb k10 > docs/validation/clustering_m13_k10.tsv
bash validation/slurm/m13_report.sh m13-9ce0adb bench > docs/validation/clustering_m13_bench.tsv
```

The committed fixtures K10's tests read come from M8's portable instrumented runs: `bash validation/slurm/fixtures_m13.sh`.

One run by hand, on a committed fixture (seconds, under env.sh's Haswell kernel):

```bash
target/release/scrust cluster --counts tests/fixtures/clustering/soupx_tenx_test/input --preset workflow --export --out /tmp/soupx_tenx_test_export
target/release/scrust validate manifest --dir /tmp/soupx_tenx_test_export
```

`README.txt` lists every file and how to read it; `sct/corrected_counts.mtx.gz`'s comment lines say what it is; `summary.json`'s `files` section is the manifest.

## M14: the loaders

M14 builds the loaders of plan C3g: two small packages in `loaders/`, outside the Cargo workspace, that read any scrust output through its manifest ([`docs/output_manifest.md`](../output_manifest.md)) and build the objects R and Python users work with, each matrix only in a slot meant for its scale. They are checked against the saved workflow objects (K11), through the workflow's own marker steps (K13), and by their own tests.

### What M14 built

| Where | M14 contents |
|---|---|
| `loaders/R/scrustio` (new) | the R package ([its README](../../loaders/R/scrustio/README.md)): `read_manifest()`, `manifest_entries()`, `check_manifest()`, `read_entry()`; the slots and their guard (`slot_table()`, `read_for_slot()`); `load_seurat()`, `load_sce()`, `add_scrust()`; `save_rds()`, the one-line `.rds` command; `seurat_feature_names()`; its tests (`tests/scrustio_tests.R`). Imports Matrix and jsonlite; Seurat, SeuratObject and SingleCellExperiment when their objects are built |
| `loaders/python/scrustio` (new) | the Python package ([its README](../../loaders/python/scrustio/README.md)): `read_manifest()`, `slot_table()`, `read_for_slot()`, `load_anndata()`, `write_h5ad()`, `python -m scrustio h5ad` (and `manifest`); its tests (pytest, `tests/`). Needs numpy, scipy, pandas and anndata |
| `validation/R/loader_compare.R`, `lib/loaders.R` (new) | K11 for the R loaders, and M17's comparison: the Seurat object and the SingleCellExperiment against the saved `.rds`, component by component, graded as M11–M13 grade; every value read against the files' text; the `.rds` round trips |
| `validation/python/loader_compare.py` (new) | K11 for the Python loader: the AnnData against M8's export of the saved object; every value read against the files, parsed again independently; the `.h5ad` round trip |
| `validation/R/markers_from_loader.R` (new) | K13: the workflow's `find_markers.R` and `combine_markers.R` steps (copied, with their source cited) on the Seurat loader's object |
| `validation/slurm` | `m14_loaders.sbatch` (K11 for the three loaders, with the 10x readers on each export; `M14_STEPS` runs a subset of its steps), `m14_markers.sbatch` (K13, an array), `m14_m17.sbatch` (M17's comparison, timed), `m14_tests.sbatch` (the packages' own tests), `m14_cargo_tests.sbatch`, `m14_report.sh`, `m14_manifest.sh` |
| `envs/loaders-py.yml`, `.linux-64.pin.txt` | pytest 9.2.0 added for the Python loader's tests; re-solved and the pin recorded |
| `crates/scrust` | the README of `scrust cluster --export` and of `scrust sctransform --corrected-counts` now names the loaders (it said they were coming); no file, role or value changes |
| `third_party/fetch_sources.sh`, `SHA256SUMS` | SingleCellExperiment 1.32.0's source, the version whose classes the R loader fills |
| `.gitignore` | `.rlib/`, the repository's own R library the loader is installed into for its tests |

### The R loader

**The manifest.** `read_manifest(dir)` reads `summary.json`'s `files` section and refuses it if its `schema` is not `scrust-files` or its `version` not 1. Every listed file must exist with its listed size (with `check_sha256 = TRUE`, its SHA-256 too), every path must stay inside the output, and no path may be listed twice. Files are found by their role (`manifest_entries(m, "sct_counts")`), never by name: a test renames the corrected counts and their row file, edits only the manifest, and loads an identical object. Each matrix's row and column names come from the files its `rows` and `columns` name, and their number is checked against its `shape` (and Matrix Market's stored entries against `nnz`). The one file name a loader reads is a levels file's: `cluster_levels` files are tied to their column only by their documented name (`<column>_levels.txt`), as the schema has no field for it ([findings](#m14-findings)).

**The slots.** `slot_table()` lists the 18 slots the R loaders fill from matrix files (9 of Seurat's, 9 of the SingleCellExperiment's), each with the one role, scale and format it takes; the raw-count slots (Seurat's RNA `counts`, the SingleCellExperiment's `counts`) take only `raw_counts`. `read_for_slot()` refuses any other role, scale or format, and reads the Matrix Market file's comment lines: a file marked "NOT raw counts" (as scrust marks every normalised matrix) is refused for a raw-count slot, and an unmarked file for a slot of SCTransform's layers.

**The Seurat object** is built with the calls the workflow's chain makes it with, so that every incidental slot (keys, `misc`, `meta.features`, the computed metadata columns) is Seurat's:

1. the raw counts (`counts/`) through `CreateSeuratObject()`: an `Assay5` named RNA, `orig.ident`, `nCount_RNA` and `nFeature_RNA`;
2. the SCT assay as SCTransform.Assay assembles it: `CreateAssayObject()` on the corrected counts, the variable features, `SetAssayData()` for the `data` layer (read from its file, and checked to be `log1p` of the counts) and for `scale.data` (which puts its rows into the assay's gene order, as Seurat keeps them; the file has them in `VariableFeatures()` order); then the `SCTAssay` class with one `SCTModel` named `model1` and `assay.orig` RNA, added with `[[<-`, which computes `nCount_SCT` and `nFeature_SCT` (checked against the export's);
3. the metadata, typed as Seurat holds them (doubles, R integers, factors in their levels' order), the identities from `seurat_clusters`, and SCT the default assay;
4. the `pca` reduction with `CreateDimReducObject()` as RunPCA.default makes it (key `PC_`, assay SCT, `misc$total.variance` from the summary), and both graphs with `as.Graph()`, assay SCT.

**The SCTModel** is rebuilt from the model files as Seurat's `PrepVSTResults()` assembles it from vst()'s result, step for step: the feature attributes from `gene_attr.tsv` (vst's five columns), `model_pars_fit.tsv` (theta, intercept, slope) and `model_pars.tsv` (the step-1 genes and their fits); the cell attributes from `cell_attr.tsv`; the clip ranges, `set_min_var` and the parameters from the summary; vst()'s 31 recorded arguments with Seurat's values; the median UMI. Two values are not in the files:

- **`cells_step1`**, which cells vst() fitted on: all of them up to `ncells` (5,000), else `sample(colnames(umi), 5000)`, the first draw after SCTransform's `set.seed()`. The export does not list them; it records the seed and `ncells`. As the coordinating session decided ([notes](#m14-deviations-and-notes)), the loader draws them again rather than reading them:
  - `n_cells = min(ncells, cells)` and the seed come from the manifest's parameters. The summary's `n_cells_step1` and the manifest's `seeds` must agree with them.
  - The cells are taken in `counts/` order, which is vst's `colnames(umi)`, and the draw happens only when `n_cells` is below the number of cells.
  - It uses R's own `sample()` under R's default generator kinds, set for the draw. The caller's kinds and `.Random.seed` are restored afterwards, or `.Random.seed` removed again if there was none.
  - **The self-check.** When the output holds the raw counts, the loader also makes vst's next draw from the same stream: the step-1 genes, chosen among the model's genes detected in at least `min_cells` of the drawn cells, which the export lists in sampling order (`model_pars.tsv`). It uses sctransform's own `row_gmean()` and `row_var()`, and refuses the output if the list differs. That needs no new file and adds a few seconds on the largest objects.

  K11 compares `cells_step1` on its own with every saved object, and M17's comparison includes it.
- **`verbosity`**, which records the R call's `verbose` (0 for the workflow's `verbose = FALSE`): the loader writes the workflow's 0.

**Feature names** follow Seurat's rules: `Read10X()`'s `make.unique()`, then `_` and `|` replaced by `-` (`CheckFeaturesNames()`); every list of features of the output (the SCT genes, the variable features, the step-1 genes, the PCA's genes) is mapped onto the counts' rows by name, or by id where names repeat. The workflow's objects need neither rule; the tests' fixtures do not either.

**The SingleCellExperiment** holds the raw counts as `counts` (rows named by the features' names, made unique; `rowData` the ids, names and types), the SCT layers as `altExp(sce, "SCT")` (assays `corrected_counts` and `logcounts`; the SCTModel's feature attributes as its `rowData`, its cell attributes as its `colData`, its other slots and the variable features in its `metadata()`), the residuals as `altExp(sce, "SCT_scale_data")` (assay `scaledata`, in `VariableFeatures()` order), the metadata as `colData`, the PCA as `reducedDim(sce, "PCA")` (attributes `rotation`, `stdev`, `total.variance`) and both graphs as `colPairs`.

**Other outputs.** `scrust sctransform --corrected-counts` loads as the SCT assay alone (it has no raw counts); `scrust doubletfinder` as the singlets' raw counts with DoubletFinder's calls as metadata columns, named as DoubletFinder names them (`pANN_<pN>_<pK>_<nExp>`, `DF.classifications_<pN>_<pK>_<nExp>`, from the summary). An output with neither raw counts nor SCT layers (`scrust cluster` without `--export`, `scrust pca`, `scrust neighbors`) is refused by the object builders; `add_scrust()` adds any output's per-cell columns (the calls, the cluster labels with their levels, the scores) to an object the user already has, by barcode.

**Writing.** `save_rds(dir, file, object = "seurat" | "sce")` builds the object and saves it: `Rscript -e 'scrustio::save_rds("export", "object.rds")'`. **Versions:** the loader checks Seurat (5.5.1) and SeuratObject (5.4.0) and warns, carrying on, with others. **Install:** `R CMD INSTALL --library=.rlib loaders/R/scrustio` in the repository's own library, never into the R env.

### The Python loader

The same manifest rules, slots (9: `layers["counts"]`, `X`, `layers["sct_corrected_counts"]`, `obsm["sct_scale_data"]`, `obsm["X_pca"]`, `varm["PCs"]`, `uns["pca"]["stdev"]`, `obsp["sct_nn"]`, `obsp["sct_snn"]`) and content rule. `load_anndata(dir)` builds:

- `layers["counts"]`: the raw counts, cells x genes;
- `X`: the log1p corrected counts (Seurat's `data`), with scanpy's `uns["log1p"] = {"base": None}`, and `layers["sct_corrected_counts"]`, both over every gene, 0 outside SCTransform's model (`var["sct_modelled"]`);
- `obsm["sct_scale_data"]`: the residuals as a data frame (cells x variable features, in `VariableFeatures()` order), with `var["sct_variable"]`;
- `uns["sct_model"]`: the SCTModel's tables and values, as the R loader rebuilds them, but for `cells_step1` above 5,000 cells (R's generator is not reproduced);
- `obs`: the metadata, R integers as int32 and factors as categoricals in their levels' order;
- `obsm["X_pca"]`, `varm["PCs"]` (0 for genes outside the PCA), `uns["pca"]` (Seurat's `stdev`, `variance`, `variance_ratio`, `total_variance`, and `features`, the loadings' genes in order);
- `obsp["sct_nn"]`, `obsp["sct_snn"]`, and `uns["sct_snn"] = {"connectivities_key": "sct_snn", "params": ...}` for scanpy's `neighbors_key`.

Text columns with repeated values are made categorical as anndata makes them when it writes `.h5ad` (`strings_to_categoricals()`), so a written and reread object equals the loader's. `python -m scrustio h5ad <dir> <file.h5ad>` is the one-line command. Matrix Market files are read with scipy's `mmread` (fast_matrix_market, correctly rounded) and tables with pandas' `float_precision="round_trip"`.

**numpy's `log1p` is not R's.** Seurat's `data` layer holds the C library's `log1p` of the corrected counts; numpy 2.5.5's `np.log1p` (vectorised code) differs from it in the last bit for 13,457 of the 65,546 values of a tier-T export, while `math.log1p` (the C library's) matches all of them. So the loader reads `X` from the file and checks it with `math.log1p` (per distinct value), and Python users should not recompute Seurat's `data` with numpy.

### How K11 compares

The exports loaded are M13's validation runs (label `m13-9ce0adb`: `scrust cluster --counts ... --preset workflow --export`, 8 threads, with stage dumps), 34 objects in both profiles: 68 exports, made by the binary of commit 9ce0adb, whose code main still has. `validation/slurm/m14_loaders.sbatch` runs, for each:

- **the R loaders** (`loader_compare.R`): the Seurat object and the SingleCellExperiment against the saved `.rds`, read directly. Compared, with `identical()` and statistics: the RNA counts (every entry and dimname; the RNA assay as a whole where the saved one is an `Assay5`); `orig.ident` and the chain's metadata columns (`nCount_RNA`, `nFeature_RNA`, `nCount_SCT`, `nFeature_SCT`, `SCT_snn_res.0.8` and `seurat_clusters` with their levels, `silhouette_width`, `neighborhood_purity`); the identities; the default assay; the SCT layers (`counts`, `data`, `scale.data` with its row order), the variable features and the assay's other slots; the SCTModel slot by slot and column by column (and as a whole), with `cells_step1` also on its own (the loaders draw it again); the PCA (embeddings, dims 1:30 and 31:50, the loadings and their genes, the standard deviations, the total variance, the other slots); both graphs (as `Graph` objects, `assay.used` included); the project name and SeuratObject version. The SingleCellExperiment's counterparts likewise (its `colPairs` as sparse matrices);
- **the Python loader** (`loader_compare.py`): the AnnData against M8's export of the saved object (`data/clustering/<object>/saved/`, and `input/` for the RNA counts): the same components; where M8 kept only R's digest (scale.data at the standard and light levels, the log1p layer at those levels) by that digest;
- **the files**: the loaders load exactly what the files hold. In R every value read from text is formatted again with `%.17g` and compared with the file's text (one-to-one on doubles), entry by entry for the Matrix Market files (row, column and value, in stored order), column by column for the tables; `scale_data.npy` as bytes; the summary's numbers by R's parse against jsonlite's. In Python every value is parsed again from the files with Python's own `float()` and compared bit for bit, and `.npy` as bytes;
- **the round trips**: `save_rds()` (Seurat and SingleCellExperiment) and `python -m scrustio h5ad`'s `write_h5ad()`, reread with `readRDS()` and `anndata.read_h5ad()`, against the loader's object (`identical()` in R; every slot, dtype and category in Python);
- **the 10x readers**: Seurat's `Read10X` and scanpy's `read_10x_mtx` (with every prefix) on every directory of the export, as M13's guard (`guard_read10x.R`, `guard_read10x.py`): only `counts/` may be read, as the input.

**Grading** is M11–M13's. Exact mode when the export's recorded OpenBLAS kernel is the saved objects' (SapphireRapids): every component must be identical. Floor mode otherwise (the portable profile's Haswell): the BLAS-free components (the counts, `nCount_RNA`, `nFeature_RNA`, the step-1 statistics `detection_rate`, `gmean`, `variance`, the model's cell attributes, the slots' layout) stay exact, and the BLAS-dependent ones may differ within `[clustering.floor.kernel]` of `validation/tolerances.toml`, as measured by M8 and M11: scale.data (`S2d_scale_data`), the model's `residual_mean` and `residual_variance` together (`S2d_residual_stats`), the embeddings (`S3_embeddings_dims_1_30`, `_31_50`, each component's sign aligned first), the loadings (`S3_loadings`, aligned), the standard deviations (`S3_sdev`), the silhouette widths (`K8_silhouette_width`) and the purities (`K9_purity`). What the floor does not list stays exact: the corrected counts and their log1p, `nCount_SCT`, `nFeature_SCT`, the labels, levels and identities, both graphs, the total variance, theta and the regularised parameters. A whole object (the SCT assay, the SCTModel, the PCA) is graded by its parts.

**DoubletFinder.** The same scripts with `--kind doubletfinder` load the `scrust doubletfinder` outputs of M13's Step 1 check (label `m13-9ce0adb`, the five tier-T and five tier-P inputs, both profiles) and compare the singlets' raw counts and DoubletFinder's two metadata columns with the saved DoubletFinder objects (R: the `.rds`; Python: M6b's and M7's export, `data/doubletfinder/<input>/saved/`).

### K11 results

Label `m13-9ce0adb` (M13's exports, binary `scrust-m13-9ce0adb`); [`m14_loaders.sbatch`](../../validation/slurm/m14_loaders.sbatch), on `shared` and `serial_requeue` nodes (Cascade Lake, Xeon Platinum 8268: no step uses BLAS), the R loaders installed in `.rlib/` (their comparisons run a second time with the final loader, `M14_STEPS=r`: [notes](#m14-deviations-and-notes)), the Python loader in an env solved from `envs/loaders-py.yml` identically to the recorded pin. One row per object, profile and loader: [`clustering_m14.tsv`](clustering_m14.tsv) ([`m14_report.sh`](../../validation/slurm/m14_report.sh)).

| outputs | tier | profile | loader | objects | passing | components identical | within the floor | beyond | differing | not rebuilt | added | values read, checked against the files | differing | round trips identical |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| exports | P | portable | seurat | 7 | 7 | 296 | 81 | 0 | 0 | 22 | 7 | 1,597,635,968 | 0 | 7 of 7 |
| exports | P | portable | sce | 7 | 7 | 241 | 60 | 0 | 0 | 0 | 0 | (Seurat's) |  | 7 of 7 |
| exports | P | portable | anndata | 7 | 7 | 263 | 57 | 0 | 0 | 20 | 0 | 1,597,593,954 | 0 | 7 of 7 |
| exports | P | production | seurat | 7 | 7 | 377 | 0 | 0 | 0 | 22 | 7 | 1,597,635,968 | 0 | 7 of 7 |
| exports | P | production | sce | 7 | 7 | 301 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 7 of 7 |
| exports | P | production | anndata | 7 | 7 | 323 | 0 | 0 | 0 | 20 | 0 | 1,597,593,954 | 0 | 7 of 7 |
| exports | T | portable | seurat | 27 | 27 | 1,289 | 168 | 0 | 0 | 82 | 27 | 11,044,969 | 0 | 27 of 27 |
| exports | T | portable | sce | 27 | 27 | 1,020 | 141 | 0 | 0 | 0 | 0 | (Seurat's) |  | 27 of 27 |
| exports | T | portable | anndata | 27 | 27 | 1,155 | 141 | 0 | 0 | 54 | 0 | 10,993,321 | 0 | 27 of 27 |
| exports | T | production | seurat | 27 | 27 | 1,457 | 0 | 0 | 0 | 82 | 27 | 11,044,969 | 0 | 27 of 27 |
| exports | T | production | sce | 27 | 27 | 1,161 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 27 of 27 |
| exports | T | production | anndata | 27 | 27 | 1,296 | 0 | 0 | 0 | 54 | 0 | 10,993,321 | 0 | 27 of 27 |
| doubletfinder | P | portable | seurat | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 204,823,477 | 0 | 5 of 5 |
| doubletfinder | P | portable | sce | 5 | 5 | 15 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 5 of 5 |
| doubletfinder | P | portable | anndata | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 204,665,311 | 0 | 5 of 5 |
| doubletfinder | P | production | seurat | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 204,823,477 | 0 | 5 of 5 |
| doubletfinder | P | production | sce | 5 | 5 | 15 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 5 of 5 |
| doubletfinder | P | production | anndata | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 204,665,311 | 0 | 5 of 5 |
| doubletfinder | T | portable | seurat | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 273,430 | 0 | 5 of 5 |
| doubletfinder | T | portable | sce | 5 | 5 | 15 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 5 of 5 |
| doubletfinder | T | portable | anndata | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 272,160 | 0 | 5 of 5 |
| doubletfinder | T | production | seurat | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 273,430 | 0 | 5 of 5 |
| doubletfinder | T | production | sce | 5 | 5 | 15 | 0 | 0 | 0 | 0 | 0 | (Seurat's) |  | 5 of 5 |
| doubletfinder | T | production | anndata | 5 | 5 | 25 | 0 | 0 | 0 | 0 | 0 | 272,160 | 0 | 5 of 5 |

*Components*: the rows of `k11.tsv` (Seurat 58, SingleCellExperiment 43, AnnData 49 to 50 per export; 5, 3 and 5 per DoubletFinder output). *Values read*: every value the Seurat loader (R) or the AnnData loader (Python) read from the export's text files and `.npy` files, checked against the files; the SingleCellExperiment reads the same files through the same readers. *Not rebuilt* and *added*: [below](#what-the-loaders-do-not-rebuild).

**The production profile: identical objects.** In all 102 comparisons of the 34 exports, every component is identical to the saved object: the RNA counts (explicit zeros included), the SCT layers, scale.data in the SCT assay's gene order, the variable features, the SCTModel column by column and as a whole object, with `cells_step1` also on its own (every cell on the 28 objects with at most 5,000 cells, the redraw on the six above), the SCT assay as a whole object, `orig.ident`, the chain's metadata columns with their types and levels, the identities, the PCA (every slot, and as a whole object), both graphs as `Graph` objects; the SingleCellExperiment's counterparts and the AnnData's against M8's exports. The RNA assay is identical as a whole object wherever the saved one is an `Assay5` (32 objects).

**The portable profile: the files exactly, the saved objects within the floor.** Every value the loaders read is the files' (the files check), and the objects differ from the saved ones only in the BLAS-dependent components M8 measured, each within `[clustering.floor.kernel]`; every other component is identical, the corrected counts, their log1p, `nCount_SCT`, `nFeature_SCT`, the labels, levels, identities, graphs, theta, the regularised parameters, the total variance and the step-1 statistics included. Over the 34 objects, for the Seurat loader (the SingleCellExperiment and the AnnData give the same statistics, the AnnData's scale.data on the four objects whose export M8 kept in full):

| floor entry | tier | objects within it | beyond | values differing (of those compared in these objects) | largest absolute difference | largest scaled difference | the floor's bounds |
|---|---|---|---|---|---|---|---|
| `K8_silhouette_width` | P | 7 | 0 | 269,054 (271,509) | 3.26e-13 | 3.48e-13 | max_abs 3.3e-13, fraction_differing 0.99 |
| `K8_silhouette_width` | T | 27 | 0 | 5,871 (5,904) | 2.41e-14 | 9.58e-14 | max_abs 3.3e-13, fraction_differing 0.99 |
| `K9_purity` | P | 4 | 0 | 8 (229,972) | 0.0187 | 0.0187 | cells_differing 19, objects 10, max_abs 0.019 |
| `K9_purity` | T | 6 | 0 | 11 (1,062) | 0.0142 | 0.0142 | cells_differing 19, objects 10, max_abs 0.019 |
| `S2d_residual_stats` | P | 7 | 0 | 1,430 (360,212) | 1.42e-14 | 3.2e-17 | max_abs 1.5e-14, max_scaled 3.21e-17 |
| `S2d_scale_data` | P | 7 | 0 | 3,559,615 (814,527,000) | 4.26e-14 | 9.54e-16 | max_abs 4.3e-14, max_fraction_differing 0.0072 |
| `S3_embeddings_dims_1_30` | P | 7 | 0 | 8,108,432 (8,145,270) | 1.53e-10 | 2.02e-12 | max_abs 1.6e-10, max_scaled 2.4e-12 |
| `S3_embeddings_dims_1_30` | T | 27 | 0 | 176,267 (177,120) | 3.66e-11 | 2.34e-12 | max_abs 1.6e-10, max_scaled 2.4e-12 |
| `S3_embeddings_dims_31_50` | P | 7 | 0 | 5,428,912 (5,430,180) | 1.36e-10 | 2.31e-12 | max_abs 4.3e-06 |
| `S3_embeddings_dims_31_50` | T | 27 | 0 | 118,067 (118,080) | 4.27e-06 | 5.69e-07 | max_abs 4.3e-06 |
| `S3_loadings` | P | 7 | 0 | 1,045,743 (1,050,000) | 1.56e-12 | 2.96e-12 | max_abs 1.2e-07, max_scaled 4.5e-07 |
| `S3_loadings` | T | 27 | 0 | 1,054,462 (1,057,600) | 1.19e-07 | 4.49e-07 | max_abs 1.2e-07, max_scaled 4.5e-07 |
| `S3_sdev` | P | 7 | 0 | 333 (350) | 4.44e-14 | 1.8e-15 | max_abs 1.9e-13, max_scaled 3.6e-14 |
| `S3_sdev` | T | 27 | 0 | 1,302 (1,350) | 1.85e-13 | 3.56e-14 | max_abs 1.9e-13, max_scaled 3.6e-14 |

These are M8's measurements of the portable R runs against the saved objects, object by object ([clustering_m8.md](clustering_m8.md#portable-runs-against-the-saved-objects)): the 19 purities in 10 objects, the 3,559,615 scale.data values of tier P, the 1,430 residual statistics, the 269,054 tier-P silhouette widths. The portable exports equal R's portable runs (M13), so loading them reproduces that floor and nothing else. On the three light-level objects the AnnData's scale.data can be compared only with M8's digest of the saved one, which differs (reported, not graded, as M11-M13 report digests in tolerance mode); the R loaders compare those three element by element with the saved `.rds`, within the floor.

**The 10x readers** fail, as required, on every directory of all 68 exports but `counts/`, which both read as the input:

| tier | profile | exports | Read10X: directories tried, as required | read_10x_mtx: attempts (every prefix), exports all as required |
|---|---|---|---|---|
| P | portable | 7 | 35, 35 | 5,362, 7 of 7 |
| P | production | 7 | 35, 35 | 5,362, 7 of 7 |
| T | portable | 27 | 135, 135 | 20,682, 27 of 27 |
| T | production | 27 | 135, 135 | 20,682, 27 of 27 |

Each export has five directories (stage dumps aside): the top directory, `counts/`, `sct/`, `pca/` and `neighbors/`; scanpy's attempts cover every beginning of every file's path as `prefix=` (M13's 963 attempts per object also covered a `scrust sctransform` output).

**The round trips:** all 264 `.rds` and `.h5ad` files written by the one-line commands reread into objects identical to the loaders' (`identical()` in R; every slot, dtype and category in Python), the largest `.rds` 2.1 GB and `.h5ad` 3.8 GB.

**DoubletFinder.** The ten `scrust doubletfinder` outputs (five tier-T, five tier-P inputs) load in both profiles into objects whose singlets' raw counts and DoubletFinder columns (`pANN_<pN>_<pK>_<nExp>`, `DF.classifications_<pN>_<pK>_<nExp>`) equal the saved DoubletFinder objects', in all 60 comparisons (pANN stays identical across kernels, as M11 found).

### What the loaders do not rebuild

The export does not hold these, so the loaders do not rebuild them; every comparison lists them (status "not rebuilt"):

| part of the saved object | objects | why |
|---|---|---|
| the `umap` reduction | 34 | RunUMAP is not part of the chain scrust reproduces; it changes nothing downstream (M8) |
| `commands` (SCTransform.RNA, RunPCA.SCT, RunUMAP.SCT.pca, FindNeighbors.SCT.pca, FindClusters) | 34 | the R calls' records; scrust records its run (command line, parameters, versions, inputs' checksums) in `summary.json` |
| metadata columns of other workflow steps | 34 | `percent.mt` (all 34); the emptyDrops objects' `nCount_originalexp`, `nFeature_originalexp`, `Sample`, `Barcode`, `sizeFactor` (`as.Seurat()`, 2); `scDblFinder.class` (the 5 scDblFinder outputs and their 5 post-hoc objects); DoubletFinder's `pANN_*` and `DF.classifications_*` (the 5 DoubletFinder outputs and their 5 post-hoc objects; the `scrust doubletfinder` outputs hold them, and `add_scrust()` adds them) |
| the RNA assay's class on the emptyDrops path | 2 | the saved RNA assay is a v3 `Assay` made by `as.Seurat()` of a SingleCellExperiment, with a `data` layer and `meta.features` of its own; the loader makes the `Assay5` `CreateSeuratObject()` makes, with the identical counts |
| `cells_step1` of the SCTModel, in the AnnData | 6 | above 5,000 cells vst() drew the step-1 cells with R's `sample()`; the Python loader does not reproduce R's generator (the R loader does) |

The export also adds two columns the workflow drops: bluster's `silhouette_other` and `purity_maximum` (status "added"). The saved objects' `tools`, `misc`, `neighbors` and `images` are empty, as the loader's.

### K13: the markers

[`m14_markers.sbatch`](../../validation/slurm/m14_markers.sbatch), arrays over the 34 objects in each profile (`shared` and `serial_requeue`, Cascade Lake), run a second time with the final loader. One row per object and profile: [`clustering_m14_k13.tsv`](clustering_m14_k13.tsv).

**68 of 68 combined tables equal the saved `<prefix>_markergenes.csv` byte for byte**, from the portable and the production exports alike (FindMarkers reads the SCT `data` layer and the identities, which do not depend on the kernel). Tier T: 27 objects, 82 to 327 rows, under a minute and 0.43 GiB each. Tier P:

| object | clusters | rows | identical (portable, production) | minutes | loading the export, s | peak RSS, GiB |
|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 15 | 22,943 | TRUE, TRUE | 0.4, 0.4 | 13, 13 | 0.93, 0.93 |
| soupx_tenx_10xmousebrain | 21 | 81,448 | TRUE, TRUE | 4.5, 4.8 | 114, 113 | 5.94, 5.94 |
| filtered_seurat_tenx_10xmousebrain | 21 | 81,317 | TRUE, TRUE | 3.8, 4.0 | 80, 80 | 5.99, 5.96 |
| soupx_emptydrops_4500871lane | 31 | 35,901 | TRUE, TRUE | 1.6, 1.5 | 53, 48 | 3.88, 3.86 |
| soupx_emptydrops_sscovelli | 36 | 57,493 | TRUE, TRUE | 4.8, 4.5 | 110, 100 | 8.48, 8.48 |
| filtered_seurat_emptydrops_sscovelli | 37 | 58,776 | TRUE, TRUE | 5.0, 4.6 | 105, 98 | 8.51, 8.51 |
| cellbender_fromraw_ghirsutum2 | 22 | 44,222 | TRUE, TRUE | 5.7, 5.4 | 139, 135 | 10.13, 10.13 |

The rows and the saved cluster orders are M8's marker references'.

### The R loader for M17

M17 compares each object's Rust export with its saved `.rds` through the Seurat loader: `loader_compare.R --object <name> --export <dir> --rds <saved .rds> --out <dir> --loaders seurat --files no --round-trip no`, which writes `k11.tsv` (one row per component) and `compare.json` (the mode, the counts by status, each step's time and peak RSS), after which the export's large files can go. It handles every kind of object in the panel: both assay paths (v5 and the emptyDrops objects' v3), objects with singletons (the six tier-P objects with 1 to 40), and the 99,452-feature object. Timed with [`m14_m17.sbatch`](../../validation/slurm/m14_m17.sbatch) on the production exports (one Cascade Lake node, Xeon Platinum 8268, one object after another, with the final loader), every comparison passing ([`clustering_m14_m17.tsv`](clustering_m14_m17.tsv)):

| object | cells | features | wall time | peak RSS | loading the export | reading the `.rds` | comparing |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 32,285 | 0:16 | 1.5 GiB | 8 s | 1 s | 3 s |
| soupx_tenx_10xmousebrain | 11,831 | 32,285 | 1:41 | 8.3 GiB | 73 s | 8 s | 16 s |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 32,285 | 1:39 | 8.2 GiB | 70 s | 8 s | 16 s |
| soupx_emptydrops_4500871lane | 26,598 | 32,285 | 1:07 | 9.0 GiB | 41 s | 7 s | 15 s |
| soupx_emptydrops_sscovelli | 79,409 | 28,754 | 2:42 | 25.8 GiB | 88 s | 19 s | 51 s |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 28,754 | 2:38 | 25.7 GiB | 86 s | 23 s | 45 s |
| cellbender_fromraw_ghirsutum2 | 59,323 | 99,452 | 3:03 | 19.0 GiB | 123 s | 18 s | 37 s |
| filtered_seurat_tenx_test | 300 | 2,000 | 0:06 | 0.4 GiB | 0 s | 0 s | 1 s |
| filtered_seurat_emptydrops_test | 267 | 2,000 | 0:06 | 0.4 GiB | 0 s | 0 s | 1 s |

Each comparison includes `SCTModel cells_step1`. The Seurat loader's self-check adds 2 to 3 s to loading on the largest objects (measured alone: 2.9 s on the 99,452-feature object, 2.1 s on a 79,409-cell one and on an 11,831-cell one) and up to 0.35 GiB to the peak. On the 99,452-feature object the peak went from 18.7 to 19.0 GiB here, and from 9.8 to 10.1 GiB in K13, where loading is the peak; the 79,409-cell objects' peak here comes later and is unchanged. The first run of this timing, before the self-check and on another node of the same type, took 3:57 on the 99,452-feature object. So M17's array tasks need about 4 minutes and 32 GB for the largest objects (25.8 GiB used), and seconds and 1 GB for the small ones; the comparison is about the size of R's own `readRDS()` of the saved object plus the loaded one.

### M14 tests

- **cargo test --workspace** passes in debug and in release, 281 tests in 41 suites each, none failing or ignored, with no compiler warnings (job 61888985, `data/clustering/m14/cargo_tests/61888985/`) ([`m14_cargo_tests.sbatch`](../../validation/slurm/m14_cargo_tests.sbatch), env.sh's Haswell kernel). M14 changes no Rust code but two README sentences, so the suite is M13's; it runs on the final tree, `validation/tolerances.toml` (which the validators embed) included.
- **The R loader's tests** (`loaders/R/scrustio/tests/scrustio_tests.R`, 72 tests, in the production R env through [`m14_tests.sbatch`](../../validation/slurm/m14_tests.sbatch), job 61878765, with this tree's release binary): every slot against every role, scale and format of schema 1 (18 slots, 8,496 wrong combinations, every one refused, the 18 right ones taken); the comment-line rule both ways; four tampered manifests refused by both builders (the corrected counts declared raw counts; their scale alone changed; the residuals declared raw counts; the counts and log1p layers swapped); the manifest's rules (another schema or version, a missing file, a changed size, a same-size change caught only with the checksum, a shape that does not match its identifiers, a path outside the output); files found by role (renamed files load the identical object); the Seurat object's layout and FindMarkers on it; the SingleCellExperiment's; `scrust sctransform --corrected-counts`', `scrust doubletfinder`'s and `scrust cluster`'s (without `--export`, refused) outputs; `add_scrust()` with DoubletFinder's calls and with the cluster labels; `save_rds()`'s `.rds` reloaded identical, both objects; the version warning. `cells_step1`, 13 tests, on an export made with `--ncells 100 --dump-stages light`:
  - every cell when there are at most `ncells`;
  - the 100 cells scrust drew, by its own record (`stages/S2a_cells_step1.npy`);
  - the same draw under R's default kinds, under L'Ecuyer-CMRG with Box-Muller and Rounding, and with no `.Random.seed`. Each time the caller's kinds and `.Random.seed` (or its absence) come back unchanged, with no warning;
  - the SingleCellExperiment's `cells_step1` the same, and the other cell attributes unchanged;
  - refused: two listed step-1 genes swapped (the self-check, in both builders); another seed in both of the manifest's records (the self-check); the two records disagreeing; and an `n_cells_step1` other than `min(ncells, cells)`.
- **The Python loader's tests** (`loaders/python/scrustio/tests`, 12 tests, pytest 9.2.0 in `.envs/loaders-py` rebuilt from the new pin, the same job): the same checks for its 9 slots (4,248 wrong combinations, every one refused), the tampered manifests, the manifest's rules, files found by role, the AnnData's layout, the `.h5ad` written by the one-line command reread equal, the other commands' outputs, and C's `log1p` against numpy's.

### M14 findings

- **numpy's `log1p` is not the C library's.** numpy 2.5.5's `np.log1p` differs in the last bit from glibc's `log1p` (which R, and so Seurat's SCT `data` layer, use) for 13,457 of the 65,546 values of one tier-T export; `math.log1p` matches every one. Seurat's `data` layer must therefore be read, not recomputed with numpy; the Python loader reads it and checks it with `math.log1p`.
- **`cells_step1` can be drawn again from the export, and checked against it.** The step-1 cells are `sample(colnames(umi), 5000)` right after SCTransform's `set.seed(12345)`. R's own `sample()` with the recorded seed gives the saved objects' cells on all six tier-P objects above 5,000 cells, in both profiles and through both R loaders. vst's next draw, the 2,000 step-1 genes, depends on which cells were drawn, and the export already lists its result in sampling order. So every load checks its own draw without a new file: 2 to 3 s on the largest objects, with sctransform's own functions making vst's draw exactly.
- **The SCT assay's incidental slots** (its key, an empty `misc`, an empty `meta.features`, `nCount_SCT` and `nFeature_SCT` as integer and double) come out as Seurat's only when the assay is built and added with Seurat's own calls in SCTransform's order; the loader does that rather than filling slots by hand, and every saved object's SCT assay, SCTModel and PCA compare identical as whole objects in the production profile.
- **A levels file is tied to its column only by its name.** Schema 1 gives `cluster_levels` files no field naming their column; the loaders read `<column>_levels.txt`'s name for it, the one name they read ([docs/output_manifest.md](../output_manifest.md), "For loaders"). A field could be added within version 1.
- **M17's comparison is affordable:** the Seurat loader's comparison of one export with its saved `.rds` takes at most 3 to 4 minutes (the 99,452-feature object: 3.1 in the final run, 4.0 in the first) and 25.8 GiB (the 79,409-cell objects), most of it loading the export's text files (one R thread parsing up to 73.6 million entries per layer); M17's arrays can take 32 GB per task for the largest objects and far less for the median ones.

### M14 corrections to the plan

Corrected in [docs/PLAN.md](../PLAN.md) on this branch:

- **C3g, the rule:** it said no loader puts corrected counts into a slot named `counts`, which contradicted C3g's own table and its acceptance test: in Seurat they go into the SCT assay's `counts` layer, where SCTransform keeps them and FindMarkers reads them. It now says never into a raw-count slot (Seurat's RNA `counts`, the SingleCellExperiment's `counts`, AnnData's `layers["counts"]`), as `docs/output_manifest.md`'s "For loaders" already did, and that each slot takes one role, scale and format.
- **The user's documentation request (2026-11-18): a new milestone, M14b**, right after M14 (as Step 1's M6b): a page per command (`doubletfinder`, `sctransform`, `pca`, `neighbors`, `cluster`) with its inputs, outputs, formats and how to read each; a page on chaining (which outputs feed which tools, such as `scrust cluster`'s labels into `scrust doubletfinder --clusters` or R decontX's `z`, and where the loaders let R and Python steps join); a schematic (a Mermaid diagram, which GitHub renders). Done when every file, column and command the pages name is checked against real outputs. **The standing rule** joins A7's cycle as step 7: every later milestone or port that adds or changes a tool, a file or a loader updates those pages and the schematic. C6's order notes M14b after M14; whether Step 2's completion waits for it is left to the user (the completion criteria are unchanged).
- **C3g, the table and the notes:** the SingleCellExperiment's SCTModel (`rowData`, `colData` and `metadata()` of `altExp(sce, "SCT")`) and its reduction's attributes; the AnnData's `uns["sct_model"]`, `var["sct_variable"]`, `uns["pca"]`'s contents and `uns["sct_snn"]`'s `connectivities_key`, and an empty `X` for outputs without SCTransform's data; how the Seurat object is built (Seurat's own calls, in SCTransform's order; Seurat's feature names); `cells_step1`; reading Seurat's data layer rather than recomputing it; what the loaders do not rebuild; `add_scrust()`; DoubletFinder's column names.
- **A2:** M14 extended `envs/loaders-py.yml` with pytest and recorded the new pin; the loader is not installed into the env.
- **C4, loader tests:** the exports compared (M13's, both profiles), the criteria per profile, the files check, the 10x readers on every export loaded, the wrong-slot and tampered-manifest tests, the round trips, and `scrust doubletfinder`'s outputs against the saved DoubletFinder objects. **C4, full panel:** M17's comparison script and its time and memory.
- **C5, K11:** identical in the production profile; in the portable profile exactly what the files hold, with the differences from the saved objects within `[clustering.floor.kernel]`; the reloaded objects identical; what an export does not hold reported as not rebuilt.
- **C7, the loaders follow Seurat's internals:** the outcome.

### M14 deviations and notes

- **`cells_step1`: the coordinating session's decision.** The export does not list the step-1 cells (the SCTModel's `cells_step1` above 5,000 cells). The coordinating session was offered three options: export them (a new file and role), have the R loader draw them again from the recorded seed (no export change), or leave them out. It chose to draw them again, on these conditions:
  - everything from the export: the seed and `ncells` from the manifest's parameters and seeds, the cells in `counts/` order, `n_cells = min(ncells, cells)`, and a draw only when that is below the number of cells;
  - a draw independent of the caller's session: R's default kinds set for it, and the caller's kinds and `.Random.seed` restored (or `.Random.seed` removed again);
  - K11 comparing `cells_step1` on its own with every saved object, and M17's comparison too;
  - a self-check only if it is cheap and needs no new file;
  - the R package's README and this report saying so.

  K11's first run already had the draw, but with `n_cells` taken from the summary's `n_cells_step1` and without the rest. After the decision:
  - the R loader took `n_cells` and the seed from the manifest's parameters, checking them against the summary and the seeds;
  - it gained the self-check;
  - it keeps the caller's generator state while loading sctransform's namespace too;
  - `loader_compare.R` gained the `SCTModel cells_step1` row.

  Then everything that runs the R loader was run again with the final loader: K11's R comparisons (`M14_STEPS=r`, 16 jobs, 61880826 to 61880938), K13 (jobs 61881015, 61881036, 61881081, 61881084) and M17's timing (job 61881086). The package tests (job 61878765) and `cargo test` ran on the final tree. Every component is as in the first run, the files checks, round trips and floor grades included, and the new row is identical in all 136 R comparisons. The SingleCellExperiment carries the same `cells_step1` (in `altExp(sce, "SCT")`'s `colData`). The AnnData leaves it out above 5,000 cells, as Python does not reproduce R's generator.
- **No new R reference runs and no new Rust validation runs.** M14 loads M13's exports (binary `scrust-m13-9ce0adb`, whose code main still has) and `scrust doubletfinder` outputs, and reads the saved `.rds` files in batch jobs; K13 runs FindMarkers on the loaded objects. The packages' tests make small outputs with the release binary of this branch's tree.
- **The export's README** of `scrust cluster --export` and `scrust sctransform --corrected-counts` now names the loaders instead of saying they were coming. Only that sentence of `README.txt` changes; M13's exports keep theirs, and no loader reads `README.txt`.
- **The R package's `Maintainer` field** holds a placeholder address (`scrust-maintainers@example.org`, a reserved domain), since `R CMD INSTALL` requires one; replace it before publishing the package.
- **Grading statistics.** `S2d_residual_stats` is graded over `residual_mean` and `residual_variance` together, as the floor measured that array (`cbind(residual_mean, residual_variance)`); seven portable comparisons first graded each column alone (all within the floor) and were run again. Two runs that read `loader_compare.R` while it was being edited (`soupx_tenx_10xmousebrain`'s portable comparison and the M17 timing of `soupx_emptydrops_4500871lane`) stopped with a parse error and were run again. Every result here is from the final scripts. The R loaders' K11 comparisons, K13 and M17's timing come from their second run, with the final loader. The Python loader's comparisons and the 10x readers come from the first run, as their code did not change.
- **The Python env** was re-solved from the extended yaml: pytest 9.2.0 and its dependencies (colorama, exceptiongroup, iniconfig, pluggy, pygments, tomli) added; libaec and numba rebuilt (same versions, new builds), platformdirs 4.14.5 to 4.14.7; every pinned package (Python 3.13.18, scanpy 1.13.0, anndata 0.13.4, numpy 2.5.5, scipy 1.18.1, pandas 3.1.1) and h5py 3.16.0 unchanged. `.envs/loaders-py` was recreated from the new pin; the K11 runs used an env created from the same yaml with the identical solve.
- **The loader comparisons' outputs** (`data/clustering/<object>/loaders/`, `data/doubletfinder/<input>/loaders/`) are listed in M14's section of `validation/manifest.toml`; the exports they load are M13's, not listed, as in M9–M13.

### Reproducing M14

From the repository's top directory, after `source env.sh`, `cargo build --release -p scrust`, M13's validation runs (label `m13-9ce0adb`, their exports in `data/clustering/<object>/rust/m13-9ce0adb/<profile>/run`), and `.envs/loaders-py` created from `envs/loaders-py.linux-64.pin.txt` (as M13's report says):

```bash
sbatch validation/slurm/m14_tests.sbatch
```

(which installs the R loader into `.rlib/`, then runs both packages' tests), then:

```bash
sbatch --mem=8G validation/slurm/m14_loaders.sbatch m13-9ce0adb portable <tier-T objects> df:<tier-T DoubletFinder inputs>
sbatch --mem=128G validation/slurm/m14_loaders.sbatch m13-9ce0adb portable <a tier-P object>
sbatch --mem=8G validation/slurm/m14_loaders.sbatch m13-9ce0adb production <tier-T objects> df:<tier-T DoubletFinder inputs>
sbatch --mem=128G validation/slurm/m14_loaders.sbatch m13-9ce0adb production <a tier-P object>
sbatch --array=0-26 --mem=4G validation/slurm/m14_markers.sbatch m13-9ce0adb portable <the 27 tier-T objects>
sbatch --array=0-6 --mem=32G validation/slurm/m14_markers.sbatch m13-9ce0adb production <the 7 tier-P objects>
sbatch --mem=64G validation/slurm/m14_m17.sbatch m13-9ce0adb production <objects>
sbatch validation/slurm/m14_cargo_tests.sbatch
bash validation/slurm/m14_report.sh m13-9ce0adb > docs/validation/clustering_m14.tsv
bash validation/slurm/m14_report.sh m13-9ce0adb k13 > docs/validation/clustering_m14_k13.tsv
bash validation/slurm/m14_report.sh m13-9ce0adb m17 > docs/validation/clustering_m14_m17.tsv
```

(the markers array likewise for the other profile and tier). With `M14_STEPS=r` before `sbatch`, `m14_loaders.sbatch` runs only the R loaders' comparisons, as their second run did: 8 GB for each profile's 27 tier-T objects at once, and 8 to 48 GB for a tier-P object. One loader by hand, on an export:

```r
library(scrustio)                       # installed in .rlib/
obj <- load_seurat("export")           # or load_sce("export"); save_rds("export", "object.rds")
```

```bash
python -m scrustio h5ad export object.h5ad    # with loaders/python/scrustio on the path or installed
```
