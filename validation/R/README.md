# validation/R

R scripts that make scRust's reference data (docs/PLAN.md A4). They run in the R reference env ([`envs/r-reference/`](../../envs/r-reference/README.md)). Scripts that read `.rds` files run only through sbatch wrappers in [`validation/slurm/`](../slurm/README.md), never on a login node.

[`rcompat_vectors.R`](rcompat_vectors.R) makes the R vectors for the `rcompat` crate (stage S0): R's random-number stream and `.Random.seed` checkpoints, `sample`, `rnorm`/`qnorm`, `round`, `^`, `seq`, `order`, `sum`/`mean`/`var`/`sd` and number formatting (M1); `as.numeric()` of strings (R's own parser, `R_strtod`, including strings it does not round correctly), `fft`, `dnorm` and `approxfun` (M3); over many seeds, sizes and edge cases. It reads no `.rds`, so it may run on a login node:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
  "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
```

Three kinds of script per tool:

- `export_<tool>.R`: one-time conversion of workflow artifacts into plain files: the tool's inputs, and the fields of the workflow's saved outputs that are compared (for DoubletFinder: the doublet set, singlet pANN, pK and nExp). Counts become Matrix Market, tables TSV, and dense arrays and `.Random.seed` checkpoints `.npy` (written by a small base-R helper).
- `reference_<tool>.R`: a direct run of the original tool through a thin script with the workflow's parameters and seeds, writing complete outputs as plain files. Each run must reproduce the corresponding saved workflow output (or test reference output), which shows the script calls the tool exactly as the workflow did.
- `instrumented_<tool>.R`: the same run, also dumping the tool's internal stages (for DoubletFinder: doublet indices, SCTransform, PCA, neighbours and RNG checkpoints), so a difference can be traced to the stage that caused it.

## DoubletFinder (milestone M2)

Datasets (tier T: the five workflow test objects; tier P: five production objects) and their pipeline paths are in [`../doubletfinder_datasets.tsv`](../doubletfinder_datasets.tsv). Results: [`docs/validation/doubletfinder_m2.md`](../../docs/validation/doubletfinder_m2.md).

- [`export_doubletfinder.R`](export_doubletfinder.R): from each input object, the counts DoubletFinder 2.0.6 reads with `sct = TRUE` (the RNA assay's `counts` layer) as 10x Matrix Market, entries exactly as stored, explicit zeros included, cells and features in object order; `seurat_clusters` and its level order. From each saved workflow output: the kept barcodes with their pANN, the doublet set (input minus output barcodes), and pK and nExp parsed from the column names. It re-reads the matrix with the reference runs' reader and checks `identical(..., num.eq = FALSE)` against the object's own matrix.
- [`reference_doubletfinder.R`](reference_doubletfinder.R): the direct run, with the workflow's calls and seeds (`set.seed(12345)`, `paramSweep(PCs = 1:10, sct = TRUE)`, `summarizeSweep(GT = FALSE)`, `find.pK`, the pK via `as.numeric(as.character(...))`, the homotypic nExp, `set.seed(12345)`, `doubletFinder(pN = 0.25, ...)`). Writes plan B3's `cells.tsv` (every input cell), `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json` (a deterministic `result` section and a `run` section with timings, peak RSS and the environment fingerprint), plus `find_pK.pdf`, the plot `find.pK` draws. Its input is an export directory (the Seurat object is rebuilt with `CreateSeuratObject()` from the identical counts) or, for diagnosis, the workflow's `.rds`.
- [`instrumented_doubletfinder.R`](instrumented_doubletfinder.R): the same run, also dumping stages S1–S7 to `stages/` (`.npy`, with `index.tsv`, `rng.npy`/`rng.tsv` and `stages.json`). Mechanism: `trace()` hooks in the installed functions that only read values and write files ([`lib/doubletfinder_hooks.R`](lib/doubletfinder_hooks.R) documents the hooks, the array layout and the index conventions). Its final outputs are byte-identical to the direct run's, checked for every object that has both.
- [`compare_doubletfinder.R`](compare_doubletfinder.R): the M2 checks on plain files: a run against the saved workflow output (doublet set, singlet pANN bitwise, column names), and two runs against each other (final files, the `result` section, every stage array block by block, RNG checkpoints), driven by [`../slurm/compare_doubletfinder.sbatch`](../slurm/compare_doubletfinder.sbatch).
- [`manifest_doubletfinder.R`](manifest_doubletfinder.R): regenerates the DoubletFinder part of `validation/manifest.toml` from the files on disk.
- [`lib/io.R`](lib/io.R): the `.npy` writer and reader, exact text output (`%.17g`), JSON, the 10x Matrix Market writer and exact reader, peak RSS and the environment fingerprint (which calls `scrust fingerprint` for the CPU, OpenBLAS core and threads, and glibc). [`lib/doubletfinder.R`](lib/doubletfinder.R): the workflow's call sequence, shared by the direct and instrumented runs.

Exactness: doubles are written as text with 17 significant digits (glibc's correctly rounded printf; a correctly rounded reader, such as Rust's, gets the exact double back). R's own decimal parser is not correctly rounded, so checks compare the written text, never R-parsed non-integers. The workflow's counts are integer-valued (SoupX `roundToInt = TRUE`, decontX `round()`, CellBender), which R parses exactly.

Rules:

- R runs with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB` and `PYTHONNOUSERSITE=1`; sourcing `env.sh`, as every wrapper does, sets this up.
- OpenBLAS's kernel and thread count come from `env.sh` (portable profile by default, see the profiles in [`validation/slurm/README.md`](../slurm/README.md)).
- Seeds and parameters are the workflow's, as written in its scripts (for DoubletFinder, `workflow/scripts/doubletfinder.R` with seed 12345).
- Outputs go to `data/` (large, gitignored) or `tests/fixtures/` (small, committed), and every file is listed in `validation/manifest.toml` with its sha256.

## DoubletFinder core (milestone M3)

The Rust core runs on the instrumented runs' own embeddings and is compared with them by `scrust validate doubletfinder` ([`../slurm/m3_doubletfinder.sbatch`](../slurm/m3_doubletfinder.sbatch)); no new R run was needed beyond the `rcompat` vectors. Results: [`docs/validation/doubletfinder_m3.md`](../../docs/validation/doubletfinder_m3.md).

- [`report_doubletfinder_m3.R`](report_doubletfinder_m3.R): the report's tables from the M3 outputs (the `validate.json` of every Rust run, the direct R runs, the saved workflow outputs and the benchmark's `/usr/bin/time -v`). Reads only plain files.

## The whole tool (milestone M6)

`scrust doubletfinder` runs from the counts, the cluster labels and the seed; it is compared with every direct and instrumented reference run of M2 and with the saved workflow outputs by `scrust validate doubletfinder` ([`../slurm/m6_doubletfinder.sbatch`](../slurm/m6_doubletfinder.sbatch)), and benchmarked against M2's direct production runs ([`../slurm/m6_benchmark.sbatch`](../slurm/m6_benchmark.sbatch)). No new R run was needed. Results: [`docs/validation/doubletfinder.md`](../../docs/validation/doubletfinder.md).

- [`report_doubletfinder_m6.R`](report_doubletfinder_m6.R) `<label> <bench label>`: the report's tables from the M6 outputs (each Rust run's `vs_<R run>.json` and `vs_saved.json`, its `summary.json` and `time_v.txt`; the benchmark's runs and the direct R runs' `time_v.txt` and `summary.json`). Reads only plain files.

## The singlets' counts (milestone M6b)

`scrust doubletfinder` also writes `singlets/`: the counts without the doublets, as a new 10x directory (docs/PLAN.md A3 and B3). The workflow's equivalent is the RNA `counts` layer of its saved DoubletFinder output, since `doubletfinder.R` keeps the singlets with `subset()`.

- [`export_doubletfinder_singlets.R`](export_doubletfinder_singlets.R): exports that layer from each saved output into `data/doubletfinder/<object>/saved/singlets/`, next to the M2 export: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` (`lib/io.R`'s writer, with the input export's feature lines) and `export.json`. It checks three things:
  - the `.rds` is the one M2 exported (sha256);
  - its cells are the kept cells, in order;
  - the directory reads back `identical()`.

  It also records whether the layer equals the input export's counts restricted to the kept cells, slot for slot. Run it through [`../slurm/export_doubletfinder_singlets.sbatch`](../slurm/export_doubletfinder_singlets.sbatch).
- `scrust validate doubletfinder --saved` compares `singlets/` with that export: the barcodes in order, the feature lines, and every stored entry and value of the matrix. [`../slurm/m6_doubletfinder.sbatch`](../slurm/m6_doubletfinder.sbatch) runs this check with the others.
- [`read10x_singlets.R`](read10x_singlets.R) `<label> <profile> <object>...` reads the Rust runs' `singlets/` back with R's readers, as a downstream R step would: `lib/io.R`'s exact reader and `Seurat::Read10X()`. It compares each with the export, slot for slot, and records whether `Read10X` keeps every stored entry. It reads only plain files.

## The full panel (milestone M7)

All 150 production DoubletFinder objects (30 samples, 5 arms) against their saved workflow outputs. [`../doubletfinder_panel.tsv`](../doubletfinder_panel.tsv) lists them with their pipeline paths and the workflow job that wrote each saved output (job, Snakemake run, node, CPU type, time, peak RSS). The five tier-P objects are among them. No new R script runs DoubletFinder: the other 145 objects get the same one-time exports as the P objects, `export_doubletfinder.R` and `export_doubletfinder_singlets.R`, through [`../slurm/m7_export.sbatch`](../slurm/m7_export.sbatch), and `scrust validate doubletfinder --saved` compares the Rust runs ([`../slurm/m7_panel.sbatch`](../slurm/m7_panel.sbatch)). Results: [`docs/validation/doubletfinder_m7.md`](../../docs/validation/doubletfinder_m7.md).

- [`report_doubletfinder_m7.R`](report_doubletfinder_m7.R) `<label> [<objects.tsv>]`: the report's tables from the panel table, the exports' `export.json` and `saved.json`, and each Rust run's `vs_saved.json`, `summary.json`, `fingerprint.json`, `node_features.txt` and `time_v.txt`; with a second argument also one row per object as TSV (`docs/validation/doubletfinder_m7.tsv`). An object passes when all nine checks are made and identical. Reads only plain files.
- [`manifest_doubletfinder.R`](manifest_doubletfinder.R) lists the panel's exports too, as tier F.

## Step 2, the clustering crate: reference data (milestone M8)

The workflow's clustering chain (docs/PLAN.md C1): `SCTransform(seed.use = 12345)`, `RunPCA(seed.use = 12345)` (50 PCs), `RunUMAP(dims = 1:30, seed.use = 12345)`, `FindNeighbors(dims = 1:30)`, `FindClusters(random.seed = 12345)` and the workflow's `add_silhouette_to_metadata()` (bluster's `approxSilhouette` and `neighborPurity`). Objects (tier T: the 27 test objects; tier P: 7 production objects), their pipeline paths and dump levels are in [`../clustering_datasets.tsv`](../clustering_datasets.tsv); the full panel's 840 objects in [`../clustering_panel.tsv`](../clustering_panel.tsv). Results: [`docs/validation/clustering_m8.md`](../../docs/validation/clustering_m8.md). Everything goes to `data/clustering/` (gitignored); the committed subset is in `tests/fixtures/clustering/`.

- [`lib/clustering.R`](lib/clustering.R): the chain, copied from the workflow's scripts with their source cited (`set_workflow_seed`, `require_min_cells_for_pca` and `add_silhouette_to_metadata` from `silhouette_utils.R`, the calls and `write_cluster_metadata` from `tenx2seuratrds.R` and the eight other scripts); the input object (`CreateSeuratObject` on the exported counts, as a v5 `Assay5`, or a v3 `Assay` for the emptyDrops objects); and one plain-file layout (`cl_write_object`) for the saved objects and the runs, so they compare file by file: `metadata.tsv` (the chain's columns, doubles exact), the levels, `sct/` (the corrected counts in full; scale.data and the log1p data by level, otherwise as digests; the SCTModel's slots), `pca/`, `graphs/`.
- [`export_clustering.R`](export_clustering.R): from a saved object, its RNA counts layer as saved (or links to an earlier M2 or M6b export of the same counts, checked equal slot for slot) and what the chain left in it, in that layout, plus a copy of its marker table. Through [`../slurm/export_clustering.sbatch`](../slurm/export_clustering.sbatch).
- [`reference_clustering.R`](reference_clustering.R): the direct run on an export; `--umap no` leaves RunUMAP out, `--assay v5|v3` forces the RNA assay's class. Writes the layout above, `cluster_ids.txt`, `nclusters.txt` and `summary.json` (a deterministic `result` section, with the optimizer's printed numbers and a digest of `.Random.seed` after every step, and a `run` section with timings and peak RSS after every step).
- [`instrumented_clustering.R`](instrumented_clustering.R): the same run, plus `stages/` through passive `trace()` hooks ([`lib/clustering_hooks.R`](lib/clustering_hooks.R), which reuses M2's SCTransform and PCA hooks unchanged): SCTransform at seed 12345 (S2a–S2d, and vst's result in `scrust sctransform`'s layout), the PCA (S3), the RNG checkpoints, the float items Annoy holds (K3), the index file from a second identical build (K4), every cell's neighbours and distances (K5), the optimizer's and GroupSingletons' results (K7), approxSilhouette's (K8) and neighborPurity's (K9) intermediates, and `sctransform::correct()`'s median depth (K10).
- [`compare_clustering.R`](compare_clustering.R): a run against its saved object, or two runs against each other (every file, every stage array, RNG checkpoint and scalar; the floor statistics). Through [`../slurm/compare_clustering.sbatch`](../slurm/compare_clustering.sbatch).
- [`markers_clustering.R`](markers_clustering.R): the workflow's `find_markers.R` and `combine_markers.R` steps, copied with their source cited, on a saved `.rds`; the combined table must equal the saved `<prefix>_markergenes.csv` byte for byte, or becomes the reference. Through [`../slurm/markers_clustering.sbatch`](../slurm/markers_clustering.sbatch).
- [`optimizer_clustering.R`](optimizer_clustering.R): the instrumented modularity optimizer (Seurat's `ModularityOptimizer.cpp` standalone, built by [`../cpp/build.sh`](../cpp/build.sh)) on a saved SNN graph, against R's optimizer result, the saved labels and the workflow job's log.
- [`k0_vectors.R`](k0_vectors.R): R's part of the K0 vectors, `set.seed(1); sample(x, 1)` as GroupSingletons draws. Reads no `.rds`.
- [`floor_clustering.R`](floor_clustering.R): the R-vs-R floor against the saved objects. Every portable-profile run against its saved object, output by output: identical or not, values differing, largest absolute and relative difference, and the PCA's own largest difference (signs aligned) beside each. The saved objects were made with OpenBLAS's SapphireRapids kernel, so only production-profile runs reproduce them exactly. It writes one row per run and output, checks that an object's portable runs all differ from the saved object in the same way, and prints the report's tables. Plain files only; through [`../slurm/floor_clustering.sbatch`](../slurm/floor_clustering.sbatch).
- [`floor_scale_data.R`](floor_scale_data.R): for the objects whose runs kept scale.data only as a digest (dump level light), the portable profile's scale.data, corrected counts and feature attributes against the saved object's, element by element. It reruns the chain's SCTransform call and reads the saved `.rds`. The new scale.data's digest must equal the reference run's. Through [`../slurm/floor_scale_data.sbatch`](../slurm/floor_scale_data.sbatch).
- [`clustering_panel.R`](clustering_panel.R): writes `validation/clustering_panel.tsv` from the Snakemake logs, sacct, scontrol, the saved objects' modification times and the jobs' logs. Reads no `.rds`.
- [`report_clustering_m8.R`](report_clustering_m8.R) and [`manifest_clustering.R`](manifest_clustering.R): the report's tables and the manifest's clustering section, from plain files.

## The neighbour stage (milestone M9)

The R side of M9's checks of `scrust neighbors` (docs/PLAN.md C4, C5 K3–K6, C7; results in [`docs/validation/clustering_m9.md`](../../docs/validation/clustering_m9.md)). The Euclidean references are M8's instrumented runs; these add the rest. No `.rds` is read and no BLAS is used.

- [`neighbors_references.R`](neighbors_references.R): FindNeighbors in the pinned env, writing each reference in the layout of an instrumented run's `stages/` and `graphs/` (the index file from Seurat's own `AnnoyBuildIndex`, saved; the items; FindNeighbors' neighbour ids and distances; both graphs, also as their `dgCMatrix` slots), with checks in `stages.json` (FindNeighbors' search equals a second, unsaved build's and the saved index's). `--mode real` takes an instrumented run's embeddings (columns 1:30) for the metrics M8 did not run (cosine, manhattan); `--mode edge` builds the constructed edge cases (inputs, references for all three metrics, `cases.tsv`, and `expected.tsv` with every reference file's sha256 for cargo test); `--mode bench` times the workflow's FindNeighbors call alone, for the benchmark. Through [`../slurm/neighbors_references.sbatch`](../slurm/neighbors_references.sbatch) and [`../slurm/m9_benchmark.sbatch`](../slurm/m9_benchmark.sbatch).
- [`report_clustering_m9.R`](report_clustering_m9.R) and [`manifest_neighbors.R`](manifest_neighbors.R): the report's tables (and `docs/validation/clustering_m9.tsv`) and the manifest's neighbour section, from plain files.

## FindClusters (milestone M10)

The R side of M10's checks of `scrust cluster` (docs/PLAN.md C4, C5 K7, C7; results in [`docs/validation/clustering_m10.md`](../../docs/validation/clustering_m10.md)). The workflow-setting references are M8's instrumented runs; these add the other settings, the edge cases and R's benchmark. No `.rds` is read and no BLAS is used.

- [`cluster_references.R`](cluster_references.R): FindClusters in the pinned env on an SNN graph held as `SCT_snn` in a minimal Seurat object, with passive `trace()` hooks on `RunModularityClustering` and `GroupSingletons` (each call's optimizer ids, final ids, singletons, and `.Random.seed` at entry and exit), each reference written in the layout of an instrumented run with the instrumented optimizer's per-start results beside it (`optimizer/<res.r>/`). The saved graphs hold their values as `%.17g` text, which R's own parser does not always read back exactly, so each distinct value is corrected to the double whose text it is, and the matrix is checked to write back the same text. `--mode settings` runs the workflow's call with one argument changed (the parameter grid's clustering settings, plus multi-resolution and singleton-label variants) on an object's saved graph, and checks that the workflow setting reproduces the saved object and M8's printed numbers; `--mode edge` builds the constructed edge cases; `--mode bench` times the workflow's FindClusters call alone. Through [`../slurm/cluster_references.sbatch`](../slurm/cluster_references.sbatch) and [`../slurm/m10_benchmark.sbatch`](../slurm/m10_benchmark.sbatch).
- [`report_clustering_m10.R`](report_clustering_m10.R) and [`manifest_cluster.R`](manifest_cluster.R): the report's tables (and `docs/validation/clustering_m10.tsv`) and the manifest's FindClusters section (with a `SHA256SUMS` in each reference directory), from plain files.

The C++ programs, in [`../cpp/`](../cpp/): `ModularityOptimizer.patch` (prints each random start's iterations, modularity as `%.17g` and partition; nothing else changes), `javarandom_vectors.cpp` (the optimizer's `JavaRandom`, from Seurat's unmodified source) and `kiss64_vectors.cpp` (RcppAnnoy's `Kiss64Random`, from its installed header). `build.sh` compiles them with the R env's compiler and `R CMD config` flags into `data/clustering/bin/` (gitignored), recording compiler, flags and checksums in `build.json`; nothing is written to the R env.
