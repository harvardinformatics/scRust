# validation/slurm

sbatch wrappers for everything that runs on compute nodes (docs/PLAN.md A2 and A4):

- the export scripts and the direct and instrumented reference runs in `validation/R/`; anything that reads `.rds` files runs this way, never on a login node;
- tier-P and full-panel validation (`scrust validate ...`, from M3 on);
- memory and time benchmarks.

Start each wrapper from [`template.sbatch`](template.sbatch), which:

- sources the node guard, `/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh`. It reads one byte of a file stored on holylfs05 OST 51 and exits with code 75 if that hasn't returned within 30 s, instead of letting the job hang unkillably on a node that lost that storage target;
- excludes the nodes known to have lost it (`--exclude=holy8a32608,holy8a24102`) and runs under the `informatics` account (`-A informatics`);
- sources [`env.sh`](../../env.sh): the project env, OpenBLAS pinned to one kernel and one thread for R and Rust alike, and no user-level R, Perl or Python libraries;
- prints an environment fingerprint (node, CPU model, glibc, OpenBLAS kernel and threads, rustc, commit) to the log.

## Reference profiles (A5b)

| Profile | `OPENBLAS_CORETYPE` | Nodes | Use |
|---|---|---|---|
| portable | `Haswell` (env.sh's default) | any x86-64 node with AVX2 | development, CI, fixtures that can be regenerated on any cluster |
| production | `auto` (OpenBLAS detects `SapphireRapids`) | Sapphire Rapids: `--partition=sapphire`, or `--constraint=sapphirerapids` elsewhere | comparisons with the saved production outputs; all 150 production DoubletFinder jobs ran on Sapphire Rapids |

Use one CPU per job for anything that calls BLAS (`--cpus-per-task=1`), as the production jobs did.

The production profile must leave `OPENBLAS_CORETYPE` unset: env.sh unsets it when it is `auto`. Passing a value to OpenBLAS does not reproduce detection: on Sapphire Rapids, `OPENBLAS_CORETYPE=SapphireRapids` and a literal `OPENBLAS_CORETYPE=auto` both select `Cooperlake` (measured with `scrust fingerprint`, M2).

[`scrust fingerprint`](../../crates/scrust/src/fingerprint.rs) prints the fingerprint as JSON: node, SLURM job, CPU model and flags, the OpenBLAS core, configuration and thread count actually in use, the OpenBLAS library loaded, and glibc. Build it once with `cargo build --release -p scrust`; the R reference scripts record its output in their `summary.json`.

## DoubletFinder (M2)

Submit from the repository's top directory. Objects and their pipeline paths are in [`../doubletfinder_datasets.tsv`](../doubletfinder_datasets.tsv) (tiers T and P) and, for the full panel's other objects, [`../doubletfinder_panel.tsv`](../doubletfinder_panel.tsv) (M7); the export wrappers look in both. Outputs go to `data/doubletfinder/<object>/` and are never overwritten.

- [`export_doubletfinder.sbatch`](export_doubletfinder.sbatch) `<object>...`: exports inputs and saved-output fields (reads `.rds`).
- [`export_doubletfinder_singlets.sbatch`](export_doubletfinder_singlets.sbatch) `<object>...` (M6b): exports the saved outputs' RNA counts, the singlets' counts, to `saved/singlets/`, next to an existing M2 export (reads `.rds`). It needs no particular CPU type, so `-p shared,sapphire` starts it sooner. 32 GB covers every object.
- [`reference_doubletfinder.sbatch`](reference_doubletfinder.sbatch) `<portable|production> <direct|instrumented> <object>...`: reference runs into `runs/<mode>-<profile>-<cpu arch>/`, under `/usr/bin/time -v`. The job stops if OpenBLAS does not report the profile's core (`Haswell`, or `SapphireRapids`). Other CPU types for the portable profile: `--partition=shared` (Cascade Lake), `--partition=serial_requeue --constraint=genoa` (AMD Genoa).
- [`compare_doubletfinder.sbatch`](compare_doubletfinder.sbatch) `<object>...`: the M2 checks (cross-check against the saved workflow output, instrumented against direct, CPU types, portable against production); small objects can run it with `bash`.

Memory used for the reference runs (`--mem`): tier T 4 GB; soupx_tenx_SRR18256555 16 GB; soupx_tenx_10xmousebrain 40–48 GB; soupx_emptydrops_4500871lane 64–72 GB; soupx_emptydrops_sscovelli 200 GB; cellbender_fromraw_ghirsutum2 250 GB (the measured peaks are in docs/validation/doubletfinder_m2.md).

## DoubletFinder core (M3)

After `cargo build --release -p scrust`, from the repository's top directory. The core uses no BLAS, so these jobs take several CPUs; results do not depend on the thread count.

- [`m3_doubletfinder.sbatch`](m3_doubletfinder.sbatch) `<object>...`: for every instrumented reference run of each object, the Rust core on that run's embeddings (`scrust doubletfinder --from-stages`, with stage dumps) into `data/doubletfinder/<object>/rust/m3/<run>/`, then `scrust validate doubletfinder` against the same run (`validate.json`, `validate.md`). 16 CPUs and 16 GB cover every object (the ten objects take 30 s in one job).
- [`m3_benchmark.sbatch`](m3_benchmark.sbatch) `[<object>] [<run>]`: wall time and peak RSS under `/usr/bin/time -v` of `scrust knn` (the final run's neighbour search) and of the whole core, at 1 thread and at `SLURM_CPUS_PER_TASK` threads (default soupx_emptydrops_sscovelli, 32 CPUs), into `data/doubletfinder/<object>/rust/m3-benchmark/<job id>/`.

Logs go to `<job name>-<job id>.sbatch.log` in the submit directory, which `.gitignore` excludes.

## PCA (M5)

After `cargo build --release -p scrust`, from the repository's top directory, with a pinned copy of the binary (`cp target/release/scrust data/doubletfinder/logs/bin/scrust-<label>`, passed as `SCRUST_BIN`), so that a rebuild cannot change what a queued job runs:

- [`m5_pca.sbatch`](m5_pca.sbatch) `<label> <object>...`: for every instrumented reference run of each object, `scrust validate pca` (stage S3 in all 7 contexts, teacher-forced on R's PCA input and after Rust's own SCTransform) under that run's OpenBLAS profile, the core checked first, into `data/doubletfinder/<object>/rust/<label>/<run>/` (`validate.json`, `validate.md`, `sources.txt`, `fingerprint.json`, `time_v.txt`). OpenBLAS runs one thread; irlba's two large products are split across the job's CPUs only where that keeps R's bits, so the CPU count changes the time only. 8 CPUs and 48 GB cover every object.
- [`m5_benchmark.sbatch`](m5_benchmark.sbatch) `rust|r <object> [<sct dir>]`: `scrust pca --npcs 10` at 1 thread and at `SLURM_CPUS_PER_TASK` threads (`--cpus-per-task=32 --mem=32G`), or R's `RunPCA()` through [`runpca_reference.R`](../R/runpca_reference.R) (`--cpus-per-task=1 --mem=48G`), under `/usr/bin/time -v`, on the same SCTransform output (default: R's own, from the M4 benchmark of the object's final run), in the portable profile; the outputs are compared byte for byte (`data/doubletfinder/<object>/rust/m5-benchmark/final/`).

## The whole tool (M6)

After `cargo build --release -p scrust`, from the repository's top directory, with a pinned copy of the binary passed as `SCRUST_BIN` (as for M5):

- [`m6_doubletfinder.sbatch`](m6_doubletfinder.sbatch) `<label> <object>...`: for each object and each OpenBLAS profile its reference runs use (portable, production), `scrust doubletfinder` from the export's counts and cluster labels, with stage dumps at the level of the profile's instrumented R run, then `scrust validate doubletfinder` against every reference run of that profile (instrumented: the final files, the `result` section, every stage array S1–S7, every `.Random.seed` checkpoint and scalar; direct: the final files and the `result` section) and against the saved workflow output, into `data/doubletfinder/<object>/rust/<label>/<profile>/` (`run/`, `vs_<R run>.{json,md}`, `vs_saved.{json,md}`, `fingerprint.json`, `time_v.txt`). The saved-output comparison includes `run/singlets/` against the export in `saved/singlets/` (M6b) wherever that export exists. The OpenBLAS core is checked before each run. 8 CPUs and 48 GB cover every object; the ten objects take about 20 minutes in one job.
- [`m6_benchmark.sbatch`](m6_benchmark.sbatch) `<label> <object> <threads>...`: `scrust doubletfinder` without dumps at each thread count under `/usr/bin/time -v`, in the production profile (`--cpus-per-task=32 --mem=32G`), compared with the object's direct production R run (whose `time_v.txt` is R's wall time and peak RSS on the same input and CPU type) and the runs with each other, into `data/doubletfinder/<object>/rust/<label>/t<threads>/`.

## The full panel (M7)

All 150 production DoubletFinder objects of [`../doubletfinder_panel.tsv`](../doubletfinder_panel.tsv), one array task per object (task i: the table's data row i). The table records, for each object, its pipeline paths and the workflow job that wrote its saved output: the job, its Snakemake run, node, CPU type (the node's SLURM feature), time and peak RSS (sacct). From the repository's top directory:

- [`m7_export.sbatch`](m7_export.sbatch): the one-time exports of an object, run with `bash` inside the task: `export_doubletfinder.sbatch` (inputs and compared fields), then `export_doubletfinder_singlets.sbatch` (the saved output's RNA counts), each unless already done. Exports need no particular CPU type (`--partition=shared,sapphire`); 1 CPU and 24 GB cover every object (peak 8.6 GiB).
- [`m7_panel.sbatch`](m7_panel.sbatch) `<label>`: `scrust doubletfinder` without stage dumps on the object's export, then `scrust validate doubletfinder --saved` against the saved workflow output (the doublet set, the kept cells' pANN, the column names, so pK and nExp, and `singlets/` against the saved RNA counts), into `data/doubletfinder/<object>/rust/<label>/production/` (`run/`, `vs_saved.{json,md}`, `fingerprint.json`, `node_features.txt`, `time_v.txt`, `stderr.txt`). The production profile: the task runs on `sapphire`, stops unless the node has the CPU feature of the object's workflow job and OpenBLAS, left to detect the CPU, reports that CPU's core (`SapphireRapids`); it fails unless all nine checks are made and identical. 8 CPUs and 32 GB cover every object (peak 6.6 GiB, at most 2.4 minutes). With a pinned copy of the binary, as for M6.
- [`doubletfinder_panel.sh`](doubletfinder_panel.sh) `<snakemake log>...`: writes the panel table to stdout from the pipeline's Snakemake logs (the DoubletFinder jobs' SLURM ids), sacct (node, state, elapsed time, peak RSS of step `.0`) and scontrol (the node's CPU-type feature). It reads no `.rds`, so a login node is enough; it stops if an object had more than one job or a job did not complete.

```bash
c=$(git rev-parse --short HEAD)
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m7-$c
mkdir -p data/doubletfinder/logs/m7
E=$(sbatch --parsable --array=1-150%50 --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_export.sbatch)
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m7-$c sbatch --export=ALL --array=1-150%50 --dependency=aftercorr:$E \
    --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_panel.sbatch m7-$c
```

`--dependency=aftercorr` starts each object's Rust run once its own export task has succeeded.

## The clustering chain's reference data (M8)

Submit from the repository's top directory. Objects, dump levels, pipeline paths and the earlier exports reused are in [`../clustering_datasets.tsv`](../clustering_datasets.tsv) (tiers T and P); every production object's workflow job, node and CPU type in [`../clustering_panel.tsv`](../clustering_panel.tsv), written by `Rscript validation/R/clustering_panel.R validation/clustering_panel.tsv <the two Snakemake logs of rerun #2>` (no `.rds`; a login node is enough). Outputs go to `data/clustering/<object>/` and are never overwritten. All 840 production jobs, and the 27 test-reference jobs, ran on Sapphire Rapids.

- [`export_clustering.sbatch`](export_clustering.sbatch) `<object>...`: the one-time exports (reads `.rds`). No particular CPU type. Peak RSS at most 11.9 GiB (filtered_seurat_emptydrops_sscovelli) and at most 8 minutes per object, longest where the counts are written rather than linked; 4 GB covers tier T, 24 GB every P object.
- [`reference_clustering.sbatch`](reference_clustering.sbatch) `<portable|production> <direct|instrumented> <object>...`: reference runs into `runs/<mode>-<profile>-<cpu arch>[-noumap][-v5]/`, under `/usr/bin/time -v`. The production profile stops unless the node has the CPU feature of the object's workflow job and OpenBLAS, left to detect the CPU, reports its core (`SapphireRapids`). Variants through the environment: `CL_UMAP=no` (no RunUMAP), `CL_ASSAY=v5` (the emptyDrops objects through the v5 path), `CL_LEVEL`, `CL_PURITY_LISTS`. Measured peaks (`/usr/bin/time -v`): at most 0.93 GiB on tier T, 2.5 to 23.4 GiB on the small and medium P objects, 55.6 to 62.2 GiB on the two 79,409-cell objects (the v5 path the most) and 104.1 GiB on cellbender_fromraw_ghirsutum2; ask for 6 GB, 32 GB, 96 GB and 160 GB. Several configurations can share one job: `sbatch --wrap "CL_UMAP=no bash validation/slurm/reference_clustering.sbatch production direct <objects>; ..."`.
- [`markers_clustering.sbatch`](markers_clustering.sbatch) `<object>...`: the marker references (reads `.rds`; no BLAS, any CPU type).
- [`optimizer_clustering.sbatch`](optimizer_clustering.sbatch) `<object>...`: the instrumented modularity optimizer on the saved SNN graph (no `.rds`; tier T runs with `bash` in seconds). Build it first: `bash validation/cpp/build.sh`.
- [`compare_clustering.sbatch`](compare_clustering.sbatch) `<object>...`: every run against its saved object, and the pairs of runs M8 compares (instrumented against direct, Haswell against SapphireRapids kernels, Cascade Lake against Sapphire Rapids, without against with UMAP, v5 against v3), into `data/clustering/compare/<object>/`. Plain files only; tier T can run it with `bash`.
- [`floor_scale_data.sbatch`](floor_scale_data.sbatch) `<object>...`: the portable profile's scale.data against the saved object's, for the objects whose runs kept it only as a digest (reads `.rds`; portable profile, Sapphire Rapids), into `floor/scale_data_portable.json`. Measured peaks: 55.6 and 58.6 GiB on the 79,409-cell objects (6 minutes each), 104.2 GiB on cellbender_fromraw_ghirsutum2 (11 minutes); ask for 96 GB and 160 GB. M8 also ran it on soupx_tenx_SRR18256555 (2.5 GiB, under a minute) as a cross-check against its full scale.data.
- [`floor_clustering.sbatch`](floor_clustering.sbatch) `<out.tsv> <tables.md>`: every portable run against its saved object, output by output (the R-vs-R floor), and the report's tables. Run it after `floor_scale_data.sbatch`. Plain files only; no particular CPU type. On all 34 objects: peak 1.9 GiB, 2.5 minutes; ask for 16 GB.
- [`fixtures_clustering.sh`](fixtures_clustering.sh): copies the committed tier-T subset into `tests/fixtures/clustering/` (plain files; a login node is enough).

During M8 the `sapphire` queue was held up by another run under the same account, so most Sapphire Rapids jobs ran on the `test` partition (Sapphire Rapids nodes; at most 5 jobs and 12 hours per job), and CPU-independent jobs on `shared` and `serial_requeue`.

## The neighbour stage (M9)

`scrust neighbors` (FindNeighbors: Annoy and the SNN graph) against R, teacher-forced on R's own embeddings; results in [`docs/validation/clustering_m9.md`](../../docs/validation/clustering_m9.md). Nothing here uses BLAS or reads `.rds`, so any CPU type will do. Submit from the repository's top directory after `cargo build --release -p scrust`; outputs are never overwritten.

- [`neighbors_references.sbatch`](neighbors_references.sbatch) `edge` | `real <object>/<run>...`: the R references M8 did not make, by [`../R/neighbors_references.R`](../R/neighbors_references.R) — the constructed edge cases into `data/neighbors/edge/`, or FindNeighbors with `annoy.metric` cosine and manhattan (`METRICS` to change) on instrumented runs' embeddings into `data/clustering/<object>/neighbors/<run>/<metric>/`. M9 ran it on the 27 tier-T objects (portable and production instrumented runs), soupx_tenx_SRR18256555 (both) and soupx_emptydrops_sscovelli (portable): 12 GB is plenty for tier T and the small object (job peak 0.6 GiB, 1.3 minutes for the 56 runs' two metrics), 24 GB for the 79,409-cell one (peak 2.9 GiB, 4 minutes); the edge cases take 42 s and 0.5 GiB.
- [`m9_neighbors.sbatch`](m9_neighbors.sbatch) `<object>...`: `scrust neighbors --annoy yes` on the embeddings of every instrumented run (Euclidean) and of every reference above (cosine, manhattan), each followed by `scrust validate neighbors` against its R run, into `data/clustering/<object>/rust/m9/<run>/<metric>/`. `SCRUST_BIN` pins a copy of the binary (M9 used `data/clustering/bin/scrust-m9-db5e474`); `EUCLIDEAN_ONLY=1` skips the other metrics. All 212 comparisons of the 34 objects took 3.4 minutes on 16 CPUs (2.1 for the 98 Euclidean ones, 1.2 for the 114 others).
- [`m9_edge.sbatch`](m9_edge.sbatch): the same on every edge case, into `data/neighbors/edge-rust/` with `status.tsv`; seconds, so `bash` on a login node is enough.
- [`m9_benchmark.sbatch`](m9_benchmark.sbatch) `<object>...`: R's FindNeighbors alone against `scrust neighbors` at 1 and 32 threads on one node, and the outputs at 1 and 32 threads compared byte for byte, into `data/clustering/<object>/rust/m9-bench/`.
- [`fixtures_neighbors.sh`](fixtures_neighbors.sh): copies the committed edge-case fixture into `tests/fixtures/neighbors/edge/` (plain files; a login node is enough).

## FindClusters (M10)

`scrust cluster` (FindClusters: the modularity optimizer and GroupSingletons) against R, teacher-forced on R's own SNN graphs; results in [`docs/validation/clustering_m10.md`](../../docs/validation/clustering_m10.md). Nothing here uses BLAS or reads `.rds`, so any CPU type will do. Submit from the repository's top directory after `cargo build --release -p scrust` and `bash validation/cpp/build.sh` (the instrumented optimizer); outputs are never overwritten.

- [`cluster_references.sbatch`](cluster_references.sbatch) `settings <object>...` | `edge` | `bench <object>...`: the R references M8 did not make, by [`../R/cluster_references.R`](../R/cluster_references.R) — FindClusters with the workflow's call and one argument changed on each object's saved graph, into `data/clustering/<object>/cluster_refs/<setting>/` (the tier decides the settings; `SETTINGS=a,b` picks some), the constructed edge cases into `data/clustering/edge/`, or FindClusters alone at the workflow's settings, timed (`cluster_refs/bench.json`). M10 ran the settings on all 34 objects, one job per tier-P object and three for tier T: 13 to 20 s per tier-T object (14 settings), 0.5 to 4.2 minutes for the tier-P objects up to 26,598 cells and 10.6 to 15.5 minutes for the three larger ones (12 settings), each within 16 GB; the edge cases take under a minute.
- [`m10_cluster.sbatch`](m10_cluster.sbatch) `<object>...`: `scrust cluster --snn --keep-stages` on the graph of every M8 instrumented run (the workflow's settings) and of every reference setting above, each followed by `scrust validate cluster` against its R run and the instrumented optimizer (and, for tier P, the workflow job's logged numbers), into `data/clustering/<object>/rust/m10/<run>/` and `rust/m10/refs/<setting>/`; where R stopped, Rust must stop with R's message. The `resolutions` setting is rerun at one thread and compared byte for byte. `SCRUST_BIN` pins a copy of the binary (M10 used `data/clustering/bin/scrust-m10-5ddbd97`); `SKIP_REFS=1` runs only the instrumented runs.
- [`m10_benchmark.sbatch`](m10_benchmark.sbatch) `<object>...`: R's FindClusters alone against `scrust cluster` at the workflow's settings (one thread) and with five resolutions at 1 and 32 threads, on one node, the five-resolution outputs compared byte for byte, into `data/clustering/<object>/rust/m10/bench/`.
- [`fixtures_cluster.sh`](fixtures_cluster.sh): copies the committed edge-case fixture into `tests/fixtures/clustering/edge/`, with the digests and printed lines of R's files (plain files and Python's standard library; a login node is enough).
