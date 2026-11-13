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

- [`m7_export.sbatch`](m7_export.sbatch): the one-time exports of an object, run with `bash` inside the task: `export_doubletfinder.sbatch` (inputs and compared fields), then `export_doubletfinder_singlets.sbatch` (the saved output's RNA counts), each unless already done. Exports need no particular CPU type (`--partition=shared,sapphire`); 1 CPU and 24 GB cover every object.
- [`m7_panel.sbatch`](m7_panel.sbatch) `<label>`: `scrust doubletfinder` without stage dumps on the object's export, then `scrust validate doubletfinder --saved` against the saved workflow output (the doublet set, the kept cells' pANN, the column names, so pK and nExp, and `singlets/` against the saved RNA counts), into `data/doubletfinder/<object>/rust/<label>/production/` (`run/`, `vs_saved.{json,md}`, `fingerprint.json`, `node_features.txt`, `time_v.txt`, `stderr.txt`). The production profile: the task runs on `sapphire`, stops unless the node has the CPU feature of the object's workflow job and OpenBLAS, left to detect the CPU, reports that CPU's core (`SapphireRapids`); it fails unless all nine checks are made and identical. 8 CPUs and 32 GB cover every object. With a pinned copy of the binary, as for M6.

```bash
c=$(git rev-parse --short HEAD)
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m7-$c
mkdir -p data/doubletfinder/logs/m7
E=$(sbatch --parsable --array=1-150%50 --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_export.sbatch)
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m7-$c sbatch --export=ALL --array=1-150%50 --dependency=aftercorr:$E \
    --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_panel.sbatch m7-$c
```

`--dependency=aftercorr` starts each object's Rust run once its own export task has succeeded.
