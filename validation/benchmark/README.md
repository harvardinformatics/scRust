## Notes

The benchmark's harness (docs/benchmarks.md describes the method; this file is appended to `data/benchmarks/benchmark.md` by `report.py`). Submit the jobs from the repository's top directory.

### Findings

- **CPU time:**
  - R ÷ Rust at one thread:
    - the clustering chain: 2.6–7.1× without `--export`, 2.1–4.6× with it;
    - DoubletFinder alone: 2.9–3.6×;
    - the whole DoubletFinder rule: 2.7–3.0×.
  - At 32 threads, Rust uses 3–20% more CPU time than at one thread and finishes 11–30× sooner than R in wall time.
- **Peak memory:**
  - R ÷ Rust: 8–40× for the clustering, 9–38× for DoubletFinder.
  - R peaks at 3.1 to 153 GiB, Rust at 0.15 to 4.5 GiB.
- **R's I/O:**
  - **Clustering rules:** reading and writing take 21–64% of R's CPU time.
    - Writing takes 15–43%. It is mostly `saveRDS`, which serializes and gzips the whole Seurat object on one CPU.
    - Reading takes 2–27%. The SoupX rules also read the base object's `.rds` and the raw matrix.
  - **DoubletFinder:** I/O takes 1–2% of R's CPU time alone and 4–11% for the whole rule, because DoubletFinder is compute-bound.
  - **Memory:** I/O does not drive R's memory.
    - By the end of reading, R holds 3–36% of its peak: the input plus R's packages.
    - Writing adds nothing.
    - The peak comes from the computation: SCTransform or RunPCA in the clustering, and DoubletFinder's final run (its distance matrix) or its sweep.
- **Rust's I/O:**
  - **Clustering, CPU time at one thread:**
    - Without `--export`, reading and writing take 5–18% of the CPU time: reading the gzipped counts 3–15%, writing 2–4%.
    - With `--export`, they take 13–53%. Most of it is writing the corrected counts as gzipped text: 46% of the CPU time on the two 11,831-cell objects, whose corrected counts are the densest.
  - **DoubletFinder, CPU time:** I/O takes 3–12% alone and 6–22% for the whole rule. The whole rule writes `singlets/` and reads them back.
  - **Memory:** reading sets much of Rust's small peak.
    - At one thread, 35–81% of the peak is reached by the end of reading: the parsed counts and their decompressed text.
    - Writing adds at most 6% at one thread, and up to 32% at 32 threads, where each thread holds compression buffers.
- **Rust's parallel writer:**
  - At 32 threads, writing the neighbour graphs and the export's matrices costs 8–22 s of system CPU time on the large objects, against at most 2 s at one thread.
  - That is most of the extra CPU time its writing takes at 32 threads (`sys_write_s` in `benchmark_kinds.tsv`). It is kernel time from the threads' memory handling, not waiting for the disk.
  - This is room for improvement in the writer, not a measurement artefact.

### What each side runs (like for like)

- **R** runs the preprocessing workflow's own scripts as its jobs run them. The scripts are copies with phase-marker lines added and nothing else changed: `validation/benchmark/workflow/scripts/`, made by `make_workflow_copies.py`, which checks that removing the markers gives back each source byte for byte. The sources are the pipeline repo's `workflow/scripts/` at HEAD `ea053fa` (clean); each one's git blob id is in `validation/benchmark/workflow/SOURCES.tsv`.
  - **The workflow's settings:**
    - 1 CPU, as `profiles/slurm` sets `cpus_per_task: 1` for these rules;
    - `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and the other four `*_NUM_THREADS` variables set to 1, as Snakemake's `shell()` sets them to the rule's threads;
    - `OPENBLAS_CORETYPE` unset;
    - the production R env (R 4.5.3, Seurat 5.5.1, DoubletFinder 2.0.6);
    - `SCRNASEQ_PREPROCESS_SEED=12345`;
    - each rule's inputs from the production run, read only.
  - **The outputs** are written under `data/benchmarks/`. The large ones are deleted after their sizes are recorded.
  - **Checks** (`checks.tsv`): each clustering run's number of clusters equals the workflow log's. Each DoubletFinder run's equals the number Rust's clustering of the same singlets finds. Each output `.rds` is within 10 bytes of the size of the workflow's saved `.rds`; the objects record the time of each command.
- **Rust** runs the binary built from main's code as verified for M13 (`bin/scrust-b6f6e4e`, sha256 `6b373bdc…`), in the production profile (OpenBLAS's SapphireRapids kernel). Each run uses `--threads 1` or `--threads 32` and gets 32 CPUs either way, as in M13's benchmark.
  - **Its inputs** are the counts each chain starts from: `data/clustering/<object>/input`, and for DoubletFinder `data/doubletfinder/<object>/input` (the counts plus `seurat_clusters.tsv`). Both were exported once from the saved objects.
  - **Checks:** each run's outputs equal the reference runs' byte for byte. For the clustering that is M13's benchmark run. For DoubletFinder it is R's direct run: `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv`. The export's corrected-count matrices differ only in a comment line that names the input's path, so their data lines are compared.
- **The clustering chain:** R's side is the rule's script without the steps no port does. Those steps are SoupX, in `soupx.R` for the rules soupx and soupx_emptydrops; emptyDrops with logNormCounts and `as.Seurat`, in `emptydrops.R`; and RunUMAP, in every script. They are timed as their own phases and shown apart, in "R's scripts, whole".
  - **R's side therefore covers:** R's startup and packages; reading the rule's inputs into the object; PercentageFeatureSet, SCTransform, RunPCA, FindNeighbors, FindClusters, and silhouette with purity; `saveRDS` and the two cluster-id files (emptydrops also writes its filtered matrix directory); and R's exit.
  - **Rust's side** is `scrust cluster --counts <input> --preset workflow`, without `--export` and with it, as two rows. With `--export` is the like for like for the outputs: R's `.rds` holds the whole SCT assay, the corrected counts included. Without `--export` is the chain and its own files alone.
- **DoubletFinder, in two rows.** The workflow's `doubletfinder.R` runs DoubletFinder (readRDS, paramSweep … doubletFinder, subset to the singlets), then clusters the singlets again (SCTransform → … → purity, plus RunUMAP), then saveRDS.
  - **DoubletFinder alone.**
    - R's side: R's startup and packages, readRDS, and the DoubletFinder phases: paramSweep, summarizeSweep and find.pK, nExp and doubletFinder, and the subset to the singlets.
    - R writes nothing in this part, because its outputs are saved only after the re-clustering. R's write share is therefore shown as "none".
    - Rust's side: `scrust doubletfinder`. Its write is `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` and `singlets/`, plus its manifest and summary.json.
  - **The whole rule:** DoubletFinder plus the re-clustering of the singlets, without RunUMAP. RunUMAP is shown apart, in "R's scripts, whole".
    - R's side: the script without RunUMAP.
    - Rust's side: `scrust doubletfinder`, then `scrust cluster --counts <its singlets/> --preset workflow --export`. The export is used because R's `.rds` holds the SCT assay.
    - These are two processes. Their CPU and wall times add up, and the peak is the larger of the two.
    - The `singlets/` hand-off counts as I/O at both ends: written by `scrust doubletfinder`, read by `scrust cluster`. R keeps the singlets in memory instead. M15's in-process route (below) keeps them in memory too.
- **M15: clustering on the fly and the re-clustering in one process** (tool `doubletfinder_m15`, one binary for its four variants, so they compare like for like):
  - `d`: `scrust doubletfinder --clusters` (DoubletFinder alone, the labels given), and `r`: `scrust cluster --counts <its singlets/> --preset workflow --export` (the two-command route's second command), as above;
  - `h`: `scrust doubletfinder` without `--clusters`: it clusters its own input in process first (clustering on the fly). Against `d`, it shows the hook's cost;
  - `p`: `scrust doubletfinder --clusters ... --recluster --recluster-export`: the whole rule in one process, the singlets handed to the clustering from memory, with the export as `r` writes it. Against `d` + `r`, it shows what the in-process hand-off saves; against R's whole rule, the port's gain.
  - The in-process clusterings record their stages on the command's clock (`clustering: <stage>`, `recluster: <stage>`), which `analyze.py` classifies as `r`'s stages. From M15 on, the chain's stage files are written in its write phase, not after each stage; `analyze.py` knows the binaries that wrote them after each stage (`STAGE_WRITER`), and applies the pca writer's calibration to those only.

### Definitions

- **Input cells:** the cells entering the clustering chain, or for DoubletFinder the cells of the input object. **Genes:** the features of the count matrix. **Nonzeros:** its stored entries. All three come from the object's export record.
- **CPU time:** user + system time of the process, over all its threads and its waited-for children (getrusage: R's `proc.time()` at each mark, `wait4` at the end). **Wall time:** from the start of the process to its exit. Both are medians over the repeats.
- **Peak memory:** peak RSS (VmHWM).
  - R: the largest of the per-phase peaks. The benchmark resets VmHWM at each phase's start by writing 5 to `/proc/self/clear_refs`; this works on these nodes' 4.18 kernel, and was tested.
  - Rust: `ru_maxrss` from `wait4`, which equals VmHWM. Rust's VmHWM is never reset.
- **Phases:**
  - **read:** the inputs read and turned into the in-memory object.
    - R: from the first input read to the object the next step uses: Read10X, readRDS, read10xCounts or Read_CellBender_h5_Mat, plus CreateSeuratObject where it is built straight from the read matrix.
    - Rust: the load phase. That is a digest pass over the inputs, or with `--export` the `counts/` copy's read, then the parse. DoubletFinder reads the counts again for `singlets/` and for its input digests; that is counted as read too.
  - **compute:** the computation.
  - **write:** every output written.
    - R: `saveRDS` (serialization with gzip, on one CPU), the cluster-id files, and emptydrops' 10x matrix directory.
    - Rust: each stage's files, which it writes right after that stage's computation; the labels and metadata; summary.json; with `--export`, the `counts/` copy and `sct/`; and the manifest, which re-reads every output to hash it.
  - **startup** (R and its packages; Rust's process start) and **exit:** in the total, not in the I/O shares.
- **CPU time due to I/O:** the read phase's CPU time and the write phase's CPU time, each as a % of the total CPU time, and their sum. The table shows them as `read + write = sum`. Each is rounded on its own, so the parts can differ from the sum by a point.
- **Memory due to reading:** the peak RSS reached by the end of the read phase, as a % of the run's peak. It covers the input held in memory, the read buffers, and for R the interpreter and its packages.
- **Memory added by writing:** the write phase's peak minus the RSS when writing began, as a % of the run's peak. Rust writes after each stage, so this is the largest such increase over its write stretches.
- **R ÷ Rust:** the ratio of the medians.

### Caveats

- **Different files, the same computation.**
  - **What each side reads:**
    - R reads the rule's original inputs, as the workflow does. Depending on the rule, that is the 10x directories, the base object's `.rds`, the raw droplet matrix, the CellBender `.h5`, or the input object's `.rds`.
    - Rust reads the counts the chain starts from, exported once from the saved objects as gzipped Matrix Market, plus the cluster labels for DoubletFinder.
  - **The bytes differ.** `benchmark_table.tsv` gives them in `R_input_bytes` and `Rust*_input_bytes`.
  - **R's read includes SoupX's own inputs.** For the SoupX rules, that is the raw matrix and the base object, and "R's scripts, whole" gives their share.
  - **What is compared is the same.** Every Rust run's outputs equal the reference runs' byte for byte. Every R run reproduces the workflow's clusters, and an `.rds` of the workflow's size.
- **R's steps that no port does** (SoupX; emptyDrops with logNormCounts and `as.Seurat`; RunUMAP) are left out of R's columns and shown apart. R's process still carries their memory: the chain's phases start from whatever those steps left allocated.
- **Lustre.**
  - Inputs and outputs are all on `/n/holylfs05`.
  - Each run's inputs were evicted from the node's page cache first, using posix_fadvise DONTNEED, checked with mincore. So every run reads from the file system, as a workflow job on a fresh node would.
  - Writes land in the page cache and are flushed in the background, with no fsync except Rust's `counts/` copy.
  - I/O wall times vary with Lustre's load, and so does R's package loading (its wall time ranges widely while its CPU time stays about the same). CPU times vary much less.
- **R's measurements are exact per phase.** CPU comes from `proc.time()` at each mark, and peaks from VmHWM reset at each phase. A mark costs about 2 ms, and a run has about 17. R's getrusage `maxrss` is not used, because the resets change it.
- **Rust's measurements are sampled per phase; its totals are exact.**
  - **Per-phase CPU:** from `/proc/<pid>/stat` sampled every 10 ms, aligned with `summary.json`'s phase marks to within about 10 ms.
  - **Writing inside a stage:** found from the files open and the bytes written. Rust formats a text file in memory before it opens it, so the formatting of a stage's first file is invisible to the sampler.
  - **The pca calibration.** The largest such file is `pca/embeddings.tsv`. `scrust pca` writes the same files with the same function and times that as a phase of its own (`runs/calib/pca`), so the unseen part is moved from compute to write. The other stages' first files are small.
  - **Peaks and totals:** the run's peak and total CPU are exact (getrusage). Per-phase peaks are sampled, unless VmHWM rose within the phase.
- **Rust's write includes its bookkeeping:** the manifest re-reads every output to hash it, and it writes summary.json. R writes no checksums.
- **Compression:** R's `saveRDS` compresses on one CPU. Rust compresses its `.gz` outputs on all its threads.
- **Threads.**
  - R runs on one CPU, as the workflow's jobs do.
  - Rust at 1 thread had 32 CPUs available but used one.
  - At 32 threads, Rust's CPU time includes the cost of parallel work, so it is higher than at one thread while its wall time is lower.
- **DoubletFinder's `singlets/`:** in the whole-rule row, Rust writes the singlets to disk and reads them back for the clustering; both are counted as I/O. R keeps them in memory.
- **Memory due to reading, for R,** includes R and its packages: 0.4 to 0.7 GiB before any input is read. Rust's binary takes a few MiB.
- **The workflow's own jobs (sacct)** are given for comparison in "R's scripts, whole".
  - Their step `.0` is Snakemake's job step, which also runs Snakemake's Python process and the conda activation.
  - Their MaxRSS is the step's processes together, sampled every 30 s by Slurm's cgroup plugin (`JobAcctGatherFrequency task=30`, page cache excluded). It misses short peaks such as SCTransform's, so the exact peaks here can exceed it.
  - Those jobs ran once each, in November, on nodes shared with the rest of the workflow, so differences of ±20% are within the repeat spread seen here.

### Files

All files are under `data/benchmarks/`, which git ignores.

- **`benchmark_table.tsv`:** one row per tool × object (× Rust variant), with medians, `_min` and `_max`.
- **`benchmark_runs.tsv`:** one row per run and segment (phase), with the raw measurements: CPU, system CPU, wall, RSS at start, peak and end, and the `/proc/<pid>/io` bytes. Runs on other CPU models are marked `included = no`.
- **`benchmark_kinds.tsv`:** one row per run, with its totals by kind and its I/O shares. The runs on 8480+ nodes are in `benchmark_other_cpu.tsv` instead.
- **`checks.tsv`:** each run's output checks.
- **`sacct_workflow.tsv`:** the workflow jobs' accounting.
- **`runs/`:** every run's own records: markers, sampler, `env.txt`, `summary.json` and `outputs.tsv`. `runs/calib/pca` holds the pca writer's calibration.
- **The harness,** in the repository's `validation/benchmark/` since M15 (it was `data/benchmarks/scripts/`, which is kept as it ran for the first tables; moved, its paths are the repository's, and the binary is `$SCRUST_BIN`):
  - `procmon.py`: the sampler;
  - `workflow/scripts/`: the R copies plus `bench_markers.R`, made by `make_workflow_copies.py`, with `SOURCES.tsv`;
  - `run_r.sbatch`, `run_rust.sbatch` and `run_rust_batch.sbatch`: the jobs, with the submit scripts and the run list `r_runs.tsv`;
  - `pca_calib.sbatch`: the calibration;
  - `test_proc_features.sbatch`: the check of the `/proc` features;
  - `analyze.py`, `checks.py`, `report.py` and `sacct_workflow.sh`: the tables.
- **Regenerating,** from the repository's top directory (`BENCH_OUT=<dir>` writes the tables and `benchmark.md` there instead of `data/benchmarks/`):

  ```
  PYTHONNOUSERSITE=1 python3 -I validation/benchmark/analyze.py && PYTHONNOUSERSITE=1 python3 -I validation/benchmark/checks.py > data/benchmarks/checks.tsv && PYTHONNOUSERSITE=1 python3 -I validation/benchmark/report.py
  ```
- **Pilots:** `runs/*/*/soupx_tenx_SRR18256555/rep0*` are the pilots. They are left out of the tables.
