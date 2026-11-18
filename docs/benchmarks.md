# R vs Rust: benchmarks

This page compares each Rust port with the workflow's R step it replaces. It reports CPU time, peak memory, the share of each that goes to reading and writing files, and the size of the input.

It grows with the ports. Each new tool adds its rows, with the same method and columns, and the findings are updated (docs/PLAN.md, A7's cycle, step 6).

**Last updated:** 2026-11-18, the first benchmark, covering DoubletFinder (Step 1) and the clustering chain (M11–M13).

## Findings (2026-11-18)

- **CPU time:** at one thread, Rust uses 2.1–7.1× less CPU time than R:
  - the clustering chain: 2.6–7.1× without `--export`, 2.1–4.6× with it;
  - DoubletFinder alone: 2.9–3.6×;
  - the whole DoubletFinder rule: 2.7–3.0×.
- **Peak memory:** Rust uses 8–40× less. R peaks at 3.1–153 GiB, Rust at 0.15–4.5 GiB.
- **Elapsed time:** at 32 threads Rust finishes 11–30× sooner than R. It uses only 3–20% more CPU time than at one thread.
- **R's file I/O:**
  - In the clustering rules, reading and writing take 21–64% of R's CPU time. Most of it is `saveRDS`, which compresses the whole Seurat object on one CPU.
  - DoubletFinder is compute-bound: I/O takes 1–2% of its CPU time.
  - I/O does not drive R's memory. Its peak comes from SCTransform and RunPCA, or from DoubletFinder's final run and parameter sweep.
- **Rust's file I/O:**
  - Without `--export`, reading and writing take 5–18% of Rust's CPU time.
  - With `--export`, they take 13–53%, mostly writing the corrected counts as gzipped text.
  - Reading the counts sets much of Rust's small peak: 35–81% of it is reached by the end of reading.
- **To fix:** at 32 threads, Rust's parallel `.mtx.gz` writer spends 8–22 s of kernel (system) CPU time on the large objects, against at most 2 s at one thread.

## Method

- **R** runs copies of the workflow's own scripts, with phase-marker lines added and nothing else changed (checked byte for byte; sources at pipeline commit `ea053fa`).
  - It runs as the workflow's jobs do: 1 CPU, the six `*_NUM_THREADS` set to 1, the production R env, seed 12345, the production inputs.
  - Steps not yet ported are timed but left out of the comparison: SoupX, emptyDrops with `logNormCounts` and `as.Seurat`, and RunUMAP. Their cost is in the [last table](#workflow-steps-not-yet-ported).
- **Rust** runs `scrust` (binary built from `b6f6e4e`) in the production profile, at `--threads 1` and `--threads 32`.
- **Repeats and nodes:** each case runs 3 times, and the tables give the median. All runs are on Xeon Platinum 8480CL nodes of the `sapphire` partition.
  - Inputs are evicted from the page cache before each run, so every read is cold.
  - Most R runs shared a node with other benchmark runs, as the workflow's 1-CPU jobs do. Repeats vary by up to 25% in CPU time for R and 13% for Rust.
- **Checks:** every run was checked before it was counted.
  - Every Rust output is byte-identical to the reference runs.
  - Every R run reproduces the workflow's number of clusters, and its `.rds` is within 10 bytes of the workflow's saved one.
- **Data and scripts** are in `data/benchmarks/`, which is not in git:
  - `benchmark.md`: the full tables, with ranges, the read and write parts and the 32-thread columns;
  - `benchmark_table.tsv`, `benchmark_runs.tsv` and `benchmark_kinds.tsv`: the numbers;
  - `checks.tsv`: the checks;
  - `scripts/`: the harness, including the sampler `procmon.py`, the R script copies, and the analysis and report scripts.

## How to read the tables

| Column | Meaning |
|---|---|
| **input cells** | The cells entering the chain, or in the DoubletFinder object. |
| **CPU s** | CPU seconds: processor time, user plus system, summed over all threads and child processes. At one thread it is close to the elapsed time. |
| **peak GiB** | The most memory (resident RAM) the run held at any moment. For R, it is the largest of the per-phase peaks, which are reset at each phase through `/proc/self/clear_refs`. For Rust, it is `ru_maxrss`. |
| **I/O: % CPU** | The share of CPU time spent reading the inputs plus writing the outputs. R writes with `saveRDS` and the cluster-id files. Rust writes its output files and the manifest, which re-reads every output to compute its checksum. R's startup (loading packages) counts in the total but not in the I/O share. |
| **I/O: % of peak** | How much of the peak memory was already in use when reading finished: the memory due to loading the data. For R this includes R and its packages (0.4–0.7 GiB). Writing adds little at one thread, usually under 1% of the peak and at most 6% (Rust, the smallest object), so it is not shown. |
| **Rust with `--export`** | The same run, also writing the export the loaders read: a copy of `counts/`, the corrected counts and their log1p as gzipped text, the residuals as `scale_data.npy`, and the model tables. |
| **R ÷ Rust** | The ratio of the medians, against Rust at one thread without `--export`. |

### Clustering chain (SCTransform → PCA → neighbours → clusters → silhouette and purity)

R runs on 1 CPU, Rust on 1 thread.

| object | input cells | R: CPU s / peak GiB | R I/O: % CPU · % of peak | Rust: CPU s / peak GiB | Rust I/O: % CPU · % of peak | Rust with `--export`: CPU s / peak GiB · I/O % CPU | R ÷ Rust: CPU · peak |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 61 / 3.05 | 35 · 35 | 10.3 / 0.15 | 8 · 58 | 13.3 / 0.16 · 21 | 5.9× · 20× |
| soupx_tenx_10xmousebrain | 11,831 | 255 / 12.1 | 64 · 36 | 35.9 / 1.44 | 18 · 80 | 76.0 / 1.44 · 53 | 7.1× · 8.4× |
| filtered_seurat_tenx_10xmousebrain | 11,831 | 199 / 12.5 | 55 · 19 | 36.4 / 1.44 | 18 · 80 | 77.2 / 1.44 · 53 | 5.5× · 8.7× |
| soupx_emptydrops_4500871lane | 26,598 | 289 / 23.9 | 46 · 17 | 79.8 / 1.48 | 10 · 81 | 99.0 / 1.48 · 16 | 3.6× · 16× |
| soupx_emptydrops_sscovelli | 79,409 | 705 / 62.3 | 36 · 10 | 271 / 3.15 | 6 · 58 | 328 / 3.34 · 13 | 2.6× · 20× |
| filtered_seurat_emptydrops_sscovelli | 79,409 | 777 / 58.5 | 34 · 4 | 270 / 3.15 | 5 · 58 | 329 / 3.35 · 13 | 2.9× · 19× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 791 / 105.2 | 21 · 3 | 286 / 2.65 | 5 · 77 | 374 / 3.06 · 17 | 2.8× · 40× |

### DoubletFinder

R runs on 1 CPU, Rust on 1 thread.

- **Alone:** R's startup, `readRDS` and DoubletFinder's sweep and final run, against all of `scrust doubletfinder`. R writes nothing in this part, so its I/O share is reading alone.
- **Whole rule:** the workflow's DoubletFinder script without RunUMAP. It runs DoubletFinder, keeps the singlets, re-clusters them and saves the `.rds`.
  - Rust's side is `scrust doubletfinder` followed by `scrust cluster --export` on its `singlets/`.
  - Their CPU times are added, and the peak is the larger of the two.
  - `singlets/` counts as I/O when written and when read; R keeps the singlets in memory.
- **Peak memory** is the same in both rows. By the end of reading, R has reached 3–10% of it and Rust 35–51%.

| object | input cells | Alone, R: CPU s / peak GiB · I/O % CPU | Alone, Rust: CPU s / peak GiB · I/O % CPU | Whole rule, R: CPU s · I/O % CPU | Whole rule, Rust: CPU s · I/O % CPU | R ÷ Rust: CPU (alone · whole) · peak |
|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 3,108 | 201 / 5.13 · 1 | 62.9 / 0.24 · 3 | 228 · 4 | 74.7 · 6 | 3.2× · 3.0× · 22× |
| soupx_tenx_10xmousebrain | 11,831 | 660 / 23.3 · 1 | 212 / 2.53 · 12 | 793 · 11 | 279 · 21 | 3.1× · 2.8× · 9.2× |
| soupx_emptydrops_4500871lane | 26,598 | 708 / 37.0 · 1 | 247 / 2.39 · 10 | 879 · 7 | 330 · 12 | 2.9× · 2.7× · 15× |
| soupx_emptydrops_sscovelli | 79,409 | 1,446 / 127.7 · 2 | 404 / 4.41 · 11 | 1,973 · 8 | 682 · 12 | 3.6× · 2.9× · 29× |
| cellbender_fromraw_ghirsutum2 | 59,323 | 1,837 / 153.2 · 1 | 566 / 4.02 · 8 | 2,481 · 6 | 892 · 12 | 3.2× · 2.8× · 38× |

### Elapsed time at 1 and 32 threads

Wall-clock seconds, start to exit. R always runs on 1 CPU. At 32 threads Rust's peak is 0.18–3.41 GiB for the clustering and 0.46–4.50 GiB for DoubletFinder.

| tool | object | R | Rust, 1 thread | Rust, 32 threads | R ÷ Rust at 32 threads |
|---|---|---|---|---|---|
| clustering (without / with `--export`) | soupx_tenx_SRR18256555 | 72.7 | 10.4 / 13.5 | 2.71 / 3.20 | 27× / 23× |
| | soupx_tenx_10xmousebrain | 257 | 36.3 / 76.3 | 12.4 / 16.6 | 21× / 16× |
| | filtered_seurat_tenx_10xmousebrain | 201 | 36.7 / 77.6 | 12.9 / 17.5 | 16× / 11× |
| | soupx_emptydrops_4500871lane | 293 | 80.3 / 99.5 | 20.2 / 23.5 | 14× / 12× |
| | soupx_emptydrops_sscovelli | 709 | 272 / 329 | 54.5 / 64.0 | 13× / 11× |
| | filtered_seurat_emptydrops_sscovelli | 785 | 271 / 330 | 51.5 / 61.4 | 15× / 13× |
| | cellbender_fromraw_ghirsutum2 | 800 | 288 / 376 | 44.5 / 60.0 | 18× / 13× |
| DoubletFinder (alone / whole rule) | soupx_tenx_SRR18256555 | 203 / 229 | 63.1 / 75.1 | 8.25 / 11.2 | 25× / 21× |
| | soupx_tenx_10xmousebrain | 664 / 798 | 213 / 280 | 43.3 / 58.0 | 15× / 14× |
| | soupx_emptydrops_4500871lane | 713 / 884 | 248 / 331 | 37.7 / 56.9 | 19× / 16× |
| | soupx_emptydrops_sscovelli | 1,452 / 1,981 | 406 / 685 | 52.2 / 109 | 28× / 18× |
| | cellbender_fromraw_ghirsutum2 | 1,844 / 2,491 | 569 / 897 | 62.1 / 116 | 30× / 22× |

### Workflow steps not yet ported

R's CPU seconds for the steps inside the benchmarked scripts that no port covers yet (median of the same 3 runs). This is today's cost for the ports still to come; each port moves its rows into a table of its own above.

| R step | rule (script) | CPU s, by object |
|---|---|---|
| SoupX | `soupx`, `soupx_emptydrops` (soupx.R) | SRR18256555 13.2; 10xmousebrain 90.5; 4500871lane 96.6; sscovelli 137.6 |
| reading SoupX's own inputs (the raw matrix and the base `.rds`) | the same | 9.6–68.3 (shown in R's read share above) |
| emptyDrops, with `logNormCounts` and `as.Seurat` | `emptydrops` (emptydrops.R) | sscovelli 215.6 |
| RunUMAP | every clustering rule; the DoubletFinder rules | clustering 9.4–95.7; DoubletFinder 7.6–81.4 |

Not measured yet: decontX, scDblFinder, CellBender (GPU), the marker steps (FindMarkers) and the post-hoc filtering.

## Updating this page

- **With each port:**
  1. Benchmark the R step and the port with this method.
  2. Add the tool's rows, using these columns, and its elapsed times.
  3. Move it out of the not-yet-ported table.
  4. Update the findings and the date.
  5. Record the scrust commit and the R env.
- **Chained tools:** when M15 runs the singlets' re-clustering inside `scrust doubletfinder`, add its row next to the whole rule's two-command route. That row measures what the in-process hand-off saves.
- **The harness:** with the next benchmark, review `data/benchmarks/scripts/` and move it into the repository (`validation/benchmark/`), so that later tables are made the same way.
- **Re-measure** the parallel writer once its kernel time is addressed.
- **Re-run everything** only if the method changes, such as another node type or R env, and say so here.
