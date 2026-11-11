# DoubletFinder core on R's own embeddings (milestone M3)

M3 ports the part of DoubletFinder 2.0.6 that follows the PCA (docs/PLAN.md B1, B3 and B6): the neighbour ranking (stage S4), pANN in the sweep and the final run (S5), the sweep summary, BCmvn and the choice of pK (S6), and the homotypic nExp and the Doublet/Singlet calls (S7). It runs on R's own inputs, *teacher forcing* (plan A5b): each instrumented R run's sweep subsample (stage S1) and PCA embeddings (stage S3) feed the Rust core, and everything after them is Rust's. SCTransform and PCA are not ported here (M4 and M5).

## Outcome

| Done when (plan B6, and this milestone's brief) | Result |
|---|---|
| S4-S7 exact on all 10 objects, every instrumented run's dumps | **Yes, all 24 instrumented runs of the 10 objects.** Every block of every stage array the R runs dumped is bit-identical: S4 digests of every real cell (and, on tier T, R's full order positions, 3.1 million of them), S5 pANN at every pK of every sweep iteration and of the final run, S6 `bkde` grids and densities, the values passed to `bimodality_coefficient` and BCreal for every (pN, pK), S7 calls. S6 and S7 are also identical when recomputed from R's own S5. No near-tie needed explaining. |
| Rust outputs equal the R reference runs' `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json` result fields | **Yes.** All three files are byte-identical to the instrumented run's and to **every direct R run** of the same object (portable and production profiles, Sapphire Rapids and Cascade Lake, wherever M2 ran them); 19 of the 20 `result` fields are identical, and the 20th, `rng_end_sha256` (the digest of R's final `.Random.seed`), does not apply: no random numbers are drawn under teacher forcing. Through R's cross-check (M2), the Rust outputs also reproduce the saved workflow outputs: doublet sets, singlet pANN bit for bit and column names, checked directly as well. |
| `R_strtod` port bit-exact on R vectors | **Yes**, on 12,785 strings: every label of the pK and pN grids and their re-parsed forms, every 4-decimal label from 0.0001 to 0.5, 7,550 random decimal strings of all magnitudes and precisions, special and malformed input, and 65 strings that R does **not** round correctly (found among 1.6 million random strings). |
| kNN and core memory and time measured on sscovelli | **Yes** (below). The final run's kNN alone: 67 s and 32 MiB at one thread, 2.4 s and 142 MiB at 32 threads; the whole core: 87 s and 66 MiB at one thread, 3.2 s and 187 MiB at 32 threads. R needs 89.7 GB for the distance matrix of that final run alone. |
| `cargo test` passes | **Yes**, including the tier-T reproduction of every committed reference file and stage array, and new R-vector tests. |

Beyond the brief, S1 is settled too: rcompat redraws the sweep subsample and every doublet draw from R's own `.Random.seed` before them, identically, and ends in R's state after them (all 24 runs).

## What was built

| Crate | M3 contents |
|---|---|
| `rcompat` | `strtod`: R's decimal parser (`R_strtod5`) and `as.numeric()` of strings, emulated operation by operation on the x87 long double; `fft`: Singleton's mixed-radix transform as R's `fft.c` computes it (goto-for-goto, interleaved layout with increment 2); `dnorm` (with R's NA handling) and `approx`/`approxfun` (`regularize.values`, `approx1`). |
| `scio` | `.npy` read and write in the R scripts' layout, C's `%.17g` (the reference files' number format), TSV, `summary.json` in the R scripts' layout plus a parser, and a standard-library gzip/DEFLATE reader for the 10x files (byte-identical to the system `gzip` on `matrix.mtx.gz`). |
| `knn` | Exact neighbours with `fields::rdist`'s arithmetic and R's `order()` (ties to the lower index), an O(N) buffer per thread, parallel over cells with output independent of the thread count; the S4 digest of the R hooks. |
| `kde` | `KernSmooth::bkde` (normal kernel, its defaults; Fortran `linbin` with gfortran's `INT`; convolution through rcompat's `fft`), `approxfun` of it, and DoubletFinder's `skewness`, `kurtosis` and `bimodality_coefficient`. |
| `doubletfinder` | The grids and the values that pass through strings; the sweep and final pANN; `summarizeSweep`, `find.pK`, `modelHomotypic` and nExp; the calls; R's column names; the workflow's call sequence over an `EmbeddingSource` (R's stage dumps now, Rust's SCTransform and PCA in M6); stage dumps read and written in R's layout; the comparison with an R run. |
| `scrust` | `scrust doubletfinder --from-stages <dir>` (the core on an instrumented R run's embeddings, writing the reference runs' layout and optionally the stage arrays), `scrust validate doubletfinder` (Rust against R, JSON and markdown report), `scrust knn` (the neighbour stage alone). |

No external crate was added: the workspace still uses only the Rust standard library.

## How the comparison works

For each instrumented R run, `scrust doubletfinder` reads the run's sweep subsample (`S1_subsample.npy`, when n > 10,000) and the PCA embeddings of the six sweep iterations and of the final run (`S3_embeddings.npy`), computes S4 to S7 and writes `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv`, `summary.json` and its own stage arrays (`--dump-stages`, same names and layout as R's). `scrust validate doubletfinder` then compares, against the same R run ([`validation/slurm/m3_doubletfinder.sbatch`](../../validation/slurm/m3_doubletfinder.sbatch); reports in `data/doubletfinder/<object>/rust/m3/<run>/validate.{json,md}`):

- the three tables byte for byte, and the `result` section of `summary.json` field by field, as text;
- every block of S4_digest, S4_positions (tier T), S5_pann_sweep, S5_pann_final, S6_bkde, S6_gkde and S7_doublet that R dumped, bit for bit. The Rust chain feeds each stage from its own previous stage; since every stage is identical, this is the same as teacher-forcing each stage;
- S6 and S7 recomputed from R's own S5 (BCreal, the `bimodality_coefficient` input, the calls from R's final pANN and nExp), so that no S6 or S7 difference could hide behind an S5 one;
- S1 redrawn from R's `.Random.seed` checkpoints.

Every criterion is exact (`validation/tolerances.toml`, S1 and S4-S7 now settled).

## Results: Rust on R's embeddings, against the same R run

Columns give `=` (identical) and the number of values compared; `S1 subsample` is `none` when n <= 10,000 (no subsample is drawn, and R's state is unchanged, which is checked).

| object | run | S1 subsample | S1 doublets | S4 digests | S5 sweep | S5 final | S6 bkde | S6 gkde | S6 BCreal (R's S5) | S7 calls | S7 (R's S5) | cells / sweep_stats / bcmvn | result fields | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_test | portable-cascadelake | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| soupx_tenx_test | portable-sapphirerapids | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| soupx_tenx_test | production-sapphirerapids | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_test | portable-cascadelake | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_test | portable-sapphirerapids | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_test | production-sapphirerapids | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| decontx_raw_tenx_test | portable-cascadelake | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| decontx_raw_tenx_test | portable-sapphirerapids | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| decontx_raw_tenx_test | production-sapphirerapids | = (none) | = (1,012) | = (4,200) | = (48,600) | = (300) | = (129,924) | = (48,600) | = (162) | = (300) | = (300) | = / = / = | 19 / 19 | yes |
| decontx_raw_emptydrops_test | portable-cascadelake | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| decontx_raw_emptydrops_test | portable-sapphirerapids | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| decontx_raw_emptydrops_test | production-sapphirerapids | = (none) | = (900) | = (3,738) | = (43,254) | = (267) | = (129,924) | = (43,254) | = (162) | = (267) | = (267) | = / = / = | 19 / 19 | yes |
| cellbender_fromraw_test | portable-cascadelake | = (none) | = (940) | = (3,892) | = (45,036) | = (278) | = (129,924) | = (45,036) | = (162) | = (278) | = (278) | = / = / = | 19 / 19 | yes |
| cellbender_fromraw_test | portable-sapphirerapids | = (none) | = (940) | = (3,892) | = (45,036) | = (278) | = (129,924) | = (45,036) | = (162) | = (278) | = (278) | = / = / = | 19 / 19 | yes |
| cellbender_fromraw_test | production-sapphirerapids | = (none) | = (940) | = (3,892) | = (45,036) | = (278) | = (129,924) | = (45,036) | = (162) | = (278) | = (278) | = / = / = | 19 / 19 | yes |
| soupx_tenx_SRR18256555 | portable-cascadelake | = (none) | = (10,476) | = (43,512) | = (578,088) | = (3,108) | = (149,172) | = (578,088) | = (186) | = (3,108) | = (3,108) | = / = / = | 19 / 19 | yes |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | = (none) | = (10,476) | = (43,512) | = (578,088) | = (3,108) | = (149,172) | = (578,088) | = (186) | = (3,108) | = (3,108) | = / = / = | 19 / 19 | yes |
| soupx_tenx_SRR18256555 | production-sapphirerapids | = (none) | = (10,476) | = (43,512) | = (578,088) | = (3,108) | = (149,172) | = (578,088) | = (186) | = (3,108) | = (3,108) | = / = / = | 19 / 19 | yes |
| soupx_tenx_10xmousebrain | portable-sapphirerapids | = (10,000) | = (34,930) | = (143,662) | = (1,920,000) | = (11,831) | = (153,984) | = (1,920,000) | = (192) | = (11,831) | = (11,831) | = / = / = | 19 / 19 | yes |
| soupx_tenx_10xmousebrain | production-sapphirerapids | = (10,000) | = (34,930) | = (143,662) | = (1,920,000) | = (11,831) | = (153,984) | = (1,920,000) | = (192) | = (11,831) | = (11,831) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_4500871lane | portable-sapphirerapids | = (10,000) | = (44,774) | = (173,196) | = (1,980,000) | = (26,598) | = (158,796) | = (1,980,000) | = (198) | = (26,598) | = (26,598) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_4500871lane | production-sapphirerapids | = (10,000) | = (44,774) | = (173,196) | = (1,980,000) | = (26,598) | = (158,796) | = (1,980,000) | = (198) | = (26,598) | = (26,598) | = / = / = | 19 / 19 | yes |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | = (10,000) | = (79,982) | = (278,818) | = (1,980,000) | = (79,409) | = (158,796) | = (1,980,000) | = (198) | = (79,409) | = (79,409) | = / = / = | 19 / 19 | yes |
| cellbender_fromraw_ghirsutum2 | portable-sapphirerapids | = (10,000) | = (66,590) | = (238,646) | = (1,980,000) | = (59,323) | = (158,796) | = (1,980,000) | = (198) | = (59,323) | = (59,323) | = / = / = | 19 / 19 | yes |

S4 order positions, dumped in full by the tier-T runs (every position of every real cell, all 7 contexts):

| object | run | contexts | positions | differing |
|---|---|---|---|---|
| soupx_tenx_test | portable-cascadelake | 7 | 242,700 | 0 |
| soupx_tenx_test | portable-sapphirerapids | 7 | 242,700 | 0 |
| soupx_tenx_test | production-sapphirerapids | 7 | 242,700 | 0 |
| soupx_emptydrops_test | portable-cascadelake | 7 | 191,973 | 0 |
| soupx_emptydrops_test | portable-sapphirerapids | 7 | 191,973 | 0 |
| soupx_emptydrops_test | production-sapphirerapids | 7 | 191,973 | 0 |
| decontx_raw_tenx_test | portable-cascadelake | 7 | 218,700 | 0 |
| decontx_raw_tenx_test | portable-sapphirerapids | 7 | 218,700 | 0 |
| decontx_raw_tenx_test | production-sapphirerapids | 7 | 218,700 | 0 |
| decontx_raw_emptydrops_test | portable-cascadelake | 7 | 186,099 | 0 |
| decontx_raw_emptydrops_test | portable-sapphirerapids | 7 | 186,099 | 0 |
| decontx_raw_emptydrops_test | production-sapphirerapids | 7 | 186,099 | 0 |
| cellbender_fromraw_test | portable-cascadelake | 7 | 202,662 | 0 |
| cellbender_fromraw_test | portable-sapphirerapids | 7 | 202,662 | 0 |
| cellbender_fromraw_test | production-sapphirerapids | 7 | 202,662 | 0 |

## Against the direct R runs and the saved workflow output

The Rust outputs from every instrumented run's embeddings, against every direct R run of the same object (`cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` byte for byte) and against the saved workflow output (the doublet set, the pANN of every kept cell as text, the column names).

| object | Rust run (embeddings of) | direct R runs with identical files | saved: doublet set | saved: singlet pANN | saved: columns |
|---|---|---|---|---|---|
| soupx_tenx_test | portable-cascadelake | 3 / 3 | = (24) | = (276) | `pANN_0.25_0.29_24` |
| soupx_tenx_test | portable-sapphirerapids | 3 / 3 | = (24) | = (276) | `pANN_0.25_0.29_24` |
| soupx_tenx_test | production-sapphirerapids | 3 / 3 | = (24) | = (276) | `pANN_0.25_0.29_24` |
| soupx_emptydrops_test | portable-cascadelake | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.28_30` |
| soupx_emptydrops_test | portable-sapphirerapids | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.28_30` |
| soupx_emptydrops_test | production-sapphirerapids | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.28_30` |
| decontx_raw_tenx_test | portable-cascadelake | 3 / 3 | = (27) | = (273) | `pANN_0.25_0.09_27` |
| decontx_raw_tenx_test | portable-sapphirerapids | 3 / 3 | = (27) | = (273) | `pANN_0.25_0.09_27` |
| decontx_raw_tenx_test | production-sapphirerapids | 3 / 3 | = (27) | = (273) | `pANN_0.25_0.09_27` |
| decontx_raw_emptydrops_test | portable-cascadelake | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.22_30` |
| decontx_raw_emptydrops_test | portable-sapphirerapids | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.22_30` |
| decontx_raw_emptydrops_test | production-sapphirerapids | 3 / 3 | = (30) | = (237) | `pANN_0.25_0.22_30` |
| cellbender_fromraw_test | portable-cascadelake | 3 / 3 | = (31) | = (247) | `pANN_0.25_0.23_31` |
| cellbender_fromraw_test | portable-sapphirerapids | 3 / 3 | = (31) | = (247) | `pANN_0.25_0.23_31` |
| cellbender_fromraw_test | production-sapphirerapids | 3 / 3 | = (31) | = (247) | `pANN_0.25_0.23_31` |
| soupx_tenx_SRR18256555 | portable-cascadelake | 3 / 3 | = (420) | = (2,688) | `pANN_0.25_0.005_420` |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | 3 / 3 | = (420) | = (2,688) | `pANN_0.25_0.005_420` |
| soupx_tenx_SRR18256555 | production-sapphirerapids | 3 / 3 | = (420) | = (2,688) | `pANN_0.25_0.005_420` |
| soupx_tenx_10xmousebrain | portable-sapphirerapids | 2 / 2 | = (1,646) | = (10,185) | `pANN_0.25_0.07_1646` |
| soupx_tenx_10xmousebrain | production-sapphirerapids | 2 / 2 | = (1,646) | = (10,185) | `pANN_0.25_0.07_1646` |
| soupx_emptydrops_4500871lane | portable-sapphirerapids | 1 / 1 | = (3,359) | = (23,239) | `pANN_0.25_0.02_3359` |
| soupx_emptydrops_4500871lane | production-sapphirerapids | 1 / 1 | = (3,359) | = (23,239) | `pANN_0.25_0.02_3359` |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | 1 / 1 | = (9,909) | = (69,500) | `pANN_0.25_0.19_9909` |
| cellbender_fromraw_ghirsutum2 | portable-sapphirerapids | 1 / 1 | = (6,769) | = (52,554) | `pANN_0.25_0.03_6769` |

## Benchmark

soupx_emptydrops_sscovelli (79,409 cells; 105,879 with the artificial doublets of the final run; k = 20,117), the embeddings of `instrumented-portable-sapphirerapids`, Intel Xeon Platinum 8480+ (`sapphire`, 32 CPUs allocated, job 59731336, [`validation/slurm/m3_benchmark.sbatch`](../../validation/slurm/m3_benchmark.sbatch)); wall time, CPU time and peak RSS from `/usr/bin/time -v`. *kNN* is `scrust knn` on the final run's embeddings: every real cell's k + 1 nearest positions and their doublet count (pANN's numerator), the computation behind R's 89.7 GB `fields::rdist` matrix. *Whole core* is `scrust doubletfinder --from-stages` without stage dumps: loading, the six sweep iterations (10,000 real cells each), `summarizeSweep`, `find.pK`, nExp, the final run and the outputs. Every core run's three tables were byte-identical to R's.

| what | threads | wall | CPU time | CPU % | peak RSS | sweep (6 pN) | final run |
|---|---|---|---|---|---|---|---|
| kNN: final run (79,409 queries, 105,879 cells, k + 1 = 20,118) | 1 | 67.3 s | 66.8 s | 99% | 32 MiB |  |  |
| kNN: final run (79,409 queries, 105,879 cells, k + 1 = 20,118) | 32 | 2.4 s | 72.1 s | 2981% | 142 MiB |  |  |
| whole core (sweep, summary, pK, nExp, final run, outputs) | 1 | 86.6 s | 86.0 s | 99% | 66 MiB | 15.8 s | 69.7 s |
| whole core (sweep, summary, pK, nExp, final run, outputs) | 32 | 3.2 s | 89.7 s | 2808% | 187 MiB | 0.7 s | 2.3 s |

R for comparison, the direct production run of M2 (one thread): wall 36.3 min, peak RSS 126.9 GiB; `paramSweep` 10.7 min and `doubletFinder` 24.6 min, but both include SCTransform and PCA, which M3 does not port, so the times are not comparable stage for stage. The memory is: R's final run holds the 105,879 x 105,879 distance matrix (89.7 GB) at its peak; Rust's kNN needs 32 MiB at one thread.

- **Threads.** 32 threads run 27-28 times faster than one. The results do not depend on the thread count (tested with 1 to 16 threads in `cargo test`, and the 1- and 32-thread runs here give identical files).
- **Memory.** The points, about 4 MB of buffers per thread and the per-cell results: 32 MiB at one thread, 142 MiB at 32 threads.
- **Where the time goes.** Per query of the final run's size (105,879 distances, k + 1 = 20,118 positions; a synthetic micro-benchmark on a login node), 0.36 ms of distances and 0.46 ms of selection on one core. The distance loop is bound by memory bandwidth (8.5 MB of coordinates per query), so its AVX2 build gains only 15%; computing several queries per pass over the points would cut that traffic, a possible later optimisation. The selection brackets the k + 1-th smallest distance with a sorted sample and selects within the bracket only (exact by construction, with a full-selection fallback). These two changes brought the one-thread kNN from 99 s and the core from 120 s (first version, job 59730750) to 67 s and 87 s.
- **Validation runs.** With stage dumps (which sort every cell's final-run neighbours for the S4 digest) and 16 threads, all 24 runs of the 10 objects took 30 s in one job:

| object | run | real cells | with doublets | k | threads | dumps | wall | peak RSS |
|---|---|---|---|---|---|---|---|---|
| soupx_tenx_test | portable-cascadelake | 300 | 400 | 116 | 16 | full | 0.1 s | 17 MiB |
| soupx_tenx_test | portable-sapphirerapids | 300 | 400 | 116 | 16 | full | 0.1 s | 17 MiB |
| soupx_tenx_test | production-sapphirerapids | 300 | 400 | 116 | 16 | full | 0.1 s | 17 MiB |
| soupx_emptydrops_test | portable-cascadelake | 267 | 356 | 100 | 16 | full | 0.1 s | 16 MiB |
| soupx_emptydrops_test | portable-sapphirerapids | 267 | 356 | 100 | 16 | full | 0.1 s | 16 MiB |
| soupx_emptydrops_test | production-sapphirerapids | 267 | 356 | 100 | 16 | full | 0.1 s | 16 MiB |
| decontx_raw_tenx_test | portable-cascadelake | 300 | 400 | 36 | 16 | full | 0.1 s | 17 MiB |
| decontx_raw_tenx_test | portable-sapphirerapids | 300 | 400 | 36 | 16 | full | 0.1 s | 17 MiB |
| decontx_raw_tenx_test | production-sapphirerapids | 300 | 400 | 36 | 16 | full | 0.1 s | 17 MiB |
| decontx_raw_emptydrops_test | portable-cascadelake | 267 | 356 | 78 | 16 | full | 0.1 s | 16 MiB |
| decontx_raw_emptydrops_test | portable-sapphirerapids | 267 | 356 | 78 | 16 | full | 0.0 s | 16 MiB |
| decontx_raw_emptydrops_test | production-sapphirerapids | 267 | 356 | 78 | 16 | full | 0.1 s | 16 MiB |
| cellbender_fromraw_test | portable-cascadelake | 278 | 371 | 85 | 16 | full | 0.0 s | 17 MiB |
| cellbender_fromraw_test | portable-sapphirerapids | 278 | 371 | 85 | 16 | full | 0.0 s | 16 MiB |
| cellbender_fromraw_test | production-sapphirerapids | 278 | 371 | 85 | 16 | full | 0.0 s | 17 MiB |
| soupx_tenx_SRR18256555 | portable-cascadelake | 3,108 | 4,144 | 21 | 16 | light | 0.2 s | 43 MiB |
| soupx_tenx_SRR18256555 | portable-sapphirerapids | 3,108 | 4,144 | 21 | 16 | light | 0.2 s | 43 MiB |
| soupx_tenx_SRR18256555 | production-sapphirerapids | 3,108 | 4,144 | 21 | 16 | light | 0.2 s | 43 MiB |
| soupx_tenx_10xmousebrain | portable-sapphirerapids | 11,831 | 15,775 | 1,104 | 16 | light | 1.4 s | 119 MiB |
| soupx_tenx_10xmousebrain | production-sapphirerapids | 11,831 | 15,775 | 1,104 | 16 | light | 1.4 s | 119 MiB |
| soupx_emptydrops_4500871lane | portable-sapphirerapids | 26,598 | 35,464 | 709 | 16 | light | 1.8 s | 137 MiB |
| soupx_emptydrops_4500871lane | production-sapphirerapids | 26,598 | 35,464 | 709 | 16 | light | 1.8 s | 137 MiB |
| soupx_emptydrops_sscovelli | portable-sapphirerapids | 79,409 | 105,879 | 20,117 | 16 | light | 9.4 s | 195 MiB |
| cellbender_fromraw_ghirsutum2 | portable-sapphirerapids | 59,323 | 79,097 | 2,373 | 16 | light | 4.0 s | 173 MiB |

## Findings for the port

- **R's parser rarely misrounds, but it does.** Of 1.6 million random decimal strings, R's `as.numeric()` gave a value other than the correctly rounded one for 65: 64 in scientific notation, nearly all with decimal exponents beyond 100 in magnitude, where the long-double power of ten built by squaring is inexact; and one ordinary 14-digit number, `23361.110619735`, whose exact long-double quotient is rounded twice (to 64 bits, then to 53). Within [0, 1], none. Every pK and pN label parses to the correctly rounded double. The port reproduces R either way.
- **Values through strings.** The pK grid is `seq()`'s: eight of its values (0.06, 0.07, 0.1, 0.15, 0.18, 0.21, 0.24, 0.29) are not the doubles their labels parse to; 0.29, for instance, is 0.29000000000000004 in the sweep and 0.28999999999999998 after `as.numeric(as.character(...))`. Two objects chose such a pK (soupx_tenx_test 0.29, soupx_tenx_10xmousebrain 0.07). In both, k comes out the same from either value (116 and 1,104), but the port takes the sweep's k from the grid value and the final run's from the parsed one, as R does, and `summary.json` records the parsed pK bit for bit.
- **Position 1.** The selection by (distance, index) puts a lower-index twin before the cell itself, as R's stable `order()` does; cellbender_fromraw_ghirsutum2's digests, which change when positions 1 and 2 swap, are identical for all 59,323 cells.
- **Distances are compared after the square root**, as R compares them: distinct sums of squares can share a square root, and then the lower index comes first.
- **nExp edge cases follow R.** `order(pANN, decreasing = TRUE)[1:nExp]` gives one cell for nExp = 0 (`1:0` is `c(1, 0)`), stops for nExp <= -1, and ignores indices beyond the cells; k = 0 reads positions 2 and 1 (`2:1`). None of these occur in the 10 objects; they matter for other inputs (a one-cluster object has a homotypic proportion of 1 and nExp = 0).
- **`bkde` stops when every pANN of a (pN, pK) pair is identical** (bandwidth 0, then `if (L == 0)` on NA); the port returns that as an error.
- **The FFT's compile-time constants.** GCC folds `cos` and `sin` of fft.c's constant 72-degree angle at compile time (correctly rounded); glibc's run-time values for that argument are the same, so the port calls glibc throughout. Its R vectors cover 56 lengths in both directions: 1 to 32, then 24 more up to 1,155, among them the primes 97 and 401, powers of 3, 5, 7, 13 and 17, and bkde's 512 and 1,024.

## Deviations from the plan and the brief

- **Teacher forcing only.** `scrust doubletfinder` requires `--from-stages` until Rust's SCTransform and PCA exist (M4-M6). Its other options are the plan's (B3): `--pcs`, `--pn`, `--nexp-rate`, `--pk auto|<value>` (a given value is read with R's parser and skips the sweep, like calling `doubletFinder` directly), `--seed` (recorded), `--threads`.
- **Inputs.** Besides the stage directory, the run reads the 10x directory's barcodes and features (for `cells.tsv` and `n_genes`) and the cluster labels with their level order (`seurat_clusters_levels.txt`, as the export wrote it). Reading the gzipped 10x files needed a gzip reader; it was written with the standard library (`scio::gz`) rather than adding a crate. M6 will read `matrix.mtx.gz` with it.
- **S4 on tier P through digests.** The light and standard dumps of the P objects keep R's per-cell digests, not the positions; tier T compares the positions themselves. The digest changes when a position is replaced or two positions swap.
- **`rng_end_sha256`** is not produced under teacher forcing (no random numbers are drawn); `summary.json` has `null` there.
- **More R vectors than `R_strtod`'s.** The same script (`validation/R/rcompat_vectors.R`, no `.rds`, run on a login node) now also writes vectors for `fft`, `dnorm` and `approx`; the 13 M1 families regenerate byte-identically. The 65 misrounded strings are literals in the script, found by a one-off search described there.
- **S1** (sampling) is settled from R's checkpoints, although the M3 Rust run takes the draws from R; the draws themselves are Rust's from M6 on.
- **Earlier code versions.** The first version of the core (commit 9802545) and the version with the faster selection (05a6393) were also validated on all 24 runs, with the same result; their outputs are kept in `data/doubletfinder/<object>/rust/m3-9802545/` and `rust/m3-05a6393/`. The tables above are from commit e16bf4c (`rust/m3/`), whose core is 05a6393's plus the S1 check.
- **Benchmark node.** 32 CPUs of a shared `sapphire` node, not an exclusive one; the one-thread times may include some interference.

## Reproducing

From the repository's top directory, after `source env.sh`:

```bash
cargo build --release -p scrust
sbatch validation/slurm/m3_doubletfinder.sbatch soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2
sbatch validation/slurm/m3_benchmark.sbatch
Rscript validation/R/report_doubletfinder_m3.R <benchmark job id>
```

One run by hand, on a tier-T fixture (seconds, any machine):

```bash
f=tests/fixtures/doubletfinder/soupx_tenx_test
target/release/scrust doubletfinder --counts $f/input --clusters $f/input/seurat_clusters.tsv --from-stages $f/stages --out /path/to/out --dump-stages light
target/release/scrust validate doubletfinder --reference $f/reference --reference-stages $f/stages --candidate /path/to/out
```

The report script reads only plain files; R runs with the hermetic flags of `validation/R/README.md`. SLURM logs of every M3 job are in `data/doubletfinder/logs/`.
