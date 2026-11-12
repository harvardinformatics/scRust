# RunPCA (irlba) in Rust (milestone M5)

M5 ports Seurat 5.5.1's `RunPCA` as DoubletFinder 2.0.6 calls it (docs/PLAN.md B3b and B6): `PrepDR`'s zero-variance filter, `RunPCA.default`, and irlba 2.4.1's augmented implicitly restarted Lanczos bidiagonalisation, with every matrix product the OpenBLAS (or LAPACK) call R makes. This is stage S3 of the DoubletFinder pipeline. It was checked on R's own inputs (*teacher forcing*, plan A5b) and chained after Rust's own SCTransform, against every instrumented reference run of M2. The full `scrust doubletfinder` (M6) comes next.

## Outcome

| Done when (plan B6, and this milestone's brief) | Result |
|---|---|
| S3 bit-identical to R in all 7 contexts of all 24 instrumented runs, each under its own OpenBLAS profile: sdev, irlba's `d`, embeddings, loadings and the RNG checkpoint after RunPCA, teacher-forced and chained | **Yes, all 24 runs of the 10 objects, both ways.** Also identical: the PCA features (PrepDR's filter), irlba's iterations and products, `total.variance`, and R's `.Random.seed` at RunPCA's entry, at irlba's entry (after `set.seed(42)`) and exit, and at RunPCA's exit. Signs included: no component needed aligning. |
| Every irlba branch ported, including the near-invariant-subspace draws | **Yes**, and each checked against R on synthetic inputs: the tiny-problem `svd()`, `nv = 1`, the `work` cap, rank-deficient matrices (up to 31 extra `rnorm()` vectors), tightly clustered spectra (the working dimension grows, once or twice), `maxit` reached. The reference data reach none of the extra-draw branches: in all 168 contexts of the 24 runs, R's irlba drew `rnorm(nfeatures)` and nothing more (`m5_tier_t.rs`, from R's own checkpoints). |
| `crates/pca` and `scrust pca`; tier-T checks in `cargo test`; fixtures in the manifest | **Yes.** `scrust pca` is byte-identical to R's `RunPCA()` on a tier-T object (Seurat's default 50 components, 40 restarts) and on the final runs of the two largest objects (10 components). |
| `tolerances.toml` settles S3 | **Yes: exact** under the run's kernel (S3 and S3-chain). |
| Benchmark: Rust against R's `RunPCA` on the same input, sscovelli and ghirsutum2, 1 and 32 threads; results independent of the thread count | **Yes** (below): outputs byte-identical at 1 and 32 threads and to R's. |
| `cargo test --workspace` passes in a fresh clone | **Yes.** |

Over the 24 runs, 4,368 checks (2,184 teacher-forced, 2,184 chained) compared about 23.1 million values, among them 1,344 RNG checkpoints; none differed. The S3 tables are in [Results](#results), the R vectors in [Reference vectors](#reference-vectors).

## What was built

| Crate | M5 contents |
|---|---|
| `pca` (new) | `run_pca`: `PrepDR` and `RunPCA.default` (`set.seed(seed.use)`, `total.variance`, `npcs <- min(npcs, nrow - 1)`, sdev, loadings, `u %*% diag(d)`); `irlba`: irlba 2.4.1's R code path, statement by statement; `products`: irlba's two large products and the only thread split that keeps R's bits; `rowvar`: Seurat's `RowVar` with Eigen 3.4's summation order. |
| `rcompat` | `blas::matprod` (R's `%*%`: its NaN check, `dgemv` or `dgemm`), `blas::la_svd` (`svd()`, i.e. `La.svd` with LAPACK `dgesdd` and R's workspace query), the OpenBLAS core name and thread count; LAPACK linked by name (`-llapack`, the same OpenBLAS file). |
| `scio` | `npy::read_f64_row_major`: a column-major `.npy` matrix streamed into row-major order with about one copy in memory. |
| `doubletfinder` | `pca`: DoubletFinder's `RunPCA(npcs = 10)` and the stage-S3 comparison (teacher-forced and chained) with an instrumented R run. |
| `scrust` | `scrust pca` (RunPCA on a `scrust sctransform` output directory) and `scrust validate pca`; OpenBLAS is pinned to one thread at start. |

No external crate was added; the workspace uses only the Rust standard library and, through FFI, the OpenBLAS (with its LAPACK) that R uses (`La_library()`: `libopenblasp-r0.3.34.so`).

## How the comparison works

For every instrumented R run of M2, `scrust validate pca` ([`validation/slurm/m5_pca.sbatch`](../../validation/slurm/m5_pca.sbatch)) takes each of the seven contexts:

- **teacher-forced (S3):** RunPCA on R's own PCA input, the `scale.data` rows `RunPCA.default` received (`S2d_scale_data.<ctx>.npy`, dumped by the full and standard runs), with the generator restored from R's `.Random.seed` at RunPCA's entry. The light runs (sscovelli, ghirsutum2) did not dump that matrix; there Rust's S2d on R's own S2c parameters stands in (M4 showed S2d identical to R's in every run; S3 being identical confirms it again);
- **chained (S3-chain):** Rust's whole SCTransform from the counts and R's S1 draws, starting from R's state at SCTransform's entry, then RunPCA on Rust's `scale.data` with the generator SCTransform left: what `scrust doubletfinder` will run in M6.

Each compares, bit for bit, irlba's `d`, the standard deviations, the embeddings and the loadings (`S3_*`), the PCA features (`S2d_pca_features`, PrepDR's filter), the scalars `pca_n_features`, `irlba_iter`, `irlba_mprod` and `total_variance` (`stages.json`, written with `%.17g`), and four `.Random.seed` checkpoints: RunPCA's entry, irlba's entry and exit, RunPCA's exit (the state the next pN's doublet draws continue from). 13 checks per context. The products are OpenBLAS calls, so each run is validated under its own profile: the portable runs with `OPENBLAS_CORETYPE=Haswell`, the production runs on Sapphire Rapids nodes with the kernel OpenBLAS detects (SapphireRapids), checked from the fingerprint before each run. The job ran 8 threads; the results do not depend on it (below).

irlba's rare branches were checked alone on reference vectors that [`validation/R/pca_vectors.R`](../../validation/R/pca_vectors.R) writes from the installed irlba and Seurat, on synthetic inputs Rust rebuilds exactly (R's `rnorm()` from a seed, R's `%*%`, exact diagonals).

## Results

Commit 5e7cc8e (`rust/m5-5e7cc8e/<run>/validate.{json,md}`), job 59988999: `scrust validate pca --mode both`, 8 threads, each run under its own OpenBLAS profile, on Intel Xeon Platinum 8480+ (`sapphire`); the job took 9:56 for all 24 runs. Each stage cell gives the checks identical out of the checks made in that context (13: four RNG checkpoints, the PCA features, four scalars, `d`, sdev, embeddings, loadings) and the values compared. "cells" includes the artificial doublets; "features" are the PCA's (PrepDR kept every variable feature everywhere); "iter" and "mprod" are irlba's restarts and products, R's and Rust's alike. The tier-T objects have fewer than 3,000 genes left after SCTransform's filter, so all of them are variable features.

#### soupx_tenx_test, `instrumented-portable-cascadelake`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 4 | 118 | 13/13, 12,252 | 13/13, 12,252 |
| sweep_pN0.1 | 333 | 836 | 3 | 104 | 13/13, 12,554 | 13/13, 12,554 |
| sweep_pN0.15 | 353 | 849 | 3 | 104 | 13/13, 12,897 | 13/13, 12,897 |
| sweep_pN0.2 | 375 | 858 | 4 | 118 | 13/13, 13,216 | 13/13, 13,216 |
| sweep_pN0.25 | 400 | 873 | 4 | 116 | 13/13, 13,631 | 13/13, 13,631 |
| sweep_pN0.3 | 429 | 896 | 3 | 102 | 13/13, 14,174 | 13/13, 14,174 |
| final | 400 | 879 | 4 | 118 | 13/13, 13,697 | 13/13, 13,697 |

#### soupx_tenx_test, `instrumented-portable-sapphirerapids`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 4 | 118 | 13/13, 12,252 | 13/13, 12,252 |
| sweep_pN0.1 | 333 | 836 | 3 | 104 | 13/13, 12,554 | 13/13, 12,554 |
| sweep_pN0.15 | 353 | 849 | 3 | 104 | 13/13, 12,897 | 13/13, 12,897 |
| sweep_pN0.2 | 375 | 858 | 4 | 118 | 13/13, 13,216 | 13/13, 13,216 |
| sweep_pN0.25 | 400 | 873 | 4 | 116 | 13/13, 13,631 | 13/13, 13,631 |
| sweep_pN0.3 | 429 | 896 | 3 | 102 | 13/13, 14,174 | 13/13, 14,174 |
| final | 400 | 879 | 4 | 118 | 13/13, 13,697 | 13/13, 13,697 |

#### soupx_tenx_test, `instrumented-production-sapphirerapids`

Validation: 0.6 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 4 | 118 | 13/13, 12,252 | 13/13, 12,252 |
| sweep_pN0.1 | 333 | 836 | 3 | 104 | 13/13, 12,554 | 13/13, 12,554 |
| sweep_pN0.15 | 353 | 849 | 3 | 104 | 13/13, 12,897 | 13/13, 12,897 |
| sweep_pN0.2 | 375 | 858 | 4 | 118 | 13/13, 13,216 | 13/13, 13,216 |
| sweep_pN0.25 | 400 | 873 | 4 | 116 | 13/13, 13,631 | 13/13, 13,631 |
| sweep_pN0.3 | 429 | 896 | 3 | 102 | 13/13, 14,174 | 13/13, 14,174 |
| final | 400 | 879 | 4 | 118 | 13/13, 13,697 | 13/13, 13,697 |

#### soupx_emptydrops_test, `instrumented-portable-cascadelake`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 4 | 118 | 13/13, 12,012 | 13/13, 12,012 |
| sweep_pN0.1 | 297 | 831 | 4 | 118 | 13/13, 12,139 | 13/13, 12,139 |
| sweep_pN0.15 | 314 | 843 | 4 | 118 | 13/13, 12,441 | 13/13, 12,441 |
| sweep_pN0.2 | 334 | 855 | 4 | 118 | 13/13, 12,773 | 13/13, 12,773 |
| sweep_pN0.25 | 356 | 872 | 4 | 120 | 13/13, 13,180 | 13/13, 13,180 |
| sweep_pN0.3 | 381 | 880 | 4 | 120 | 13/13, 13,518 | 13/13, 13,518 |
| final | 356 | 879 | 4 | 120 | 13/13, 13,257 | 13/13, 13,257 |

#### soupx_emptydrops_test, `instrumented-portable-sapphirerapids`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 4 | 118 | 13/13, 12,012 | 13/13, 12,012 |
| sweep_pN0.1 | 297 | 831 | 4 | 118 | 13/13, 12,139 | 13/13, 12,139 |
| sweep_pN0.15 | 314 | 843 | 4 | 118 | 13/13, 12,441 | 13/13, 12,441 |
| sweep_pN0.2 | 334 | 855 | 4 | 118 | 13/13, 12,773 | 13/13, 12,773 |
| sweep_pN0.25 | 356 | 872 | 4 | 120 | 13/13, 13,180 | 13/13, 13,180 |
| sweep_pN0.3 | 381 | 880 | 4 | 120 | 13/13, 13,518 | 13/13, 13,518 |
| final | 356 | 879 | 4 | 120 | 13/13, 13,257 | 13/13, 13,257 |

#### soupx_emptydrops_test, `instrumented-production-sapphirerapids`

Validation: 0.9 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 4 | 118 | 13/13, 12,012 | 13/13, 12,012 |
| sweep_pN0.1 | 297 | 831 | 4 | 118 | 13/13, 12,139 | 13/13, 12,139 |
| sweep_pN0.15 | 314 | 843 | 4 | 118 | 13/13, 12,441 | 13/13, 12,441 |
| sweep_pN0.2 | 334 | 855 | 4 | 118 | 13/13, 12,773 | 13/13, 12,773 |
| sweep_pN0.25 | 356 | 872 | 4 | 120 | 13/13, 13,180 | 13/13, 13,180 |
| sweep_pN0.3 | 381 | 880 | 4 | 120 | 13/13, 13,518 | 13/13, 13,518 |
| final | 356 | 879 | 4 | 120 | 13/13, 13,257 | 13/13, 13,257 |

#### decontx_raw_tenx_test, `instrumented-portable-cascadelake`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 4 | 120 | 13/13, 11,911 | 13/13, 11,911 |
| sweep_pN0.1 | 333 | 800 | 4 | 118 | 13/13, 12,158 | 13/13, 12,158 |
| sweep_pN0.15 | 353 | 811 | 4 | 118 | 13/13, 12,479 | 13/13, 12,479 |
| sweep_pN0.2 | 375 | 823 | 3 | 104 | 13/13, 12,831 | 13/13, 12,831 |
| sweep_pN0.25 | 400 | 841 | 3 | 104 | 13/13, 13,279 | 13/13, 13,279 |
| sweep_pN0.3 | 429 | 854 | 3 | 104 | 13/13, 13,712 | 13/13, 13,712 |
| final | 400 | 850 | 4 | 118 | 13/13, 13,378 | 13/13, 13,378 |

#### decontx_raw_tenx_test, `instrumented-portable-sapphirerapids`

Validation: 0.7 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 4 | 120 | 13/13, 11,911 | 13/13, 11,911 |
| sweep_pN0.1 | 333 | 800 | 4 | 118 | 13/13, 12,158 | 13/13, 12,158 |
| sweep_pN0.15 | 353 | 811 | 4 | 118 | 13/13, 12,479 | 13/13, 12,479 |
| sweep_pN0.2 | 375 | 823 | 3 | 104 | 13/13, 12,831 | 13/13, 12,831 |
| sweep_pN0.25 | 400 | 841 | 3 | 104 | 13/13, 13,279 | 13/13, 13,279 |
| sweep_pN0.3 | 429 | 854 | 3 | 104 | 13/13, 13,712 | 13/13, 13,712 |
| final | 400 | 850 | 4 | 118 | 13/13, 13,378 | 13/13, 13,378 |

#### decontx_raw_tenx_test, `instrumented-production-sapphirerapids`

Validation: 0.9 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 4 | 120 | 13/13, 11,911 | 13/13, 11,911 |
| sweep_pN0.1 | 333 | 800 | 4 | 118 | 13/13, 12,158 | 13/13, 12,158 |
| sweep_pN0.15 | 353 | 811 | 4 | 118 | 13/13, 12,479 | 13/13, 12,479 |
| sweep_pN0.2 | 375 | 823 | 3 | 104 | 13/13, 12,831 | 13/13, 12,831 |
| sweep_pN0.25 | 400 | 841 | 3 | 104 | 13/13, 13,279 | 13/13, 13,279 |
| sweep_pN0.3 | 429 | 854 | 3 | 104 | 13/13, 13,712 | 13/13, 13,712 |
| final | 400 | 850 | 4 | 118 | 13/13, 13,378 | 13/13, 13,378 |

#### decontx_raw_emptydrops_test, `instrumented-portable-cascadelake`

Validation: 0.7 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 4 | 118 | 13/13, 11,748 | 13/13, 11,748 |
| sweep_pN0.1 | 297 | 812 | 4 | 118 | 13/13, 11,930 | 13/13, 11,930 |
| sweep_pN0.15 | 314 | 825 | 4 | 118 | 13/13, 12,243 | 13/13, 12,243 |
| sweep_pN0.2 | 334 | 828 | 4 | 120 | 13/13, 12,476 | 13/13, 12,476 |
| sweep_pN0.25 | 356 | 849 | 5 | 132 | 13/13, 12,927 | 13/13, 12,927 |
| sweep_pN0.3 | 381 | 855 | 4 | 118 | 13/13, 13,243 | 13/13, 13,243 |
| final | 356 | 857 | 4 | 118 | 13/13, 13,015 | 13/13, 13,015 |

#### decontx_raw_emptydrops_test, `instrumented-portable-sapphirerapids`

Validation: 0.7 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 4 | 118 | 13/13, 11,748 | 13/13, 11,748 |
| sweep_pN0.1 | 297 | 812 | 4 | 118 | 13/13, 11,930 | 13/13, 11,930 |
| sweep_pN0.15 | 314 | 825 | 4 | 118 | 13/13, 12,243 | 13/13, 12,243 |
| sweep_pN0.2 | 334 | 828 | 4 | 120 | 13/13, 12,476 | 13/13, 12,476 |
| sweep_pN0.25 | 356 | 849 | 5 | 132 | 13/13, 12,927 | 13/13, 12,927 |
| sweep_pN0.3 | 381 | 855 | 4 | 118 | 13/13, 13,243 | 13/13, 13,243 |
| final | 356 | 857 | 4 | 118 | 13/13, 13,015 | 13/13, 13,015 |

#### decontx_raw_emptydrops_test, `instrumented-production-sapphirerapids`

Validation: 0.9 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 4 | 118 | 13/13, 11,748 | 13/13, 11,748 |
| sweep_pN0.1 | 297 | 812 | 4 | 118 | 13/13, 11,930 | 13/13, 11,930 |
| sweep_pN0.15 | 314 | 825 | 4 | 118 | 13/13, 12,243 | 13/13, 12,243 |
| sweep_pN0.2 | 334 | 828 | 4 | 120 | 13/13, 12,476 | 13/13, 12,476 |
| sweep_pN0.25 | 356 | 849 | 5 | 132 | 13/13, 12,927 | 13/13, 12,927 |
| sweep_pN0.3 | 381 | 855 | 4 | 118 | 13/13, 13,243 | 13/13, 13,243 |
| final | 356 | 857 | 4 | 118 | 13/13, 13,015 | 13/13, 13,015 |

#### cellbender_fromraw_test, `instrumented-portable-cascadelake`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 4 | 116 | 13/13, 12,044 | 13/13, 12,044 |
| sweep_pN0.1 | 309 | 840 | 3 | 102 | 13/13, 12,358 | 13/13, 12,358 |
| sweep_pN0.15 | 327 | 851 | 3 | 106 | 13/13, 12,659 | 13/13, 12,659 |
| sweep_pN0.2 | 348 | 873 | 4 | 118 | 13/13, 13,111 | 13/13, 13,111 |
| sweep_pN0.25 | 371 | 892 | 4 | 118 | 13/13, 13,550 | 13/13, 13,550 |
| sweep_pN0.3 | 397 | 898 | 3 | 106 | 13/13, 13,876 | 13/13, 13,876 |
| final | 371 | 885 | 4 | 118 | 13/13, 13,473 | 13/13, 13,473 |

#### cellbender_fromraw_test, `instrumented-portable-sapphirerapids`

Validation: 0.8 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 4 | 116 | 13/13, 12,044 | 13/13, 12,044 |
| sweep_pN0.1 | 309 | 840 | 3 | 102 | 13/13, 12,358 | 13/13, 12,358 |
| sweep_pN0.15 | 327 | 851 | 3 | 106 | 13/13, 12,659 | 13/13, 12,659 |
| sweep_pN0.2 | 348 | 873 | 4 | 118 | 13/13, 13,111 | 13/13, 13,111 |
| sweep_pN0.25 | 371 | 892 | 4 | 118 | 13/13, 13,550 | 13/13, 13,550 |
| sweep_pN0.3 | 397 | 898 | 3 | 106 | 13/13, 13,876 | 13/13, 13,876 |
| final | 371 | 885 | 4 | 118 | 13/13, 13,473 | 13/13, 13,473 |

#### cellbender_fromraw_test, `instrumented-production-sapphirerapids`

Validation: 0.9 s wall, 0.0 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 4 | 116 | 13/13, 12,044 | 13/13, 12,044 |
| sweep_pN0.1 | 309 | 840 | 3 | 102 | 13/13, 12,358 | 13/13, 12,358 |
| sweep_pN0.15 | 327 | 851 | 3 | 106 | 13/13, 12,659 | 13/13, 12,659 |
| sweep_pN0.2 | 348 | 873 | 4 | 118 | 13/13, 13,111 | 13/13, 13,111 |
| sweep_pN0.25 | 371 | 892 | 4 | 118 | 13/13, 13,550 | 13/13, 13,550 |
| sweep_pN0.3 | 397 | 898 | 3 | 106 | 13/13, 13,876 | 13/13, 13,876 |
| final | 371 | 885 | 4 | 118 | 13/13, 13,473 | 13/13, 13,473 |

#### soupx_tenx_SRR18256555, `instrumented-portable-cascadelake`

Validation: 16.6 s wall, 0.6 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 3,000 | 2 | 84 | 13/13, 65,748 | 13/13, 65,748 |
| sweep_pN0.1 | 3,453 | 3,000 | 2 | 84 | 13/13, 67,558 | 13/13, 67,558 |
| sweep_pN0.15 | 3,656 | 3,000 | 2 | 84 | 13/13, 69,588 | 13/13, 69,588 |
| sweep_pN0.2 | 3,885 | 3,000 | 2 | 84 | 13/13, 71,878 | 13/13, 71,878 |
| sweep_pN0.25 | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |
| sweep_pN0.3 | 4,440 | 3,000 | 2 | 84 | 13/13, 77,428 | 13/13, 77,428 |
| final | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |

#### soupx_tenx_SRR18256555, `instrumented-portable-sapphirerapids`

Validation: 16.6 s wall, 0.6 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 3,000 | 2 | 84 | 13/13, 65,748 | 13/13, 65,748 |
| sweep_pN0.1 | 3,453 | 3,000 | 2 | 84 | 13/13, 67,558 | 13/13, 67,558 |
| sweep_pN0.15 | 3,656 | 3,000 | 2 | 84 | 13/13, 69,588 | 13/13, 69,588 |
| sweep_pN0.2 | 3,885 | 3,000 | 2 | 84 | 13/13, 71,878 | 13/13, 71,878 |
| sweep_pN0.25 | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |
| sweep_pN0.3 | 4,440 | 3,000 | 2 | 84 | 13/13, 77,428 | 13/13, 77,428 |
| final | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |

#### soupx_tenx_SRR18256555, `instrumented-production-sapphirerapids`

Validation: 16.3 s wall, 0.6 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 3,000 | 2 | 84 | 13/13, 65,748 | 13/13, 65,748 |
| sweep_pN0.1 | 3,453 | 3,000 | 2 | 84 | 13/13, 67,558 | 13/13, 67,558 |
| sweep_pN0.15 | 3,656 | 3,000 | 2 | 84 | 13/13, 69,588 | 13/13, 69,588 |
| sweep_pN0.2 | 3,885 | 3,000 | 2 | 84 | 13/13, 71,878 | 13/13, 71,878 |
| sweep_pN0.25 | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |
| sweep_pN0.3 | 4,440 | 3,000 | 2 | 84 | 13/13, 77,428 | 13/13, 77,428 |
| final | 4,144 | 3,000 | 2 | 84 | 13/13, 74,468 | 13/13, 74,468 |

#### soupx_tenx_10xmousebrain, `instrumented-portable-sapphirerapids`

Validation: 62.8 s wall, 2.8 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 98 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 100 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 98 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 98 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 98 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 2 | 84 | 13/13, 175,888 | 13/13, 175,888 |
| final | 15,775 | 3,000 | 3 | 102 | 13/13, 190,778 | 13/13, 190,778 |

#### soupx_tenx_10xmousebrain, `instrumented-production-sapphirerapids`

Validation: 61.9 s wall, 2.7 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 98 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 100 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 98 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 98 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 98 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 2 | 84 | 13/13, 175,888 | 13/13, 175,888 |
| final | 15,775 | 3,000 | 3 | 102 | 13/13, 190,778 | 13/13, 190,778 |

#### soupx_emptydrops_4500871lane, `instrumented-portable-sapphirerapids`

Validation: 61.7 s wall, 3.7 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 102 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 106 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 104 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 102 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 102 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 3 | 102 | 13/13, 175,888 | 13/13, 175,888 |
| final | 35,464 | 3,000 | 3 | 102 | 13/13, 387,668 | 13/13, 387,668 |

#### soupx_emptydrops_4500871lane, `instrumented-production-sapphirerapids`

Validation: 60.6 s wall, 3.7 GiB peak RSS (8 threads); teacher-forced on R's PCA input.

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 102 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 106 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 104 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 102 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 102 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 3 | 102 | 13/13, 175,888 | 13/13, 175,888 |
| final | 35,464 | 3,000 | 3 | 102 | 13/13, 387,668 | 13/13, 387,668 |

#### soupx_emptydrops_sscovelli, `instrumented-portable-sapphirerapids`

Validation: 107.8 s wall, 7.9 GiB peak RSS (8 threads); teacher-forced on Rust's S2d (R dumped no PCA input).

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 104 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 106 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 100 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 102 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 100 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 3 | 100 | 13/13, 175,888 | 13/13, 175,888 |
| final | 105,879 | 3,000 | 2 | 84 | 13/13, 1,091,818 | 13/13, 1,091,818 |

#### cellbender_fromraw_ghirsutum2, `instrumented-portable-sapphirerapids`

Validation: 2:47.02 wall, 7.0 GiB peak RSS (8 threads); teacher-forced on Rust's S2d (R dumped no PCA input).

| context | cells | features | iter | mprod | S3 | S3-chain |
|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 3,000 | 3 | 98 | 13/13, 138,288 | 13/13, 138,288 |
| sweep_pN0.1 | 11,111 | 3,000 | 3 | 102 | 13/13, 144,138 | 13/13, 144,138 |
| sweep_pN0.15 | 11,765 | 3,000 | 3 | 98 | 13/13, 150,678 | 13/13, 150,678 |
| sweep_pN0.2 | 12,500 | 3,000 | 3 | 100 | 13/13, 158,028 | 13/13, 158,028 |
| sweep_pN0.25 | 13,333 | 3,000 | 3 | 102 | 13/13, 166,358 | 13/13, 166,358 |
| sweep_pN0.3 | 14,286 | 3,000 | 3 | 104 | 13/13, 175,888 | 13/13, 175,888 |
| final | 79,097 | 3,000 | 3 | 102 | 13/13, 823,998 | 13/13, 823,998 |

### Totals over the 24 runs

| stage | checks identical | elements compared | elements differing | RNG checkpoints identical |
|---|---|---|---|---|
| S3 | 2,184 of 2,184 | 11,526,810 | 0 | 672 of 672 |
| S3-chain | 2,184 of 2,184 | 11,526,810 | 0 | 672 of 672 |

All 24 runs pass: every check of every context is identical, teacher-forced and chained.

## Reference vectors

| family | cases | what they cover |
|---|---|---|
| `irlba` | 19 | `irlba(A, nv)` after `set.seed()`: `d`, `u`, `v`, `iter`, `mprod`, convergence and `.Random.seed` afterwards, at 1 and 3 threads. Random 60 x 40, 30 x 70, 300 x 120 (Seurat-like); `nv = 1` (R drops `Bsvd$v[, 1:k]` to a vector: the restart's product is `dgemv`); the tiny-problem path (`min(m, n) < 6`: `svd()`, no draw); the `work` cap (`work >= min(m, n)`) and a small `work`; rank-1, 3 and 6 matrices (the near-invariant-subspace branches: 28 to 31 extra `rnorm()` vectors on both bases); exactly repeated singular values (17 extra vectors) and near-repeats; four tight clusters of the largest singular values, gaps 1e-9 to 1e-11 (the clustered-spectrum test expands `work` by 5, once or twice); `maxit = 1` and `2` (not converged) |
| `rowvar` | 17 | Seurat's `RowVar` (and R's `sum()` of it): rows of every length from 1 to 9 (the SSE2 packet loop and its tails), long rows with an offset of 1e8 that makes the summation order visible, constant rows (0.1 x 316 has a positive variance in Eigen's order, 2^-10 x 316 a zero one) |
| `cli_soupx_tenx_test` | 4 files | `scrust pca`'s outputs from Seurat's `RunPCA()` with its defaults (50 components) on R's SCTransform of a tier-T object, through DoubletFinder's path ([`runpca_reference.R`](../../validation/R/runpca_reference.R)): `stdev.tsv`, and the SHA-256 of `embeddings.tsv`, `loadings.tsv` and `random_seed.txt` |

## Benchmark

RunPCA (`npcs = 10`, DoubletFinder's) on the final run of the two largest objects, Rust and R on the same input: R's own SCTransform output of that run from the M4 benchmark (`rust/m4-benchmark/final/r/out/`: `scale_data.npy`, `variable_features.tsv`, `cell_attr.tsv`), read by `scrust pca --sct` and by [`validation/R/runpca_reference.R`](../../validation/R/runpca_reference.R), which builds the SCT assay and calls `RunPCA()` as DoubletFinder does. Portable profile (OPENBLAS_CORETYPE=Haswell, the instrumented runs'), Intel Xeon Platinum 8480CL (`sapphire`, shared nodes), [`validation/slurm/m5_benchmark.sbatch`](../../validation/slurm/m5_benchmark.sbatch): Rust jobs 59989006 and 59989015 (32 CPUs allocated, binary `scrust-m5-5e7cc8e`), R jobs 59989010 and 59989017 (one CPU). *RunPCA* is R's `RunPCA()` call alone (`proc.time`) and Rust's `pca` phase (`summary.json`: PrepDR, irlba, the embeddings). Process wall, CPU time and peak RSS come from `/usr/bin/time -v`; they include reading the 1.9-2.5 GB `scale_data.npy` and building the input (R also builds a Seurat object), and writing the outputs.

| object (final run) | cells x features | side | threads | irlba iter, products | RunPCA (s) | process wall (s) | process CPU (s) | peak RSS (GiB) | RSS after loading (GiB) |
|---|---|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | 105,879 x 3,000 | R | 1 | 2, 84 | 37.8 | 53.4 | 52.1 | 8.6 | 5.2 |
| soupx_emptydrops_sscovelli | 105,879 x 3,000 | Rust | 1 | 2, 84 | 23.8 | 26.9 | 25.9 | 2.4 | 2.4 |
| soupx_emptydrops_sscovelli | 105,879 x 3,000 | Rust | 32 | 2, 84 | 3.6 | 6.2 | 54.9 | 2.5 | 2.4 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 3,000 | R | 1 | 3, 102 | 32.0 | 44.2 | 43.9 | 6.8 | 4.0 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 3,000 | Rust | 1 | 3, 102 | 22.2 | 24.6 | 23.7 | 1.8 | 1.8 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 3,000 | Rust | 32 | 3, 102 | 3.0 | 4.7 | 39.9 | 1.9 | 1.8 |

**Every output was byte-identical:** `embeddings.tsv`, `loadings.tsv`, `stdev.tsv` and `random_seed.txt` of R and of Rust at 1 and at 32 threads, for both objects.

- **Time.** On one thread Rust's RunPCA takes 0.63 (sscovelli) and 0.69 (ghirsutum2) of R's: both spend nearly all of it in the same `dgemm` calls, which exactness requires (below), and Rust saves R's copies (`data.use[features, ]` twice, `t(object)`) and its second `RowVar` pass. On 32 threads it is 10.5 and 10.8 times faster than R.
- **Memory.** Rust holds the matrix once (2.5 and 1.9 GB as doubles) plus irlba's bases (`W`: 22 MB for sscovelli): 2.4-2.5 and 1.8-1.9 GiB at either thread count, 3.4-3.7 times less than R. R peaks at 8.6 and 6.8 GiB here, 3.4 and 2.7 GiB above its input (the matrix and the Seurat object built from it); in DoubletFinder its RunPCA adds those copies to the Seurat object it already holds.
- **Where the time goes.** irlba converged in 2 restarts (84 products) on sscovelli and 3 (102) on ghirsutum2. Each product reads the whole matrix through OpenBLAS's level-3 driver (it packs blocks of A even for a one-column product), about 9 GB/s on one core. On 32 threads the two products run in 32 pieces each; the reorthogonalisation and the restarts stay single calls, as R makes them.

## Threads: splitting irlba's products

irlba's two large products, `A %*% v` (`dgemm("N", "N", m, 1, n)`) and `t(w) %*% A` (`dgemm("T", "N", 1, n, m)`), dominate the time. The brief allows a split across threads only if the result is R's single call bit for bit, at every thread count. OpenBLAS 0.3.34's source ([`third_party/`](../../third_party/README.md); the conda-forge build, `DYNAMIC_ARCH=1`) settles where that holds ([`crates/pca/src/products.rs`](../../crates/pca/src/products.rs)):

- **No forwarding to `gemv`** on x86-64 (`GEMM_GEMV_FORWARD` is set for ARM64, POWER and RISC-V only): both products go through `dgemm`.
- **Small-matrix kernels.** On the SkylakeX kernels (also Cooperlake's and SapphireRapids'), a call with `M N K <= 10^6` (for `"T", "N"` also `M N <= 1200` and `K >= 32`) goes to separate small-matrix kernels, which round differently; Haswell has none. Every other call takes the level-3 driver.
- **The level-3 driver** blocks K by `GEMM_Q` and `GEMM_UNROLL_M` only, so an element's sum is the same chain of K blocks however M and N are cut; the micro-kernel's tile decides the order within a block. Haswell's 4-row tiles of a one-column product keep four partial sums (`KERNEL4x1`, k modulo 4), its 2- and 1-row tails one; every one-row tile keeps one sum per column. SkylakeX's one-column tiles all keep one sum per row; its 12-column tile of a one-row product keeps four, but only the driver's SkylakeX branch hands it 12 columns, and a `DYNAMIC_ARCH` build compiles the driver once for the base target (only `kernel/` per core), so the driver passes 6 columns at a time and that tile is never reached.

So the pieces start at multiples of 16 rows (Haswell needs 4) and 12 columns (a margin; any column cut keeps the bits), the last piece holding the call's own tail; on the SkylakeX family neither the call nor any piece may qualify for the small-matrix kernels; on any other kernel nothing is split. Each piece is a `dgemm` call with the same K and leading dimensions on a block of `A`. Tests ([`products.rs`](../../crates/pca/src/products.rs)) check, under Haswell and under SapphireRapids, that every plan from 2 to 13 pieces gives the single call's bits on shapes that cover every tile and tail, and that cuts breaking a rule do change bits: a row cut at 1,002 changes 4 rows under Haswell (the two moved to a 2-row tail and the two at the end), and pieces small enough for the small-matrix kernels change 5 to 225 elements under SapphireRapids. Products below 2^20 elements of `A` are not split (for speed only), which leaves out the tier-T objects and the R vectors; both large products of every P context were split in 8 pieces in the validation and in 32 in the benchmark, and gave R's bits.

## Findings for the port

- **`W[, 1:j]` drops to a vector.** In the first Lanczos step of the first restart, irlba reorthogonalises `W[, 2]` against `W[, 1:1]`, which R drops to a plain vector. irlba's `cross()` then takes its generic branch, `t(Y) %*% X`: R's `matprod` with a 1 x m matrix, `dgemv("N", 1, m)`, instead of `direct_dgemm_c`'s `dgemm("T", "N", 1, 1, m)`; and `X %*% doty` becomes `dgemv("N", m, 1)`. The kernels round these differently, so the port follows R's operands' shapes, not only its functions.
- **Two ways to multiply.** irlba's own `gemm()` (`direct_dgemm_c`) checks nothing; R's `%*%` (`orthog`'s `X %*% doty`, the restarts, the final `U` and `V`, the embeddings' `u %*% diag(d)`) first runs R's imprecise NaN/Inf test and falls back to its own double loop, then picks `dgemv` or `dgemm` by shape. The port does both, as R does.
- **`svd(B)`** is `La.svd` with `jobu = "S"` (B is square): `dgesdd` on a copy, after a workspace query whose size selects `dgesdd`'s code path; the port makes the same two calls.
- **The clustered-spectrum test** compares the gap between Ritz values with `sqrt(eps)`, and the change of the `k`-th one with `svtol`; `lastsv` is updated just before, so that change is always 0 and only the gap decides (as the plan noted). Two consecutive clustered restarts add 5 to `work`.
- **`invs > 3` never holds**: irlba resets `invs` at every Lanczos step and adds at most 2 in one. The expansion it guards is dead code, ported as written.
- **PrepDR's filter uses Seurat's `RowVar`**, an Eigen loop with SSE2 packets (Seurat was built with `-march=nocona` against RcppEigen 0.3.4.0.2, Eigen 3.4.0): its sums run in two packet accumulators, then the odd element. The order matters at the margin: a constant row of 0.1 over 316 cells has a positive variance (2.3e-32), since its mean rounds; the port follows Eigen's order (checked on R vectors). In the data no feature was dropped.
- **`RunPCA` with one component** stops in R when `weight.by.var = TRUE` (`diag()` of a single `d` is an identity matrix of size `d`), except when `1 <= d < 2`; the port refuses one component instead.
- **OpenBLAS's threads.** R's reference runs used one OpenBLAS thread; with more, OpenBLAS splits large products itself. `scrust` now sets one OpenBLAS thread at start (`openblas_set_num_threads(1)`, a no-op under `env.sh`).

## Plan corrections

- **B3b, "every product is OpenBLAS `dgemm` through `direct_dgemm_c`; its `orthog` step also uses R's `%*%`":** `direct_dgemm_c` makes the two large products, `norm2` and `orthog`'s coefficients; R's `%*%` makes `orthog`'s `X %*% doty` (`dgemv`, after R's NaN test), the restarts and the final vectors (`dgemm`, or `dgemv` for one column), and, in the first Lanczos step, `cross()` itself (`dgemv` with one row), because `W[, 1:1]` drops to a vector.
- **B3b, `work = 26`:** irlba caps `work` at `min(m, n) - 2` when `nv + 16` reaches `min(m, n)`; DoubletFinder's inputs (at least 793 features and 281 cells) always have 26.
- **B5, S3 "embeddings (up to sign per component)":** under the matching kernel the embeddings are identical, signs included.
- **B3b, the zero-variance filter** is Seurat's `RowVar` with Eigen's rounding, not R's `var()`; "zero variance" means a variance that rounds to zero in that order.
- **Plan A5 and A5b say OpenBLAS runs one thread and work is parallel over independent units.** M5 adds a second, narrow mechanism: one BLAS call split into BLAS calls on blocks, only on kernels whose tiles were read and tested to give the same bits, and only along the dimension the tiles allow.

## Deviations

- **Only RunPCA's path for DoubletFinder.** Not ported: irlba's options Seurat does not pass (`center`, `scale`, `shift`, restarting from a previous `v`, `smallest`, `right_only`, complex matrices) and RunPCA's `rev.pca = TRUE`, `approx = FALSE` (`prcomp`) and sparse or BPCells inputs; one component with `weight.by.var` (above).
- **The products are split only on verified kernels** (Haswell, SkylakeX, Cooperlake, SapphireRapids); on any other kernel they run as single calls, so the results stay R's there too but the time does not scale.
- **The light runs' teacher-forced input** is Rust's S2d (R did not dump it); the chained run covers them too, and the benchmark compares the final runs with R's own RunPCA.
- **The benchmark runs the portable profile**, that of the instrumented runs and of the M4 benchmark that made its input; the production profile's split (SapphireRapids) is exercised by the validation of the production runs and by the tests.
- **Validation label.** All tables come from `rust/m5-5e7cc8e/` (job 59988999, binary `data/doubletfinder/logs/bin/scrust-m5-5e7cc8e`, sha256 `5365860b74ea170d`), built from commit 5e7cc8e; the commits after it change tests and documentation only.

## Reproducing

From the repository's top directory, after `source env.sh`:

```bash
cargo build --release -p scrust
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m5-$(git rev-parse --short HEAD)
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m5-$(git rev-parse --short HEAD) sbatch --export=ALL validation/slurm/m5_pca.sbatch m5-$(git rev-parse --short HEAD) soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2
SCRUST_BIN=... sbatch --export=ALL --cpus-per-task=32 --mem=32G validation/slurm/m5_benchmark.sbatch rust soupx_emptydrops_sscovelli
sbatch --cpus-per-task=1 --mem=48G validation/slurm/m5_benchmark.sbatch r soupx_emptydrops_sscovelli
Rscript validation/R/report_doubletfinder_m5.R m5-$(git rev-parse --short HEAD)
```

The R reference vectors and the CLI fixture (no `.rds`; a login node is enough), with the production R env and the hermetic flags:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
R="env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 $R_ENV/bin/Rscript"
$R validation/R/pca_vectors.R tests/fixtures/pca
$R validation/R/sctransform_cli_reference.R tests/fixtures/doubletfinder/soupx_tenx_test/input <sct dir>
$R validation/R/runpca_reference.R <sct dir> <out dir>
```

One context by hand, on a tier-T fixture (seconds):

```bash
target/release/scrust validate pca --input tests/fixtures/doubletfinder/soupx_tenx_test/input --stages tests/fixtures/doubletfinder/soupx_tenx_test/stages --mode chained --contexts final
target/release/scrust sctransform --counts tests/fixtures/doubletfinder/soupx_tenx_test/input --out /tmp/sct
target/release/scrust pca --sct /tmp/sct --out /tmp/pca --npcs 10
```
