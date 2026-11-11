# SCTransform v2 in Rust (milestone M4)

M4 ports SCTransform v2 as DoubletFinder 2.0.6 calls it (docs/PLAN.md B3a and B6): Seurat 5.5.1's `SCTransform()` with its defaults, that is sctransform 0.4.3's `vst(vst.flavor = "v2")` with glmGamPoi 1.22.0's fits, the Pearson residuals, the 3,000 variable features and `scale.data`. This is stage S2 of the DoubletFinder pipeline, in four sub-stages: S2a setup, S2b fits, S2c regularisation, S2d residuals. Each sub-stage was checked on R's own inputs (*teacher forcing*, plan A5b) and as one chain on Rust's own results, against every instrumented reference run of M2. The PCA (M5) and the end-to-end tool (M6) come next.

## Outcome

| Done when (plan B6, and this milestone's brief) | Result |
|---|---|
| S2a-S2d exact everywhere (S2d where the BLAS kernel matches the run's profile) | **Yes, in all 7 contexts (six sweep pN values and the final run) of all 24 instrumented runs of the 10 objects**, each run under its own OpenBLAS profile. Every array R dumped is bit-identical: cell attributes, the gene filter and gene statistics, the fitting cells and genes and the density-weighted sampling (S2a); every intermediate of every gene's glmGamPoi fit and the model parameters (S2b); the regularised parameters of every kept gene, the outliers, the Poisson and smoothing genes and the bandwidth (S2c); the residual mean and variance of every kept gene, the variable features in order, and the full residual matrix and `scale.data` wherever R dumped them (S2d). |
| The chained run exact | **Yes, in all 24 runs:** the whole `SCTransform()` from the counts and R's S1 draws, each sub-stage on the previous one's Rust results, matches every S2 array and scalar, including the states of R's random-number generator at `vst`'s and `SCTransform`'s exits. |
| The RNG checkpoints match | **Yes:** `.Random.seed` at `vst`'s entry, after the cell sampling, at the fits' entry, at each bin of the fits, at the fits' exit, and at `vst`'s and `SCTransform`'s exits, in every context of every run. |
| `nlminb` (PORT) bit-exact; the decision point | **`nlminb` is bit-exact**, with and without the Hessian, on 42 R test cases and in every fit of the data. The decision point's fallback was not needed. |
| The new rcompat pieces exact on R vectors | **Yes:** `nlminb`, `density`/`bw.nrd`, `quantile`/`median`/`mad`/`IQR`, glibc's and R's gamma functions, matrixStats (`rowMeans2`, `colMeans2`, `rowVars`, `weightedMedian`), `dnbinom` (mean form), `optimize`, `uniroot`, `bw.SJ`, `ksmooth`, `.bincode`; and a std-only OpenBLAS binding (`dgemm`, `dgemv`, `ddot`, `dnrm2`). |
| `scrust sctransform` | **Yes**, and its outputs are byte-identical to R's `SCTransform()` on the same input: on a tier-T object (a `cargo test`) and on the final runs of the two largest objects (the benchmark). |
| Benchmark measured | **Yes** (below). |
| `cargo test` passes | **Yes.** |

Over the 24 runs, 19,244 checks compared about 3.04 billion values; none differed. The S2 tables are in [Results](#results); the R vectors in [Reference vectors](#reference-vectors).

## What was built

| Crate | M4 contents |
|---|---|
| `sctransform` (new) | `SctParams` (Seurat's and vst's defaults); the four sub-stages as separate functions, each taking the previous stage's results or R's: `setup` (S2a), `fit` (S2b; `glm`: glmGamPoi's `glm_gp` for the intercept-only model), `regularize` (S2c), `residuals` (S2d); `sctransform()`, the whole chain from a counts matrix and a seed. |
| `rcompat` | `port`: R's `nlminb`, the PORT routines DRMNHB (with Hessian) and DRMNGB (gradient only) with bounds, reverse communication as R drives them; `gamma`: glibc's `lgamma`, R's `lgammafn`, `digamma`, `trigamma`; `matrixstats`: matrixStats 1.5.0's `rowMeans2`, `colMeans2`, `rowVars`, `weightedMedian` (both interpolation modes) and R's `R_qsort_I`; `nbinom`: `dnbinom(mu =, log = TRUE)` with `dbinom_raw`, `stirlerr`, `bd0`, `lgamma1p`; `optimize`: `Brent_fmin`; `zeroin`: `uniroot` over `R_zeroin2`; `bandwidth`: `bw.SJ`; `smooth`: `ksmooth` and `.bincode`; `density`: `density(bw = "nrd")` and `bw.nrd`; `quantile`: type-7 quantiles, `median`, `mad`, `IQR`; `blas`: a std-only FFI binding of OpenBLAS's `dgemm`, `dgemv`, `ddot` and `dnrm2`, and R's `tcrossprod`. |
| `scio` | CSC and CSR sparse matrices; a Matrix Market reader that keeps explicit zeros. |
| `doubletfinder` | DoubletFinder's SCT input (the real cells and their doublets, with explicit zeros where Matrix keeps them); the S2 harness (teacher-forced checks of S2a-S2d, the chained run, the RNG checkpoints and stage scalars). |
| `scrust` | `scrust sctransform` (a 10x directory, or with `--from-stages` and `--context` one context of a DoubletFinder run); `scrust validate sctransform` (`--mode` both, teacher or chained). |

No external crate was added; the workspace uses only the Rust standard library, plus the system OpenBLAS that R itself uses, through FFI.

## How the comparison works

For every instrumented R run of M2, `scrust validate sctransform` rebuilds the SCT input of each of the seven contexts from the object's counts (the export, `input/`) and R's S1 draws (the sweep subsample, when there is one, and the doublet pairs), then:

- **teacher-forced checks:** S2a from the SCT input and R's `.Random.seed` at `SCTransform`'s entry; S2b from R's S2a dumps (fitting genes and cells, cell attributes, gene moments) and R's state at the fits' entry; S2c from R's S2b `model_pars` and S2a statistics; S2d from R's S2c `model_pars_fit` and the cells' `log_umi`. A difference would show up in the sub-stage that caused it;
- **the chained run:** `sctransform::sctransform()` from the SCT input alone, each stage on Rust's own results, compared with all the same dumps, scalars and checkpoints, plus the RNG states at `vst`'s and `SCTransform`'s exits.

Every comparison is bitwise: doubles by their bit patterns, indices and counts by value. Stage scalars (`stages.json`) were written with `%.17g` and are compared exactly. The residuals' products are OpenBLAS calls, so each run is validated under its own profile: the portable runs with `OPENBLAS_CORETYPE=Haswell`, the production runs on Sapphire Rapids nodes with the kernel OpenBLAS detects (SapphireRapids), checked from the fingerprint before each run ([`validation/slurm/m4_sctransform.sbatch`](../../validation/slurm/m4_sctransform.sbatch)). The "full" runs (tier T) dumped the whole residual matrix and `scale.data` of every context; the "standard" runs dumped `scale.data`; the "light" runs (sscovelli, ghirsutum2) only the residual statistics and the features, whose `scale.data` the benchmark compares instead.

The rcompat pieces were checked alone first, on reference vectors that `validation/R/rcompat_vectors.R` writes from R itself (base R, stats, matrixStats). Some sctransform and glmGamPoi functions were also checked directly against the installed packages on synthetic inputs, through `validation/R/sctransform_vectors.R`: the SCT input, glmGamPoi's `fitBeta_one_group` with its fallback, and sctransform's `is_outlier`.

## Results

Commit 7a273be (`rust/m4d-7a273be/<run>/validate.{json,md}`), job 59770277: `scrust validate sctransform --mode both`, 8 threads, each run under its own OpenBLAS profile, on Intel Xeon Platinum 8480+ (`sapphire`). Each stage cell gives the checks identical out of the checks made, and the values compared; a check is one array, scalar or RNG state of one context. The S2-chain column repeats every S2a-S2d check on the chained run's results and adds the mean cell sum and the RNG states at `vst`'s and `SCTransform`'s exits. "Fitting genes" is the number of genes in the fits (2,000 when density sampling applies), "theta Inf" the fitted thetas set to infinity, "outliers" the fitting genes left out of the smoothing (outliers in a parameter, or infinite theta). The cells include the artificial doublets.

#### soupx_tenx_test, `instrumented-portable-cascadelake`

Validation: 0:02.03 wall, 0.1 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 707 | 0 | 117 | 13 | 694 | 15/15, 8,016 | 18/18, 10,608 | 11/11, 3,995 | 11/11, 524,892 | 57/57, 547,513 |
| sweep_pN0.1 | 333 | 836 | 700 | 0 | 136 | 13 | 687 | 15/15, 8,118 | 18/18, 10,503 | 11/11, 4,036 | 11/11, 560,960 | 57/57, 583,619 |
| sweep_pN0.15 | 353 | 849 | 669 | 0 | 180 | 8 | 661 | 15/15, 8,163 | 18/18, 10,038 | 11/11, 4,062 | 11/11, 603,643 | 57/57, 625,908 |
| sweep_pN0.2 | 375 | 858 | 639 | 0 | 219 | 6 | 633 | 15/15, 8,193 | 18/18, 9,588 | 11/11, 4,070 | 11/11, 647,794 | 57/57, 669,647 |
| sweep_pN0.25 | 400 | 873 | 602 | 0 | 271 | 9 | 593 | 15/15, 8,247 | 18/18, 9,033 | 11/11, 4,090 | 11/11, 702,769 | 57/57, 724,141 |
| sweep_pN0.3 | 429 | 896 | 573 | 0 | 323 | 5 | 568 | 15/15, 8,385 | 18/18, 8,598 | 11/11, 4,157 | 11/11, 773,252 | 57/57, 794,394 |
| final | 400 | 879 | 607 | 0 | 272 | 5 | 602 | 15/15, 8,298 | 18/18, 9,108 | 11/11, 4,123 | 11/11, 707,599 | 57/57, 729,130 |

#### soupx_tenx_test, `instrumented-portable-sapphirerapids`

Validation: 0:02.07 wall, 0.1 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 707 | 0 | 117 | 13 | 694 | 15/15, 8,016 | 18/18, 10,608 | 11/11, 3,995 | 11/11, 524,892 | 57/57, 547,513 |
| sweep_pN0.1 | 333 | 836 | 700 | 0 | 136 | 13 | 687 | 15/15, 8,118 | 18/18, 10,503 | 11/11, 4,036 | 11/11, 560,960 | 57/57, 583,619 |
| sweep_pN0.15 | 353 | 849 | 669 | 0 | 180 | 8 | 661 | 15/15, 8,163 | 18/18, 10,038 | 11/11, 4,062 | 11/11, 603,643 | 57/57, 625,908 |
| sweep_pN0.2 | 375 | 858 | 639 | 0 | 219 | 6 | 633 | 15/15, 8,193 | 18/18, 9,588 | 11/11, 4,070 | 11/11, 647,794 | 57/57, 669,647 |
| sweep_pN0.25 | 400 | 873 | 602 | 0 | 271 | 9 | 593 | 15/15, 8,247 | 18/18, 9,033 | 11/11, 4,090 | 11/11, 702,769 | 57/57, 724,141 |
| sweep_pN0.3 | 429 | 896 | 573 | 0 | 323 | 5 | 568 | 15/15, 8,385 | 18/18, 8,598 | 11/11, 4,157 | 11/11, 773,252 | 57/57, 794,394 |
| final | 400 | 879 | 607 | 0 | 272 | 5 | 602 | 15/15, 8,298 | 18/18, 9,108 | 11/11, 4,123 | 11/11, 707,599 | 57/57, 729,130 |

#### soupx_tenx_test, `instrumented-production-sapphirerapids`

Validation: 0:02.05 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 824 | 707 | 0 | 117 | 13 | 694 | 15/15, 8,016 | 18/18, 10,608 | 11/11, 3,995 | 11/11, 524,892 | 57/57, 547,513 |
| sweep_pN0.1 | 333 | 836 | 700 | 0 | 136 | 13 | 687 | 15/15, 8,118 | 18/18, 10,503 | 11/11, 4,036 | 11/11, 560,960 | 57/57, 583,619 |
| sweep_pN0.15 | 353 | 849 | 669 | 0 | 180 | 8 | 661 | 15/15, 8,163 | 18/18, 10,038 | 11/11, 4,062 | 11/11, 603,643 | 57/57, 625,908 |
| sweep_pN0.2 | 375 | 858 | 639 | 0 | 219 | 6 | 633 | 15/15, 8,193 | 18/18, 9,588 | 11/11, 4,070 | 11/11, 647,794 | 57/57, 669,647 |
| sweep_pN0.25 | 400 | 873 | 602 | 0 | 271 | 9 | 593 | 15/15, 8,247 | 18/18, 9,033 | 11/11, 4,090 | 11/11, 702,769 | 57/57, 724,141 |
| sweep_pN0.3 | 429 | 896 | 573 | 0 | 323 | 5 | 568 | 15/15, 8,385 | 18/18, 8,598 | 11/11, 4,157 | 11/11, 773,252 | 57/57, 794,394 |
| final | 400 | 879 | 607 | 0 | 272 | 5 | 602 | 15/15, 8,298 | 18/18, 9,108 | 11/11, 4,123 | 11/11, 707,599 | 57/57, 729,130 |

#### soupx_emptydrops_test, `instrumented-portable-cascadelake`

Validation: 0:01.96 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 686 | 0 | 148 | 11 | 675 | 15/15, 7,908 | 18/18, 10,293 | 11/11, 4,016 | 11/11, 472,882 | 57/57, 495,101 |
| sweep_pN0.1 | 297 | 831 | 659 | 0 | 172 | 11 | 648 | 15/15, 7,857 | 18/18, 9,888 | 11/11, 3,977 | 11/11, 497,773 | 57/57, 519,497 |
| sweep_pN0.15 | 314 | 843 | 633 | 0 | 210 | 11 | 622 | 15/15, 7,902 | 18/18, 9,498 | 11/11, 3,999 | 11/11, 533,623 | 57/57, 555,024 |
| sweep_pN0.2 | 334 | 855 | 605 | 0 | 250 | 10 | 595 | 15/15, 7,950 | 18/18, 9,078 | 11/11, 4,020 | 11/11, 575,419 | 57/57, 596,469 |
| sweep_pN0.25 | 356 | 872 | 565 | 0 | 307 | 9 | 556 | 15/15, 7,998 | 18/18, 8,478 | 11/11, 4,049 | 11/11, 625,228 | 57/57, 645,755 |
| sweep_pN0.3 | 381 | 880 | 543 | 0 | 337 | 6 | 537 | 15/15, 8,055 | 18/18, 8,148 | 11/11, 4,062 | 11/11, 674,964 | 57/57, 695,231 |
| final | 356 | 879 | 570 | 0 | 309 | 9 | 561 | 15/15, 8,055 | 18/18, 8,553 | 11/11, 4,082 | 11/11, 630,247 | 57/57, 650,939 |

#### soupx_emptydrops_test, `instrumented-portable-sapphirerapids`

Validation: 0:01.88 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 686 | 0 | 148 | 11 | 675 | 15/15, 7,908 | 18/18, 10,293 | 11/11, 4,016 | 11/11, 472,882 | 57/57, 495,101 |
| sweep_pN0.1 | 297 | 831 | 659 | 0 | 172 | 11 | 648 | 15/15, 7,857 | 18/18, 9,888 | 11/11, 3,977 | 11/11, 497,773 | 57/57, 519,497 |
| sweep_pN0.15 | 314 | 843 | 633 | 0 | 210 | 11 | 622 | 15/15, 7,902 | 18/18, 9,498 | 11/11, 3,999 | 11/11, 533,623 | 57/57, 555,024 |
| sweep_pN0.2 | 334 | 855 | 605 | 0 | 250 | 10 | 595 | 15/15, 7,950 | 18/18, 9,078 | 11/11, 4,020 | 11/11, 575,419 | 57/57, 596,469 |
| sweep_pN0.25 | 356 | 872 | 565 | 0 | 307 | 9 | 556 | 15/15, 7,998 | 18/18, 8,478 | 11/11, 4,049 | 11/11, 625,228 | 57/57, 645,755 |
| sweep_pN0.3 | 381 | 880 | 543 | 0 | 337 | 6 | 537 | 15/15, 8,055 | 18/18, 8,148 | 11/11, 4,062 | 11/11, 674,964 | 57/57, 695,231 |
| final | 356 | 879 | 570 | 0 | 309 | 9 | 561 | 15/15, 8,055 | 18/18, 8,553 | 11/11, 4,082 | 11/11, 630,247 | 57/57, 650,939 |

#### soupx_emptydrops_test, `instrumented-production-sapphirerapids`

Validation: 0:01.97 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 834 | 686 | 0 | 148 | 11 | 675 | 15/15, 7,908 | 18/18, 10,293 | 11/11, 4,016 | 11/11, 472,882 | 57/57, 495,101 |
| sweep_pN0.1 | 297 | 831 | 659 | 0 | 172 | 11 | 648 | 15/15, 7,857 | 18/18, 9,888 | 11/11, 3,977 | 11/11, 497,773 | 57/57, 519,497 |
| sweep_pN0.15 | 314 | 843 | 633 | 0 | 210 | 11 | 622 | 15/15, 7,902 | 18/18, 9,498 | 11/11, 3,999 | 11/11, 533,623 | 57/57, 555,024 |
| sweep_pN0.2 | 334 | 855 | 605 | 0 | 250 | 10 | 595 | 15/15, 7,950 | 18/18, 9,078 | 11/11, 4,020 | 11/11, 575,419 | 57/57, 596,469 |
| sweep_pN0.25 | 356 | 872 | 565 | 0 | 307 | 9 | 556 | 15/15, 7,998 | 18/18, 8,478 | 11/11, 4,049 | 11/11, 625,228 | 57/57, 645,755 |
| sweep_pN0.3 | 381 | 880 | 543 | 0 | 337 | 6 | 537 | 15/15, 8,055 | 18/18, 8,148 | 11/11, 4,062 | 11/11, 674,964 | 57/57, 695,231 |
| final | 356 | 879 | 570 | 0 | 309 | 9 | 561 | 15/15, 8,055 | 18/18, 8,553 | 11/11, 4,082 | 11/11, 630,247 | 57/57, 650,939 |

#### decontx_raw_tenx_test, `instrumented-portable-cascadelake`

Validation: 0:01.79 wall, 0.1 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 572 | 0 | 221 | 11 | 561 | 15/15, 7,425 | 18/18, 8,583 | 11/11, 3,738 | 11/11, 505,145 | 57/57, 524,893 |
| sweep_pN0.1 | 333 | 800 | 548 | 0 | 252 | 9 | 539 | 15/15, 7,446 | 18/18, 8,223 | 11/11, 3,744 | 11/11, 536,804 | 57/57, 556,219 |
| sweep_pN0.15 | 353 | 811 | 518 | 0 | 293 | 8 | 510 | 15/15, 7,482 | 18/18, 7,773 | 11/11, 3,759 | 11/11, 576,625 | 57/57, 595,641 |
| sweep_pN0.2 | 375 | 823 | 498 | 0 | 325 | 6 | 492 | 15/15, 7,560 | 17/17, 7,472 | 11/11, 3,789 | 11/11, 621,369 | 56/56, 640,192 |
| sweep_pN0.25 | 400 | 841 | 467 | 0 | 374 | 6 | 461 | 15/15, 7,650 | 17/17, 7,007 | 11/11, 3,830 | 11/11, 677,009 | 56/56, 695,498 |
| sweep_pN0.3 | 429 | 854 | 440 | 0 | 414 | 6 | 434 | 15/15, 7,734 | 17/17, 6,602 | 11/11, 3,855 | 11/11, 737,006 | 56/56, 755,199 |
| final | 400 | 850 | 471 | 0 | 379 | 6 | 465 | 15/15, 7,716 | 17/17, 7,067 | 11/11, 3,870 | 11/11, 684,254 | 56/56, 702,909 |

#### decontx_raw_tenx_test, `instrumented-portable-sapphirerapids`

Validation: 0:01.89 wall, 0.1 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 572 | 0 | 221 | 11 | 561 | 15/15, 7,425 | 18/18, 8,583 | 11/11, 3,738 | 11/11, 505,145 | 57/57, 524,893 |
| sweep_pN0.1 | 333 | 800 | 548 | 0 | 252 | 9 | 539 | 15/15, 7,446 | 18/18, 8,223 | 11/11, 3,744 | 11/11, 536,804 | 57/57, 556,219 |
| sweep_pN0.15 | 353 | 811 | 518 | 0 | 293 | 8 | 510 | 15/15, 7,482 | 18/18, 7,773 | 11/11, 3,759 | 11/11, 576,625 | 57/57, 595,641 |
| sweep_pN0.2 | 375 | 823 | 498 | 0 | 325 | 6 | 492 | 15/15, 7,560 | 17/17, 7,472 | 11/11, 3,789 | 11/11, 621,369 | 56/56, 640,192 |
| sweep_pN0.25 | 400 | 841 | 467 | 0 | 374 | 6 | 461 | 15/15, 7,650 | 17/17, 7,007 | 11/11, 3,830 | 11/11, 677,009 | 56/56, 695,498 |
| sweep_pN0.3 | 429 | 854 | 440 | 0 | 414 | 6 | 434 | 15/15, 7,734 | 17/17, 6,602 | 11/11, 3,855 | 11/11, 737,006 | 56/56, 755,199 |
| final | 400 | 850 | 471 | 0 | 379 | 6 | 465 | 15/15, 7,716 | 17/17, 7,067 | 11/11, 3,870 | 11/11, 684,254 | 56/56, 702,909 |

#### decontx_raw_tenx_test, `instrumented-production-sapphirerapids`

Validation: 0:01.73 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 316 | 793 | 572 | 0 | 221 | 11 | 561 | 15/15, 7,425 | 18/18, 8,583 | 11/11, 3,738 | 11/11, 505,145 | 57/57, 524,893 |
| sweep_pN0.1 | 333 | 800 | 548 | 0 | 252 | 9 | 539 | 15/15, 7,446 | 18/18, 8,223 | 11/11, 3,744 | 11/11, 536,804 | 57/57, 556,219 |
| sweep_pN0.15 | 353 | 811 | 518 | 0 | 293 | 8 | 510 | 15/15, 7,482 | 18/18, 7,773 | 11/11, 3,759 | 11/11, 576,625 | 57/57, 595,641 |
| sweep_pN0.2 | 375 | 823 | 498 | 0 | 325 | 6 | 492 | 15/15, 7,560 | 17/17, 7,472 | 11/11, 3,789 | 11/11, 621,369 | 56/56, 640,192 |
| sweep_pN0.25 | 400 | 841 | 467 | 0 | 374 | 6 | 461 | 15/15, 7,650 | 17/17, 7,007 | 11/11, 3,830 | 11/11, 677,009 | 56/56, 695,498 |
| sweep_pN0.3 | 429 | 854 | 440 | 0 | 414 | 6 | 434 | 15/15, 7,734 | 17/17, 6,602 | 11/11, 3,855 | 11/11, 737,006 | 56/56, 755,199 |
| final | 400 | 850 | 471 | 0 | 379 | 6 | 465 | 15/15, 7,716 | 17/17, 7,067 | 11/11, 3,870 | 11/11, 684,254 | 56/56, 702,909 |

#### decontx_raw_emptydrops_test, `instrumented-portable-cascadelake`

Validation: 0:01.83 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 587 | 0 | 223 | 11 | 576 | 15/15, 7,467 | 18/18, 8,808 | 11/11, 3,821 | 11/11, 459,274 | 57/57, 479,372 |
| sweep_pN0.1 | 297 | 812 | 550 | 0 | 262 | 7 | 543 | 15/15, 7,416 | 18/18, 8,253 | 11/11, 3,796 | 11/11, 486,392 | 57/57, 505,859 |
| sweep_pN0.15 | 314 | 825 | 517 | 0 | 308 | 9 | 508 | 15/15, 7,446 | 18/18, 7,758 | 11/11, 3,813 | 11/11, 522,229 | 57/57, 541,248 |
| sweep_pN0.2 | 334 | 828 | 486 | 0 | 342 | 9 | 477 | 15/15, 7,431 | 17/17, 7,292 | 11/11, 3,794 | 11/11, 557,248 | 56/56, 575,767 |
| sweep_pN0.25 | 356 | 849 | 450 | 0 | 399 | 6 | 444 | 15/15, 7,515 | 17/17, 6,752 | 11/11, 3,845 | 11/11, 608,737 | 56/56, 626,851 |
| sweep_pN0.3 | 381 | 855 | 419 | 0 | 436 | 4 | 415 | 15/15, 7,533 | 17/17, 6,287 | 11/11, 3,840 | 11/11, 655,789 | 56/56, 673,451 |
| final | 356 | 857 | 466 | 0 | 391 | 5 | 461 | 15/15, 7,611 | 17/17, 6,992 | 11/11, 3,894 | 11/11, 614,473 | 56/56, 632,972 |

#### decontx_raw_emptydrops_test, `instrumented-portable-sapphirerapids`

Validation: 0:01.86 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 587 | 0 | 223 | 11 | 576 | 15/15, 7,467 | 18/18, 8,808 | 11/11, 3,821 | 11/11, 459,274 | 57/57, 479,372 |
| sweep_pN0.1 | 297 | 812 | 550 | 0 | 262 | 7 | 543 | 15/15, 7,416 | 18/18, 8,253 | 11/11, 3,796 | 11/11, 486,392 | 57/57, 505,859 |
| sweep_pN0.15 | 314 | 825 | 517 | 0 | 308 | 9 | 508 | 15/15, 7,446 | 18/18, 7,758 | 11/11, 3,813 | 11/11, 522,229 | 57/57, 541,248 |
| sweep_pN0.2 | 334 | 828 | 486 | 0 | 342 | 9 | 477 | 15/15, 7,431 | 17/17, 7,292 | 11/11, 3,794 | 11/11, 557,248 | 56/56, 575,767 |
| sweep_pN0.25 | 356 | 849 | 450 | 0 | 399 | 6 | 444 | 15/15, 7,515 | 17/17, 6,752 | 11/11, 3,845 | 11/11, 608,737 | 56/56, 626,851 |
| sweep_pN0.3 | 381 | 855 | 419 | 0 | 436 | 4 | 415 | 15/15, 7,533 | 17/17, 6,287 | 11/11, 3,840 | 11/11, 655,789 | 56/56, 673,451 |
| final | 356 | 857 | 466 | 0 | 391 | 5 | 461 | 15/15, 7,611 | 17/17, 6,992 | 11/11, 3,894 | 11/11, 614,473 | 56/56, 632,972 |

#### decontx_raw_emptydrops_test, `instrumented-production-sapphirerapids`

Validation: 0:01.65 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 281 | 810 | 587 | 0 | 223 | 11 | 576 | 15/15, 7,467 | 18/18, 8,808 | 11/11, 3,821 | 11/11, 459,274 | 57/57, 479,372 |
| sweep_pN0.1 | 297 | 812 | 550 | 0 | 262 | 7 | 543 | 15/15, 7,416 | 18/18, 8,253 | 11/11, 3,796 | 11/11, 486,392 | 57/57, 505,859 |
| sweep_pN0.15 | 314 | 825 | 517 | 0 | 308 | 9 | 508 | 15/15, 7,446 | 18/18, 7,758 | 11/11, 3,813 | 11/11, 522,229 | 57/57, 541,248 |
| sweep_pN0.2 | 334 | 828 | 486 | 0 | 342 | 9 | 477 | 15/15, 7,431 | 17/17, 7,292 | 11/11, 3,794 | 11/11, 557,248 | 56/56, 575,767 |
| sweep_pN0.25 | 356 | 849 | 450 | 0 | 399 | 6 | 444 | 15/15, 7,515 | 17/17, 6,752 | 11/11, 3,845 | 11/11, 608,737 | 56/56, 626,851 |
| sweep_pN0.3 | 381 | 855 | 419 | 0 | 436 | 4 | 415 | 15/15, 7,533 | 17/17, 6,287 | 11/11, 3,840 | 11/11, 655,789 | 56/56, 673,451 |
| final | 356 | 857 | 466 | 0 | 391 | 5 | 461 | 15/15, 7,611 | 17/17, 6,992 | 11/11, 3,894 | 11/11, 614,473 | 56/56, 632,972 |

#### cellbender_fromraw_test, `instrumented-portable-cascadelake`

Validation: 0:01.91 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 687 | 0 | 139 | 14 | 673 | 15/15, 7,899 | 18/18, 10,308 | 11/11, 3,982 | 11/11, 488,170 | 57/57, 510,361 |
| sweep_pN0.1 | 309 | 840 | 668 | 0 | 172 | 12 | 656 | 15/15, 7,974 | 18/18, 10,023 | 11/11, 4,021 | 11/11, 523,324 | 57/57, 545,344 |
| sweep_pN0.15 | 327 | 851 | 640 | 0 | 211 | 9 | 631 | 15/15, 8,010 | 18/18, 9,603 | 11/11, 4,040 | 11/11, 560,813 | 57/57, 582,468 |
| sweep_pN0.2 | 348 | 873 | 614 | 0 | 259 | 9 | 605 | 15/15, 8,127 | 18/18, 9,213 | 11/11, 4,102 | 11/11, 611,977 | 57/57, 633,421 |
| sweep_pN0.25 | 371 | 892 | 592 | 0 | 300 | 11 | 581 | 15/15, 8,244 | 18/18, 8,883 | 11/11, 4,154 | 11/11, 666,328 | 57/57, 687,611 |
| sweep_pN0.3 | 397 | 898 | 551 | 0 | 347 | 8 | 543 | 15/15, 8,235 | 18/18, 8,268 | 11/11, 4,140 | 11/11, 717,506 | 57/57, 738,151 |
| final | 371 | 885 | 573 | 0 | 312 | 10 | 563 | 15/15, 8,145 | 18/18, 8,598 | 11/11, 4,108 | 11/11, 661,099 | 57/57, 681,952 |

#### cellbender_fromraw_test, `instrumented-portable-sapphirerapids`

Validation: 0:01.87 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 687 | 0 | 139 | 14 | 673 | 15/15, 7,899 | 18/18, 10,308 | 11/11, 3,982 | 11/11, 488,170 | 57/57, 510,361 |
| sweep_pN0.1 | 309 | 840 | 668 | 0 | 172 | 12 | 656 | 15/15, 7,974 | 18/18, 10,023 | 11/11, 4,021 | 11/11, 523,324 | 57/57, 545,344 |
| sweep_pN0.15 | 327 | 851 | 640 | 0 | 211 | 9 | 631 | 15/15, 8,010 | 18/18, 9,603 | 11/11, 4,040 | 11/11, 560,813 | 57/57, 582,468 |
| sweep_pN0.2 | 348 | 873 | 614 | 0 | 259 | 9 | 605 | 15/15, 8,127 | 18/18, 9,213 | 11/11, 4,102 | 11/11, 611,977 | 57/57, 633,421 |
| sweep_pN0.25 | 371 | 892 | 592 | 0 | 300 | 11 | 581 | 15/15, 8,244 | 18/18, 8,883 | 11/11, 4,154 | 11/11, 666,328 | 57/57, 687,611 |
| sweep_pN0.3 | 397 | 898 | 551 | 0 | 347 | 8 | 543 | 15/15, 8,235 | 18/18, 8,268 | 11/11, 4,140 | 11/11, 717,506 | 57/57, 738,151 |
| final | 371 | 885 | 573 | 0 | 312 | 10 | 563 | 15/15, 8,145 | 18/18, 8,598 | 11/11, 4,108 | 11/11, 661,099 | 57/57, 681,952 |

#### cellbender_fromraw_test, `instrumented-production-sapphirerapids`

Validation: 0:01.93 wall, 0.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 293 | 826 | 687 | 0 | 139 | 14 | 673 | 15/15, 7,899 | 18/18, 10,308 | 11/11, 3,982 | 11/11, 488,170 | 57/57, 510,361 |
| sweep_pN0.1 | 309 | 840 | 668 | 0 | 172 | 12 | 656 | 15/15, 7,974 | 18/18, 10,023 | 11/11, 4,021 | 11/11, 523,324 | 57/57, 545,344 |
| sweep_pN0.15 | 327 | 851 | 640 | 0 | 211 | 9 | 631 | 15/15, 8,010 | 18/18, 9,603 | 11/11, 4,040 | 11/11, 560,813 | 57/57, 582,468 |
| sweep_pN0.2 | 348 | 873 | 614 | 0 | 259 | 9 | 605 | 15/15, 8,127 | 18/18, 9,213 | 11/11, 4,102 | 11/11, 611,977 | 57/57, 633,421 |
| sweep_pN0.25 | 371 | 892 | 592 | 0 | 300 | 11 | 581 | 15/15, 8,244 | 18/18, 8,883 | 11/11, 4,154 | 11/11, 666,328 | 57/57, 687,611 |
| sweep_pN0.3 | 397 | 898 | 551 | 0 | 347 | 8 | 543 | 15/15, 8,235 | 18/18, 8,268 | 11/11, 4,140 | 11/11, 717,506 | 57/57, 738,151 |
| final | 371 | 885 | 573 | 0 | 312 | 10 | 563 | 15/15, 8,145 | 18/18, 8,598 | 11/11, 4,108 | 11/11, 661,099 | 57/57, 681,952 |

#### soupx_tenx_SRR18256555, `instrumented-portable-cascadelake`

Validation: 0:34.22 wall, 0.9 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 15,452 | 2,000 | 0 | 3,783 | 35 | 1,965 | 18/18, 130,893 | 20/20, 30,005 | 11/11, 54,109 | 10/10, 9,855,908 | 61/61, 10,070,917 |
| sweep_pN0.1 | 3,453 | 15,621 | 2,000 | 0 | 4,897 | 19 | 1,981 | 18/18, 130,560 | 20/20, 30,005 | 11/11, 55,746 | 10/10, 10,399,246 | 61/61, 10,615,559 |
| sweep_pN0.15 | 3,656 | 15,829 | 2,000 | 0 | 6,270 | 29 | 1,971 | 18/18, 130,087 | 20/20, 30,005 | 11/11, 57,733 | 10/10, 11,008,662 | 61/61, 11,226,489 |
| sweep_pN0.2 | 3,885 | 16,019 | 2,000 | 0 | 7,694 | 19 | 1,981 | 18/18, 129,446 | 20/20, 30,005 | 11/11, 59,737 | 10/10, 11,696,042 | 61/61, 11,915,232 |
| sweep_pN0.25 | 4,144 | 16,197 | 2,000 | 0 | 9,016 | 16 | 1,984 | 18/18, 129,003 | 20/20, 30,005 | 11/11, 61,596 | 10/10, 12,473,398 | 61/61, 12,694,004 |
| sweep_pN0.3 | 4,440 | 16,447 | 2,000 | 0 | 10,177 | 15 | 1,985 | 18/18, 129,569 | 20/20, 30,005 | 11/11, 63,508 | 10/10, 13,361,898 | 61/61, 13,584,982 |
| final | 4,144 | 16,207 | 2,000 | 0 | 9,139 | 22 | 1,978 | 18/18, 128,837 | 20/20, 30,005 | 11/11, 61,743 | 10/10, 12,473,418 | 61/61, 12,694,005 |

#### soupx_tenx_SRR18256555, `instrumented-portable-sapphirerapids`

Validation: 0:33.97 wall, 1.0 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 15,452 | 2,000 | 0 | 3,783 | 35 | 1,965 | 18/18, 130,893 | 20/20, 30,005 | 11/11, 54,109 | 10/10, 9,855,908 | 61/61, 10,070,917 |
| sweep_pN0.1 | 3,453 | 15,621 | 2,000 | 0 | 4,897 | 19 | 1,981 | 18/18, 130,560 | 20/20, 30,005 | 11/11, 55,746 | 10/10, 10,399,246 | 61/61, 10,615,559 |
| sweep_pN0.15 | 3,656 | 15,829 | 2,000 | 0 | 6,270 | 29 | 1,971 | 18/18, 130,087 | 20/20, 30,005 | 11/11, 57,733 | 10/10, 11,008,662 | 61/61, 11,226,489 |
| sweep_pN0.2 | 3,885 | 16,019 | 2,000 | 0 | 7,694 | 19 | 1,981 | 18/18, 129,446 | 20/20, 30,005 | 11/11, 59,737 | 10/10, 11,696,042 | 61/61, 11,915,232 |
| sweep_pN0.25 | 4,144 | 16,197 | 2,000 | 0 | 9,016 | 16 | 1,984 | 18/18, 129,003 | 20/20, 30,005 | 11/11, 61,596 | 10/10, 12,473,398 | 61/61, 12,694,004 |
| sweep_pN0.3 | 4,440 | 16,447 | 2,000 | 0 | 10,177 | 15 | 1,985 | 18/18, 129,569 | 20/20, 30,005 | 11/11, 63,508 | 10/10, 13,361,898 | 61/61, 13,584,982 |
| final | 4,144 | 16,207 | 2,000 | 0 | 9,139 | 22 | 1,978 | 18/18, 128,837 | 20/20, 30,005 | 11/11, 61,743 | 10/10, 12,473,418 | 61/61, 12,694,005 |

#### soupx_tenx_SRR18256555, `instrumented-production-sapphirerapids`

Validation: 0:33.48 wall, 0.9 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 3,272 | 15,452 | 2,000 | 0 | 3,783 | 35 | 1,965 | 18/18, 130,893 | 20/20, 30,005 | 11/11, 54,109 | 10/10, 9,855,908 | 61/61, 10,070,917 |
| sweep_pN0.1 | 3,453 | 15,621 | 2,000 | 0 | 4,897 | 19 | 1,981 | 18/18, 130,560 | 20/20, 30,005 | 11/11, 55,746 | 10/10, 10,399,246 | 61/61, 10,615,559 |
| sweep_pN0.15 | 3,656 | 15,829 | 2,000 | 0 | 6,270 | 29 | 1,971 | 18/18, 130,087 | 20/20, 30,005 | 11/11, 57,733 | 10/10, 11,008,662 | 61/61, 11,226,489 |
| sweep_pN0.2 | 3,885 | 16,019 | 2,000 | 0 | 7,694 | 19 | 1,981 | 18/18, 129,446 | 20/20, 30,005 | 11/11, 59,737 | 10/10, 11,696,042 | 61/61, 11,915,232 |
| sweep_pN0.25 | 4,144 | 16,197 | 2,000 | 0 | 9,016 | 16 | 1,984 | 18/18, 129,003 | 20/20, 30,005 | 11/11, 61,596 | 10/10, 12,473,398 | 61/61, 12,694,004 |
| sweep_pN0.3 | 4,440 | 16,447 | 2,000 | 0 | 10,177 | 15 | 1,985 | 18/18, 129,569 | 20/20, 30,005 | 11/11, 63,508 | 10/10, 13,361,898 | 61/61, 13,584,982 |
| final | 4,144 | 16,207 | 2,000 | 0 | 9,139 | 22 | 1,978 | 18/18, 128,837 | 20/20, 30,005 | 11/11, 61,743 | 10/10, 12,473,418 | 61/61, 12,694,005 |

#### soupx_tenx_10xmousebrain, `instrumented-portable-sapphirerapids`

Validation: 2:10.06 wall, 3.4 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 20,816 | 2,000 | 39 | 3,923 | 75 | 1,899 | 18/18, 189,539 | 20/20, 30,005 | 11/11, 70,275 | 10/10, 31,628,636 | 61/61, 31,918,457 |
| sweep_pN0.1 | 11,111 | 21,029 | 2,000 | 18 | 5,032 | 50 | 1,918 | 18/18, 190,259 | 20/20, 30,005 | 11/11, 72,042 | 10/10, 33,384,062 | 61/61, 33,676,370 |
| sweep_pN0.15 | 11,765 | 21,265 | 2,000 | 11 | 6,213 | 42 | 1,918 | 18/18, 191,101 | 20/20, 30,005 | 11/11, 73,931 | 10/10, 35,346,534 | 61/61, 35,641,573 |
| sweep_pN0.2 | 12,500 | 21,501 | 2,000 | 4 | 7,430 | 28 | 1,908 | 18/18, 192,131 | 20/20, 30,005 | 11/11, 75,846 | 10/10, 37,552,006 | 61/61, 37,849,990 |
| sweep_pN0.25 | 13,333 | 21,746 | 2,000 | 3 | 8,645 | 25 | 1,911 | 18/18, 193,361 | 20/20, 30,005 | 11/11, 77,799 | 10/10, 40,051,496 | 61/61, 40,352,663 |
| sweep_pN0.3 | 14,286 | 21,937 | 2,000 | 5 | 9,812 | 23 | 1,915 | 18/18, 194,479 | 20/20, 30,005 | 11/11, 79,543 | 10/10, 42,910,878 | 61/61, 43,214,907 |
| final | 15,775 | 22,064 | 2,000 | 3 | 8,945 | 25 | 1,912 | 18/18, 200,201 | 20/20, 30,005 | 11/11, 79,054 | 10/10, 47,378,132 | 61/61, 47,687,394 |

#### soupx_tenx_10xmousebrain, `instrumented-production-sapphirerapids`

Validation: 2:07.93 wall, 3.4 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 20,816 | 2,000 | 39 | 3,923 | 75 | 1,899 | 18/18, 189,539 | 20/20, 30,005 | 11/11, 70,275 | 10/10, 31,628,636 | 61/61, 31,918,457 |
| sweep_pN0.1 | 11,111 | 21,029 | 2,000 | 18 | 5,032 | 50 | 1,918 | 18/18, 190,259 | 20/20, 30,005 | 11/11, 72,042 | 10/10, 33,384,062 | 61/61, 33,676,370 |
| sweep_pN0.15 | 11,765 | 21,265 | 2,000 | 11 | 6,213 | 42 | 1,918 | 18/18, 191,101 | 20/20, 30,005 | 11/11, 73,931 | 10/10, 35,346,534 | 61/61, 35,641,573 |
| sweep_pN0.2 | 12,500 | 21,501 | 2,000 | 4 | 7,430 | 28 | 1,908 | 18/18, 192,131 | 20/20, 30,005 | 11/11, 75,846 | 10/10, 37,552,006 | 61/61, 37,849,990 |
| sweep_pN0.25 | 13,333 | 21,746 | 2,000 | 3 | 8,645 | 25 | 1,911 | 18/18, 193,361 | 20/20, 30,005 | 11/11, 77,799 | 10/10, 40,051,496 | 61/61, 40,352,663 |
| sweep_pN0.3 | 14,286 | 21,937 | 2,000 | 5 | 9,812 | 23 | 1,915 | 18/18, 194,479 | 20/20, 30,005 | 11/11, 79,543 | 10/10, 42,910,878 | 61/61, 43,214,907 |
| final | 15,775 | 22,064 | 2,000 | 3 | 8,945 | 25 | 1,912 | 18/18, 200,201 | 20/20, 30,005 | 11/11, 79,054 | 10/10, 47,378,132 | 61/61, 47,687,394 |

#### soupx_emptydrops_4500871lane, `instrumented-portable-sapphirerapids`

Validation: 1:55.95 wall, 5.3 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 20,805 | 2,000 | 27 | 3,569 | 86 | 1,862 | 18/18, 190,237 | 20/20, 30,005 | 11/11, 69,851 | 10/10, 31,628,614 | 61/61, 31,918,709 |
| sweep_pN0.1 | 11,111 | 21,030 | 2,000 | 12 | 4,273 | 58 | 1,898 | 18/18, 191,815 | 20/20, 30,005 | 11/11, 71,266 | 10/10, 33,384,064 | 61/61, 33,677,152 |
| sweep_pN0.15 | 11,765 | 21,267 | 2,000 | 6 | 5,064 | 29 | 1,917 | 18/18, 193,527 | 20/20, 30,005 | 11/11, 72,787 | 10/10, 35,346,538 | 61/61, 35,642,859 |
| sweep_pN0.2 | 12,500 | 21,490 | 2,000 | 7 | 5,980 | 30 | 1,906 | 18/18, 194,967 | 20/20, 30,005 | 11/11, 74,361 | 10/10, 37,551,984 | 61/61, 37,851,319 |
| sweep_pN0.25 | 13,333 | 21,710 | 2,000 | 4 | 6,936 | 30 | 1,903 | 18/18, 196,503 | 20/20, 30,005 | 11/11, 75,974 | 10/10, 40,051,424 | 61/61, 40,353,908 |
| sweep_pN0.3 | 14,286 | 21,947 | 2,000 | 2 | 7,857 | 20 | 1,906 | 18/18, 198,477 | 20/20, 30,005 | 11/11, 77,609 | 10/10, 42,910,898 | 61/61, 43,216,991 |
| final | 35,464 | 23,860 | 2,000 | 2 | 8,895 | 22 | 1,890 | 18/18, 253,969 | 20/20, 30,005 | 11/11, 84,370 | 10/10, 106,448,724 | 61/61, 106,817,070 |

#### soupx_emptydrops_4500871lane, `instrumented-production-sapphirerapids`

Validation: 1:53.71 wall, 5.3 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 20,805 | 2,000 | 27 | 3,569 | 86 | 1,862 | 18/18, 190,237 | 20/20, 30,005 | 11/11, 69,851 | 10/10, 31,628,614 | 61/61, 31,918,709 |
| sweep_pN0.1 | 11,111 | 21,030 | 2,000 | 12 | 4,273 | 58 | 1,898 | 18/18, 191,815 | 20/20, 30,005 | 11/11, 71,266 | 10/10, 33,384,064 | 61/61, 33,677,152 |
| sweep_pN0.15 | 11,765 | 21,267 | 2,000 | 6 | 5,064 | 29 | 1,917 | 18/18, 193,527 | 20/20, 30,005 | 11/11, 72,787 | 10/10, 35,346,538 | 61/61, 35,642,859 |
| sweep_pN0.2 | 12,500 | 21,490 | 2,000 | 7 | 5,980 | 30 | 1,906 | 18/18, 194,967 | 20/20, 30,005 | 11/11, 74,361 | 10/10, 37,551,984 | 61/61, 37,851,319 |
| sweep_pN0.25 | 13,333 | 21,710 | 2,000 | 4 | 6,936 | 30 | 1,903 | 18/18, 196,503 | 20/20, 30,005 | 11/11, 75,974 | 10/10, 40,051,424 | 61/61, 40,353,908 |
| sweep_pN0.3 | 14,286 | 21,947 | 2,000 | 2 | 7,857 | 20 | 1,906 | 18/18, 198,477 | 20/20, 30,005 | 11/11, 77,609 | 10/10, 42,910,898 | 61/61, 43,216,991 |
| final | 35,464 | 23,860 | 2,000 | 2 | 8,895 | 22 | 1,890 | 18/18, 253,969 | 20/20, 30,005 | 11/11, 84,370 | 10/10, 106,448,724 | 61/61, 106,817,070 |

#### soupx_emptydrops_sscovelli, `instrumented-portable-sapphirerapids`

Validation: 2:06.96 wall, 10.8 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 20,282 | 2,000 | 14 | 1,665 | 22 | 1,939 | 18/18, 189,689 | 20/20, 30,005 | 11/11, 66,455 | 9/9, 49,568 | 60/60, 335,719 |
| sweep_pN0.1 | 11,111 | 20,380 | 2,000 | 6 | 1,882 | 9 | 1,947 | 18/18, 191,125 | 20/20, 30,005 | 11/11, 66,974 | 9/9, 49,764 | 60/60, 337,870 |
| sweep_pN0.15 | 11,765 | 20,502 | 2,000 | 6 | 2,282 | 10 | 1,943 | 18/18, 192,707 | 20/20, 30,005 | 11/11, 67,736 | 9/9, 50,008 | 60/60, 340,458 |
| sweep_pN0.2 | 12,500 | 20,621 | 2,000 | 2 | 2,703 | 6 | 1,961 | 18/18, 194,109 | 20/20, 30,005 | 11/11, 68,532 | 9/9, 50,246 | 60/60, 342,894 |
| sweep_pN0.25 | 13,333 | 20,734 | 2,000 | 1 | 3,223 | 3 | 1,946 | 18/18, 195,799 | 20/20, 30,005 | 11/11, 69,376 | 9/9, 50,472 | 60/60, 345,654 |
| sweep_pN0.3 | 14,286 | 20,829 | 2,000 | 0 | 3,754 | 3 | 1,941 | 18/18, 197,455 | 20/20, 30,005 | 11/11, 70,187 | 9/9, 50,662 | 60/60, 348,311 |
| final | 105,879 | 22,247 | 2,000 | 0 | 3,795 | 6 | 1,912 | 18/18, 391,847 | 20/20, 30,005 | 11/11, 74,453 | 9/9, 53,498 | 60/60, 549,805 |

#### cellbender_fromraw_ghirsutum2, `instrumented-portable-sapphirerapids`

Validation: 3:12.75 wall, 9.3 GiB peak RSS (8 threads).

| context | cells | genes kept | fitting genes | theta Inf | Poisson genes | outliers | smoothing genes | S2a | S2b | S2c | S2d | S2-chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sweep_pN0.05 | 10,526 | 48,747 | 2,000 | 38 | 10,367 | 47 | 1,876 | 18/18, 400,289 | 20/20, 30,005 | 11/11, 160,489 | 9/9, 106,498 | 60/60, 697,283 |
| sweep_pN0.1 | 11,111 | 49,191 | 2,000 | 18 | 12,446 | 33 | 1,875 | 18/18, 400,861 | 20/20, 30,005 | 11/11, 163,899 | 9/9, 107,386 | 60/60, 702,153 |
| sweep_pN0.15 | 11,765 | 49,730 | 2,000 | 20 | 15,702 | 23 | 1,892 | 18/18, 399,817 | 20/20, 30,005 | 11/11, 168,789 | 9/9, 108,464 | 60/60, 707,077 |
| sweep_pN0.2 | 12,500 | 50,219 | 2,000 | 6 | 19,146 | 14 | 1,885 | 18/18, 398,467 | 20/20, 30,005 | 11/11, 173,693 | 9/9, 109,442 | 60/60, 711,609 |
| sweep_pN0.25 | 13,333 | 50,630 | 2,000 | 12 | 22,510 | 19 | 1,863 | 18/18, 397,019 | 20/20, 30,005 | 11/11, 178,268 | 9/9, 110,264 | 60/60, 715,558 |
| sweep_pN0.3 | 14,286 | 51,113 | 2,000 | 3 | 25,680 | 11 | 1,863 | 18/18, 396,625 | 20/20, 30,005 | 11/11, 182,887 | 9/9, 111,230 | 60/60, 720,749 |
| final | 79,097 | 57,899 | 2,000 | 5 | 25,439 | 21 | 1,837 | 18/18, 580,699 | 20/20, 30,005 | 11/11, 202,978 | 9/9, 124,802 | 60/60, 938,486 |

### Totals over the 24 runs

| stage | checks identical | elements compared | elements differing |
|---|---|---|---|
| S2a | 2,709 of 2,709 | 13,618,458 | 0 |
| S2b | 3,126 of 3,126 | 2,793,846 | 0 |
| S2c | 1,848 of 1,848 | 5,482,667 | 0 |
| S2d | 1,771 of 1,771 | 1,499,093,848 | 0 |
| S2-chain | 9,790 of 9,790 | 1,520,989,155 | 0 |

All 24 runs pass: every check of every context is identical, teacher-forced and chained.

## Reference vectors

| family | cases | what they cover |
|---|---|---|
| `orderstat` (rcompat) | 11 | `quantile(type = 7)`, `median` (the long-double mean of the middle pair), `mad`, `IQR`, with ties and even lengths |
| `density` (rcompat) | 6 | `bw.nrd` and `density(bw = "nrd")`: `BinDist`, R >= 4.4's grid, the FFT convolution, `approx` |
| `gamma` (rcompat) | 401 | R's `lgamma` (`lgammafn`), `digamma`, `trigamma` from 1e-310 to 1e300, negative non-integers, and the arguments glmGamPoi's likelihood reaches (counts and half counts plus 1/theta, theta from 1e-16 to 1e16); glibc's `lgamma` is called directly |
| `matrixstats` (rcompat) | 5 | `rowMeans2`, `colMeans2`, `rowVars` (double accumulators, `refine = TRUE`), checked equal to DelayedMatrixStats' |
| `weighted_median` (rcompat) | 30 | `weightedMedian` with and without interpolation, with `dnorm` weights as `loc_median_fit` uses them: ties (`R_qsort_I` is not stable), zeros and a huge value |
| `nlminb` (rcompat) | 42 | 21 problems with and without the Hessian, with bounds hit and not: every evaluation point, the result, iterations, evaluations and the PORT return code |
| `dnbinom` (rcompat) | 720 (11,520 values) | `dnbinom(x, size, mu = mu, log = TRUE)` on both sides of every `stirlerr` branch and series cut-off, x = 0, x < 1e-10 size, non-integer x |
| `optimize` (rcompat) | 16 | `optimize()` minimum and maximum: every evaluation point; plateaus, kinks, Inf, -Inf, NaN |
| `uniroot` (rcompat) | 9 | `uniroot()`: every evaluation point, root, iterations, precision; -Inf at an end point, a step function, `maxiter` reached, a root at an end point |
| `bw_sj` (rcompat) | 19 | `bw.SJ()` from n = 2 to 3,000: unbinned and binned pair counts, ties, bimodal and heavy-tailed samples, one error case |
| `ksmooth` (rcompat) | 12 | `ksmooth()`, box and normal kernels, ties, points outside the data (NA) |
| `bincode` (rcompat) | 32 | `.bincode()`, right- and left-closed, `include.lowest`, points on the breaks, NaN |
| `sct_input` (sctransform) | 5 | DoubletFinder's SCT input as Matrix builds it (explicit zeros), the gene filter, `median(umi@x)`, `row_gmean`, `row_var` |
| `fit_beta` (sctransform) | 24 | glmGamPoi's `fitBeta_one_group`: 14 genes on which Newton reaches NaN and `optimize()` takes over (two with half-integer counts), ordinary genes, an all-zero gene, theta = 0 |
| `is_outlier` (sctransform) | 10 | sctransform's `is_outlier`, n = 30 to 2,000, planted outliers, ties |
| `cli_soupx_tenx_test` | 6 files | R's `SCTransform()` on a tier-T object's counts in `scrust sctransform`'s layout (scale.data by its SHA-256) |

## Benchmark

SCTransform on the final run of the two largest objects, Rust and R on the same input: the counts with that run's artificial doublets, built from the export and the instrumented run's S1 draws (`scrust sctransform --from-stages ... --context final`, and [`validation/R/sctransform_benchmark.R`](../../validation/R/sctransform_benchmark.R), which builds it as DoubletFinder 2.0.6 does). Portable profile (OPENBLAS_CORETYPE=Haswell, the instrumented runs'), Intel Xeon Platinum 8480+ and 8480CL (`sapphire`, shared nodes), [`validation/slurm/m4_benchmark.sbatch`](../../validation/slurm/m4_benchmark.sbatch): Rust jobs 59770299 and 59770303 (commit 7a273be, 32 CPUs allocated) and 59769186 and 59769206 (b70bb14, before rows were dropped early), R jobs 59769194 and 59769218 (one CPU). *SCTransform* is R's `SCTransform()` call alone (`proc.time`) and Rust's `sctransform` phase (`summary.json`); process wall, CPU time and peak RSS come from `/usr/bin/time -v`, and include reading `matrix.mtx.gz` and building the input (R: `mtx_read_10x`, 35-40 s; Rust: 15 s).

| object (final run) | cells x genes kept | label | side | threads | SCTransform (s) | process wall (s) | process CPU (s) | peak RSS (GiB) | RSS after loading (GiB) |
|---|---|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | 105,879 x 22,247 | m4-benchmark-7a273be | Rust | 1 | 134.8 | 149.3 | 148.8 | 9.5 | 3.1 |
| soupx_emptydrops_sscovelli | 105,879 x 22,247 | m4-benchmark-7a273be | Rust | 32 | 15.2 | 29.7 | 193.0 | 27.4 | 3.1 |
| soupx_emptydrops_sscovelli | 105,879 x 22,247 | m4-benchmark | R | 1 | 353.9 | 414.1 | 410.7 | 88.8 | 4.3 |
| soupx_emptydrops_sscovelli | 105,879 x 22,247 | m4-benchmark | Rust | 1 | 133.6 | 148.3 | 150.6 | 9.5 | 3.1 |
| soupx_emptydrops_sscovelli | 105,879 x 22,247 | m4-benchmark | Rust | 32 | 13.9 | 28.6 | 175.2 | 27.6 | 3.1 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 57,899 | m4-benchmark-7a273be | Rust | 1 | 240.8 | 256.7 | 255.9 | 8.0 | 3.5 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 57,899 | m4-benchmark-7a273be | Rust | 32 | 18.7 | 34.6 | 295.9 | 21.6 | 3.5 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 57,899 | m4-benchmark | R | 1 | 641.9 | 699.9 | 695.0 | 144.1 | 4.7 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 57,899 | m4-benchmark | Rust | 1 | 255.3 | 270.5 | 272.4 | 8.0 | 3.5 |
| cellbender_fromraw_ghirsutum2 | 79,097 x 57,899 | m4-benchmark | Rust | 32 | 19.0 | 34.0 | 302.2 | 31.0 | 3.5 |

**Every output was byte-identical to R's:** `cell_attr.tsv`, `gene_attr.tsv`, `model_pars.tsv`, `model_pars_fit.tsv`, `variable_features.tsv` and `scale_data.npy` (2.5 GB for sscovelli, 1.9 GB for ghirsutum2), for every Rust run, at 1 and 32 threads.

- **Time.** On one thread Rust's SCTransform is 2.6 (sscovelli) and 2.7 (ghirsutum2) times faster than R's; on 32 threads 23 and 34 times. R's time includes `sctransform::correct()`, the corrected counts that `do.correct.umi = TRUE` adds to the SCT assay and DoubletFinder never reads; Rust does not compute them.
- **Memory.** R holds the full residual matrix: 88.8 and 144.1 GiB at its peak. Rust needs 9.5 and 8.0 GiB on one thread, 9 and 18 times less: the counts and the SCT input (3.1-3.5 GiB after loading), the rows of the 3,000 variable features (2.5 and 1.9 GB), and per thread one bin's product (500 genes x cells, 0.4 and 0.3 GB) while it is turned into residuals. On 32 threads the bins in flight dominate (27.4 and 21.6 GiB): workers drop a row once its variance is below the 3,000th best so far, but sscovelli has only 45 bins, and the first 32 start before 3,000 genes are ranked. Computing each bin's rows on all threads, one bin at a time, would make the memory independent of the thread count (a possible follow-up; the bins' products must keep R's `dgemm` shapes).
- **R's SCTransform share, from M2.** M2's runs time DoubletFinder's phases, not SCTransform itself. In the instrumented portable runs, the final run's `doubletFinder()` (SCTransform, PCA, the distance matrix, pANN) took 1,712 s (sscovelli) and 1,527 s (ghirsutum2), with peaks of 127.5 and 153.1 GiB. SCTransform on the same input measured here is 354 s (21%) and 642 s (42%) of that, at 88.8 and 144.1 GiB. In the direct production runs (SapphireRapids kernels) the phase took 1,477 and 1,258 s, of which the same SCTransform would be 24% and 51%. The sweep's six SCTransforms (on 10,000 real cells and their doublets each) are inside `paramSweep`'s 458 and 753 s, which M2 does not break down further.

## Findings for the port

- **glmGamPoi falls back to `optimize()` more quietly than it looks.** `fitBeta_one_group` calls `optimize()` over the negative binomial log-likelihood when a Newton step makes beta NaN, not only when the 100 iterations run out, and it still reports the Newton iteration count. So R's dumped iterations (at most 9 in the data) did not reveal that one gene of cellbender_fromraw_ghirsutum2 (fitting gene 1,910 in the first sweep context) takes the fallback: its first Newton fit overshoots from -9.84 to 5.8e8 by the fourth step, and R's answer, -12.517085358459251, is `optimize()`'s. The port runs the fallback the same way. For a half-integer count (a doublet), `dnbinom` is -Inf for every beta, so `optimize()` minimises a constant `DBL_MAX` and returns where its golden-section steps end; the port does too (checked on R vectors).
- **R's compiled PORT code differs from its source.** DRMNGB calls `DL7TVM(N, V(X01), V(L), V(X01))` with the same array as input and output. gfortran assumes no aliasing and stores `X(I) = 0` before loading `Y(I)` (disassembly of `dl7tvm_` in R's `stats.so`), so the result is all zeros, not L'y. That zero vector then makes R accept steps that the source's arithmetic would reject. The port reproduces the compiled order. DRMNGB (the gradient-only branch) is reached in the data: glmGamPoi falls back to it when the Hessian version stops with "false convergence" (soupx_emptydrops_test).
- **`nlminb` calls BLAS.** PORT's `dd7tpr` and `dv2nrm` are `ddot` and `dnrm2` in R's build, so even the overdispersion MLE goes through OpenBLAS; at these sizes (one or two parameters) the kernels agree, and S2b is identical across profiles.
- **Two gamma functions.** glmGamPoi's count-table path uses the C library's `lgamma` (glibc), its Rcpp sugar path R's own `lgammafn`, `digamma` and `trigamma`; and inside R's math library `stirlerr` calls glibc's `lgamma` too. Each call site uses the one R uses.
- **matrixStats computes in double.** Its `LDOUBLE` is `long double` only under `HAVE_LONG_DOUBLE`, which R defines in its private `config.h`, not in the headers packages are built with. Base R's `mean`, `var`, `sum` and `rowMeans` do use the x87 long double (emulated by `rcompat::F80`).
- **The count table's order.** glmGamPoi's `make_table_if_small` builds a libstdc++ `unordered_map<long, size_t>` with `reserve()`; its iteration order (and so the Cox-Reid sums' order) follows the bucket count from libstdc++'s prime list (304 primes, read from the R env's `libstdc++.so.6.0.36`) and the singly linked node list's insertion rules. Keys are `(long)` casts, so half-integer counts truncate.
- **`future.seed` costs two draws per bin.** future 1.82.0's `make_rng_seeds` and future.apply 1.21.0's `future_xapply` each call `next_random_seed()` (one uniform), on R's main stream; the L'Ecuyer kind switch in between is undone.
- **Explicit zeros matter.** `dgCMatrix + dgCMatrix` keeps the union of the patterns with explicit zeros, and sctransform's `row_var_dgcmatrix` adds stored entries' squared deviations in storage order and the implicit zeros at the end, so an explicit zero changes the bits; `median(umi@x)` counts them too.
- **The last bin's product.** R's `tcrossprod` switches to `dgemv` when a bin has one gene; the port follows R's dispatch (and its own triple loop when NaN or Inf may be present).
- **Memory.** R holds the full genes x cells residual matrix (37 GB for ghirsutum2's final run) and copies it while clipping; the port computes bins in parallel and keeps only the running best 3,000 genes' rows.

## Plan corrections

The sources and the dumps settle these; plan B3a (and the brief's hazard list) said otherwise:

- **Two RNG draws per fitting bin, not one** (B3a step 4): `make_rng_seeds()` and `future_xapply()` each call `next_random_seed()`.
- **Seurat 5.5.1's default path computes the full residuals inside `vst`** (`sct.method = "default"`, `residual_type = "pearson"`), not through `get_residual_var`/`get_residuals`; the residual mean is base `rowMeans` (long double) and the residual variance `matrixStats::rowVars` (double, three passes), not a two-pass variance.
- **`scale.data` is centred in plain R**: `ScaleData` calls `FastRowScale`, `mat - matrixStats::rowMeans2(mat, na.rm = TRUE)`, not an Eigen kernel.
- **matrixStats accumulates in double**, so "long-double mean/var" applies to base R only (glmGamPoi's `mean()` of `colMeans2`, the MLE's start, the residual mean).
- **BLAS is used in three places, not two:** the residual products (`dgemm`, `dgemv` for a one-gene bin), the Cox-Reid term (Armadillo's `dgemv`), and `nlminb` itself (`ddot`, `dnrm2`).
- **glmGamPoi was built against RcppArmadillo 15.2.3-1** (Armadillo 15.2), not the 15.6.2-1 the env now holds; the port follows 15.2's code (the 1 x 1 `lu`, `inv_sympd` and `trace` of the Cox-Reid term).
- **`fitBeta_one_group` has an `optimize()` fallback** that the data reach (above); the plan's list of fits did not include it.
- **The exclude-Poisson rule has a one-row quirk:** `rownames(model_pars[diff_theta < 1e-3, ])` is NULL when exactly one gene qualifies, so that gene keeps its theta. In `reg_model_pars` the same construct only feeds messages.
- **Infinite thetas also come from MLE = 0** (an overdispersion estimate of exactly zero), and the data reach both `nlminb` bounds (log 1e-16 and log 1e16).
- **The Poisson genes' intercept uses `mean(colSums(umi))` of the kept genes** (Matrix's double `colSums` of the filtered matrix, then R's long-double `mean`), not the cell attribute `umi`, which sums all genes.

## Deviations

- **Not ported, never reached:** the Poisson limit of `dnbinom` (`size = Inf`, which needs nmath's `dpois_raw` and `ebd0`), reached only if glmGamPoi's `optimize()` fallback ran with theta = 0; the port returns an error there instead of a number. glmGamPoi's `optimize()` fallback is otherwise ported and checked.
- **Only DoubletFinder's settings.** `SctParams` rejects what the port does not implement: `exclude_poisson = FALSE`, `gmean_eps` other than 1, parallel `future` plans; batch variables, non-regularised latent variables, deviance residuals, `min_variance` "model_mean"/"model_median", `conserve.memory`, `residual.features` and reference models are not offered.
- **Corrected counts are not computed.** `do.correct.umi = TRUE` makes R also run `sctransform::correct()` for the SCT assay's counts; DoubletFinder never reads them, and `scrust sctransform` does not write them (the benchmark's R time includes them).
- **Implemented, checked on R vectors, not reached in the data:** `weightedMedian` with interpolation (used by `loc_median_fit` when a bin has fewer genes than its window), glmGamPoi's narrow-data trend branch, `uniroot` failures, and the threshold selection of variable features (`variable_features_n = None`).
- **`scale.data` of the light runs** (sscovelli, ghirsutum2) is not in R's dumps; their S2d check covers the residual mean and variance of every kept gene and the features, and the benchmark compares their final runs' `scale.data` with R's byte for byte.
- **Validation labels.** The first full S2b run (`rust/m4b/`, M4b's code before it was committed) stopped on the ghirsutum2 gene above and is kept as a record. The passing runs of each sub-milestone are kept too: `rust/m4b-34cce6f/` (S2a-S2b, job 59759589), `rust/m4c-f9b645b/` (S2a-S2c, job 59763140), `rust/m4d-9e607ac/` (S2a-S2d teacher-forced, 16 of 24 runs: the job was cancelled because a rebuild during it would have changed the binary of its later runs; the wrapper now runs a pinned copy, `SCRUST_BIN`), `rust/m4d-3e810f1/` (teacher-forced and chained, all 24 runs, job 59767410). The tables above are from `rust/m4d-7a273be/`, the final code. A cancelled job left `soupx_tenx_10xmousebrain/rust/m4d-9e607ac/instrumented-portable-sapphirerapids.partial.59765712/`, an incomplete directory.
- **The benchmark runs the portable profile** (Haswell kernels on Sapphire Rapids), the profile of the instrumented runs it compares with; production runs detect SapphireRapids kernels.

## Reproducing

From the repository's top directory, after `source env.sh`:

```bash
cargo build --release -p scrust
cp target/release/scrust data/doubletfinder/logs/bin/scrust-$(git rev-parse --short HEAD)
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-$(git rev-parse --short HEAD) sbatch --cpus-per-task=8 --export=ALL validation/slurm/m4_sctransform.sbatch m4d-$(git rev-parse --short HEAD) soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2
SCRUST_BIN=... sbatch --cpus-per-task=32 --mem=64G --export=ALL validation/slurm/m4_benchmark.sbatch rust soupx_emptydrops_sscovelli
sbatch --cpus-per-task=1 --mem=250G validation/slurm/m4_benchmark.sbatch r soupx_emptydrops_sscovelli
```

The R reference vectors (no `.rds`; a login node is enough), with the production R env and the hermetic flags:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 "$R_ENV/bin/Rscript" validation/R/sctransform_vectors.R tests/fixtures/sctransform
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 "$R_ENV/bin/Rscript" validation/R/sctransform_cli_reference.R tests/fixtures/doubletfinder/soupx_tenx_test/input <out dir>
```

One context by hand, on a tier-T fixture (seconds):

```bash
target/release/scrust validate sctransform --input tests/fixtures/doubletfinder/soupx_tenx_test/input --stages tests/fixtures/doubletfinder/soupx_tenx_test/stages --contexts final
target/release/scrust sctransform --counts tests/fixtures/doubletfinder/soupx_tenx_test/input --out /tmp/sct
```
