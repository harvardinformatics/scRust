# scRust plan: Rust reimplementations of scRNA-seq tools, starting with DoubletFinder

## Context

scRust (`harvardinformatics/scRust`, checked out at `/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRust`) is a standalone home for Rust reimplementations of single-cell RNA-seq tools. The aim is for each to use less memory, run faster, or ideally both.

**In scope now:** emptyDrops (DropletUtils), SoupX, CellBender, decontX, DoubletFinder and scDblFinder. Downsampling and FindMarkers are out of scope for the time being.

**Acceptance:** a Rust implementation is accepted when it reproduces the original tool's outputs on the reference data, exactly or near-exactly. Only then may it replace the original in the Snakemake preprocessing workflow. That substitution is a separate, later change in the pipeline repo, not part of scRust.

**Why these tools:** memory. In the production run, the workflow jobs around these tools peaked at:
- DoubletFinder 125–137 GiB (it builds a full cell-by-cell distance matrix);
- scDblFinder 86 GiB;
- SoupX 61 GiB, decontX 60 GiB, emptyDrops 57 GiB.

Those peaks include the Seurat processing done in the same jobs, so each tool's own share is measured directly, as part of the reference runs (A4).

**First tool: DoubletFinder.** It has the largest memory problem, and it brings in shared pieces that later tools reuse: R's random-number stream, SCTransform, PCA and nearest neighbours.

**Next: a clustering crate** (A7b) that reproduces Seurat's clustering without Seurat. All of its parameters can be tuned, for future cluster-parameter optimisation.

**Decisions so far:**
- one repository, organised as a Cargo workspace;
- GPL-3 license, replacing the MIT license the repo was created with, so the GPL-3 R packages can be ported faithfully;
- no Seurat objects: tools read sparse matrices and small tables, and write sparse matrices and tables;
- validation uses the workflow's test data plus selected production samples.

---

## Part A: General framework

### A1. Repository layout

```
scRust/
  Cargo.toml                 workspace
  crates/
    rcompat/                 R-compatible random numbers and R arithmetic semantics (A5)
    scio/                    sparse matrix I/O (10x Matrix Market and HDF5); TSV/JSON/.npy
    sctransform/             SCTransform v2 port (with glmGamPoi's fitting)
    pca/                     port of irlba as Seurat calls it
    knn/                     exact nearest neighbours, memory-bounded, parallel over cells
    kde/                     KernSmooth::bkde, approxfun, moments
    doubletfinder/           first tool (library)
    neighbors/               Annoy (C++ through FFI) and Seurat's shared-nearest-neighbour graph
    cluster/                 Seurat's Louvain optimizer (Leiden later); the clustering chain (A7b)
    (later) emptydrops/, decontx/, scdblfinder/, cellbender/, soupx/ (separate binary, A6)
    scrust/                  single command-line binary: `scrust <tool> ...`
  validation/
    R/, python/              export scripts and direct reference runs of the original tools
    slurm/                   sbatch wrappers (node guard, hermetic env flags, kernel pinning)
    manifest.toml            reference datasets: paths, sha256, sizes, provenance, roles
    tolerances.toml          match criterion for each tool, stage and output
  tests/fixtures/            small committed fixtures (test data, RNG vectors); ≤ ~20 MB
  data/                      large fixtures, gitignored, listed in manifest.toml
  envs/                      rust-dev.yml; r-reference/ (pinned copy of the workflow's R env spec)
  third_party/               fetch_sources.sh; src/ (gitignored): exact upstream sources used as porting references
  docs/                      PLAN.md (this document), validation reports
```

- **Tools:** each tool is a library crate plus a subcommand of one binary, so the repo ships as one conda package.
- **Shared stages are library crates linked into every tool that needs them, never copied.**
  - SCTransform, for example, is `crates/sctransform` with a typed API such as `sctransform::vst(&counts, &SctParams, &mut RRng) -> SctFit`.
  - DoubletFinder calls it 7 times per object (6 sweep values plus the final run) on matrices that exist only inside the tool, so it can't be a separate step that passes files between processes.
  - One implementation, validated once, is inherited by every tool. `SctParams` mirrors R's arguments for callers that use different settings.
  - Each shared stage is also a standalone subcommand (`scrust sctransform`, `scrust pca`) for validation and direct use.
- **One R-compatible random-number generator is passed explicitly** (`&mut RRng`) through every call, like R's single global stream. A function that calls `set.seed` in R reseeds it the same way in Rust.

### A2. Environments and HPC rules

- **Rust toolchain:** comes from `envs/rust-dev.yml`, created as a project-local env. `env.sh` sets `CARGO_HOME` inside the repo (`.cargo-home/`, gitignored), so nothing is written to the home directory or `~/.local`. The env pins:
  - conda-forge `rust` 1.99.x;
  - the exact BLAS builds of the R env, so FFI calls hit the same kernels: `libopenblas 0.3.34 pthreads_hf13c14d_2`, `libblas`/`liblapack 3.11.0-12 *_openblas`.
- **R reference env:** the production env (`/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_`): R 4.5.3, Seurat 5.5.1, sctransform 0.4.3, glmGamPoi 1.22.0, irlba 2.4.1, fields 17.3, DoubletFinder 2.0.6 (commit `1b244d8`), DropletUtils 1.30.0, SoupX 1.6.2, decontX 1.8.0, scDblFinder 1.24.10, xgboost 3.2.0. Copy its spec into `envs/r-reference/` so it can be rebuilt. The CellBender env is added when that tool starts.
- **Running R:** always with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB` and `PYTHONNOUSERSITE=1`. Anything that reads `.rds` runs through sbatch or srun, never on a login node. sbatch wrappers source the node guard and `--exclude` the two nodes known to lose a Lustre OST.
- **Pipeline repo:** read-only from scRust. Never write there and never run Snakemake there.

### A3. Data formats: no Seurat objects

Each tool's real inputs and outputs are sparse matrices and small per-cell tables:

| Tool | Inputs | Outputs |
|---|---|---|
| emptyDrops | raw matrix | per-barcode statistics, called barcodes |
| SoupX | raw and filtered matrices, cluster labels | contamination estimate, corrected matrix |
| CellBender | raw matrix | cell probabilities, ambient profile, denoised matrix |
| decontX | filtered matrix, cluster labels, raw matrix (background) | contamination per cell, decontaminated matrix |
| DoubletFinder | counts matrix, cluster labels (for the homotypic adjustment) | pANN, classification, pK, nExp, sweep table |
| Clustering (Seurat chain) | counts matrix, parameters | cluster labels; optionally the intermediates (SCTransform, PCA, neighbour graph) for reuse |
| scDblFinder | counts matrix | score and class per cell |

- **Formats:** matrices in 10x Matrix Market (and 10x HDF5 later); tables as TSV; run metadata (parameters, seed, versions, timings) in `summary.json`.
- **Where Seurat appears:**
  - In the workflow, Seurat objects are only containers, plus R's processing after each tool (SCTransform, clustering), which is not part of these tools.
  - Some tool inputs (cluster labels) were made by Seurat upstream; they become plain input files.
  - DoubletFinder internally runs Seurat's SCTransform and PCA. The Rust version reimplements those calculations, but no Seurat object is involved.
  - `.rds` files are read only by the one-time export scripts (A4).

### A4. Reference data and validation

**Reference outputs** are what the original tool produces on a given input, captured completely in plain files. They are made by running the tool directly, through a thin script with the workflow's parameters and seeds, in the pinned production env.

**Sources:**

| Tier | What | Committed? | Use |
|---|---|---|---|
| T | The workflow's test data: one small mouse sample (raw 2,000 genes × 20,000 barcodes; 300 cells) plus the test reference outputs in `tests/reference_outputs/` | yes, after export | fast `cargo test`, CI |
| P | A handful of production samples per tool, chosen to cover the size range and size-dependent code paths | no: `data/` plus manifest | scale paths, memory and time benchmarks |
| Full panel | All production samples, compared against the outputs the production run already saved | no | final gate before workflow substitution |

- **Test data is necessary but not sufficient.**
  - It is too small for code paths that switch on only at scale: DoubletFinder's 10,000-cell sweep subsample, SCTransform's 5,000-cell and 2,000-gene sampling (the test data has only 2,000 genes, so gene sampling never runs), and scDblFinder's behaviour on real data sizes (its empty-training case appeared only on the tiny test object).
  - It cannot show memory or time gains.
- **Not every sample needs a new R run.** Production already ran every tool on every sample. Those saved outputs serve as the full-panel references after a one-time export, so the full panel costs only the Rust runs.
- **Direct reference runs are still needed for P and T**, for two reasons:
  - The workflow's saved outputs are incomplete. DoubletFinder and scDblFinder `.rds` files keep only singlets, so the scores of removed cells are gone. Tools' internal statistics are not saved either.
  - Workflow jobs mix the tool with Seurat processing, so their memory and time are not a baseline for the tool itself.
- **Cross-check:** each direct reference run must reproduce the corresponding production (or test-reference) output. That proves the thin script calls the tool exactly as the workflow did.
- **Export:** `validation/R/export_*.R` (run via sbatch) converts workflow artifacts to plain formats:
  - counts as Matrix Market, tables as TSV;
  - dense arrays as `.npy`, written by a ~20-line base-R helper;
  - RNG state checkpoints (`.Random.seed`) as `.npy`.
- **Instrumented runs** dump a tool's internal stages (for DoubletFinder: doublet indices, SCTransform, PCA, neighbours), so differences can be traced to the stage that caused them.
- **Match criteria** (`tolerances.toml`), set per tool, stage and output:
  - `exact`: bitwise or integer equality;
  - `numeric`: relative/absolute tolerance;
  - `agreement`: fraction identical, correlation.

  The R-vs-R floor sets realistic tolerances: the same tool re-run on another CPU type. Earlier work found R bit-identical for the same CPU type and thread count, and about 1e-11 apart in PCA across CPU types.
- **CPU type:** OpenBLAS chooses its kernels by CPU, so each production job's node CPU type is recorded (sacct `NodeList` → node features). Production-profile runs use `--constraint` for the same CPU type. All 150 production DoubletFinder jobs ran on Sapphire Rapids nodes (38 nodes, AVX-512), so validation jobs use `--constraint=sapphirerapids`.
- **Harness:**
  - Tier-T checks run in `cargo test`.
  - Tier-P and full-panel checks run through `scrust validate <tool> --reference <dir> --candidate <dir>` (sbatch), which writes a markdown and JSON report to `docs/validation/`.
- **Benchmarks:** wall time and peak RSS (from `/usr/bin/time -v` and sacct) for Rust against the direct reference run of the original tool on the same input.

### A5. R-semantics fidelity rules (the `rcompat` crate)

Exact matching depends on copying R's behaviour in these places. Each is unit-tested against R-generated vectors:

- **Random numbers:**
  - Mersenne-Twister, including `set.seed`'s initial scrambling and `unif_rand`'s fixup;
  - `R_unif_index` rejection sampling (sample.kind "Rejection");
  - `sample()` with and without replacement, and weighted sampling without replacement;
  - `rnorm` by inversion (`qnorm5`).
- **Draw counts:** the number of draws each R function consumes (for example irlba's start vector), because later samples continue the same stream.
- **`round()` rounds half to even** (`f64::round_ties_even`). Production logs confirm it: n = 1,122 at pN = 0.2 gave 280 artificial doublets, not 281.
- **`seq(from, to, by)` arithmetic:** `from + i*by`, then `pmin(x, to)`.
- **`order()`:** radix, stable; ties go to the lower index, also with `decreasing = TRUE`.
- **`sum()` and `mean()`:** accumulate in 80-bit long double. Emulate this exactly (for example with the `rustc_apfloat` x87 type). `var` and `sd` are two-pass.
- **Number to string:** R formats with 15 significant digits and chooses between fixed and scientific notation (`5e-04` vs `0.001`). This matters wherever values pass through strings, such as DoubletFinder's pK.
- **Summation order:** follow R's loop order for distances and other reductions; parallelise only across independent units (cells).
- **Powers:** `x^2` is `x*x`. Every other power goes through glibc `pow`, which the Rust binary links at runtime.
- **FFT:** R's `fft` (Singleton's mixed-radix algorithm) is ported as-is wherever an R function uses it, e.g. `bkde` and `density`.
- **BLAS and LAPACK:** R's matrix products go through OpenBLAS 0.3.34, built with `DYNAMIC_ARCH`, so the kernel (and its rounding) is chosen per CPU at runtime. Ports of BLAS-dependent code call that same conda OpenBLAS through FFI (`dgemm_`, `dgemv_`, `dgesdd_`), with the same routines and arguments R uses.
- **No FMA:** R and its packages are compiled with `-march=nocona -O2` and use no FMA. Rust never fuses operations implicitly, so plain `f64` code matches.

### A5b. Reproducibility across machines and clusters

**What varies by machine:**
- **OpenBLAS kernel.** OpenBLAS picks its kernel per CPU at runtime, so BLAS-based results change in the last bits from one CPU to another. Measured on a Sapphire Rapids login node: auto-detection chose `SapphireRapids`; forcing `Haswell` changed a matrix-vector product, and forcing `Sandybridge` changed a matrix-matrix product.
- **Thread count.**
- **Host libm.** glibc picks the same FMA variants on any modern x86-64, so this is a minor factor.

**What doesn't vary:** Rust's own arithmetic. There is no implicit FMA and no fast-math, and the parallel work is designed to be deterministic, so it gives identical results on any machine.

**Mechanisms:**
1. **Pin the BLAS kernel and threads.** `env.sh` and every sbatch wrapper set `OPENBLAS_CORETYPE` and `OPENBLAS_NUM_THREADS=1`, for both reference generation and Rust runs. Both envs use the identical OpenBLAS build, so the same machine code runs and results are the same on any CPU that supports that kernel's instruction set. M2 checks this across CPU types.
2. **Two reference profiles.**
   - `portable`, for development and CI: references generated with `OPENBLAS_CORETYPE=Haswell` (AVX2 + FMA, present on essentially all Intel and AMD x86-64 since 2013). They can be regenerated on any cluster, and Rust runs with the same pin.
   - `production`: the kernel production actually used (auto-detected `SapphireRapids`), for comparisons with the saved production outputs, run on Cannon Sapphire Rapids nodes. Forcing that kernel by name on other AVX-512 machines remains to be verified: in the test, `OPENBLAS_CORETYPE=SapphireRapids` reported `Cooperlake`.
3. **Environment fingerprint.** Every reference and every Rust run records:
   - CPU model and OpenBLAS core name (via `openblas_get_corename`);
   - thread count and glibc version;
   - library and package versions.

   The harness runs exact checks only when the BLAS-relevant fields match. Otherwise it switches to tolerance mode (the R-vs-R floor) and says so in the report. A different machine therefore produces neither false failures nor false "exact" claims.
4. **Teacher forcing keeps most tests machine-independent.** Each stage can take the reference's own intermediate as input. Stages without BLAS are then exact on any machine: RNG, sampling, filters, the glmGamPoi Newton fits and `nlminb`, regularisation, distances, ranks, pANN, BC and pK. For DoubletFinder, only the Cox–Reid term in S2b, S2d and S3 depend on BLAS.
5. **Fixtures are recipes.** The manifest stores a sha256 for each portable fixture. Developers on another cluster regenerate the fixtures from the pinned env instead of copying gigabytes; matching checksums prove the environment is identical before any Rust test runs.
6. **Non-x86 machines (ARM, macOS):** stages that use BLAS or libm run in tolerance mode; BLAS-free stages stay exact.
7. **CI:** a GitHub Actions x86-64 runner uses the portable profile with the small committed fixtures. Production-profile tests are `#[ignore]`d and run through sbatch on Cannon.

### A6. Licensing and attribution

- Replace the MIT `LICENSE` with GPL-3.0 (SPDX `GPL-3.0-only`, matching the GPL-3 packages being ported).
- Add a `NOTICE` that credits the original authors and licenses:
  - DoubletFinder (CC0);
  - sctransform, glmGamPoi, irlba, DropletUtils, scDblFinder (GPL-3);
  - decontX/celda and Seurat (MIT);
  - xgboost (Apache-2.0);
  - CellBender (BSD-3-Clause);
  - R (GPL ≥ 2).
- SoupX is GPL-2-only, which is incompatible with GPL-3. It will get its own binary (it can stay in this repo), or an implementation written from the paper.

### A7. Tools and order

| Order | Tool | Original | Notes |
|---|---|---|---|
| 1 | DoubletFinder | R, DoubletFinder 2.0.6 | Part B. Brings `rcompat`, `sctransform`, `pca`, `knn`, `kde` |
| 2 | Clustering | R, Seurat 5.5.1 (SCTransform → PCA → FindNeighbors → FindClusters) | A7b. Reuses `sctransform` and `pca`; adds `neighbors` and `cluster`. Gives SoupX, decontX and DoubletFinder a Seurat-free source of clusters, and supports cluster-parameter optimisation |
| 3 | scDblFinder | R, scDblFinder 1.24.10 | Reuses PCA and kNN. Adds xgboost (its C++ library through FFI; R's xgboost draws from R's RNG) |
| 4 | emptyDrops | R, DropletUtils 1.30.0 | Ambient profile, multinomial log-probabilities, Monte Carlo p-values (RNG to be mapped), knee detection |
| 5 | decontX | R, decontX 1.8.0 | Variational EM. The workflow passes clusters, which avoids decontX's UMAP+DBSCAN initialisation |
| 6 | SoupX | R, SoupX 1.6.2 | Separate GPL-2 binary. `autoEstCont`, `adjustCounts` with random rounding (`rbinom`) |
| 7 | CellBender | Python/PyTorch | Different in kind: a neural network trained stochastically on GPU, which doesn't reproduce itself exactly across runs. The criterion has to be agreement within its own run-to-run variation. Whether Rust can beat GPU PyTorch on speed needs checking before committing to it |

Tools that take cluster labels (SoupX, decontX, DoubletFinder) accept them as an input file from any source. The clustering crate is one such source, not a hidden dependency.

- Each tool follows the same cycle:
  1. Read the exact code path in the pinned env.
  2. Write the export, direct reference and instrumented scripts.
  3. Port it.
  4. Validate stage by stage.
  5. Validate on T, P and the full panel.
  6. Benchmark.
- The details of each later tool (parameters as the workflow uses them, RNG use, exactness hazards) are researched when that tool starts, as was done for DoubletFinder.

### A7b. Clustering crate (after DoubletFinder)

**Independent of DoubletFinder.** DoubletFinder runs only SCTransform and PCA; it never builds a neighbour graph or clusters. The clustering crate calls the same general-purpose `sctransform` and `pca` crates with its own settings. Neither inherits the other's.

**What it reproduces:** the workflow's base clustering (`tenx2seuratrds.R` and the other conversion scripts):
- `SCTransform(seed.use = 12345)`;
- `RunPCA(seed.use = 12345)`, 50 PCs;
- `FindNeighbors(dims = 1:30)`: k = 20, Annoy with 50 trees and Euclidean distance, shared-nearest-neighbour pruning 1/15;
- `FindClusters(random.seed = 12345)`: original Louvain, resolution 0.8, 10 starts, 10 iterations.

UMAP, which the workflow runs between PCA and the neighbour search, is not needed. Whether it affects anything downstream will be confirmed when this crate starts.

**New components:**
- `neighbors`: Annoy through FFI to the same C++ version RcppAnnoy 0.0.23 bundles, with Seurat's insertion order and float conversion; Seurat's shared-nearest-neighbour (Jaccard) graph and pruning;
- `cluster`: a port of Seurat's C++ modularity optimizer (Louvain algorithms 1–3) with its Java-style RNG, and Seurat's ordering of cluster IDs. Leiden can follow.

**Built for parameter optimisation:**
- **Every Seurat argument is a parameter, with Seurat's default:**
  - SCTransform: seed, fitting cells, number of variable features, clipping;
  - PCA: number of PCs, seed, weighting by variance;
  - neighbours: dims, k, Annoy trees and metric, pruning;
  - clustering: algorithm, resolution, starts, iterations, seed, singleton grouping.

  They are set by CLI flags or a TOML parameter file, so sweeps need no code changes.
- **Staged subcommands** read and write their intermediates: `scrust sctransform` → `scrust pca` → `scrust neighbors` → `scrust cluster`. A sweep over k or resolution therefore reuses one SCTransform and PCA.
- **Several resolutions per run,** as `FindClusters` allows.

**Validation:**
- **Workflow settings:** compared against `seurat_clusters` in the workflow's base objects, for every sample, after export.
- **Beyond the defaults:** a parameter grid (for example k ∈ {10, 20, 30, 50}, resolution ∈ {0.2, 0.4, 0.8, 1.2, 2.0}, dims ∈ {1:10, 1:30}) on the test data and two production samples, against Seurat run directly with the same settings. Rust can then stand in for Seurat across the range you'd optimise over.
- **Criterion:** identical labels (exact target), with a floor of ARI ≥ 0.99 where every difference is traced to a near-tie.
- **Milestones** are written when the crate starts, following Part B's pattern.

---

## Part B: Step 1, DoubletFinder

### B1. Reference behaviour

The reference is DoubletFinder 2.0.6 (commit `1b244d8`) as `workflow/scripts/doubletfinder.R` calls it, with seed 12345 (from `SCRNASEQ_PREPROCESS_SEED`). The workflow's choices of pK and nExp are part of that behaviour, and they are parameters of the Rust tool.

1. `set.seed(12345)`, then `paramSweep(seu, PCs = 1:10, sct = TRUE)`:
   - **pK grid:** `c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))`, keeping only values with `round(pK*min.cells) >= 1`, where `min.cells = round(n/0.95 - n)`.
   - **pN grid:** `seq(0.05, 0.3, by = 0.05)`.
   - **Subsample:** if n > 10,000 cells, `sample(1:n, 10000)` without replacement.
   - **For each pN, in order:**
     - `n_doublets = round(n.real/(1-pN) - n.real)`;
     - two `sample(real.cells, n_doublets, replace = TRUE)` draws;
     - doublets = (counts₁ + counts₂)/2, appended after the real cells;
     - SCTransform with defaults, then `RunPCA(npcs = 10)`;
     - `fields::rdist(pca[, 1:10])[, 1:n.real]`;
     - for each real cell, `order()` its distances and keep the top `round(nCells*max(pK)) + 5` rows;
     - for each pK, with `pk.temp = round(nCells*pK)`: pANN = (number of indices > n.real among order positions 2…pk.temp+1) / pk.temp.
2. `summarizeSweep(GT = FALSE)`: pN and pK are parsed back from the list names (strings). For each (pN, pK):
   - `approxfun(bkde(pANN, kernel = "normal"))`, evaluated on `seq(min, max, length.out = n.real)`;
   - bimodality coefficient from DoubletFinder's own `skewness` and `kurtosis`.
3. `find.pK`: for each pK, BCmetric = mean(BC over pN) / var(BC). The chosen pK is `as.numeric(as.character(pK[which.max(BCmetric)]))`, so the value passes through R's string formatting.
4. `nExp = round(round(0.15*n) * (1 - sum((table(seurat_clusters)/n)^2)))`.
5. `set.seed(12345)`, then `doubletFinder(pN = 0.25, pK, nExp, sct = TRUE)`:
   - same doublet construction, SCTransform and PCA, on all n cells;
   - full `rdist`, which is O(N²) memory;
   - k = `round(nCells*pK)`, and pANN uses order positions 2…k+1. Position 1 is skipped whatever cell it is; it is not necessarily the cell itself;
   - the top nExp cells by `order(pANN, decreasing = TRUE)` are labelled "Doublet".
   - Columns written: `pANN_0.25_<pK>_<nExp>` and `DF.classifications_0.25_<pK>_<nExp>`.

What the wrapper does afterwards (keeping singlets, then Seurat re-processing) is not part of DoubletFinder, so it is neither reimplemented nor compared. DoubletFinder changes no counts, so matching the doublet set is the same as matching the filtered matrix.

### B2. Map of random-number use

| Point | Draws | Notes |
|---|---|---|
| paramSweep start (seed 12345) | n > 10k: 10,000 without replacement | `sample(1:n, 10000)` |
| per pN: doublet pairs | 2 × n_doublets with replacement | continues from the previous pN's PCA |
| inside SCTransform | resets to seed 1448145 (three times along the call chain) | N > 5,000: 5,000 cells without replacement; more than 2,000 eligible genes: weighted draw of 2,000 genes; then 1 draw per 500-gene bin |
| inside RunPCA | resets to seed 42, then irlba's start vector | `rnorm(nfeat)`, i.e. 2·nfeat uniforms; more `rnorm` only in irlba's near-invariant-subspace branch (rare, but ported) |
| doubletFinder (seed 12345) | 2 × n_doublets with replacement | then SCTransform and PCA as above |

Consequence: the doublets for pN 2–6 do not come from seed 12345. They come from `set.seed(42)` plus whatever the previous pN's irlba consumed. The port has to reproduce that count, and the RNG checkpoints check it.

The underlying primitives were all checked bit-exactly against R by emulation:
- `R_unif_index` rejection sampling, including the quirk that 16-bit draws cost two uniforms;
- swap-remove sampling without replacement;
- weighted sampling without replacement via `revsort`;
- `rnorm` by inversion, two uniforms per normal.

### B3. Rust design

**Command:**
```
scrust doubletfinder --counts <10x mtx dir> --clusters <tsv> --out <dir>
                     [--seed 12345] [--pcs 10] [--pn 0.25] [--nexp-rate 0.15] [--pk auto|<value>] [--threads N]
```

**Outputs:**
- `cells.tsv`: barcode, pANN, classification, for every input cell;
- `sweep_stats.tsv`: pN, pK, BCreal;
- `bcmvn.tsv`;
- `summary.json`: pK, pK label exactly as R formats it, nExp_poi, homotypic proportion, nExp, n_doublets, k, the column names R would write, versions, seed, timings, peak RSS.

**Components:**
- **`rcompat`:** all of A5.
- **`sctransform`:** see B3a.
- **`pca`:** irlba port, see B3b.
- **`knn` (memory design):** never build the N×N matrix.
  - For each real cell, compute distances to all N cells into an O(N) buffer per thread, with the same arithmetic as `fields::rdist`.
  - Then stable-select the first `ind` (sweep) or k+1 (final) by (distance, index), which is exactly R's stable `order`.
  - The sweep stores only the cumulative doublet-neighbour count at each pK.
  - Work is parallel over cells, and each result is independent of scheduling, so output is deterministic.
  - Memory is O(N·10 + threads·N + n·#pK), versus about 84 GB for N = 102k in R. Time is about 1e11 multiply-adds on the largest object: seconds to minutes with threads.
- **`kde`:** port of `KernSmooth::bkde` (binning and convolution), `approxfun` (linear, rule 1), and DoubletFinder's skewness, kurtosis and bimodality coefficient, using R's `sum`/`mean` semantics.
- **`doubletfinder`:** the sweep, summary, pK choice (first maximum, NaN ignored), homotypic adjustment, classification (stable descending order), and R-identical column names.

**B3a. SCTransform v2 (defaults as DoubletFinder calls it).** This is the hardest part of the port.

- **Call chain:** `SCTransform.Seurat` → `.StdAssay` → `.Assay` → `.default`. Each calls `set.seed(1448145)`; the last one governs `vst`. Then `sctransform::vst(vst.flavor = "v2")`.
- **Settings:**
  - method `glmGamPoi_offset`, `exclude_poisson = TRUE`, `min_variance = (median(umi@x)/5)^2`;
  - fitting uses min(5,000, N) cells and 2,000 genes; `min_cells = 5`; `bin_size = 500`;
  - theta regularisation `od_factor`, `bw_adjust = 3`;
  - residuals clipped at ±sqrt(N) for the variance, ±sqrt(N/30) for `scale.data`;
  - 3,000 variable features; sequential `future` plan.
- **Steps:**
  1. `log_umi = log10(colSums)` over all genes. Keep genes with `rowSums(umi >= 0.01) >= 5`. Geometric means come from sctransform's C++.
  2. **Fitting cells:** if N > 5,000, `sample(colnames, 5000)`, a partial Fisher–Yates. The draw order becomes the column order.
  3. **Fitting genes:** Poisson exclusion (var − mean > 0, over all cells). If more than 2,000 remain, weighted sampling without replacement, with probability 1/density at the gene's log geometric mean (`density(bw = "nrd")`, then `approx`). It sorts with R's `revsort`, a heapsort that is not stable, so tie order matters. The draw order defines the 500-gene bins.
  4. **Per-gene fits, bin by bin.** Each bin also advances the main RNG by one draw (from `future.seed`). The fit is glmGamPoi `glm_gp(~1, offset = log(10^log_umi))`:
     - starting moments, then Newton `fitBeta_one_group`;
     - overdispersion MLE via `nlminb` (R's PORT library) with a Cox–Reid adjustment. **Its count lookup table truncates the doublets' half-integer counts (1.5 → 1).** It is a libstdc++ `unordered_map`, so its iteration order also sets the summation order;
     - running-median dispersion trend (`loc_median_fit`, which couples genes within a bin);
     - a second Newton fit.
     - theta = 1/overdispersion, then the rules for infinite theta and the 1e-7 floor.
  5. **Regularisation** (`reg_model_pars`):
     - flag Poisson genes;
     - remove outliers (bins sized by `bw.SJ`, MAD scores);
     - `ksmooth` (Normal kernel, bandwidth 3·`bw.SJ`) of the intercept, the log_umi coefficient and log10(1 + gmean/θ) against log10 gmean;
     - back-transform; Poisson genes are overwritten.
  6. **Residuals** for all genes, in bins of 500:
     - μ = exp(B·[1, log_umi]) via OpenBLAS `dgemm`;
     - var = μ + μ²/θ, floored at `min_variance`;
     - Pearson residual, clipped at ±sqrt(N);
     - two-pass residual variance; the top 3,000 by a stable descending sort.
  7. **`scale.data`** = those residuals clipped at ±sqrt(N/30) and row-centred. `VariableFeatures` are ordered by decreasing residual variance.
- **What exactness requires:**
  - faithful ports of glmGamPoi's C++: the Newton steps, the score and Hessian, the Cox–Reid term, and the lookup table;
  - R's `nlminb` (PORT `drmnhb`/`drmngb`), `density` (Singleton FFT, `bw.nrd`), `bw.SJ` with `uniroot`/`zeroin`, `ksmooth`, matrixStats `weightedMedian`, and long-double `mean`/`var`;
  - OpenBLAS through FFI wherever R uses it: the residual products and the Cox–Reid term.
  - None of these C/C++ sources are installed (see M1b).
- **Validation by teacher forcing:** each sub-stage (S2a–S2d) can take R's own intermediate as input. A difference therefore shows up in the sub-stage that caused it, not three stages later.

**B3b. RunPCA → irlba (npcs = 10).**
- **Call path:** `RunPCA.Seurat` → `RunPCA.Assay` → `PrepDR` → `RunPCA.default`.
- **Input matrix:** the SCT `scale.data` rows for `VariableFeatures` (the top 3,000 by residual variance, stable descending sort), dropping zero-variance rows. Columns are real cells, then `X1…Xn`. Values are clipped at ±sqrt(ncells/30) and centred, not scaled.
- **irlba call:** `set.seed(42)`, then `irlba(t(x), nv = 10)` with defaults:
  - `work = 26`, `maxit = 1000`, `tol = svtol = 1e-7`, `reorth = TRUE`;
  - irlba 2.4.1 is pure R, and every product is OpenBLAS `dgemm` through `direct_dgemm_c`;
  - its `orthog` step also uses R's `%*%`, and the small `svd(B)` is LAPACK `dgesdd`.
- **Outputs:** embeddings are `u %*% diag(d)`; `sdev = d/sqrt(ncells-1)`; there is no sign normalisation (signs don't affect distances).
- **What to port:**
  - irlba's loop as written: Lanczos steps, `convtests`, restarts;
  - the "clustered spectrum" quirk, where `rel_change_k` is always 0;
  - every product through the same OpenBLAS routine with the same arguments.
- **Expected agreement:**
  - bit-identical to R on the same CPU type with single-threaded OpenBLAS (production jobs ran with 1 CPU, and Snakemake sets the `OMP_*`/`OPENBLAS_*` thread counts to 1);
  - otherwise about 1e-13 relative.

**`fields::rdist`** is Fortran. For each element: `s = (x1[i,1]-x2[j,1])^2`, then for c = 2…10, `s = (x1[i,c]-x2[j,c])^2 + s`, then `sqrt(s)`. It uses SSE2 without FMA, so a scalar Rust loop is exact. The diagonal is exactly 0 and the matrix is symmetric.

**`order`** is a stable radix sort. Ties go to the lower original index, and that holds with `decreasing = TRUE` too: use a stable sort with a descending comparator, not a reversed ascending sort. +0 and −0 compare equal.

**`bkde`** (KernSmooth 2.23.27):
- **Bandwidth and grid:** `h = (1/(4π))^(1/10) · (243/(35n))^(1/5) · sd(x)`; grid of 401 points over `[min − 4h, max + 4h]` built with `seq(length.out = 401)`.
- **Algorithm:** Fortran linear binning; a Normal kernel truncated at 4; convolution through R's `fft` at the next power of two, P = 512 or 1024.
- **Edge case:** if every pANN is identical, h = 0 and `bkde` errors. Copy that as an error.
- **`approxfun`:** linear, `rule = 1`, `ties = mean`. The evaluation formula is tested bit-exact.

### B4. Validation data for DoubletFinder

**Tier T, the workflow's test data.**
- **Inputs:** the five test objects that DoubletFinder reads, `tests/reference_outputs/testdata/results/<arm>/seurat_<arm>_test.rds`. Their counts and `seurat_clusters` are exported once to Matrix Market and TSV.
- **References:** direct R DoubletFinder runs on those inputs, giving complete outputs plus intermediates. They are cross-checked against the stored test references: kept barcodes, singlet pANN and column names.
- **Limits:** the objects are small (267–300 cells, 2,000 genes), so none of the cell or gene subsampling paths run. They use a 27-value pK grid (0.04–0.30).

| Object | Cells | Doublets made | Singlets kept |
|---|---|---|---|
| soupx_tenx_test | 300 | 100 | 276 |
| decontx_raw_tenx_test | 300 | 100 | 273 |
| soupx_emptydrops_test | 267 | 89 | 237 |
| decontx_raw_emptydrops_test | 267 | 89 | 237 |
| cellbender_fromraw_test | 278 | 93 | 247 |

**Tier P, selected production objects.**
- **Inputs:** counts and clusters exported from the production input objects under `results/`.
- **References:** direct R DoubletFinder runs on Sapphire Rapids, giving complete outputs and intermediates (residual matrices only for small and medium). They are cross-checked against `results/doubletfinder/seurat_doubletfinder_<object>.rds`.

| Role | Object | Cells | Workflow job peak / time | Why |
|---|---|---|---|---|
| small | soupx_tenx_SRR18256555 | 3,108 | 5.2 GiB / 4.6 min | no cell subsampling; SCTransform's 2,000-gene subsample active |
| medium | soupx_tenx_10xmousebrain | 11,831 | 19.7 GiB / 14 min | 10k sweep subsample and 5k SCTransform cell subsample both active |
| medium 2 | soupx_emptydrops_4500871lane | 26,598 | 38.6 GiB / 16 min | high homotypic proportion (0.158) |
| large | soupx_emptydrops_sscovelli | 79,409 | 125.4 GiB / 41 min | memory benchmark; pK = 0.19, read off the find.pK plot that R saved (`Rplots.pdf`) |
| memory outlier | cellbender_fromraw_ghirsutum2 | 59,323 | 137.1 GiB / 44 min | 57,899 genes |

**Full panel:** all 150 production DoubletFinder runs (1,120–79,409 cells), compared against their saved outputs:
- the doublet set (input barcodes minus output barcodes);
- singlet pANN;
- pK and nExp, from the column names.

This needs only a one-time export of the inputs and those fields; R is not re-run.

### B5. Acceptance criteria (in `tolerances.toml`)

| Stage | Check | Criterion |
|---|---|---|
| S0 RNG and R semantics | vectors and checkpoints against R | exact |
| S1 sampling | subsample and doublet index pairs, RNG state after each step | exact |
| S2a SCT setup | log_umi, filtered genes, gmeans, fitting cells and genes, bins, RNG state | exact |
| S2b SCT fits | per-gene Beta and theta, given R's fitting cells and genes | target exact; fallback relative ≤ 1e-6 if `nlminb` can't be made bit-exact (decision point, see M4b) |
| S2c SCT regularisation | regularised parameters, given R's fit outputs | exact |
| S2d SCT residuals | residual variance, the 3,000 variable features (ordered), `scale.data`, given R's parameters | exact when the BLAS kernel matches (OpenBLAS through FFI, A5b); otherwise numeric tolerance |
| S3 PCA | standard deviations; embeddings (up to sign per component); RNG state afterwards | identical when the BLAS kernel matches the fixture (A5b); otherwise relative ≤ 1e-10; RNG state always exact |
| S4 neighbours | first `ind` order positions per real cell | identical for ≥ 99.99% of cells; each difference explained by a near-tie |
| S5 pANN (sweep and final) | per cell | identical for ≥ 99.9%; otherwise |Δ| ≤ 2/k |
| S6 sweep | BCreal; chosen pK | relative ≤ 1e-6; pK identical |
| S7 calls | nExp; Doublet/Singlet set | nExp identical; target identical set, floor ≥ 99.9% agreement with every difference at the pANN cutoff |
| Against saved workflow outputs (T, P, full panel) | doublet set, singlet pANN, pK and nExp as R writes them | set as in S7; pANN as in S5; pK and nExp identical |
| Performance | peak RSS and wall time against the direct R DoubletFinder run on the same input | lower memory and time; sscovelli ≤ 10 GiB |

**How the comparisons work:** everything is compared as plain files. The harness compares:
- Rust's `cells.tsv`, sweep tables and `summary.json` against the direct reference run (every cell, every intermediate);
- the same outputs against the fields exported from the saved workflow outputs.

`.rds` files are read only by the export scripts.

### B6. Milestones

| ID | Content | Done when |
|---|---|---|
| M0 | GPL-3 `LICENSE` and `NOTICE`, README, `.gitignore`, workspace skeleton, `envs/rust-dev.yml`, `env.sh`, `docs/PLAN.md`, validation skeleton | `cargo build` and `cargo test` pass in the conda env; nothing written under `~` |
| M1 | `rcompat` (A5) and R-generated test vectors | all S0 tests are bit-exact |
| M1b | `third_party/fetch_sources.sh`: download the exact source tarballs into a gitignored `third_party/src/`, with checksums. Covers R 4.5.3, sctransform 0.4.3, glmGamPoi 1.22.0, matrixStats 1.5.0, Seurat 5.5.1 and SeuratObject 5.4.0, irlba 2.4.1, KernSmooth 2.23.27, fields 17.3, DoubletFinder `1b244d8`. These are the porting references | sources present; the internals described here from knowledge (`revsort`, `rbits`, `ProbSampleNoReplace`, `fft`, `qnorm5`, radix sort) are checked against R 4.5.3's code |
| M2 | export scripts; direct and instrumented reference runs for T and P in both profiles (A5b); environment fingerprinting | each reference run reproduces its saved workflow output; portable fixtures have identical checksums on at least 2 CPU types on Cannon; production core name recorded; R-vs-R floor and the tool's own memory and time measured |
| M3 | DoubletFinder core on R's own inputs: R's doublet indices and PCA embeddings feed Rust kNN, pANN, sweep, BC, pK, nExp and calls | S4–S7 exact on all T and P objects; kNN memory and time measured on sscovelli |
| M4a | SCT setup: filters, gmeans, fitting cells and genes (`density`, `revsort` weighted sampling), bins | S2a exact on T, small, medium and large |
| M4b | SCT fits: glmGamPoi Newton, `nlminb` (PORT), Cox–Reid, lookup table, `loc_median_fit` | S2b. **Decision point:** if bit-exact `nlminb` proves impractical, report the measured downstream effect (S5–S7) and agree the fallback with you |
| M4c | SCT regularisation: `bw.SJ`, outlier bins, `ksmooth` | S2c exact (teacher-forced) |
| M4d | SCT residuals, variable features, `scale.data` | S2d exact when the BLAS kernel matches (teacher-forced) |
| M5 | `pca` (irlba) port | S3 met; RNG checkpoint after PCA matches |
| M6 | full `scrust doubletfinder` | S0–S7 and the saved-output checks met on T and the P objects; benchmarks in `docs/validation/doubletfinder.md` |
| M7 | full panel: all 150 production objects (sbatch array) | doublet sets, singlet pANN, pK and nExp meet the criteria for all 150, with a summary table |

M3 comes first on purpose: it solves the memory problem and checks all the bookkeeping before the hard numerical ports.

### B7. Risks and mitigations

- **SCTransform/glmGamPoi numerics are the main risk.**
  - `nlminb` is a large Fortran optimiser.
  - `bw.SJ` uses a root-finding tolerance of about 1%, so tiny input differences can move its result.
  - glmGamPoi truncates half-integer counts and iterates over an `unordered_map`, which both have to be copied.
  - The dispersion trend couples genes within a bin.
  - Mitigation: teacher-forced sub-stages S2a–S2d, sources from M1b, an explicit decision point at M4b, and tolerances based on the R-vs-R floor.
- **PCA floating point:** R itself differs by node, because OpenBLAS picks its kernel at runtime.
  - Mitigation: call the same OpenBLAS through FFI, and pin the kernel (A5b).
  - Fallback: about 1e-13 differences, whose effects S4 and S5 confine to near-ties.
- **RNG draw counts inside SCTransform and irlba:** check with `.Random.seed` checkpoints after every call.
- **Ties at the nExp cutoff:** pANN values are multiples of 1/k, so ties at the cutoff are common. Copy R's tie order exactly and test it on purpose.
- **Values that pass through strings** (the pK label): the formatting function is tested over the whole pK grid.

---

## Verification

- **This document:** check it against the source facts: DoubletFinder 2.0.6 function bodies, `doubletfinder.R`, and the production logs and sacct numbers.
- **Each milestone:** its "Done when" column. `cargo test` covers tier T on every change; sbatch validation runs produce tier-P and full-panel reports in `docs/validation/`.
- **Step 1 complete:** M6 met (with M7 as the full-panel confirmation), and the benchmark table shows the memory and time gains over the direct R run.
