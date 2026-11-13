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

**Next: a clustering crate** (A7b) that reproduces Seurat's clustering without Seurat, with the workflow's silhouette and purity scores. All of its parameters can be tuned, for future cluster-parameter optimisation. It also exports SCTransform's normalised data with the counts and per-cell metadata, in a form that can't be mistaken for raw counts, and comes with loaders for Seurat, SingleCellExperiment and AnnData.

**Decisions so far:**
- one repository, organised as a Cargo workspace;
- GPL-3 license, replacing the MIT license the repo was created with, so the GPL-3 R packages can be ported faithfully;
- no Seurat objects: tools read sparse matrices and small tables, and write sparse matrices and tables. R and Python loaders build Seurat, SingleCellExperiment and AnnData objects from those files (C3g);
- validation uses the workflow's test data plus selected production samples;
- Step 2: Annoy is ported to Rust; the full panel is all 840 clustered production objects; normalised data is exported as plain files, with loaders that are part of Step 2. h5ad and `.rds` files come from the loaders (C3f, C3g).

---

## Part A: General framework

### A1. Repository layout

```
scRust/
  Cargo.toml                 workspace
  crates/
    rcompat/                 R-compatible random numbers and R arithmetic semantics (A5)
    scio/                    sparse matrix I/O (10x Matrix Market and HDF5); TSV/JSON/.npy; the output manifest
    sctransform/             SCTransform v2 port (with glmGamPoi's fitting), and its corrected counts
    pca/                     port of irlba as Seurat calls it
    knn/                     exact nearest neighbours, memory-bounded, parallel over cells (DoubletFinder; purity)
    kde/                     KernSmooth::bkde, approxfun, moments
    doubletfinder/           first tool (library)
    neighbors/               Annoy (RcppAnnoy's Annoy 1.17.3, ported) and Seurat's shared-nearest-neighbour graph
    cluster/                 Seurat's Louvain optimizer (Leiden later); silhouette and purity; the clustering chain (A7b)
    (later) emptydrops/, decontx/, scdblfinder/, cellbender/, soupx/ (separate binary, A6)
    scrust/                  single command-line binary: `scrust <tool> ...`
  loaders/                   R/ and python/: packages that load scrust's outputs into Seurat, SingleCellExperiment, AnnData (C3g)
  validation/
    R/, python/              export scripts and direct reference runs of the original tools
    slurm/                   sbatch wrappers (node guard, hermetic env flags, kernel pinning)
    manifest.toml            reference datasets: paths, sha256, sizes, provenance, roles
    tolerances.toml          match criterion for each tool, stage and output
  tests/fixtures/            small committed fixtures (test data, RNG vectors); ≤ ~20 MB
  data/                      large fixtures, gitignored, listed in manifest.toml
  envs/                      rust-dev.yml; r-reference/ (pinned copy of the workflow's R env spec); loaders-py.yml
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
- **Loader test envs (C3g):** the R loader is tested in the R reference env, installed into a library inside the repo (`.rlib/`, gitignored), never into the env itself. The Python loader is tested in `envs/loaders-py.yml`, a project-local env.
- **Running R:** always with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB` and `PYTHONNOUSERSITE=1`. Anything that reads `.rds` runs through sbatch or srun, never on a login node. sbatch wrappers source the node guard and `--exclude` the two nodes known to lose a Lustre OST.
- **Pipeline repo:** read-only from scRust. Never write there and never run Snakemake there.

### A3. Data formats: no Seurat objects

Each tool's real inputs and outputs are sparse matrices and small per-cell tables:

| Tool | Inputs | Outputs |
|---|---|---|
| emptyDrops | raw matrix | per-barcode statistics, called barcodes, **the called cells' matrix** |
| SoupX | raw and filtered matrices, cluster labels | contamination estimate, **corrected matrix** |
| CellBender | raw matrix | cell probabilities, ambient profile, **denoised matrix** |
| decontX | filtered matrix, cluster labels, raw matrix (background) | contamination per cell, **decontaminated matrix** |
| DoubletFinder | counts matrix, cluster labels (for the homotypic adjustment) | pANN, classification, pK, nExp, sweep table, **the singlets' matrix** |
| Clustering (Seurat chain) | counts matrix, parameters | cluster labels, silhouette width and neighbourhood purity per cell; optionally the intermediates (SCTransform, PCA, neighbour graph) for reuse, and an export of SCTransform's normalised data with the counts and metadata (C3f) |
| scDblFinder | counts matrix | score and class per cell, **the singlets' matrix** |

- **Every output of the original, and more.** A port keeps all of the original tool's outputs, such as DoubletFinder's per-barcode pANN and classification table and its sweep tables. The 10x count data below is written in addition to them, never instead. Ports expand the original tools' functionality and remove none of it.
- **Count data out (the matrices in bold):** every tool that processes the counts writes its result as new 10x-style count data, whether it keeps cells (the filters: emptyDrops, DoubletFinder, scDblFinder) or changes values (the ambient correctors: SoupX, CellBender, decontX). Their results stay on the raw count scale. The tools then chain on plain files, each reading the last one's directory, with no Seurat object between them. The directory is gzipped 10x Matrix Market (Cell Ranger v3 layout, which Seurat's `Read10X`, `DropletUtils` and scanpy read). It holds the entries exactly as stored, explicit zeros included, with values as `%.17g` writes them, so a correctly rounded reader recovers every double. The input's feature lines are carried over unchanged. DoubletFinder writes `singlets/` from M6b on (B3).
- **Normalised data is the exception.** Data on another scale is never written as 10x count data. That covers SCTransform's corrected counts, their log1p and the residuals.
  - The corrected counts are integers (`correct_counts` rounds them). A tool given them as a 10x directory would take them for raw counts without complaint.
  - Normalised matrices go into the export of C3f instead, next to the original counts and a per-cell metadata table. They are Matrix Market or `.npy` files under names that no 10x reader accepts, each marked with what it holds and how it was made.
  - The loaders (C3g) put them only into slots meant for normalised data.
- **Formats:** count matrices in 10x Matrix Market (and 10x HDF5 later); normalised matrices in Matrix Market or `.npy` under other names; tables as TSV; run metadata (parameters, seed, versions, timings) in `summary.json`. From M13 on, `summary.json` also lists every output file with its role (C3f).
- **Where Seurat appears:**
  - In the workflow, Seurat objects are only containers, plus R's processing after each tool (SCTransform, clustering), which is not part of these tools.
  - Some tool inputs (cluster labels) were made by Seurat upstream; they become plain input files.
  - DoubletFinder internally runs Seurat's SCTransform and PCA. The Rust version reimplements those calculations, but no Seurat object is involved.
  - `.rds` files are read only by the one-time export scripts (A4), and by the comparisons of the loaders' objects with the saved ones (C4).
  - For users, the R loader (C3g) builds Seurat objects from scrust's files. scrust itself never writes one.

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
- **CPU type:** OpenBLAS chooses its kernels by CPU, so each production job's node CPU type is recorded (sacct `NodeList` → node features). Production-profile runs use `--constraint` for the same CPU type. All 150 production DoubletFinder jobs ran on Sapphire Rapids nodes (38 nodes, AVX-512), so validation jobs use `--constraint=sapphirerapids`. M7 records each saved output's job, node and CPU type in `validation/doubletfinder_panel.tsv`.
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
| 2 | Clustering | R, Seurat 5.5.1 (SCTransform → PCA → FindNeighbors → FindClusters); bluster 1.20.0 (silhouette, purity) | A7b. Reuses `sctransform`, `pca` and `knn`; adds `neighbors` and `cluster`. Gives SoupX, decontX and DoubletFinder a Seurat-free source of clusters, supports cluster-parameter optimisation, and exports normalised data with loaders for Seurat, SingleCellExperiment and AnnData |
| 3 | scDblFinder | R, scDblFinder 1.24.10 | Reuses PCA and kNN. Adds xgboost (its C++ library through FFI; R's xgboost draws from R's RNG) |
| 4 | emptyDrops | R, DropletUtils 1.30.0 | Ambient profile, multinomial log-probabilities, Monte Carlo p-values (RNG to be mapped), knee detection |
| 5 | decontX | R, decontX 1.8.0 | Variational EM. The workflow passes the base object's clusters as `z`, which avoids decontX's own clustering (UMAP, then DBSCAN), and falls back to that clustering when the base object has fewer than two clusters. The port takes clusters from a file, from clustering on the fly (C3h), or from decontX's own clustering (A7b) |
| 6 | SoupX | R, SoupX 1.6.2 | Separate GPL-2 binary. `autoEstCont`, `adjustCounts` with random rounding (`rbinom`) |
| 7 | CellBender | Python/PyTorch | Different in kind: a neural network trained stochastically on GPU, which doesn't reproduce itself exactly across runs. The criterion has to be agreement within its own run-to-run variation. Whether Rust can beat GPU PyTorch on speed needs checking before committing to it |

Tools that take cluster labels (SoupX, decontX, DoubletFinder) accept them as an input file from any source. The clustering crate is one such source, not a hidden dependency. Given no labels, a tool can also cluster its own input in-process, through the same library and parameters as `scrust cluster` (C3h). DoubletFinder gets this in Step 2, the other tools in their own steps.

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
- `FindClusters(random.seed = 12345)`: original Louvain, resolution 0.8, 10 starts, 10 iterations;
- the workflow's own `add_silhouette_to_metadata`: bluster's `approxSilhouette` and `neighborPurity` on PCA dims 1:30;
- the SCT assay that SCTransform leaves: the corrected counts, their log1p and the residuals.

UMAP, which the workflow runs between PCA and the neighbour search, is not needed: it affects nothing downstream (C1).

**New components:**
- `neighbors`: Annoy as RcppAnnoy 0.0.23 builds it (Annoy 1.17.3), ported to Rust (C3a), with Seurat's insertion order and float conversion; Seurat's shared-nearest-neighbour (Jaccard) graph and pruning;
- `cluster`: a port of Seurat's C++ modularity optimizer (Louvain algorithms 1–3) with its Java-style RNG, and Seurat's ordering of cluster IDs; bluster's silhouette and purity. Leiden can follow;
- `sctransform` gains `correct_counts`, and `knn` gains the exact search that purity makes (C3e);
- an export of the normalised data with the counts and metadata, built so it can't be mistaken for raw counts (C3f), and loaders that read it into Seurat, SingleCellExperiment and AnnData (C3g);
- clustering on the fly: a tool that needs clusters and is given none computes them in-process (C3h).

**decontX's clusters, a requirement of the future decontX port.** decontX's own clustering (UMAP, then DBSCAN) can under-cluster severely, so the port must be able to take predetermined clusters from the clustering crate. It offers three sources:
- (a) `--clusters <file>`, with labels from `scrust cluster` or any other source;
- (b) clustering on the fly, through C3h's library entry point;
- (c) decontX's own clustering, kept because a port never removes an original feature (A3). It comes with the workflow's fallback to it when fewer than two clusters are given, since decontX needs at least two populations when `z` is supplied.

Which source is the default is left to decontX's planning. Today, `decontx.R` already passes the base Seurat object's clusters as `z`, aligned to the filtered matrix's barcodes, and falls back to `z = NULL` when the base object has a single cluster. After Step 2, the labels from `scrust cluster` can already be given to the original R decontX as `z`, read through the R loader or straight from `seurat_clusters.tsv`.

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
- **Workflow settings:** compared against `seurat_clusters` in every object the workflow clustered, the base objects included, after export (C4).
- **Beyond the defaults:** a parameter grid (for example k ∈ {10, 20, 30, 50}, resolution ∈ {0.2, 0.4, 0.8, 1.2, 2.0}, dims ∈ {1:10, 1:30}) on the test data and two production samples, against Seurat run directly with the same settings. Rust can then stand in for Seurat across the range you'd optimise over.
- **Criterion:** identical labels (exact target), with a floor of ARI ≥ 0.99 where every difference is traced to a near-tie.
- **Normalised data and loaders:** the corrected counts, their log1p, the residuals, the PCA and the graphs, read through each loader, against the saved objects (C4). On T and P, FindMarkers on the loaded Seurat object must reproduce the workflow's marker tables exactly.
- **On the fly:** `scrust doubletfinder` given no labels reproduces M7's 150 saved outputs.
- **Plan and milestones:** Part C (milestones M8–M17).

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
- `summary.json`: pK, pK label exactly as R formats it, nExp_poi, homotypic proportion, nExp, n_doublets, k, the column names R would write, versions, seed, timings, peak RSS;
- `singlets/` (M6b; A3): the counts without the doublets as a new 10x directory (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`). It is the workflow step's result at the level of the counts: `doubletfinder.R` keeps the singlets with `subset()` and passes the object on. Each singlet's column is written exactly as stored, in input order, with the input's feature lines. The counts are read again from the input directory once the calls are made, so the final run's memory peak does not grow.

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

**Full panel:** all 150 production DoubletFinder runs (30 samples × 5 arms; 1,120–79,409 cells), compared against their saved outputs:
- the doublet set (input barcodes minus output barcodes);
- singlet pANN;
- pK and nExp, from the column names;
- the singlets' counts (from M6b on): `singlets/` against the saved object's RNA counts, the barcodes, the feature lines and every stored entry and value.

This needs only a one-time export of the inputs and those fields (`export_doubletfinder.R`), and of the saved objects' RNA counts (`export_doubletfinder_singlets.R`); R is not re-run. `validation/doubletfinder_panel.tsv` lists the 150 objects with their pipeline paths and the workflow job that wrote each saved output: its node and CPU type (from sacct and the node's SLURM features), time and peak RSS. Each Rust run uses that CPU type, with the kernel OpenBLAS detects there (M7).

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
| S5 pANN (sweep and final) | per cell | identical for ≥ 99.9%; otherwise \|Δ\| ≤ 2/k |
| S6 sweep | BCreal; chosen pK | relative ≤ 1e-6; pK identical |
| S7 calls | nExp; Doublet/Singlet set | nExp identical; target identical set, floor ≥ 99.9% agreement with every difference at the pANN cutoff |
| Against saved workflow outputs (T, P, full panel) | doublet set, singlet pANN, pK and nExp as R writes them | set as in S7; pANN as in S5; pK and nExp identical |
| Singlets' counts (M6b: T and P; M7: full panel) | `singlets/` against the saved object's RNA counts: barcodes, feature lines, every stored entry and value | identical |
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
| M6b | count data out (A3): `singlets/`, a new 10x directory, through a standard-library gzip writer; a one-time export of the saved objects' RNA counts | the singlets' counts identical to the saved objects' on T and P in both profiles; independent of the thread count; the memory peak unchanged |
| M7 | full panel: all 150 production objects (sbatch arrays): a one-time export of each object's inputs, compared fields and RNA counts; the Rust run on the CPU type of the object's workflow job | doublet sets, singlet pANN, pK and nExp meet the criteria, and the singlets' counts are identical to the saved objects' (M6b), for all 150, with a summary table |

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

## Part C: Step 2, the clustering crate

Step 2 reproduces the workflow's base clustering (A7b) without Seurat: from a counts matrix to the cluster labels and the workflow's two per-cell scores, silhouette width and neighbourhood purity. It serves two goals:
- **a fast, standalone SCTransform and clustering**, whose outputs can be passed to tools that expect normalised data and cluster labels;
- **annotating clustered cells** with silhouette and purity, to explore the output of any stage of a workflow.

For both, it exports SCTransform's normalised data with the original counts and the per-cell metadata, in a form that can't be mistaken for raw counts (C3f), and provides loaders for Seurat, SingleCellExperiment and AnnData (C3g). A tool that needs clusters can also compute them in-process (C3h).

It reuses `sctransform`, `pca` and `knn` and adds two crates: `neighbors` (Annoy and the shared-nearest-neighbour graph) and `cluster` (Seurat's modularity optimizer, silhouette and purity). C1 and C2 record the code path as researched in the pinned env, the first step of A7's cycle. Milestones M8–M17 (C6) carry out the rest.

### C1. Reference behaviour

The reference is Seurat 5.5.1 as the workflow's scripts call it, with seed 12345 (`set_workflow_seed()`, from `SCRNASEQ_PREPROCESS_SEED`). Around it: SeuratObject 5.4.0, sctransform 0.4.3, glmGamPoi 1.22.0, irlba 2.4.1, RcppAnnoy 0.0.23 (which bundles Annoy 1.17.3) and RcppEigen 0.3.4.0.2 (which bundles Eigen 3.4.0). conda-forge built Seurat and RcppAnnoy with GCC 14.3 and `-O2 -march=nocona -mtune=haswell`. Neither shared library contains FMA or AVX instructions (checked with `objdump`).

**Where the chain runs.** Nine scripts end with the same calls; eight of them ran in production (`posthocfilter_mad.R` is not configured). All 840 Seurat objects saved by the production run went through the chain: rerun #2, the Snakemake runs of 8 and 9 November, in the pinned env. Each object has its job's log under `results/logs/`.

| Script (rule) | The chain's input | Saved objects | Their clusters are read by |
|---|---|---|---|
| `tenx2seuratrds.R` (`tenx2seuratrds`) | Cell Ranger's filtered matrix (`Read10X`, `CreateSeuratObject`) | 30, `seurat_filtered/` | SoupX and decontX (tenx arms) |
| `emptydrops.R` | emptyDrops' called cells (`as.Seurat` of a SingleCellExperiment, which makes a v3 `Assay`) | 30, `emptydrops/` | SoupX and decontX (emptydrops arms) |
| `tenx2seuratrds.R` (`cellbendercells2seurat`) | CellBender's cell calls on the raw counts | 30, `cellbendercells/` | markers and downsampling only |
| `soupx.R`, `decontx.R`, `cellbender2seurat.R` | the corrected counts | 60, 60, 30 | DoubletFinder (the homotypic adjustment) |
| `doubletfinder.R`, `scdblfinder.R` | the singlets (`subset()`) | 150, 150 | markers and downsampling only |
| `posthocfilter_threshold.R` | the cells that pass the thresholds (`subset()`, `DietSeurat()`) | 300 | markers and downsampling only |

- The 60 objects of the first two rows are the base objects that A7b names.
- Every object's clusters also drive its marker genes and the cluster downsampling (FindMarkers and downsampling are out of scope).
- In every object, the chain's input is the RNA `counts` layer as saved. Each script builds or subsets the object first and changes no counts after the chain.

**The calls,** the same in all nine scripts. The post-hoc scripts run the guard before SCTransform; it only stops a run, so the order changes nothing:

```r
obj <- SCTransform(obj, seed.use = 12345, verbose = FALSE)
require_min_cells_for_pca(obj)          # the workflow's guard: stops when there are 50 cells or fewer
obj <- RunPCA(obj, seed.use = 12345, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 12345)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, random.seed = 12345)
obj <- add_silhouette_to_metadata(obj)  # the workflow's own, through bluster
```

1. **SCTransform** follows B3a's code path, with `seed.use = 12345` in every `set.seed()` call.
   - The counts are integers, so glmGamPoi's truncation of half-integers never applies.
   - Objects built with `CreateSeuratObject` hold a v5 `Assay5`. They reach `SCTransform.Assay` through `SCTransform.StdAssay`, which wraps each counts layer in a v3 `Assay`.
   - The emptyDrops objects hold a v3 `Assay`, so they dispatch to `SCTransform.Assay` directly.
   - **The SCT assay** (`SCTransform.default`, then `SCTransform.Assay`):
     - `counts` holds the corrected counts, `data` their `log1p`, and `scale.data` the variable features' residuals, which RunPCA reads;
     - the model goes into the assay's `SCTModel`: feature and cell attributes, the regularised parameters, `median_umi`, the clip ranges and vst's arguments;
     - adding the assay to the object computes `nCount_SCT` (`Matrix::colSums`) and `nFeature_SCT` (the entries > 0).
   - **The corrected counts** come from sctransform 0.4.3's `correct_counts` (`do.correct.umi = TRUE`, the default; C3f). It covers every gene in the model (those detected in at least 5 cells) and every cell. It rounds its results, so they are integers. Both assay paths call it the same way, on all of the cells.
2. **RunPCA** follows B3b's path with `npcs = 50` (Seurat's default), and `set.seed(12345)` before irlba's start vector.
   - irlba runs with `nv = 50` and `work = 66`; for small objects `work` is capped at min(m, n) − 2.
   - The variable features are SCTransform's 3,000, or every kept gene where fewer remain (tier T).
3. **RunUMAP** calls `set.seed(12345)`, then uwot. It adds the `umap` reduction and nothing else, and it advances R's random stream.
   - **It affects nothing downstream.** FindNeighbors reads only the PCA embeddings. Annoy and the modularity optimizer have their own generators. GroupSingletons, the only later use of R's stream, calls `set.seed(1)` first (C2).
   - M8 confirms this with a run without RunUMAP.
4. **FindNeighbors** is called with `dims = 1:30` and its defaults: `k.param = 20`, `prune.SNN = 1/15`, `nn.method = "annoy"`, `n.trees = 50`, `annoy.metric = "euclidean"`, `l2.norm = FALSE`.
   - **The input** is the first 30 columns of the PCA embeddings (`u %*% diag(d)`): cells × 30, doubles. With fewer than k cells, k becomes the number of cells − 1.
   - **`AnnoyBuildIndex`:** `new(AnnoyEuclidean, 30)`; then `addItem(i − 1, row i)` for every cell in order, each double rounded to float; then `build(50)`.
   - **`AnnoySearch`:** for each cell, `getNNsByVectorList(row, 20, search_k = −1, TRUE)`. With `search_k = −1`, Annoy collects candidates until it holds at least 20 × 50 = 1,000. The query is rounded to float again. The `future_lapply` around the loop has `future.seed = FALSE`, so it draws no random numbers.
   - **The `SCT_nn` graph** is `sparseMatrix(i, j = neighbours, x = 1)`.
   - **`ComputeSNN`** (C++): with A the 0/1 kNN matrix, S = A·Aᵀ counts the neighbours s that two cells share. Each stored entry becomes s/(k + (k − s)), and is removed when it is `<` 1/15 (for k = 20, when s ≤ 2). The diagonal stays, with value 1.
   - **The `SCT_snn` graph** is that matrix. Eigen computes the product into row-major storage and converts it, so its row indices are sorted within each column.
5. **FindClusters** is called with `graph.name = "SCT_snn"`, `resolution = 0.8`, `algorithm = 1` (original Louvain), `n.start = 10`, `n.iter = 10`, `modularity.fxn = 1`, `group.singletons = TRUE` and `verbose = TRUE`.
   - **`RunModularityClusteringCpp`** is Seurat's C++ translation of Waltman and van Eck's ModularityOptimizer 1.3.0 (C3c):
     - the network is the strictly lower triangle of the SNN graph, so the diagonal is dropped. Neighbour lists are in ascending order, and each node's weight is the sum of its edge weights;
     - the resolution is rescaled to `0.8 / (2 × total edge weight)`;
     - one Java-style generator, seeded once with 12345, serves all 10 starts. Each start begins from singletons and runs up to 10 Louvain iterations, stopping after an iteration that changes nothing;
     - the start with the highest modularity wins. The comparison is a strict `>`, so the first of equal starts wins;
     - the clusters are renumbered by decreasing size with `stable_sort`, so clusters of equal size keep their order.
   - **`GroupSingletons`:** each cluster of one cell joins the cluster with the highest mean SNN weight to that cell (C3d).
   - **Assembly:**
     - the labels go into the column `SCT_snn_res.0.8`, named by `paste0("res.", r)`. Seurat 5.5.1 sorts its levels numerically;
     - `seurat_clusters` and the identities take the last resolution's labels, with levels sorted numerically.
   - **The job's log** records the number of nodes and edges, the maximum modularity of the 10 starts (`%.4f`) and the number of communities, printed by the optimizer. GroupSingletons adds the number of singletons. 36 of the 840 production logs report singletons, up to 40 (cellbender_fromraw_ghirsutum2).
6. **The workflow's additions,** which are not Seurat:
   - **`add_silhouette_to_metadata`** (`silhouette_utils.R`) takes PCA dims 1:30 and `seurat_clusters`:
     - `bluster::approxSilhouette` gives the `silhouette_width` column (C3d);
     - `bluster::neighborPurity(k = 50, weighted = TRUE, BNPARAM = NULL, BPPARAM = NULL)` gives `neighborhood_purity` (C3e). `BNPARAM = NULL` means BiocNeighbors 2.4.0's default, `KmknnParam()`: an exact Euclidean search by knncolle's KMKNN, on one thread;
     - `neighborPurity` first takes each cell's distance to its 50th nearest other cell (`findDistance`). The median of those distances becomes one radius for all cells, and `findNeighbors` collects every other cell within it;
     - each cell's purity is the weighted share of its own cluster among itself and those neighbours. Each cell weighs 1/(its cluster's size);
   - **`write_cluster_metadata`:** writes the identities' levels into `<object>_cluster_ids.txt` and their count into `<object>_nclusters.txt`, temporary files that drive the marker checkpoint.

**What the step leaves:**
- in the metadata: `nCount_SCT`, `nFeature_SCT`, `SCT_snn_res.0.8`, `seurat_clusters`, `silhouette_width` and `neighborhood_purity`;
- the SCT assay (the corrected counts, their log1p, the residuals and the model), the `pca` and `umap` reductions, and the `SCT_nn` and `SCT_snn` graphs.

The labels with their level order are the result that later steps read.

**The marker step reads the SCT assay.** SCTransform makes SCT the default assay, and no script resets it before saving: the post-hoc scripts set RNA only before rebuilding everything. `find_markers.R` sets `Idents` to `seurat_clusters` and calls `FindMarkers(ident.1 = <cluster>)` with Seurat's defaults, which dispatch to `FindMarkers.SCTAssay`:
- the Wilcoxon test on SCT `data` (through presto, which the pinned env has);
- `avg_log2FC` from log2((rowSums(expm1(data)) + 1)/n) per group, so from the corrected counts;
- `pct.1` and `pct.2`: the fractions of cells with `data` > 0;
- rows sorted by p-value, then by |pct.1 − pct.2|, with exact ties kept in feature order; `p_val_adj` is Bonferroni over all SCT features;
- `recorrect_umi` acts only when the assay holds more than one SCT model, and each object holds one.

The script keeps the rows with `p_val_adj` ≤ 0.05 and writes one file per cluster, an empty one for clusters under 3 cells. `combine_markers.R` joins them into `<prefix>_markergenes.csv`, which the production run kept for all 840 objects. FindMarkers stays out of scope as a port, but it serves as the Seurat loader's acceptance test (C3g).

### C2. Map of random-number use

There are four generators. R's is reseeded before each stage that draws from it, so no stage's draws depend on an earlier stage's:

| Point | Generator | Seeding | Draws |
|---|---|---|---|
| SCTransform | R's Mersenne-Twister | `set.seed(12345)` in each method along the call chain | as in B2 and M4: above 5,000 cells, 5,000 cells without replacement; above 2,000 fitting genes, the density-weighted draw of 2,000; two per 500-gene fitting bin |
| RunPCA | R's | `set.seed(12345)` | irlba's start vector, `rnorm(nfeatures)`; more only in irlba's rare branches (M5) |
| RunUMAP | R's (uwot seeds its own PCG streams from it) | `set.seed(12345)` | not ported; nothing later reads R's stream before `set.seed(1)` |
| Annoy build | `Kiss64Random`, one per index | the default seed, 1234567890987654321, because Seurat never calls `setSeed`. RcppAnnoy builds single-threaded, so one stream runs through all 50 trees | per split attempt: two `index()` draws for the starting pair, then 200 for the two-means iterations, then one `flip()` for each point whose margin is exactly 0; up to three attempts per node; in the fallback, one `flip()` per point per round |
| Annoy search, SNN, `correct_counts`, silhouette, purity's sums | none | | |
| Modularity optimizer | `JavaRandom`, java.util.Random's 48-bit LCG | created once per resolution, with `random.seed` (12345) | each local-moving pass draws a permutation of its network's nodes: n calls to `nextInt(n)`. This happens at every level of every iteration of every start, in one continuing stream |
| GroupSingletons | R's | `set.seed(1)` before each candidate cluster's connectivity | one `sample(<tied clusters>, 1)` per singleton |
| Purity's search index (KMKNN) | `std::mt19937_64`, in knncolle's k-means++ | the library's fixed seed, 6523 | k-means++'s draws for the index's ⌈√n⌉ centres. They decide only where the search looks, not what it finds (C3e) |

Consequences:
- **R's stream needs no tracking past PCA.** RunUMAP's draws reach nothing that is compared, and GroupSingletons reseeds before its draw. R's `.Random.seed` is compared after SCTransform and after RunPCA, as in M6, and after GroupSingletons for objects with singletons.
- **Annoy's draws depend on the data.** Margins of exactly 0 and failed splits consume extra draws, so one miscounted draw changes every later tree. The saved index (C4) checks the whole stream.
- **`nextInt` rejects as Java does.** Its loop `while (bits - val + (n - 1) < 0)` relies on 32-bit signed wrap-around, which is undefined behaviour in C++. The compiled `Seurat.so` keeps the check (an `add` and a `js` on 32-bit registers), so the port implements Java's semantics.
  - When n is a power of two, `nextInt` takes Java's multiply-and-shift branch instead. Reduced networks can be that small.
- **Each resolution starts the optimizer afresh** from `random.seed`. Resolutions are therefore independent and can run in parallel. Within a resolution, the 10 starts are sequential.
- **KMKNN's generator needs no port,** unless C3e's certificate finds a cell at a rounding boundary (M12's decision point).

### C3. Rust design

**Commands.** The four staged subcommands that A7b asks for, each also a library call:

```
scrust sctransform --counts <10x dir> --out <dir> [--seed 1448145] [--corrected-counts] ...  (exists, M4)
scrust pca         --sct <dir> --out <dir> [--npcs 50] [--seed 42]                         (exists, M5)
scrust neighbors   --pca <dir> --out <dir> [--dims 1:10] [--k 20] [--trees 50]
                   [--metric euclidean|cosine|manhattan] [--prune 1/15] [--search-k -1]
scrust cluster     --snn <dir> [--pca <dir>] --out <dir> [--resolution 0.8[,1.2,...]] [--algorithm 1|2|3]
                   [--starts 10] [--iterations 10] [--seed 0] [--modularity 1|2] [--no-group-singletons]
                   [--stats silhouette,purity|none] [--purity-k 50] [--purity-weighted true|false]
```

The whole chain also runs in one process, for the workflow:

```
scrust cluster --counts <10x dir> --out <dir> [--preset workflow] [--params <file.toml>]
               [--export] [--keep-stages] [--metadata <tsv>] [--threads N]
```

`scrust cluster` can start from any stage's output (`--sct`, `--pca`, `--snn`), so a sweep over k or the resolution reuses one SCTransform and PCA. The scores need the PCA, so starting from `--snn` they need `--pca` too.

**Parameters.**
- **Every Seurat argument listed in A7b** is a parameter, with Seurat's default. It is set by a flag or in a TOML file, read by a small parser in `scio` (standard library only).
- **`--preset workflow`** sets the workflow's calls: seed 12345 in SCTransform, RunPCA and FindClusters, and `dims = 1:30`.
- **Other code paths are refused.** Options that would take Seurat down another code path, such as `l2.norm`, `nn.method = "rann"` or `"hnsw"`, Leiden and `vars.to.regress`, are rejected with an error, as `SctParams::check` does now.
- **Fractions are divided.** A value written as `1/15` is parsed and divided, so it is R's double.
- **The scores** (`--stats`, both by default, as in the workflow) use the neighbour search's PCA dims, so 1:30 under the preset. Purity takes bluster's defaults, k = 50 and weights.

**Outputs** (plain files, A3):
- `seurat_clusters.tsv` (barcode and label, in the counts' column order) and `seurat_clusters_levels.txt`: the layout `scrust doubletfinder --clusters` reads, so the tools chain without R;
- `clusters.tsv`: one column per resolution, named as Seurat names it (`SCT_snn_res.0.8`), with each column's levels;
- `cluster_ids.txt` and `nclusters.txt`, as `write_cluster_metadata` writes them;
- `metadata.tsv`: one row per cell, in the counts' column order, with the columns the workflow's objects hold:
  - `nCount_RNA`, `nFeature_RNA`, each `res.` column, `seurat_clusters`, `silhouette_width` and `neighborhood_purity`;
  - `nCount_SCT` and `nFeature_SCT` with `--export`, since they come from the corrected counts;
  - bluster's other outputs, which the workflow drops: `silhouette_other` (the closest other cluster) and `purity_maximum` (the cluster with the largest share);
  - the columns of `--metadata`, carried over by barcode;
- `summary.json`, in three sections:
  - `result`: the network's nodes and edges, every start's modularity, the best start and its modularity as a double, the communities, the singletons and the levels;
  - `run`: parameters, the versions reproduced, timings, peak RSS and the environment fingerprint;
  - `files`: the manifest of C3f;
- stage outputs, in subdirectories:
  - `pca/` and `neighbors/`, always: they are small, and the loaders need them for the reductions and the graphs. `neighbors/` holds `knn_index.tsv` and `knn_distance.tsv` (each cell's k neighbours, 1-based as in R, and their float distances), and `nn.mtx.gz` and `snn.mtx.gz` (Matrix Market, entries in stored order, values as `%.17g` writes them);
  - `sct/`, in `scrust sctransform`'s layout, with `--keep-stages` or `--export`. Its residuals are large: 1.9 GB on the 79,409-cell objects;
  - on request, `annoy.ann`: the index in Annoy's own file format;
- with `--export`, the rest of the export (C3f): the input counts, the corrected counts and their log1p, and a README.

The clustering changes no counts, so it writes no new count data (A3). The export copies the input counts unchanged.

**Components:**
- **`neighbors`** (new): Annoy as RcppAnnoy builds and searches it (C3a), the kNN graph and `ComputeSNN` (C3b). It is ported to Rust rather than called through FFI, as decided (C3a).
- **`cluster`** (new):
  - the modularity optimizer (C3c), FindClusters' result (C3d), `approxSilhouette`, and `neighborPurity`'s weights and sums (C3e);
  - the whole chain as one library entry point, for `scrust cluster` and for clustering on the fly: in DoubletFinder now, in the decontX and SoupX ports later (C3h);
  - Leiden (`algorithm = 4`, which goes through leidenbase) can follow later; it is not in Step 2.
- **`sctransform`:** used as it is, with the workflow's parameters. New: `correct_counts` (C3f).
- **`pca`:** used as it is.
- **`knn`:** new: the search `neighborPurity` makes, with knncolle's arithmetic (C3e).
- **`rcompat`:** `set.seed` and `sample` for GroupSingletons; long-double `mean` and `sum`, whose accumulator R's `colMeans`, `colSums` and `rowSums` also use (small additions); R's `median`; R's number formatting for the `res.` names; and `order`.
- **`scio`:** 10x reading, TSV, `.npy` and JSON. New: the TOML reader; Matrix Market writers for real and integer matrices, with comment lines; the manifest; copying a directory with its checksums.
- **`scrust`:** the commands above, `scrust validate cluster`, and the new options of `scrust doubletfinder` (C3h).
- **`loaders/`** (new, outside the Cargo workspace): the R and Python packages (C3g).

**What the existing crates cover with the workflow's settings:**

| Stage | The workflow's setting | Validated in Step 1 | Left to do |
|---|---|---|---|
| SCTransform | `seed.use = 12345`; integer counts; v5 and v3 assays | `seed.use = 1448145`; counts with half-integer doublets; v5 assays (168 calls, M4) | the seed is already a parameter (`SctParams::seed_use`, `--seed`); validate seed 12345 and the v3 path of the emptyDrops objects |
| RunPCA | `npcs = 50`, `seed.use = 12345` | `npcs = 10`, seed 42 (168 calls, M5); `scrust pca` with 50 components on one tier-T object (M5) | both are already parameters (`PcaParams`, `--npcs`, `--seed`); validate 50 PCs at production scale (`work = 66`, more restarts) |
| Neighbours | Annoy, k = 20, 30 dims | DoubletFinder's exact distances (`knn`, `fields::rdist`), a different algorithm | everything: `neighbors` |
| Clustering | Louvain | none | everything: `cluster` |
| Silhouette and purity | `approxSilhouette`; `neighborPurity` through KMKNN | `knn`'s exact distances, the same sums as knncolle's (`fields::rdist`'s order), with DoubletFinder's selection | both scores: `cluster`; knncolle's squared comparisons, the 50th distance and the range query: `knn` (C3d, C3e) |
| Corrected counts | `correct_counts`, as `do.correct.umi = TRUE` runs it | none: DoubletFinder reads only the residuals | the port (C3f): the residual code's products, plus the same at the median depth |
| R semantics | `set.seed(1)` and `sample`, long-double `mean`, `round`, number formatting | all exact (M1) | `median`, `colSums`, `rowSums` (small additions) |

The `knn` crate doesn't replace Annoy: Annoy's search is approximate, and Seurat's result depends on its exact approximations. `knn` serves purity, whose search is exact.

**C3a. Annoy (FindNeighbors).** This is a port of Annoy 1.17.3's `AnnoyIndex<int32_t, float, Euclidean, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy>`, the type RcppAnnoy 0.0.23 instantiates (as `RcppAnnoy.so`'s symbols show).
- **Why a port, not FFI:**
  - the arithmetic is plain IEEE single precision with no BLAS. There is no vector code (Annoy's AVX paths are compiled out under `-march=nocona`), no FMA, and every reduction is a sequential loop. Scalar Rust reproduces it exactly, on any machine;
  - RcppAnnoy can save R's index to a file, so the port's forest is compared byte for byte (C4);
  - FFI would need a C++ compiler in `envs/rust-dev.yml`, a build step and vendored C++ sources. The workspace is now the standard library plus OpenBLAS's C interface;
  - FFI stays the fallback, should the port fail to be byte-identical.
- **Items:** each embedding row is rounded to `f32`, to nearest even, as C++'s conversion and Rcpp's `as<std::vector<float>>` do. Items are added in cell order.
- **Nodes** have `16 + 4f` bytes: `n_descendants`, the offset `a`, two children and the f-vector. A leaf holds up to K = f + 2 item ids, 32 for 30 dims.
- **Node ids:**
  - the items take 0…n−1;
  - then, tree by tree, each node takes the next id when it is created. A leaf is created when it is reached, and a split node after both its children. The smaller child is built first. A subset of one item is that item, not a new node;
  - the last 50 ids are copies of the roots.
- **The split** (`two_means`, Euclidean):
  - it picks two distinct random items, then 200 random items. Each is assigned to the nearer centroid, distances weighted by the centroid's count (`ic * dist`), and that centroid becomes `(p*ic + x/1)/(ic + 1)` in float. On a tie, neither moves;
  - the normal is p − q, normalised. The offset accumulates `a += -v*(p + q)/2` in float;
  - each point goes left or right by the sign of `a + v·x`, with a `flip()` at exactly 0;
  - there are up to three attempts while the larger side holds 95% or more. `_split_imbalance` mixes float and double, as written;
  - if a split stays above 99%, the sides are drawn at random until it is not. The node keeps its last `a`, with a zero normal.
- **The search:**
  - a max-priority queue holds pairs of (float priority, node id), with the roots at +∞. Children enter at `min(d, ±margin)`;
  - leaves add all their items, until at least 1,000 candidates have been collected. The candidates are sorted by id and deduplicated, and their distances computed in float;
  - `partial_sort` by (distance, id) keeps the 20 smallest. Distances come back as `sqrt(max(d, 0))`, computed in float.
- **Metrics:** Euclidean for the workflow. Angular (`annoy.metric = "cosine"`) and Manhattan serve the parameter grid.
  - Seurat means to convert Angular distances to cosine distances, but `res$dist <- 0.5 * (res$dist * res$dist)` adds a new list element, and the distances it returns are the Angular ones. The port returns those too.
  - The distances feed no graph, so this changes no label.
- **File format:** `annoy.ann` is the node array, as Annoy's `save()` writes it.
- **Threads:** the build is sequential, because one generator runs through all the trees. Queries run in parallel over cells, each result placed by its index.

**C3b. The SNN graph.**
- For each cell and each of its neighbours, the number of neighbours they share is an exact integer.
- The Jaccard value s/(k + (k − s)) is computed in that order, in double, and removed when it is `< prune`.
- The graph is stored as sorted compressed columns, diagonal included.
- The `SCT_nn` graph is the 0/1 kNN matrix.

**C3c. The modularity optimizer** is ported statement by statement from Seurat 5.5.1's `ModularityOptimizer.cpp` and `RModularityOptimizer.cpp`:
- **The network** is built from the SNN graph's strictly lower triangle (`matrixToNetwork`). Modularity function 1 weights each node by its total edge weight; function 2 weights each node by 1 and leaves the resolution unscaled.
- **Local moving:**
  - the nodes are visited in a random permutation, cycling until every node is stable;
  - each node joins the neighbouring cluster with the largest gain, `w(j→l) − k_j·K_l·γ`, evaluated left to right. Equal gains go to the lowest cluster index;
  - a node whose best gain is 0 moves to the last unused cluster;
  - at the end of the pass, the empty clusters are dropped and the others renumbered in increasing order of their ids.
- **The algorithms:**
  - Louvain (algorithm 1): local moving, then the same on the reduced network while it shrinks (`createReducedNetwork`), with the result merged back;
  - multilevel refinement (2);
  - smart local moving (3), through `createSubnetworks`. As written, its loop always runs all its iterations.
- **The rest:** `calcQualityFunction` with the same summation order; the best start, by strict `>`; and `orderClustersByNNodes`.

**C3d. FindClusters' result.**
- **GroupSingletons:**
  - singletons are handled in order of their first appearance among the cells (`intersect(unique(ids), ...)`), and the candidate clusters likewise (`setdiff`);
  - a candidate's connectivity is the long-double `mean()` of the singleton's SNN row over the candidate's current cells, so singletons merged earlier count;
  - the highest connectivity wins. Ties go to `sample(names, 1)` after `set.seed(1)`;
  - a cell with no SNN edges ties at 0 with every cluster, so that draw places it;
  - with `group.singletons = FALSE`, the singletons are labelled `"singleton"` instead.
- **Labels and levels:**
  - the labels are the optimizer's ids, as character once a singleton has been merged;
  - in each `res.` column, numeric labels are sorted numerically and come first, then any others alphabetically (Seurat 5.5.1);
  - `seurat_clusters` takes the last resolution's labels, with levels sorted numerically.
- **Column names:** the graph name, then `res.`, then the resolution as R formats numbers (15 significant digits).
- **Numbering:** clusters are numbered by size before the singletons are merged. The final labels are therefore not always in decreasing size, and the port keeps that.
- **`approxSilhouette`:**
  - centroids and within-cluster variances come from R's `colMeans` and `sum` in long double;
  - each cell's root-mean-square distance to every cluster takes `colSums` over the transposed embeddings;
  - the nearest other cluster is chosen by strict `<`, in level order;
  - the width is (other − own)/max(other, own);
  - its `other` column, the nearest other cluster, is written as `silhouette_other`.

**C3e. Neighbourhood purity,** as the workflow calls `bluster::neighborPurity` (1.20.0) on BiocNeighbors 2.4.0.
- **The search is exact, so only its arithmetic matters, not its algorithm.**
  - knncolle computes each Euclidean distance squared, in double, as `out += delta * delta` over the dims in order, and compares the squares.
  - Any exact search that computes and compares distances in that order finds the same neighbours. `knn`'s distance kernel already sums in that order (`fields::rdist`'s); only its selection changes.
- **Steps:**
  1. each cell's distance to its 50th nearest other cell (`findDistance`): the 50th smallest squared distance to the other cells, then its square root. The cell itself, at distance 0, is excluded;
  2. the radius: R's `median` of those n distances. With an even n, it is the long-double mean of the middle two;
  3. each cell's neighbours (`findNeighbors`): every other cell whose squared distance is `<=` the radius times itself, in double (knncolle's `denormalize`). Squaring a square root doesn't always give back the starting value, so the comparison is computed exactly as written;
  4. the weights: 1/(cluster size) per cell (`1/table(clusters)`);
  5. per cell, one sum per cluster (`sum_neighbor_weights`): the cell's own weight, then each neighbour's weight added to its cluster's sum. Within a cluster every added weight is the same, so the order of the neighbours doesn't change the sums;
  6. the total: R's `rowSums` over the clusters in level order, in long double. The purity is the own cluster's sum divided by the total. A cell always counts itself, so the total is never 0.
- **Memory and time:** there is no n × n matrix. As in DoubletFinder's search (B3), each cell's distances to all cells go into a buffer per thread. Step 1 takes the 50th from it, and step 3 counts the cells within the radius by cluster. That makes two passes over all pairs: about 4·10¹¹ multiply-adds on the 79,409-cell objects, spread over threads.
- **Where KMKNN could differ.** KMKNN skips groups of cells by the triangle inequality, with bounds computed in floating point.
  - In exact arithmetic those bounds never skip a cell that belongs in the result. With rounding, they can skip one only if its distance is within a few units in the last place of the radius (step 3) or of the cell's 50th distance (step 1).
  - The Rust search reports every pair that close, as a certificate. If none is reported, the result is identical to KMKNN's, whatever KMKNN's k-means did.
  - If some are reported, M12's decision point applies: port knncolle's KMKNN as BiocNeighbors 2.4.0 builds it, or agree with you how to treat the flagged cells. The KMKNN port means ⌈√n⌉ centres, k-means++ seeded with 6523 through `std::mt19937_64`, Hartigan–Wong refinement, and the bounds as written; the upper bound is compiled out (`KNNCOLLE_KMKNN_USE_UPPER` is undefined).
- **`neighborPurity`'s other output,** `maximum`, is the cluster with the largest share, ties to the first (`max.col`). It is written as `purity_maximum`.

**C3f. Corrected counts and the export.**

*The corrected counts* port sctransform 0.4.3's `correct_counts`, as `SCTransform.default` calls it, with no `scale_factor`:
- **The regressors:** each cell's `[1, log_umi]`, and the same with `log_umi` replaced by its median over all cells (R's `median`).
- **For each bin of 500 model genes,** in the counts' gene order:
  - μ is exp of the regularised parameters times the regressors (`tcrossprod`), through the same OpenBLAS `dgemm` call R makes, with the same shapes;
  - the variance is μ + μ²/θ, which is μ where θ = ∞;
  - the Pearson residual is (y − μ)/√variance. Unlike vst's residuals, it is neither clipped nor floored at `min_variance`;
  - μ and the variance are computed again at the median depth, giving μ + residual·√variance.
- **Rounding:** `round` (half to even), negatives set to 0, zeros dropped. The results are integer counts, stored sparse.
- **The rest of the assay:** `data` is their `log1p` (glibc's), and `nCount_SCT` and `nFeature_SCT` are their column sums and nonzero counts.
- They are computed only for the export, because clustering doesn't need them.

*The export* keeps the original counts, the normalised data and the per-cell metadata together, as plain files. `scrust cluster --export` writes all of it; `scrust sctransform --corrected-counts` writes the SCT part alone:

```
<out>/
  README.txt                       what each file holds and how to load it (C3g)
  summary.json                     the manifest: every file's role, scale, shape, and row and column identifiers;
                                   provenance (scrust version, command line, the input's path and sha256,
                                   the R functions reproduced, parameters, seeds)
  metadata.tsv                     one row per cell (C3's outputs)
  counts/                          the input counts, copied byte for byte: the export's only 10x directory
  sct/                             the SCT assay, never in 10x layout:
    corrected_counts.mtx.gz        Seurat's SCT `counts` (integer Matrix Market)
    log1p_corrected_counts.mtx.gz  Seurat's SCT `data`
    scale_data.npy                 Seurat's SCT `scale.data`: the variable features' residuals
    sct_features.tsv.gz            the model genes (both matrices' rows), in the counts' order
    sct_cells.tsv.gz               the cells (their columns)
    cell_attr.tsv, gene_attr.tsv, model_pars.tsv, model_pars_fit.tsv, variable_features.tsv
                                   the model, in `scrust sctransform`'s layout
  pca/, neighbors/                 as in C3's outputs
  seurat_clusters.tsv, clusters.tsv, ...   as without --export
```

*How it guards against being taken for counts:*
- **No 10x names for normalised data.** No directory holding normalised data has a `matrix.mtx`, `barcodes.tsv`, `features.tsv` or `genes.tsv` file, and no normalised file's name ends in `matrix.mtx(.gz)`. Seurat's `Read10X` and scanpy's `read_10x_mtx` (even with its `prefix=` argument) therefore fail on `sct/` and on the export's top level. Tests check both.
- **Each normalised file says what it holds.** Comment lines right after the Matrix Market banner state what the matrix is, that it is not raw counts, how it was made (the scrust version and command) and the input's sha256. Matrix Market readers skip comments, and Cell Ranger's own matrices carry one, so the files stay readable.
- **The manifest gives every file a role and a scale** (`raw_counts`, `corrected_counts`, `log1p`, `residuals`, …). The loaders put a matrix only into a slot meant for its scale, and refuse a file whose role doesn't fit.
- **The only 10x directory is `counts/`,** which holds the raw counts, as A3 requires.

*One manifest for every command.* From M13 on, every scrust command's `summary.json` gains a `files` section in this schema, `scrust doubletfinder`'s included (`cells.tsv`, the sweep tables, `singlets/`). One loader per language therefore reads any scrust output.

*Why plain files rather than h5ad written by Rust:*
- **Plain files keep A3's formats and the workspace's dependencies** (the standard library and OpenBLAS). Writing h5ad needs an HDF5 writer: libhdf5 linked through FFI, or our own writer for a complex format. Either way, anndata's on-disk spec has to be tracked as it changes.
- **A Seurat object can't come from Rust anyway.** Rebuilding a real SCT assay needs R and Seurat's classes, so an R loader is needed whatever the container.
- **h5ad and `.rds` files come from the loaders when wanted** (`write_h5ad`, `saveRDS`), each also as a one-line command.
- Reading 10x HDF5 (A3, later) is a separate, smaller matter.

**C3g. Loaders.** Small packages in the repo, `loaders/R/scrustio` and `loaders/python/scrustio`, read any scrust output through its manifest:

| Object | Raw counts | SCT corrected counts and their log1p | Residuals (variable features only) | Metadata, PCA, graphs |
|---|---|---|---|---|
| Seurat (R) | `RNA` assay, `counts` layer | an `SCTAssay` named `SCT`: `counts` and `data`, with the `SCTModel` rebuilt from the model files as SCTransform makes it | the SCT assay's `scale.data` | `meta.data`, with factors in their levels' order; the `pca` reduction; the `SCT_nn` and `SCT_snn` graphs; identities from `seurat_clusters`; default assay `SCT` |
| SingleCellExperiment (R) | the `counts` assay | `altExp(sce, "SCT")`, with assays `corrected_counts` and `logcounts` (the model genes only) | `altExp(sce, "SCT_scale_data")`, assay `scaledata` | `colData`; `reducedDim(sce, "PCA")`, with the loadings as its `rotation`; the `colPairs` `SCT_nn` and `SCT_snn` |
| AnnData (Python) | `layers["counts"]` | `X` = their log1p (Seurat's `data`), with scanpy's `uns["log1p"]`; `layers["sct_corrected_counts"]`. Genes outside the model are 0 there, flagged by `var["sct_modelled"]` | `obsm["sct_scale_data"]`, a data frame with the genes as columns | `obs`, with categoricals in their levels' order; `obsm["X_pca"]`, `varm["PCs"]`, `uns["pca"]`; `obsp["sct_nn"]` and `obsp["sct_snn"]`, with `uns["sct_snn"]` for scanpy's `neighbors_key` |

- **The rule:** normalised data goes only into slots meant for normalised data. No loader puts corrected counts into a slot named `counts`, and the raw counts always come from `counts/`.
- **The Seurat loader's acceptance test:** FindMarkers on a loaded object must reproduce the workflow's marker tables exactly (C1, K13). That holds only if the loader:
  - puts the corrected counts in SCT `counts` and their log1p in `data`;
  - keeps the model genes in the counts' order, since exact ties in FindMarkers' sorted output keep the feature order;
  - builds exactly one `SCTModel`;
  - gives `seurat_clusters` its saved level order, which the combined table follows.
- **The residuals** cover only the variable features, so they go where a different set of genes is allowed: Seurat's `scale.data`, an alternative experiment, an `obsm` table.
- **Other commands' outputs:** `scrust doubletfinder`'s calls become metadata columns, and `singlets/` loads as raw counts.
- **Packaging:** an R package (Matrix, SeuratObject and jsonlite, plus Seurat or SingleCellExperiment for those objects) and a Python package (numpy, scipy, pandas, anndata). Each has a one-line command that writes `.rds` or `.h5ad`. Both are tested in hermetic envs (A2).
- **Versions:** the Seurat loader builds Seurat 5.5.1's and SeuratObject 5.4.0's classes, and checks the installed versions. With other versions it carries on and warns that they are untested.

**C3h. Clustering on the fly.** The crates are libraries: `doubletfinder` already calls `sctransform` and `pca` in-process. The chain becomes a library entry point in `cluster` too. It returns the labels with their levels and, if asked, silhouette and purity.
- **An entry point, not tool code:** `doubletfinder` calls it in Step 2, and the decontX and SoupX ports will call it in their own steps. `scrust doubletfinder`'s options are a thin layer over it.
- **The hook:** `scrust doubletfinder` given no `--clusters` clusters its own input. It takes the same parameters and presets as `scrust cluster`, through `--cluster-preset` or `--cluster-params <file.toml>`. The preset defaults to `workflow`, as the tool's other defaults are the workflow's.
- **`--clusters` stays,** for three reasons: the original tools take labels; users may have curated annotations; and in the workflow, SoupX and decontX share one base object's labels.
- **Memory:** only the labels are kept. Each stage is freed once the next has what it needs, and the corrected counts are never computed. The chain runs before DoubletFinder's own work, and DoubletFinder's SCTransform on more cells (the doublets added) should still set the peak. M15 measures it.
- **Outputs:** the labels it computed, in `--clusters`' layout (`clusters/seurat_clusters.tsv` and the levels); silhouette and purity if asked (`--cluster-stats`); the clustering's numbers in `summary.json`.
- **The check:** the workflow clustered DoubletFinder's 150 inputs with this same chain. So `scrust doubletfinder` clustering on the fly must reproduce M7's 150 saved outputs exactly, with no new R runs.
- **Later tools:** the SoupX and decontX ports call the same entry point in their own steps. For decontX it is a requirement (A7b): its port takes `--clusters`, clustering on the fly, or decontX's own clustering (UMAP, then DBSCAN) with the workflow's fallback when there are fewer than two clusters. Porting decontX's own clustering, and choosing the default, belong to decontX's planning, not to Step 2.

**C3i. What exactness requires.** Each of these hazards is tested on purpose:
- **Float conversion:** R rounds the embeddings to float twice, at `addItem` and at the query. The port rounds once and reuses the result, which is the same value. Float distances convert back to double exactly.
- **Float arithmetic in Annoy:** the operation order of `two_means`, `create_split`, `normalize`, the margins and the distances, including the conversions from integer to float, and `sqrt` in float.
- **Tie orders:**
  - the queue pops equal priorities by the larger node id; +0 and −0 compare equal;
  - `partial_sort` keeps equal distances by the lower item id. Duplicated cells (identical embeddings) produce such ties, and a cell can then be missing from its own neighbour list, as in R;
  - the optimizer breaks equal gains by the lowest cluster index, and equal modularities by the first start;
  - `stable_sort` numbers equal-sized clusters in their previous order;
  - GroupSingletons breaks ties by R's random draw, not by order.
- **Summation orders in the optimizer:** the node weights; the cluster weights, updated by `+=` and `−=` in visit order; the reduced network's edge weights; the modularity. A difference in the last bit can flip a near-tie, so each sum is copied as written, never rewritten.
- **Cluster-ID ordering:** numbering by size before the singletons are merged; numeric level sorting; character labels after a merge.
- **PCA's last bits** depend on the OpenBLAS kernel (A5b). Rounding to float absorbs most differences, but not all. End-to-end runs therefore use the reference's profile, and the full panel uses each job's CPU type.
- **The emptyDrops objects' v3 assay** takes `SCTransform.Assay` directly. M8's instrumented runs check that it gives the same RNG checkpoints and outputs as the v5 path.
- **Purity:** the squared distances' summation order; the exact 50th distance, never an approximation; R's `median`; the radius squared in double and compared with `<=`; the long-double `rowSums`; and C3e's certificate.
- **Corrected counts:**
  - the same `dgemm` calls as R, with the same bin shapes. The kernel's last bits can move a value across a .5 boundary, so exact checks use the reference's kernel;
  - residuals that are neither clipped nor floored; θ = ∞; `round` half to even; glibc's `exp` and `log1p`.

### C4. Validation data

**Reference runs** (A4) use the workflow's calls and seed, in the pinned env:
- **`export_clustering.R`:**
  - from each saved object: the RNA counts layer as stored (the chain's input); `seurat_clusters` and `SCT_snn_res.0.8` with their levels; `nCount_SCT`, `nFeature_SCT`, and the silhouette and purity columns;
  - for T and P also: the SCT assay's three layers and its `SCTModel`'s slots; the PCA (embeddings, loadings, stdev); the variable features; and both graphs' slots as stored;
  - existing exports are reused. M2 and M7 already hold the counts and `seurat_clusters` of the 150 DoubletFinder inputs, and M6b and M7 the counts of the 150 DoubletFinder outputs.
- **`reference_clustering.R`:** the direct run, on an export, with the workflow's calls, UMAP included. It builds the object as the workflow's script did (a v3 `Assay` for the emptyDrops objects) and writes complete outputs as plain files. Each run must reproduce its object's saved labels and levels, silhouette widths, purities and SCT layers.
- **`instrumented_clustering.R`:** the same run, also dumping every stage:
  - SCTransform, as M2's S2a–S2d;
  - the PCA and the RNG checkpoints;
  - the float items (`getItemsVector`);
  - the index file, saved with `save()` from a second, identical build: Annoy's `save()` reloads the file and re-derives the roots;
  - every cell's neighbours and float distances;
  - both graphs;
  - purity's intermediates: each cell's 50th distance, the radius, and each cell's neighbour count; for T and the small P object, the neighbour lists too;
  - `correct_counts`' median depth.

  A variant without RunUMAP shows that UMAP changes nothing downstream.
- **An instrumented optimizer.** R returns only the final labels, but Seurat's `ModularityOptimizer.cpp` has a standalone mode with its own `main()`.
  - It is compiled in the R env, with the env's compiler and `Makeconf` flags, and patched only to print each start's modularity (`%.17g`) and partition.
  - It reads the SNN graph as an edge list written with `%.17g`, and must first reproduce R's labels.
  - It also supplies the `JavaRandom` test vectors. A small harness around RcppAnnoy's header supplies the `Kiss64Random` vectors.
- **The R-vs-R floor:** the same kernel on two CPU types, and the Haswell and SapphireRapids kernels on one CPU type, on T and P. It answers how often PCA's last bits change the float items, the neighbours and the labels. That decides what portable-profile results can claim for production objects.
- **CPU types:** `validation/clustering_panel.tsv` is written from the Snakemake logs, sacct and scontrol, as M7's table was. For every object it lists:
  - the input and saved paths;
  - the job that wrote the object: node, CPU type, time and peak RSS;
  - the optimizer's numbers from the job's log.

  Production-profile runs use each job's CPU type.

**Tier T, the workflow's test data:** the 27 test objects in `tests/reference_outputs/`, each with at most 300 cells and 2,000 genes, and each clustered by the chain:
- 2 base objects (tenx and emptydrops);
- the 5 DoubletFinder inputs, whose counts and labels are already M2 fixtures;
- 5 DoubletFinder outputs and 5 scDblFinder outputs;
- 10 post-hoc objects.

`cargo test` runs the 5 existing fixtures plus the 2 base objects. The other 20 live in `data/`, listed in the manifest, which keeps the committed fixtures near A1's 20 MB. Limits of tier T: SCTransform subsamples neither cells nor genes, every kept gene is a variable feature, and the trees are shallow.

**Tier P:**

| Role | Object | Cells | Why |
|---|---|---|---|
| small | soupx_tenx_SRR18256555 | 3,108 | Step 1's P object (export exists); SCTransform's gene subsample; grid sample 1 |
| medium | soupx_tenx_10xmousebrain | 11,831 | Step 1's; SCTransform's 5,000-cell subsample; grid sample 2 |
| base, tenx | filtered_seurat_tenx_10xmousebrain | 11,831 | a base object, from Cell Ranger's matrix |
| medium 2 | soupx_emptydrops_4500871lane | 26,598 | Step 1's; one singleton |
| large | soupx_emptydrops_sscovelli | 79,409 | Step 1's; 18 singletons; the memory and time benchmark |
| base, emptyDrops | filtered_seurat_emptydrops_sscovelli | 79,409 | the v3 assay path; 22 singletons |
| many singletons | cellbender_fromraw_ghirsutum2 | 59,323 | Step 1's; 40 singletons, the most in the panel; 99,452 features |

**Loader tests:** each loader's object is compared with the saved object: the counts, the SCT layers and model, the chain's metadata columns, the PCA and the graphs.
- The R loaders (Seurat, SingleCellExperiment) are compared with the `.rds` directly, in batch jobs.
- The Python loader is compared with T's and P's exports.
- **Markers (the Seurat loader):** `markers_from_loader.R` runs the steps of the workflow's `find_markers.R` and `combine_markers.R` on the loaded object, cluster by cluster in level order, in the pinned env.
  - The steps are copied into `validation/R/`, with their source cited; the pipeline repo's scripts are not run from there.
  - The combined table must equal the saved `<prefix>_markergenes.csv` byte for byte. The references are the production run's tables for P, and the test data's 27 in `tests/reference_outputs/` for T.
  - M8 first checks the references: FindMarkers on each saved `.rds` must reproduce its saved table. Where one doesn't, FindMarkers on the saved `.rds` becomes the reference.
  - It runs on T and P, as sbatch arrays. On the full panel, the loader's comparison of the SCT layers, the model and the identities already covers everything FindMarkers reads, so FindMarkers is not rerun on all 840.

**On the fly:** `scrust doubletfinder` given no `--clusters` runs on the 150 DoubletFinder inputs (M2's and M7's exports). R is not run.
- The outputs are compared with the saved outputs, as in M7, and file by file with M7's own Rust outputs (label `m7-b8c8a8c`).
- The labels it computed are compared with the inputs' saved `seurat_clusters`.

**Full panel:** all 840 objects the production run clustered (C1's table). Each is compared with its saved labels and levels, its silhouette widths and purities, and the optimizer's numbers in its log.
- Full exports are needed for 540 objects: the 90 base objects, the 150 scDblFinder outputs and the 300 post-hoc objects.
- The 300 DoubletFinder inputs and outputs already have their counts exported (M2, M6b and M7). A small export adds the outputs' labels, and the silhouette widths and purities of all 300.
- **The SCT layers, the PCA and the graphs** are compared through the R loader, in the same batch job. The Rust run writes its export (`--export`), the loader reads it, and the comparison reads the saved `.rds`. The export's large files are then deleted; the comparison tables are kept.
- R is not re-run.

**The parameter grid** (A7b) compares the Rust labels with Seurat run directly with the same settings. `grid_clustering.R` runs one SCTransform and PCA per sample and setting, FindNeighbors per neighbour setting, and all resolutions in one FindClusters call. The samples are the base tenx test object and grid samples 1 and 2:
- k ∈ {10, 20, 30, 50} × resolution ∈ {0.2, 0.4, 0.8, 1.2, 2.0} × dims ∈ {1:10, 1:30}: 40 label sets per sample;
- one parameter at a time, around the workflow's settings:
  - neighbours: trees {10, 100}; metric cosine and manhattan; prune {0, 0.2};
  - clustering: algorithms 2 and 3; starts {1, 20}; iterations {1, 50}; FindClusters seeds {0, 1}; modularity function 2; `group.singletons = FALSE`;
  - SCTransform and PCA: npcs {30, 100}; 2,000 and 5,000 variable features; 2,000 fitting cells; Seurat's default seeds.

### C5. Acceptance criteria (in `tolerances.toml`)

| Stage | Check | Criterion |
|---|---|---|
| K0 generators | `Kiss64Random` and `JavaRandom` vectors from the original code; `set.seed(1)` and `sample` | exact |
| K1 SCTransform | S2a–S2d at seed 12345, v5 and v3 assays | as B5's S2a–S2d: exact, with S2d under the reference's kernel |
| K2 RunPCA | 50 PCs: sdev, embeddings, loadings; `.Random.seed` afterwards | identical under the reference's kernel; RNG state always exact |
| K3 Annoy items | the first `dims` embedding columns as float, given R's embeddings | exact |
| K4 Annoy index | `annoy.ann` against R's saved index, given R's embeddings | byte-identical |
| K5 neighbours | every cell's k ids and float distances, given R's embeddings | identical |
| K6 graphs | `SCT_nn` and `SCT_snn`: dimensions, stored entries in order, values bitwise | identical |
| K7 clustering | given R's SNN graph: the labels, the levels and every `res.` column; the nodes, edges, communities and singletons; the maximum modularity (`%.4f`, as logged) and every start's modularity as a double (against the instrumented optimizer) | identical |
| K8 end to end | from the counts: K1–K7, the labels and the silhouette widths against the instrumented runs; the files against the direct runs | identical under the reference's kernel |
| K9 purity | given R's embeddings and labels: each cell's 50th distance, the radius, each cell's neighbour count and purity; the certificate | identical; no flagged pair (otherwise M12's decision point) |
| K10 corrected counts | given R's model parameters: every corrected count, the log1p data, `nCount_SCT` and `nFeature_SCT`; the same from the counts, end to end | identical under the reference's kernel |
| K11 export and loaders | each loader's object against the saved object: the counts, the SCT layers and model, the metadata, the PCA and the graphs; `Read10X` and `read_10x_mtx` on `sct/` and on the export's top level; a file offered to the wrong slot | identical; both readers fail; the loader refuses |
| K12 on the fly | `scrust doubletfinder` given no `--clusters`, on the 150 DoubletFinder inputs: the labels it computed against the inputs' `seurat_clusters`; the outputs against the saved outputs (M7's criteria) and against M7's Rust outputs | identical; files byte-identical apart from timings and the clustering's own outputs |
| K13 markers (Seurat loader) | the workflow's marker steps on the loaded object, for every cluster: the combined table against the saved `<prefix>_markergenes.csv` (T and P) | byte-identical |
| Saved outputs (T, P, full panel) | `seurat_clusters` and its levels, `silhouette_width`, `neighborhood_purity` (from M12), `nCount_SCT` and `nFeature_SCT` (from M13), the logged numbers; `SCT_snn_res.0.8` for T and P | identical |
| Parameter grid | every grid cell's labels against Seurat run directly | target identical; floor ARI ≥ 0.99, with every difference traced to a near-tie (A7b) |
| Threads | the outputs at 1 and at many threads | byte-identical |
| Performance | peak RSS and wall time against R's direct run of the chain without UMAP (with UMAP, for context), with purity's and the export's shares reported | lower memory and time; at most 10 GiB on the 79,409-cell objects |

K3–K7 and K9 use no BLAS, so under teacher forcing they are exact on any machine (A5b). K10's μ is a BLAS product, so its exact check needs the reference's kernel even when teacher-forced.

### C6. Milestones

| ID | Content | Done when |
|---|---|---|
| M8 | Sources and reference data:<br>- `fetch_sources.sh` adds RcppAnnoy 0.0.23, RcppEigen 0.3.4.0.2, bluster 1.20.0, BiocNeighbors 2.4.0 and assorthead 1.4.0 (knncolle's headers);<br>- the panel table;<br>- exports for T and P, the SCT assay included;<br>- direct and instrumented runs in both profiles (A5b), with purity's and `correct_counts`' intermediates, and the run without UMAP;<br>- the instrumented optimizer and the K0 vectors;<br>- the marker references checked on T and P (C4);<br>- the R-vs-R floor | every direct run reproduces its saved labels and levels, silhouette widths, purities and SCT layers; the run without UMAP gives identical graphs and labels; the marker reference of every T and P object settled; the instrumented optimizer reproduces R's labels; the floor, and R's memory and time for the chain, measured; fixtures in the manifest |
| M9 | `neighbors`: `Kiss64Random`; the Annoy build (Euclidean, Angular, Manhattan) and its file format; the search; the nn and SNN graphs; `scrust neighbors` | K0 (Kiss64) and K3–K6 exact on T and P from R's own embeddings; time and memory measured on the largest object |
| M10 | `cluster`: `JavaRandom`; the network; algorithms 1–3 and modularity functions 1 and 2; the best start and the numbering; GroupSingletons; FindClusters' columns and levels; several resolutions; `scrust cluster --snn` | K0 (Java) and K7 exact on T and P from R's own SNN graphs, the objects with singletons included |
| M11 | The chain at the workflow's settings:<br>- SCTransform with seed 12345 (v5 and v3 assays) and 50 PCs;<br>- `approxSilhouette`;<br>- `scrust cluster --counts` and the staged starts, with the parameters, the preset and the TOML file;<br>- `scrust validate cluster`;<br>- the benchmark | K1–K8 met on T and P in both profiles; the direct runs' files and the saved labels identical; independent of the thread count; the benchmark in `docs/validation/clustering.md` |
| M12 | Neighbourhood purity: `knn`'s exact search with knncolle's arithmetic (the 50th distance, the radius, the range query) and the certificate; `neighborPurity`'s weights and sums; purity in `scrust cluster` | K9 on T and P, teacher-forced and end to end; the saved purities identical on T and P; time measured on the largest objects. **Decision point:** if the certificate flags any pair, port KMKNN (C3e) or agree with you how to treat the flagged cells |
| M13 | Corrected counts and the export:<br>- the `correct_counts` port, the log1p, `nCount_SCT` and `nFeature_SCT`;<br>- `scrust sctransform --corrected-counts` and `scrust cluster --export`;<br>- the export's layout, comment lines and README;<br>- the manifest for every command, `scrust doubletfinder`'s included | K10 on T and P in both profiles; every file the manifest lists present and described; `Read10X` and `read_10x_mtx` fail on the export as designed; the export's memory and time measured |
| M14 | Loaders: the R package (Seurat, SingleCellExperiment) and the Python package (AnnData), with their commands, in `loaders/`; their tests in hermetic envs (A2) | K11 and K13 on T and P; the R loader ready for M17's comparisons |
| M15 | Clustering on the fly: the chain as a library entry point in `cluster`, for later ports too (C3h); `scrust doubletfinder` given no `--clusters` (`--cluster-preset`, `--cluster-params`, `--cluster-stats`) | K12 on all 150 DoubletFinder inputs; the memory peak measured against M7's runs |
| M16 | The parameter grid (C4) | every grid cell identical, or A7b's floor met with every difference traced |
| M17 | Full panel: the exports (sbatch arrays); the Rust runs with `--export`, on each job's CPU type; the comparison with the saved labels, silhouette widths, purities and logged numbers; the R loader's comparison of the SCT layers, the PCA and the graphs | all 840 objects identical; a summary table and TSV |

M9 and M10 come before M11 for the reason M3 came first in Step 1: the new code is checked on R's own inputs, on any machine, before it is chained after the stages that depend on BLAS.

**Order after M11:** M12, M13, M15 and M16 each need only M11. M14 needs M13, and M17 needs M11–M14. Apart from M13 coming before M14, M12–M16 can be taken in any order.

**What belongs in Step 2, and what can follow:**
- **In Step 2:** M8–M17. Step 2 is complete when M11–M15 are met, with M16 and M17 as confirmation (Verification).
- **Can follow, outside Step 2:**
  - a port of KMKNN for purity, if M12's certificate never fires. It would only make the search faster on very large objects;
  - in decontX's port, its three sources of clusters (A7b): `--clusters`, clustering on the fly through C3h's entry point, and decontX's own clustering with the workflow's fallback when there are fewer than two clusters. Which is the default is decided in decontX's planning;
  - in SoupX's port, the same entry point;
  - publishing the loaders (Bioconductor or CRAN, PyPI, conda-forge);
  - writing h5ad directly from Rust, if a need appears (C3f);
  - Leiden (`algorithm = 4`).

### C7. Risks and mitigations

- **Annoy is float code with random draws that depend on the data.** One wrong operation order or draw count changes every later tree.
  - Mitigation: the byte comparison of the whole index (K4) under teacher forcing, and unit tests on constructed edge cases (duplicated points, margins of exactly 0, failed splits, n ≤ K).
  - Fallback: FFI to Annoy's C++.
- **Near-ties in the optimizer.** Gains and modularities are compared exactly, so a last-bit difference in a sum flips a decision, and every decision after it.
  - Mitigation: a statement-by-statement port of Seurat's C++ (not of the Java original), the instrumented optimizer's per-start values, and K7 under teacher forcing.
- **PCA's bits across kernels:** labels may move between the portable and production profiles.
  - Mitigation: profiles and CPU types matched (A5b, M7); the floor measures how often (M8).
- **irlba with 50 PCs at scale:** a larger `work` and more restarts than in DoubletFinder's runs, and M5 checked 50 components on one tier-T object only.
  - Mitigation: K2 on every P object.
- **SCTransform at seed 12345:** new RNG states, and the v3 path for the emptyDrops objects.
  - Mitigation: K1, with the RNG checkpoints.
- **GroupSingletons** depends on order and draws at random among ties; 36 production objects have singletons.
  - Mitigation: three P objects with 18 to 40 singletons, and isolated cells constructed in tests.
- **The panel's scale:** 540 new exports (about 40 GB, judging by M7's), 840 Rust runs with `--export`, and 840 comparisons in R.
  - Mitigation: sbatch arrays, as in M7; label-only exports where the counts already exist; each export's large files deleted after its comparison.
- **Purity at a rounding boundary:** KMKNN's floating-point bounds could skip a cell that an exhaustive search keeps.
  - Mitigation: C3e's certificate on every object (T, P and the full panel), and M12's decision point.
  - Fallback: the KMKNN port.
- **Corrected counts near .5:** μ's last bits depend on the OpenBLAS kernel, and rounding can then move a count by one.
  - Mitigation: exact checks under the reference's kernel (A5b) and K10 teacher-forced; the floor (M8) measures how often kernels disagree.
- **The loaders follow Seurat's internals** (`SCTAssay`, `SCTModel`), which change between versions.
  - Mitigation: built and tested against Seurat 5.5.1 and SeuratObject 5.4.0 (the pinned env); version checks; comparisons with the saved objects (K11).
- **The export's size:** the residuals are dense, 1.9 GB as `.npy` for 3,000 × 79,409 doubles.
  - Mitigation: written only with `--export` or `--keep-stages`.
- **Not reproduced:** UMAP (A7b).

---

## Verification

- **This document:** check it against the source facts: DoubletFinder 2.0.6 function bodies, `doubletfinder.R`, and the production logs and sacct numbers.
- **Each milestone:** its "Done when" column. `cargo test` covers tier T on every change; sbatch validation runs produce tier-P and full-panel reports in `docs/validation/`.
- **Step 1 complete:** M6 met (with M7 as the full-panel confirmation), and the benchmark table shows the memory and time gains over the direct R run.
- **Part C:** check it against the sources (Seurat 5.5.1, SeuratObject 5.4.0, sctransform 0.4.3, RcppAnnoy 0.0.23, bluster 1.20.0, BiocNeighbors 2.4.0 and assorthead 1.4.0's knncolle), the workflow's scripts, and the production logs.
- **Step 2 complete:** M11–M15 met, with M16 (the parameter grid) and M17 (the full panel) as confirmation, and the benchmark showing the memory and time gains over the direct R run of the chain.
