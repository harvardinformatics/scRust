# scRust

Rust reimplementations of single-cell RNA-seq tools, aimed at using less memory and running faster.

Each reimplementation is validated against the original tool's outputs, exactly or near-exactly, before it is offered as a replacement in the scRNA-seq preprocessing workflow.

**Planned tools, in order:**
1. DoubletFinder
2. a Seurat-compatible clustering chain
3. scDblFinder
4. emptyDrops
5. decontX
6. SoupX
7. CellBender

**Status:** Step 1, DoubletFinder, is complete through milestone M7. `scrust doubletfinder` runs DoubletFinder 2.0.6 as the workflow calls it, from a 10x counts directory and the cluster labels, every stage in Rust. It writes all of DoubletFinder's outputs. In addition, it writes the singlets' counts as a new 10x directory. It reproduces R bit for bit, under R's OpenBLAS kernel ([docs/validation/doubletfinder.md](docs/validation/doubletfinder.md)):

- **Instrumented R runs (24):** every stage, every random-number checkpoint and every scalar identical.
- **Direct R runs (23):** every output file identical.
- **Saved workflow outputs:** identical for all ten objects (the five test objects and five production objects). That covers the doublet set, the kept cells' pANN, pK and nExp, and, from M6b, the singlets' counts: every stored entry of the saved object's RNA counts.
- **The full panel (M7):** identical for all 150 production objects of the workflow, the same checks against each saved output, on the CPU type of each workflow job ([docs/validation/doubletfinder_m7.md](docs/validation/doubletfinder_m7.md)). That is 1,350 of 1,350 checks and 3.3 billion values, none differing, with a peak of at most 6.6 GiB on 8 threads, where the workflow's jobs peaked at up to 137 GiB.

On the largest object, soupx_emptydrops_sscovelli (79,409 cells), it takes 7 minutes 12 s on one thread and 54 s on 32, in 4.4–4.5 GiB. Writing `singlets/` is 34 s and 10 s of that. R's direct run takes 36 minutes and 127 GiB.

The steps:

- **M0:** the Cargo workspace, environments and validation skeleton.
- **M1b:** the upstream sources ([`third_party/`](third_party/README.md)).
- **M1:** `rcompat`, R's random numbers and arithmetic semantics, bit for bit.
- **M2:** the R reference data: exports, direct and instrumented reference runs, the R-vs-R floor ([doubletfinder_m2.md](docs/validation/doubletfinder_m2.md)).
- **M3:** the DoubletFinder core after PCA, without the N x N distance matrix ([doubletfinder_m3.md](docs/validation/doubletfinder_m3.md)).
- **M4:** SCTransform v2 ([doubletfinder_m4.md](docs/validation/doubletfinder_m4.md)).
- **M5:** RunPCA with irlba ([doubletfinder_m5.md](docs/validation/doubletfinder_m5.md)).
- **M6:** the whole tool ([doubletfinder.md](docs/validation/doubletfinder.md)).
- **M6b:** the singlets' counts as a new 10x directory, `singlets/` ([doubletfinder.md](docs/validation/doubletfinder.md#m6b-the-singlets-counts)).
- **M7:** the full panel, all 150 production objects against their saved workflow outputs ([doubletfinder_m7.md](docs/validation/doubletfinder_m7.md)).

Next, by the plan: the clustering crate ([docs/PLAN.md](docs/PLAN.md) A7b). Substituting `scrust` into the preprocessing workflow is a separate, later change.

- Plan and validation framework: [docs/PLAN.md](docs/PLAN.md)
- License: GPL-3.0 ([LICENSE](LICENSE)). Acknowledgements and upstream licenses: [NOTICE](NOTICE)

## Running DoubletFinder

```bash
source env.sh                     # portable profile: OPENBLAS_CORETYPE=Haswell
cargo build --release -p scrust
target/release/scrust doubletfinder --counts <10x dir> --clusters <seurat_clusters.tsv> --out <dir> --threads 8
```

- **Inputs.** `--counts` is a 10x Matrix Market directory (Cell Ranger v3 layout, gzipped or not), such as the one the workflow's emptyDrops step writes. `--clusters` is a TSV with a header and each cell's barcode and cluster label, in the matrix's column order. The factor's level order is read from `--cluster-levels`, or `seurat_clusters_levels.txt` next to the clusters file.
- **Outputs.**
  - `cells.tsv`: the barcode, pANN and classification of every cell.
  - `sweep_stats.tsv` and `bcmvn.tsv`.
  - `summary.json`: pK, nExp, the column names R writes, timings and peak memory.
  - In addition, `singlets/`: the counts without the doublets, as a new 10x directory (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) that the next tool can read. Each singlet's column is written exactly as stored, and the feature lines are the input's.

  `scrust doubletfinder --help` lists the options: seed, PCs, pN, nExp rate, a fixed pK, threads and stage dumps.
- **Every output of the original, plus count data out.** Each port keeps all of the original tool's outputs and adds none at their expense. Every tool here that processes counts also writes its result as a new 10x directory, so the tools chain on plain files ([docs/PLAN.md](docs/PLAN.md) A3).
- **Thread count.** The results do not depend on it.
- **The OpenBLAS kernel matters.** To reproduce the workflow's saved outputs, which were made on Sapphire Rapids nodes with the kernel OpenBLAS detects, run there with `OPENBLAS_CORETYPE=auto source env.sh` (which leaves the variable unset).

## Development

Everything stays inside the repository: the conda env in `.envs/rust-dev/` and cargo's caches in `.cargo-home/`, both gitignored. Nothing is installed in the home directory, and there is no rustup.

On Cannon, conda comes from `module load python`. Keep conda's package cache off the home directory: check `conda config --show pkgs_dirs`, and if needed set `CONDA_PKGS_DIRS` to lab storage before creating the env.

```bash
module load python
conda env create -p .envs/rust-dev -f envs/rust-dev.yml
source env.sh
cargo build
cargo test
```

To rebuild the env exactly as solved, use `conda create -p .envs/rust-dev --file envs/rust-dev.linux-64.pin.txt` instead of the `conda env create` line.

- [`envs/rust-dev.yml`](envs/rust-dev.yml): Rust 1.99, the conda GCC as linker, and the exact OpenBLAS and LAPACK builds of the R reference env.
- [`env.sh`](env.sh): activates the env, keeps cargo's caches in the repository, pins OpenBLAS to one thread and to the `Haswell` kernel (override with `OPENBLAS_CORETYPE=<core>`, or `auto` to let OpenBLAS detect the CPU), and keeps user-level R, Perl and Python libraries out.
- [`envs/r-reference/`](envs/r-reference/README.md): the R environment the original tools run in.
- [`validation/`](validation/): scripts that make reference data, sbatch wrappers, the dataset manifest and the match criteria. Anything that reads `.rds` files runs through sbatch, never on a login node.

| Crate | Contents |
|---|---|
| `rcompat` | R-compatible random numbers and R arithmetic semantics (sampling, `rnorm`, long-double `sum`/`mean`/`var`, `round`, `seq`, `order`, number formatting, R's decimal parser, `fft`, `dnorm`, `approx`), and R's numerical routines that the ports reach (`nlminb`, `optimize`, `uniroot`, `density`, `bw.SJ`, `ksmooth`, quantiles, gamma functions, `dnbinom`, matrixStats), plus a std-only OpenBLAS binding (BLAS, R's `%*%`, LAPACK `dgesdd` as R's `svd()` calls it) |
| `scio` | plain-file I/O: `.npy`, `%.17g` text, TSV, `summary.json`, gzip, 10x Matrix Market directories and plain Matrix Market files, sparse matrices |
| `sctransform` | SCTransform v2 (sctransform 0.4.3 with glmGamPoi 1.22.0's fitting), as Seurat 5.5.1 runs it |
| `pca` | Seurat's `RunPCA` with irlba 2.4.1 (PrepDR's filter with Seurat's `RowVar`, the Lanczos bidiagonalisation with its restarts and rare branches), its two large products split across threads only where OpenBLAS's kernels keep R's bits |
| `knn` | exact nearest neighbours, memory-bounded |
| `kde` | `KernSmooth::bkde`, `approxfun` and moments |
| `neighbors` | Seurat's `FindNeighbors`: Annoy 1.17.3 as RcppAnnoy 0.0.23 builds and searches it (`Kiss64Random`, the Euclidean, Angular and Manhattan forests, the index file, libstdc++'s heap order), the `SCT_nn` graph and `ComputeSNN`; byte-identical to R's index and graphs |
| `doubletfinder` | DoubletFinder 2.0.6: the whole tool from the counts (the doublets drawn from R's random-number stream, SCTransform, PCA), the sweep, `summarizeSweep`, `find.pK`, nExp and the calls; stage dumps in the instrumented R runs' layout and the comparison with R runs and saved workflow outputs |
| `scrust` | the command-line binary, `scrust <tool> ...`: `doubletfinder` (the whole tool; `--from-stages` runs it on an instrumented R run's embeddings), `sctransform`, `pca`, `neighbors`, `knn`, `validate doubletfinder`, `validate sctransform`, `validate pca` and `validate neighbors` (Rust against R), and `fingerprint` (CPU, OpenBLAS kernel and threads, glibc) |
