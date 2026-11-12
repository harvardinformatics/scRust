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

**Status:** milestone M4. The Cargo workspace, environments and validation skeleton are in place (M0); the upstream sources can be fetched (M1b, [`third_party/`](third_party/README.md)); `rcompat`, R's random numbers and arithmetic semantics, matches R 4.5.3 bit for bit on its reference vectors (M1); DoubletFinder's R reference data are built (M2: exports of the five workflow test objects and five production objects, direct and instrumented reference runs that reproduce every saved workflow output, stage dumps S1–S7, the R-vs-R floor; [docs/validation/doubletfinder_m2.md](docs/validation/doubletfinder_m2.md)); the DoubletFinder core after PCA is ported (M3): on R's own embeddings, `scrust doubletfinder --from-stages` reproduces every instrumented reference run of all ten objects exactly, stage by stage and file for file, using 32 MiB where R's distance matrix needs 89.7 GB ([docs/validation/doubletfinder_m3.md](docs/validation/doubletfinder_m3.md)); and SCTransform v2 is ported (M4): Seurat 5.5.1's `SCTransform()` with sctransform 0.4.3 and glmGamPoi 1.22.0, every sub-stage (setup, fits with R's `nlminb`, regularisation, residuals) bit-identical to R in all seven SCTransform calls of all 24 instrumented runs, teacher-forced and chained, and `scrust sctransform` byte-identical to R's outputs on the largest objects, 2.3-2.5 times faster than R on one thread and 29-35 times on 32, in a ninth to an eighteenth of its memory at either thread count ([docs/validation/doubletfinder_m4.md](docs/validation/doubletfinder_m4.md)). PCA (M5) comes next; no tool runs end to end yet.

- Plan and validation framework: [docs/PLAN.md](docs/PLAN.md)
- License: GPL-3.0 ([LICENSE](LICENSE)). Acknowledgements and upstream licenses: [NOTICE](NOTICE)

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
| `rcompat` | R-compatible random numbers and R arithmetic semantics (sampling, `rnorm`, long-double `sum`/`mean`/`var`, `round`, `seq`, `order`, number formatting, R's decimal parser, `fft`, `dnorm`, `approx`), and R's numerical routines that the ports reach (`nlminb`, `optimize`, `uniroot`, `density`, `bw.SJ`, `ksmooth`, quantiles, gamma functions, `dnbinom`, matrixStats), plus a std-only OpenBLAS binding |
| `scio` | plain-file I/O: `.npy`, `%.17g` text, TSV, `summary.json`, gzip and 10x Matrix Market directories, sparse matrices |
| `sctransform` | SCTransform v2 (sctransform 0.4.3 with glmGamPoi 1.22.0's fitting), as Seurat 5.5.1 runs it |
| `pca` | irlba as Seurat's `RunPCA` calls it |
| `knn` | exact nearest neighbours, memory-bounded |
| `kde` | `KernSmooth::bkde`, `approxfun` and moments |
| `doubletfinder` | DoubletFinder |
| `scrust` | the command-line binary, `scrust <tool> ...`: `doubletfinder` (so far on an instrumented R run's embeddings, `--from-stages`), `sctransform`, `knn`, `validate doubletfinder` and `validate sctransform` (Rust against R), and `fingerprint` (CPU, OpenBLAS kernel and threads, glibc) |
