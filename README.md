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

**Status:** milestone M3. The Cargo workspace, environments and validation skeleton are in place (M0); the upstream sources can be fetched (M1b, [`third_party/`](third_party/README.md)); `rcompat`, R's random numbers and arithmetic semantics, matches R 4.5.3 bit for bit on its reference vectors (M1); DoubletFinder's R reference data are built (M2: exports of the five workflow test objects and five production objects, direct and instrumented reference runs that reproduce every saved workflow output, stage dumps S1–S7, the R-vs-R floor; [docs/validation/doubletfinder_m2.md](docs/validation/doubletfinder_m2.md)); and the DoubletFinder core after PCA is ported (M3): on R's own embeddings, `scrust doubletfinder --from-stages` reproduces every instrumented reference run of all ten objects exactly, stage by stage (neighbours, pANN, the sweep summary and pK, nExp, the calls) and file for file, using 32 MiB where R's distance matrix needs 89.7 GB ([docs/validation/doubletfinder_m3.md](docs/validation/doubletfinder_m3.md)). SCTransform and PCA (M4, M5) come next; no tool runs end to end yet.

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
| `rcompat` | R-compatible random numbers and R arithmetic semantics (sampling, `rnorm`, long-double `sum`/`mean`/`var`, `round`, `seq`, `order`, number formatting, R's decimal parser, `fft`, `dnorm`, `approx`) |
| `scio` | plain-file I/O: `.npy`, `%.17g` text, TSV, `summary.json`, gzip and 10x Matrix Market directories |
| `sctransform` | SCTransform v2 (sctransform with glmGamPoi's fitting) |
| `pca` | irlba as Seurat's `RunPCA` calls it |
| `knn` | exact nearest neighbours, memory-bounded |
| `kde` | `KernSmooth::bkde`, `approxfun` and moments |
| `doubletfinder` | DoubletFinder |
| `scrust` | the command-line binary, `scrust <tool> ...`: `doubletfinder` (so far on an instrumented R run's embeddings, `--from-stages`), `knn`, `validate doubletfinder` (Rust against R), and `fingerprint` (CPU, OpenBLAS kernel and threads, glibc) |
