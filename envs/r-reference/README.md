# R reference environment

The R environment in which the original tools run when scRust makes its reference outputs (docs/PLAN.md A2 and A4). It is the preprocessing workflow's single R environment. Both files here are read-only copies from the pipeline repo (`scRNAseq-preprocessing`, `workflow/envs/`), where they were last changed in commit `57a2ab8` (2026-11-09):

- `scrnaseq_r.yml`: the requested package set;
- `scrnaseq_r.linux-64.pin.txt`: the exact solve (explicit package URLs), from which Snakemake builds the environment.

Update these copies only together with the workflow's environment, and keep the BLAS and LAPACK builds pinned in `envs/rust-dev.yml` identical to the ones in the pin file.

## Production environment on Cannon

`/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_`

Snakemake built it from this pin file: the `.pin.txt` it saved next to the environment is identical to the copy here (checked 2026-11-10).

Main versions: R 4.5.3, Seurat 5.5.1, SeuratObject 5.4.0, sctransform 0.4.3, glmGamPoi 1.22.0, irlba 2.4.1, fields 17.3, KernSmooth 2.23-27, matrixStats 1.5.0, Matrix 1.7-6, RcppAnnoy 0.0.23, DoubletFinder 2.0.6, DropletUtils 1.30.0, SoupX 1.6.2, decontX 1.8.0, celda 1.26.0, scDblFinder 1.24.10, xgboost 3.2.0. BLAS: OpenBLAS 0.3.34, build `pthreads_hf13c14d_2`, with `libblas`, `libcblas` and `liblapack` 3.11.0 build 12 (openblas).

## DoubletFinder

DoubletFinder is not on conda. One fixed commit (v2.0.6) is installed into the environment from GitHub:

```r
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd", dependencies = FALSE, upgrade = "never")
```

The workflow does this in its rule `install_doubletfinder` (`workflow/rules/install_doubletfinder.smk`), which also passes `force = TRUE`. The production environment's DoubletFinder records `RemoteSha: 1b244d8f0d54b4b1cb4365639931bbb16f01e1cd`.

## Rebuilding it elsewhere

On another cluster, with the conda package cache kept off the home directory:

```bash
conda create -p <path> --file envs/r-reference/scrnaseq_r.linux-64.pin.txt
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
  <path>/bin/Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd", dependencies = FALSE, upgrade = "never")'
```
