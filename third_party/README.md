# third_party

The exact upstream sources that scRust's ports are read from (docs/PLAN.md A1 and B6, milestone M1b). They are reading references only: nothing in them is built, sourced or run. Their authors and licenses are credited in [NOTICE](../NOTICE).

```bash
third_party/fetch_sources.sh
```

This downloads each tarball into `third_party/src/tarballs/` (gitignored), checks it against [`SHA256SUMS`](SHA256SUMS), and unpacks it into its own directory, `third_party/src/<name>-<version>/`. Run it again at any time: it downloads only what is missing and re-verifies everything. `--record` records the sha256 of a newly added source; it never replaces a recorded one.

| Source | Version | From |
|---|---|---|
| R | 4.5.3 | CRAN `src/base/R-4/` |
| sctransform | 0.4.3 | CRAN |
| glmGamPoi | 1.22.0 | Bioconductor 3.22 |
| matrixStats | 1.5.0 | CRAN |
| Seurat | 5.5.1 | CRAN archive |
| SeuratObject | 5.4.0 | CRAN |
| irlba | 2.4.1 | CRAN |
| KernSmooth | 2.23-27 | CRAN |
| fields | 17.3 | CRAN archive |
| DoubletFinder | commit `1b244d8f0d54b4b1cb4365639931bbb16f01e1cd` (2.0.6) | GitHub archive of chris-mcginnis-ucsf/DoubletFinder |
| Matrix | 1.7-6 | CRAN (M4: `dgCMatrix` arithmetic, margins, subsetting) |
| RcppArmadillo | 15.2.3-1 | CRAN archive (M4: the headers glmGamPoi 1.22.0 and sctransform 0.4.3 were compiled with; the env has 15.6.2-1 installed) |
| Rcpp | 1.1.1 | CRAN archive (M4: the headers glmGamPoi was compiled with; sugar `lgamma`, `digamma`, `sum`) |
| future | 1.82.0 | CRAN (M4: `future.seed` and the main RNG stream) |
| future.apply | 1.21.0 | CRAN (M4: `future_lapply`) |
| parallelly | 1.50.0 | CRAN (M4: `supportsMulticore`, number of workers) |
| MatrixGenerics | 1.22.0 | Bioconductor 3.22 (M4: dispatch of `rowMeans2` etc.) |
| DelayedMatrixStats | 1.32.0 | Bioconductor 3.22 (M4: glmGamPoi's `rowMeans2`, `rowVars`, `colMeans2`) |
| OpenBLAS | 0.3.34 | GitHub archive of OpenMathLib/OpenBLAS, tag `v0.3.34` (M5: the kernels and drivers whose rounding decides how irlba's products may be split across threads; nothing is ported from it) |

CRAN and Bioconductor move a version to their archive once a newer one is released, so the script tries the current index first and then the archive; the sha256 guarantees the same bytes either way.

**These are the sources of the production R env** (`envs/r-reference/`), except where a package's compiled code includes another package's headers: then it is the version the binary was built against (M4: glmGamPoi 1.22.0 and sctransform 0.4.3 were built with RcppArmadillo 15.2.3-1 and Rcpp 1.1.1, per the host requirements in their conda recipes, `info/recipe/meta.yaml`). Checked on 2026-11-10 (M1b) and 2026-11-11 (M4 additions): each CRAN sha256 in `SHA256SUMS` equals the source checksum recorded in the conda-forge recipe (`info/recipe/meta.yaml`) of the package build installed in the production env (for RcppArmadillo and Rcpp, of the builds glmGamPoi was compiled against), and the Bioconductor tarballs (glmGamPoi, MatrixGenerics, DelayedMatrixStats) have the md5 recorded by their bioconda recipes. DoubletFinder's tarball records its commit in the tar header, which the script checks; the production env's DoubletFinder has the same `RemoteSha`.

OpenBLAS is not an R package: both envs load the same build (`libopenblasp-r0.3.34.so`, conda-forge `libopenblas 0.3.34 pthreads_hf13c14d_2`, the library R's `La_library()` reports), and the tarball's sha256 is the one in that build's recipe (`info/recipe/recipe.yaml`, source `https://github.com/OpenMathLib/OpenBLAS/archive/v0.3.34.tar.gz`). It was built with `DYNAMIC_ARCH=1 TARGET=PRESCOTT` (`info/recipe/build.sh`): the kernels are compiled once per CPU core, the drivers once.

GitHub generates archive tarballs on request. If DoubletFinder's sha256 ever stops matching, its compression may have changed rather than its content: compare the unpacked files with the old ones before recording the new checksum.
