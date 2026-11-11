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

CRAN and Bioconductor move a version to their archive once a newer one is released, so the script tries the current index first and then the archive; the sha256 guarantees the same bytes either way.

**These are the sources of the production R env** (`envs/r-reference/`). Checked on 2026-11-10: each sha256 in `SHA256SUMS` equals the source checksum recorded in the conda-forge recipe (`info/recipe/meta.yaml`) of the package build installed in the production env, and glmGamPoi's tarball has the md5 recorded by its bioconda recipe. DoubletFinder's tarball records its commit in the tar header, which the script checks; the production env's DoubletFinder has the same `RemoteSha`.

GitHub generates archive tarballs on request. If DoubletFinder's sha256 ever stops matching, its compression may have changed rather than its content: compare the unpacked files with the old ones before recording the new checksum.
