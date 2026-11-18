# scrustio (R): scrust outputs as Seurat and SingleCellExperiment objects

`scrustio` reads the plain files that scRust's command-line tool, `scrust`, writes, and builds the objects R users work with. It finds every file through the manifest in the output's `summary.json` (schema `scrust-files`, version 1: [`docs/output_manifest.md`](../../../docs/output_manifest.md)), by the file's role, never by its name. It puts each matrix only into a slot meant for its scale: the raw counts into the raw-count slots, and SCTransform's corrected counts, their log1p and the residuals only into slots for normalised data (docs/PLAN.md C3g, milestone M14).

## What it builds

| scrust output | `load_seurat()` | `load_sce()` |
|---|---|---|
| `scrust cluster --export` | the raw counts in the `RNA` assay's `counts` layer (an `Assay5`); an `SCTAssay` named `SCT`: the corrected counts in `counts`, their log1p in `data`, the residuals in `scale.data` (rows in the SCT assay's gene order, as Seurat keeps them), the variable features and exactly one `SCTModel`, rebuilt as SCTransform makes it; `meta.data` (factors in their levels' order); the `pca` reduction; the `SCT_nn` and `SCT_snn` graphs; identities from `seurat_clusters`; `SCT` the default assay | the `counts` assay; `altExp(sce, "SCT")` with the assays `corrected_counts` and `logcounts` (the model genes), the SCTModel's feature attributes as its `rowData`, its cell attributes as its `colData` and its other slots in `metadata()`; `altExp(sce, "SCT_scale_data")` with the assay `scaledata` (the variable features, in `VariableFeatures()` order); `colData`; `reducedDim(sce, "PCA")` with the loadings as its `"rotation"` attribute (and `"stdev"`, `"total.variance"`); the `colPairs` `SCT_nn` and `SCT_snn` |
| `scrust sctransform --corrected-counts` | the `SCT` assay alone (this output holds no raw counts), with the counts' per-cell columns | no rows in the main experiment; the SCT data in the altExps |
| `scrust doubletfinder` | the singlets' raw counts (`singlets/`) in the `RNA` assay; DoubletFinder's calls as metadata columns, named as DoubletFinder names them (`pANN_<pN>_<pK>_<nExp>`, `DF.classifications_<pN>_<pK>_<nExp>`) | the same as `counts` and `colData` |

`add_scrust(object, dir)` adds an output's per-cell columns to an object you already have, matching cells by barcode: `scrust doubletfinder`'s calls (for all the input's cells), or `scrust cluster`'s labels (as factors in their levels' order), silhouette widths and purities.

Seurat names features as `Read10X()` and `CreateSeuratObject()` do (`make.unique()`, then `_` and `|` replaced by `-`), and every list of features of the output is mapped onto the counts' rows; the SingleCellExperiment uses the features file's names, made unique, with the ids, names and types in `rowData`.

**What is not rebuilt** (it is not in a scrust output): the `umap` reduction (scrust does not run RunUMAP, which changes nothing downstream), the `commands` (scrust records its run in `summary.json`), metadata columns of other workflow steps (such as `percent.mt`), and, for objects whose RNA assay was a v3 `Assay` (the workflow's emptyDrops objects, made by `as.Seurat()`), that class: the loader makes the `Assay5` that `CreateSeuratObject()` makes, with the same counts.

**`cells_step1` is drawn again, not read.** The SCTModel's cell attributes say which cells vst() fitted on in its first step. SCTransform calls `set.seed(seed.use)` and passes vst `min(ncells, cells)`; vst takes every cell when there are at most `ncells` (5,000 by default), and otherwise draws them with `sample(colnames(umi), ncells)`, the first draw after that `set.seed()`. The output does not list the drawn cells; its manifest records the seed and `ncells` (the parameters and the seeds). That is enough to make the same draw, so no file was added for them. The loader therefore:

- takes `ncells` and the seed from the manifest (the summary's `n_cells_step1` and the manifest's two records of the seed must agree with them) and the cells in `counts/` order, which is vst's `colnames(umi)`;
- draws only when `min(ncells, cells)` is below the number of cells, with R's own `sample()` under R's default generator kinds (`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`), set for the draw whatever yours are;
- puts your generator kinds and `.Random.seed` back afterwards, and removes `.Random.seed` again if you had none;
- checks the draw against the output, when it holds the raw counts. vst's next draw, from the same stream, picks the step-1 genes from those detected in at least `min_cells` of the drawn cells, and the output lists those genes in vst's sampling order (`sct_model_pars`). The loader makes that draw too, with sctransform's own `row_gmean()` and `row_var()`, and refuses the output if it does not give that list. A different set of cells would not give it. This adds a few seconds on the largest objects and needs no extra file.
- runs that check only under the versions it mirrors and was validated with: R 4.5.3 and sctransform 0.4.3, on Linux. The gene draw is weighted, so a last-bit change in sctransform's internal functions, in R's `density()` (R 4.4.0 changed it) or in the C library's `log()` and `exp()` could change which genes it picks. Elsewhere, and without sctransform, the loader draws the cells without the check (the draw needs only R's `sample()`) and gives one warning saying why. Then nothing is refused.

`load_sce()` keeps the same `cells_step1` in the `colData` of `altExp(sce, "SCT")`, drawn and checked the same way.

## Install

Into a library of your choice (scRust's validation installs it into the repository's `.rlib/`, never into the R env it tests with):

```bash
R CMD INSTALL --library=<library> loaders/R/scrustio
```

It needs `Matrix` and `jsonlite`; `Seurat` (5.5.1) and `SeuratObject` (5.4.0) for `load_seurat()`, and `SingleCellExperiment`, `SummarizedExperiment` and `S4Vectors` for `load_sce()`. It builds Seurat 5.5.1's and SeuratObject 5.4.0's classes and is tested with them; with other versions it warns and carries on.

## Use

```r
library(scrustio)
obj <- load_seurat("export")               # a Seurat object
sce <- load_sce("export")                  # a SingleCellExperiment
m <- read_manifest("export")               # the manifest: files by role (manifest_entries(m, "sct_counts"))
obj <- load_seurat("export", check_sha256 = TRUE)   # check every file's SHA-256 first
```

As one line, writing `.rds`:

```bash
Rscript -e 'scrustio::save_rds("export", "object.rds")'
```

```bash
Rscript -e 'scrustio::save_rds("export", "sce.rds", object = "sce")'
```

(with `R_LIBS=<library>` if it is not on R's library path).

## The rule, enforced

`slot_table()` lists the slots the loaders fill from matrix files, each with the one role, scale and format it takes. `read_for_slot()` refuses a file of any other role, scale or format, a file whose Matrix Market comment lines say "NOT raw counts" for a raw-count slot, and an unmarked file for a slot of SCTransform's layers. `read_manifest()` refuses another schema or version, a listed file that is missing or has another size (or, with `check_sha256 = TRUE`, another SHA-256), and a path outside the output; the identifiers of each matrix's rows and columns are checked against its shape.

## Tests

`tests/scrustio_tests.R` (in the R env, through sbatch on Cannon: `validation/slurm/m14_tests.sbatch`) makes outputs with the `scrust` binary from the committed fixtures and checks: every wrong combination of role, scale and format for every slot is refused; tampered manifests are refused; the manifest's rules; the objects' layout; FindMarkers on the Seurat object; `cells_step1` on an export with `--ncells 100` (the cells scrust drew, by its own record; the caller's generator state restored under several kinds; the check refusing listed step-1 genes out of order and another seed, and skipped with one warning, the same cells and nothing refused, under other versions); the other commands' outputs; `add_scrust()`; and that `save_rds()`'s `.rds` reloads into the loader's object. Milestone M14's validation compares the objects with the saved workflow objects (docs/validation/clustering.md, "M14").
