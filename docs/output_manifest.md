# The output manifest (`summary.json`'s `files` section)

Every scrust command that writes `summary.json` (`scrust sctransform`, `pca`, `neighbors`, `cluster`, `doubletfinder`) lists there every file of its output directory: what the file is (its **role**), the **scale** of its values, its format and shape, where its row and column identifiers are, its size and SHA-256, and how the output was made (its **provenance**). This document is the schema's specification (docs/PLAN.md C3f; milestone M13). The loaders (C3g, milestone M14) are written against it: they find each file by its role, and put a matrix only into a slot meant for its scale.

The schema is **`scrust-files`, version 1**. A reader checks both and refuses other versions; a later version that changes the meaning of a field changes the version number. Fields may be added within a version; a reader ignores the fields it does not know.

The implementation is `crates/scio/src/manifest.rs`; `scrust validate manifest --dir <output dir>` checks a directory against its manifest.

## Layout

`summary.json` is one JSON object. Every command writes `result` (what was computed; deterministic), `run` (how the run went: timings, memory, the environment fingerprint) and, last, `files`, the manifest; `scrust sctransform`, `pca`, `neighbors` and `cluster` also write `tool`, `reference` and `scrust_version` at the top, and all but `cluster` (whose parameters are in `run`) `parameters`. A reader takes the manifest from `files` and ignores the other members:

```json
"files": {
  "schema": "scrust-files",
  "version": 1,
  "documentation": "docs/output_manifest.md",
  "command": "cluster",
  "provenance": { ... },
  "entries": [ { ... }, ... ]
}
```

| field | type | meaning |
|---|---|---|
| `schema` | string | always `"scrust-files"` |
| `version` | integer | `1` |
| `documentation` | string | this document's path in the scRust repository |
| `command` | string | the scrust command that wrote the directory: `sctransform`, `pca`, `neighbors`, `cluster`, `doubletfinder` |
| `provenance` | object | how the output was made (below) |
| `entries` | array of objects | one per file of the directory (below) |

## Provenance

| field | type | meaning |
|---|---|---|
| `scrust_version` | string | scrust's version (`scrust --version`) |
| `command_line` | array of strings | the command line as given: `"scrust"`, the subcommand, its arguments |
| `command` | string | the command as the normalised matrices' comment lines give it: the command line without the options that change no value of the normalised matrices (`--out`, `--threads`, and for `scrust cluster` also `--keep-stages`, `--dump-stages`, `--purity-lists`, `--annoy`, `--metadata`), so that those files do not depend on where or how fast the run went |
| `inputs` | array of objects | each input the command read: `option` (the option that named it, e.g. `counts`), `path` (as given), `files` (each file read: `name`, relative to `path`; `bytes`; `sha256`). `files` is empty only for `--from-stages` (a validation input) |
| `reproduces` | array of strings | the original software the outputs reproduce: R packages, functions, versions |
| `parameters` | object | the parameters of the run (the same object as `run.parameters` or `run`'s parameters) |
| `seeds` | object | the random seeds the run reproduces, by where R sets them |

For a 10x input, `inputs[].files` lists the matrix, the barcodes and the features, in that order; their SHA-256 also appear in the normalised matrices' comment lines, as `sha256sum` prints them.

## Entries

| field | type | meaning |
|---|---|---|
| `path` | string | relative to the output directory, `/` between directories |
| `role` | string | what the file is ([Roles](#roles)) |
| `scale` | string | the scale of its values ([Scales](#scales)); `none` for files that hold no matrix of values |
| `format` | string | `matrix_market`, `tsv`, `npy`, `json`, `text` or `annoy` ([Formats](#formats)) |
| `gzip` | boolean | the file is gzipped (its name ends in `.gz`) |
| `header` | boolean or null | `tsv`: whether the first line names the columns; null for other formats |
| `shape` | array of integers or null | `matrix_market`: rows, columns; `tsv`: data rows (without the header), columns; `npy`: the array's shape; null for other formats |
| `nnz` | integer or null | `matrix_market`: stored entries |
| `rows`, `columns` | object or null | where the identifiers of the rows and of the columns are ([Identifiers](#identifiers)); null when there are none, or (a table's columns) when the file's own header names them |
| `description` | string | what the file holds, in words |
| `bytes` | integer | the file's size; null only for `summary.json` itself |
| `sha256` | string | the file's SHA-256, 64 lowercase hex digits; null only for `summary.json` itself |

`summary.json` is listed with role `summary`, without size or checksum: it holds them.

### Completeness

- Every file in the output directory, in any subdirectory, is listed once, except the files of a nested output (below), which its own manifest lists.
- Every listed file exists, with the listed size and SHA-256.
- Every file named in `rows` or `columns` is listed.

A command checks its own manifest before writing `summary.json` and fails if any rule does not hold.

### Nested outputs (M15)

A subdirectory can itself be a scrust output, with its own `summary.json` and manifest: `scrust doubletfinder`'s `clusters/` (clustering on the fly) and `recluster/` (`--recluster`, the singlets' re-clustering), each as `scrust cluster` writes its output.

- The parent's manifest lists the nested output's `summary.json` alone, with role `nested_output` (format `json`, scale `none`), its size and its SHA-256. It lists none of the nested output's other files: the nested manifest does, under the nested directory.
- The checks recurse: a nested output's files must satisfy the completeness rules of its own manifest, and the naming rules hold within each manifest. `scrust validate manifest` reports a nested problem under the nested directory's name.
- Loaders find the parent's own files by role, one per role, as before, and ignore `nested_output`. A nested output loads on its own, from its directory. So the singlets' re-clustering, which has its own `raw_counts` with `--recluster-export`, its own metadata and graphs, never shadows the DoubletFinder output's files.

### Identifiers

```json
"rows":    {"kind": "features", "file": "sct/sct_features.tsv.gz", "column": 1},
"columns": {"kind": "cells",    "file": "metadata.tsv",           "column": "barcode"}
```

| field | type | meaning |
|---|---|---|
| `kind` | string | what the identifiers name: `cells`, `features`, `variable_features`, `components`, `clusters` |
| `file` | string | the listed file that holds them; it may be the entry itself (a table whose first column names its rows) |
| `column` | string or integer | a header name when that file has a header (`header` true), else the column's 1-based number (`1` for a file of one identifier per line) |

The identifiers are in the file's row order, one per row of the entry's dimension. Features of a 10x directory are identified by the features file's column 2 (the names, Seurat's and scanpy's default); its column 1 holds the ids.

### Formats

| `format` | what | R | Python |
|---|---|---|---|
| `matrix_market` | coordinate Matrix Market, `integer` or `real` `general`, 1-based, entries in stored (column-major) order, values as `%.17g` writes them; comment lines after the banner | `Matrix::readMM(path)` | `scipy.io.mmread(path)` |
| `tsv` | tab-separated text; `header` says whether the first line names the columns | `read.delim(path, header = ...)` | `pandas.read_csv(path, sep="\t", header=...)` |
| `npy` | NumPy format 1.0, little-endian `<f8` or `<i4`, matrices in Fortran (column-major) order | `RcppCNPy::npyLoad(path)` | `numpy.load(path)` |
| `json` | JSON | `jsonlite::fromJSON(path)` | `json.load` |
| `text` | plain text (README, one value per line) | `readLines(path)` | `open(path).read()` |
| `annoy` | Annoy's index as RcppAnnoy 0.0.23's `save()` writes it | RcppAnnoy | annoy |

Doubles in text are written with 17 significant digits, so a correctly rounded reader recovers every bit.

### Scales

The scale tells a loader which slot a matrix may go into. Only `raw_counts` is raw counts.

| `scale` | values |
|---|---|
| `raw_counts` | the input's UMI counts (or a subset of its cells), as stored |
| `corrected_counts` | SCTransform's corrected counts: sctransform's `correct()`, every cell corrected to the median sequencing depth and rounded; integers, but **not raw counts** (Seurat's SCT `counts` layer) |
| `log1p_corrected_counts` | `log1p` of the corrected counts (Seurat's SCT `data` layer) |
| `residuals` | SCTransform's Pearson residuals of the variable features, clipped at ±sqrt(cells/30) and centred (Seurat's SCT `scale.data`) |
| `pca_embeddings`, `pca_loadings`, `pca_stdev` | RunPCA's embeddings (`u %*% diag(d)`), loadings and standard deviations |
| `graph_nn`, `graph_snn` | FindNeighbors' `SCT_nn` (0/1) and `SCT_snn` (Jaccard) graphs |
| `neighbour_indices`, `neighbour_distances` | each cell's k nearest neighbours (1-based cell numbers) and their distances |
| `none` | no matrix of values: tables of mixed columns, labels, levels, text, stage dumps |

### Roles

| `role` | file (command) | format, scale | rows; columns |
|---|---|---|---|
| `summary` | `summary.json` (all) | json | |
| `readme` | `README.txt` (`cluster --export`, `sctransform --corrected-counts`) | text | |
| `raw_counts` | `counts/matrix.mtx(.gz)` (`cluster --export`: the input, copied byte for byte); `singlets/matrix.mtx.gz` (`doubletfinder`) | matrix_market, `raw_counts` | features (the features file, column 2); cells (the barcodes file) |
| `raw_counts_cells` | `counts/barcodes.tsv(.gz)`, `singlets/barcodes.tsv.gz` | tsv without header | |
| `raw_counts_features` | `counts/features.tsv(.gz)`, `singlets/features.tsv.gz`: id, name, type | tsv without header | |
| `rna_cell_metadata` | `cell_metadata.tsv` (`sctransform`, `pca`, `neighbors`): barcode, nCount_RNA, nFeature_RNA | tsv | cells (`barcode`) |
| `cell_metadata` | `metadata.tsv` (`cluster`): barcode, nCount_RNA, nFeature_RNA, nCount_SCT and nFeature_SCT (with `--export`), each `res.` column, seurat_clusters, silhouette_width, silhouette_other, neighborhood_purity, purity_maximum, then `--metadata`'s columns | tsv | cells (`barcode`) |
| `cluster_labels` | `seurat_clusters.tsv`: barcode, seurat_clusters | tsv | cells (`barcode`) |
| `cluster_levels` | `seurat_clusters_levels.txt`, `<column>_levels.txt`: a factor's levels in order, one per line | tsv without header | |
| `cluster_columns` | `clusters.tsv`: barcode, one column per resolution | tsv | cells (`barcode`) |
| `workflow_cluster_ids`, `workflow_nclusters` | `cluster_ids.txt`, `nclusters.txt`, as the workflow's `write_cluster_metadata` writes them | tsv without header; text | |
| `optimizer_ids` | `optimizer_ids.tsv`: each resolution's optimizer result before GroupSingletons | tsv | cells (`barcode`) |
| `sct_cell_attributes` | `cell_attr.tsv`: barcode, umi, log_umi (Seurat's SCTModel `cell.attributes`) | tsv | cells (`barcode`) |
| `sct_feature_attributes` | `gene_attr.tsv`: the model genes' id, name and attributes (SCTModel `feature.attributes`) | tsv | features (`name`) |
| `sct_model_pars` | `model_pars.tsv`: the fitting genes' theta, (Intercept), log_umi | tsv | features (`name`) |
| `sct_model_pars_fit` | `model_pars_fit.tsv`: the regularised parameters of every model gene (SCTModel `model.pars`) | tsv | features (`name`) |
| `sct_variable_features` | `variable_features.tsv`: id, name, highest residual variance first (`VariableFeatures()`) | tsv | variable_features (`name`) |
| `sct_scale_data` | `scale_data.npy` | npy, `residuals` | variable_features (`variable_features.tsv`, `name`); cells (`cell_attr.tsv`, `barcode`) |
| `sct_counts` | `corrected_counts.mtx.gz` | matrix_market, `corrected_counts` | features (`sct_features.tsv.gz`, column 1); cells (`sct_cells.tsv.gz`, column 1) |
| `sct_data` | `log1p_corrected_counts.mtx.gz` | matrix_market, `log1p_corrected_counts` | as `sct_counts` |
| `sct_features` | `sct_features.tsv.gz`: the model genes' names, in the counts' order | tsv without header | |
| `sct_cells` | `sct_cells.tsv.gz`: the cells' barcodes | tsv without header | |
| `pca_embeddings` | `pca/embeddings.tsv`: barcode, PC_1, ... | tsv, `pca_embeddings` | cells (`barcode`) |
| `pca_loadings` | `pca/loadings.tsv`: id, name, PC_1, ... | tsv, `pca_loadings` | variable_features (`name`) |
| `pca_stdev` | `pca/stdev.tsv`: pc, stdev | tsv, `pca_stdev` | components (`pc`) |
| `rng_state` | `pca/random_seed.txt`: R's `.Random.seed` after RunPCA | text | |
| `knn_index`, `knn_distance` | `neighbors/knn_index.tsv`, `knn_distance.tsv`: barcode, nn_1, ..., nn_k | tsv, `neighbour_indices`, `neighbour_distances` | cells (`barcode`) |
| `graph_nn`, `graph_snn` | `neighbors/nn.mtx.gz`, `snn.mtx.gz` | matrix_market, `graph_nn`, `graph_snn` | cells, both: `metadata.tsv`'s `barcode` in `scrust cluster`'s output (whose `neighbors/` may hold R's graphs without `knn_index.tsv`), `knn_index.tsv`'s in `scrust neighbors`' |
| `annoy_index` | `neighbors/annoy.ann` (`--annoy yes`) | annoy | |
| `doublet_calls` | `cells.tsv` (`doubletfinder`): barcode, pANN, classification | tsv | cells (`barcode`) |
| `doublet_sweep_stats`, `doublet_bcmvn` | `sweep_stats.tsv`, `bcmvn.tsv` (`doubletfinder`) | tsv | |
| `stage_dump`, `stage_index`, `stage_scalars`, `start_partitions` | `stages/` (`--dump-stages`, `--keep-stages`): the instrumented R runs' layout, for validation | npy, tsv, json | |
| `nested_output` | `clusters/summary.json`, `recluster/summary.json` (`doubletfinder`, M15): a nested scrust output's own summary, whose manifest lists its files ([Nested outputs](#nested-outputs-m15)) | json | |

In `scrust cluster`'s output these files sit under `sct/`, `pca/` and `neighbors/` as listed; in `scrust sctransform`'s, `scrust pca`'s and `scrust neighbors`' outputs, at the top. A `path` is always the full relative path.

## The naming rules

They keep normalised data from being taken for raw counts (C3f):

1. A file named as a 10x reader looks for (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, `genes.tsv`, gzipped or not) has role `raw_counts`, `raw_counts_cells` or `raw_counts_features`.
2. A directory holding such a file holds nothing else: no directory holding normalised data has a 10x name, and `scrust cluster --export`'s only 10x directory is `counts/`.
3. No other file's name ends in `matrix.mtx` or `matrix.mtx.gz` (scanpy's `read_10x_mtx(path, prefix=p)` reads `<p>matrix.mtx.gz`).
4. A file with role `raw_counts` has scale `raw_counts`.

So Seurat's `Read10X` and scanpy's `read_10x_mtx`, with or without `prefix=`, fail on `sct/` and on an export's top directory, and read `counts/` (checked by `validation/R/guard_read10x.R` and `validation/python/guard_read10x.py`, M13).

The normalised matrices also say what they are: right after the Matrix Market banner, comment lines state what the matrix holds, that it is **not** raw counts, scrust's version and command (`provenance.command`), and the SHA-256 of each input file.

## For loaders (M14)

The loaders, `loaders/R/scrustio` (Seurat, SingleCellExperiment) and `loaders/python/scrustio` (AnnData), are written against this document (milestone M14; docs/PLAN.md C3g):

- Find files by `role`, never by name; refuse a manifest whose `schema` or `version` is not this one.
- Put a matrix only into a slot meant for its `scale`: `raw_counts` into counts slots (Seurat's RNA `counts`, SingleCellExperiment's `counts`, AnnData's `layers["counts"]`); `corrected_counts` and `log1p_corrected_counts` only into the SCT assay's `counts` and `data` (or `altExp(sce, "SCT")`, AnnData's `layers["sct_corrected_counts"]` and `X`); `residuals` only into `scale.data` (or its equivalents). Each slot takes exactly one role, one scale and one format; the loaders refuse any other, a Matrix Market file whose comment lines say "NOT raw counts" for a raw-count slot, and a file not so marked for a slot of SCTransform's layers.
- Take row and column names from `rows` and `columns`, and check their counts against `shape`.
- `sct_scale_data`'s rows are in `VariableFeatures()` order; Seurat keeps `scale.data`'s rows in the SCT assay's gene order (C3g), which `sct_features` gives.
- Check `sha256` before reading, when the user asks for it (a file can be large).
- A file of role `cluster_levels` holds the levels of the column it is named after (`<column>_levels.txt`). Version 1 has no field that links a levels file to its column, so this is the one place the loaders read a file's name: to know which column a file found by its role belongs to, not to find it. Fields may be added within the version: a field naming the column (`levels_of`, say) would let a loader do without the name.
- Seurat's SCT `data` layer holds the C library's `log1p` of the corrected counts (R's). The loaders read it from the `sct_data` file and check it against `log1p` of `sct_counts`; Python must use `math.log1p` for that check, since numpy's `np.log1p` can differ in the last bit (M14).
