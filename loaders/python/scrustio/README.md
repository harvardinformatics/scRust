# scrustio (Python): scrust outputs as AnnData objects

`scrustio` reads the plain files that scRust's command-line tool, `scrust`, writes, and builds AnnData objects. It finds every file through the manifest in the output's `summary.json` (schema `scrust-files`, version 1: [`docs/output_manifest.md`](../../../docs/output_manifest.md)), by the file's role, never by its name. It puts each matrix only into a slot meant for its scale: the raw counts into `layers["counts"]`, and SCTransform's corrected counts, their log1p and the residuals only into slots for normalised data (docs/PLAN.md C3g, milestone M14).

## What it builds

From `scrust cluster --export`:

| slot | holds |
|---|---|
| `layers["counts"]` | the raw counts (`counts/`), cells x genes |
| `X` | SCTransform's log1p corrected counts (Seurat's SCT `data`), with scanpy's `uns["log1p"]` |
| `layers["sct_corrected_counts"]` | SCTransform's corrected counts (Seurat's SCT `counts`); `X` and this layer span every gene, and genes outside SCTransform's model are 0 there, flagged by `var["sct_modelled"]` |
| `obsm["sct_scale_data"]` | the residuals (Seurat's SCT `scale.data`), a data frame with the variable features as columns, in `VariableFeatures()` order; `var["sct_variable"]` flags them |
| `uns["sct_model"]` | SCTransform's model as Seurat's SCTModel holds it: `feature_attributes`, `cell_attributes`, `clips`, `umi_assay`, `model`, `arguments`, `median_umi`, and `variable_features` in order |
| `obs` | the per-cell columns (`metadata.tsv`), categoricals in their levels' order |
| `obsm["X_pca"]`, `varm["PCs"]`, `uns["pca"]` | the PCA's embeddings; its loadings (0 for genes outside the PCA); its standard deviations (`stdev`, Seurat's), `variance`, `variance_ratio`, `total_variance` and `features` (the loadings' genes, in order) |
| `obsp["sct_nn"]`, `obsp["sct_snn"]` | Seurat's `SCT_nn` and `SCT_snn` graphs; `uns["sct_snn"]` names `sct_snn` as the connectivities, for scanpy's `neighbors_key` (or `obsp=` in newer scanpy) |

Genes are named by the features file's names, made unique as scanpy does (`var_names_make_unique`); `var["gene_ids"]` and `var["feature_types"]` hold the rest. From `scrust sctransform --corrected-counts` (no raw counts there) `var` holds the model genes; from `scrust doubletfinder`, `layers["counts"]` holds the singlets' raw counts, `obs` DoubletFinder's calls (named as DoubletFinder names them), and `X` stays empty: only SCTransform's log1p data go there.

**What is not rebuilt:** the UMAP; metadata columns of other workflow steps (`percent.mt` is loaded, as a float column of `obs`, when the output computed it: `scrust cluster --counts` and `scrust sctransform` with `--percent-mt`, `--mt-genes` or `--vars-to-regress percent.mt`; milestone M15b); and `uns["sct_model"]["cell_attributes"]["cells_step1"]` for objects above 5,000 cells (vst() drew those cells with R's `sample()`, which this loader does not reproduce; the R loader does).

**Seurat's data layer is read, not recomputed.** Its values are the C library's `log1p` (R's). numpy's `np.log1p` may differ in the last bit (in 13,457 of 65,546 values of one test export, numpy 2.5.5), so `X` must come from the file; the loader checks it against `math.log1p` of the corrected counts.

## Install

```bash
pip install loaders/python/scrustio
```

It needs numpy, scipy, pandas and anndata (tested with anndata 0.13.4, pandas 3.1.1, numpy 2.5.5, scipy 1.18.1, Python 3.13; scRust's hermetic env is `envs/loaders-py.yml`).

## Use

```python
import scrustio
adata = scrustio.load_anndata("export")
m = scrustio.read_manifest("export")             # the manifest: files by role (m.entry("sct_counts"))
adata = scrustio.load_anndata("export", check_sha256=True)
```

As one line, writing `.h5ad`:

```bash
python -m scrustio h5ad export object.h5ad
```

(`scrustio h5ad ...` once installed; `python -m scrustio manifest export` lists the files by role.)

## The rule, enforced

`scrustio.slot_table()` lists the slots filled from matrix files, each with the one role, scale and format it takes; `scrustio.read_for_slot()` refuses any other, a Matrix Market file whose comment lines say "NOT raw counts" for `layers["counts"]`, and an unmarked file for `X` or `layers["sct_corrected_counts"]`. `scrustio.read_manifest()` refuses another schema or version, a listed file that is missing or has another size (or, with `check_sha256=True`, another SHA-256), and a path outside the output; the identifiers of each matrix's rows and columns are checked against its shape.

## Tests

```bash
cd loaders/python/scrustio && PYTHONNOUSERSITE=1 python -I -m pytest
```

`pyproject.toml` puts the package on the path for pytest, so nothing is installed. The tests make outputs with the `scrust` binary (`SCRUST_BIN`, else `target/release/scrust`) from the committed fixtures, and check every wrong combination of role, scale and format for every slot, tampered manifests, the manifest's rules, the AnnData's layout, the other commands' outputs and the `.h5ad` round trip. Milestone M14's validation compares the objects with the saved workflow objects (docs/validation/clustering.md, "M14").
