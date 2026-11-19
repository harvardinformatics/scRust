"""AnnData objects from scrust outputs, as docs/PLAN.md C3g lays them out.

- ``layers["counts"]``: the raw counts (cells x genes);
- ``X``: SCTransform's log1p corrected counts (Seurat's SCT ``data``), with scanpy's
  ``uns["log1p"]``; ``layers["sct_corrected_counts"]``: the corrected counts. Both span every gene,
  and genes outside SCTransform's model are 0 there, flagged by ``var["sct_modelled"]``;
- ``obsm["sct_scale_data"]``: the residuals (Seurat's SCT ``scale.data``), a data frame with the
  variable features as columns, in ``VariableFeatures()`` order; ``var["sct_variable"]``;
- ``uns["sct_model"]``: SCTransform's model (Seurat's SCTModel);
- ``obs``: the per-cell columns, categoricals in their levels' order;
- ``obsm["X_pca"]``, ``varm["PCs"]`` (0 for genes outside the PCA), ``uns["pca"]``;
- ``obsp["sct_nn"]`` and ``obsp["sct_snn"]``, with ``uns["sct_snn"]`` for scanpy's
  ``neighbors_key``.

Genes are named by the features file's names, made unique as scanpy does
(``var_names_make_unique``); ``var["gene_ids"]`` and ``var["feature_types"]`` hold the rest.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .manifest import Manifest, ScrustError, read_manifest
from .read import read_first_line, read_lines, read_tsv
from .sctmodel import model_parts, sct_tables
from .slots import check_slot, read_for_slot

META_DOUBLE = ("nCount_RNA", "percent.mt", "nCount_SCT", "silhouette_width", "neighborhood_purity")
META_INTEGER = ("nFeature_RNA", "nFeature_SCT")


def libm_log1p(x: np.ndarray) -> np.ndarray:
    """log1p as the C library computes it (Python's math.log1p), which R's log1p() and scrust use
    for Seurat's SCT data layer. numpy's np.log1p may use its own vectorised code, which differs in
    the last bit for some values (13,457 of 65,546 on one tier-T export with numpy 2.5.5), so
    Seurat's data layer must be read, not recomputed with numpy. Computed once per distinct value."""
    import math

    u, inv = np.unique(x, return_inverse=True)
    lu = np.fromiter((math.log1p(v) for v in u), dtype=np.float64, count=len(u))
    return lu[inv].reshape(x.shape)


def _make_unique(names) -> pd.Index:
    """Names made unique as anndata's var_names_make_unique() does ("-1", "-2", ...)."""
    import anndata

    idx = pd.Index(list(names), dtype=object)
    return anndata.utils.make_index_unique(idx) if idx.has_duplicates else idx


def counts_features(m: Manifest):
    """The raw counts' features (role raw_counts_features): ids, names, types; None if absent."""
    e = m.entry("raw_counts_features", required=False)
    if e is None:
        return None
    rows = [line.split("\t") for line in read_lines(m.path(e))]
    ids = [r[0] for r in rows]
    names = [r[1] if len(r) > 1 else r[0] for r in rows]
    types = [r[2] if len(r) > 2 else "" for r in rows]
    return {"id": ids, "name": names, "type": types, "var_names": _make_unique(names)}


def map_features(feat, names, ids=None, what="features") -> np.ndarray:
    """Positions of a list of features (by name, or by id where names repeat) among the counts'."""
    if feat is None:
        raise ScrustError("this output holds no raw counts' features to map onto")
    if len(set(feat["name"])) == len(feat["name"]):
        index, keys = pd.Index(feat["name"]), names
    elif ids is not None and len(set(feat["id"])) == len(feat["id"]):
        index, keys = pd.Index(feat["id"]), ids
    else:
        raise ScrustError(f"the counts' feature names repeat and the {what} cannot be matched to them by id")
    k = index.get_indexer(list(keys))
    if (k < 0).any():
        first = list(names)[int(np.flatnonzero(k < 0)[0])]
        raise ScrustError(f"{int((k < 0).sum())} of the {what} are not among the counts' features (first: {first})")
    return k


def cells_of(m: Manifest) -> list:
    for role in ("raw_counts", "sct_counts", "cell_metadata", "rna_cell_metadata"):
        e = m.entry(role, required=False)
        if e is None:
            continue
        ids = e["columns"] if role in ("raw_counts", "sct_counts") else e["rows"]
        return m.identifiers(ids)
    raise ScrustError(f"this {m.command} output names no cells")


def cluster_levels(m: Manifest) -> dict:
    """Factor levels from the files of role cluster_levels: each named after its column
    (``<column>_levels.txt``, docs/output_manifest.md), one level per line."""
    out = {}
    for e in m.entries_by_role("cluster_levels"):
        base = e["path"].rsplit("/", 1)[-1]
        if not base.endswith("_levels.txt"):
            raise ScrustError(f"levels file {e['path']} is not named <column>_levels.txt")
        out[base[: -len("_levels.txt")]] = read_lines(m.path(e))
    return out


def _categorical(values, levels, col):
    bad = set(pd.unique(pd.Series(values).dropna())) - set(levels)
    if bad:
        raise ScrustError(f"column {col} holds labels missing from its levels: {sorted(bad)[:5]}")
    return pd.Categorical(values, categories=levels)


def read_metadata(m: Manifest) -> pd.DataFrame | None:
    """scrust cluster's metadata.tsv (or scrust sctransform's cell_metadata.tsv), typed as Seurat
    holds it: doubles, R integers (int32), and factors in their levels' order."""
    e = m.entry("cell_metadata", required=False) or m.entry("rna_cell_metadata", required=False)
    if e is None:
        return None
    m.check_shape(e)
    path = m.path(e)
    header = read_first_line(path).split("\t")
    levels = cluster_levels(m)
    numeric = [c for c in header if c in META_DOUBLE]
    integer = [c for c in header if c in META_INTEGER]
    df = read_tsv(path, numeric=numeric, integer=integer)
    df.index = pd.Index(df["barcode"].tolist(), dtype=object)
    df = df.drop(columns="barcode")
    for col in list(df.columns):
        if col in levels:
            df[col] = _categorical(df[col].tolist(), levels[col], col)
        elif col in ("silhouette_other", "purity_maximum") and "seurat_clusters" in levels:
            vals = [None if v == "NA" else v for v in df[col].tolist()]
            df[col] = _categorical(vals, levels["seurat_clusters"], col)
        elif col not in numeric and col not in integer:
            conv = pd.to_numeric(df[col].replace("NA", np.nan), errors="coerce")
            if conv.notna().sum() == (df[col] != "NA").sum():
                df[col] = conv
            else:
                df[col] = df[col].astype(object)
    return df


def read_doublet_calls(m: Manifest) -> pd.DataFrame | None:
    """scrust doubletfinder's calls, named as DoubletFinder 2.0.6 names its metadata columns
    (pANN_<pN>_<pK>_<nExp> and DF.classifications_<pN>_<pK>_<nExp>, which the summary gives)."""
    e = m.entry("doublet_calls", required=False)
    if e is None:
        return None
    m.check_shape(e)
    df = read_tsv(m.path(e), numeric=["pANN"])
    res = m.summary.get("result") or {}
    pann = res.get("pANN_column") or "pANN"
    cls = res.get("classification_column") or "DF.classifications"
    out = pd.DataFrame({pann: df["pANN"].to_numpy(dtype=np.float64), cls: df["classification"].astype(object).to_numpy()},
                       index=pd.Index(df["barcode"].tolist(), dtype=object))
    return out


def _embed(mat_genes_by_cells: sp.csc_matrix, pos: np.ndarray, n_var: int) -> sp.csr_matrix:
    """A genes x cells matrix of some genes as cells x all genes (CSR), the genes at ``pos``."""
    t = mat_genes_by_cells.T.tocsr()
    t.sort_indices()
    out = sp.csr_matrix((t.data.astype(np.float64, copy=False), pos[t.indices].astype(np.int32), t.indptr.copy()),
                        shape=(t.shape[0], n_var))
    out.has_sorted_indices = True
    return out


def load_anndata(path, check_sha256: bool = False):
    """Build an AnnData from a scrust output.

    From ``scrust cluster --export`` (everything above), ``scrust sctransform --corrected-counts``
    (the SCT data; no raw counts there, so ``var`` holds the model genes) or ``scrust
    doubletfinder`` (the singlets' raw counts in ``layers["counts"]`` and DoubletFinder's calls in
    ``obs``; ``X`` stays empty: only SCTransform's log1p data go there). Files are found through
    the manifest, by role, and each matrix goes only into a slot meant for its scale.
    """
    import anndata

    m = path if isinstance(path, Manifest) else read_manifest(path, check_sha256=check_sha256)
    feat = counts_features(m)
    if feat is None and not m.has_role("sct_counts"):
        raise ScrustError(f"this {m.command} output holds neither raw counts nor SCTransform's layers "
                          "(scrust cluster --export holds both)")
    cells = cells_of(m)
    obs = read_metadata(m)
    calls = read_doublet_calls(m)
    if calls is not None:
        missing = pd.Index(cells).difference(calls.index)
        if len(missing):
            raise ScrustError(f"{len(missing)} cells have no DoubletFinder call")
        calls = calls.loc[cells]
        obs = calls if obs is None else pd.concat([obs, calls], axis=1)
    if obs is None:
        obs = pd.DataFrame(index=pd.Index(cells, dtype=object))
    if list(obs.index) != list(cells):
        raise ScrustError("the metadata's cells differ from the output's")
    obs.index = pd.Index(cells, dtype=object)

    layers, obsm, varm, obsp, uns = {}, {}, {}, {}, {}
    X = None
    tabs = None
    if feat is not None:
        var = pd.DataFrame({"gene_ids": pd.array(feat["id"], dtype=object),
                            "feature_types": pd.array(feat["type"], dtype=object)}, index=feat["var_names"])
        counts = read_for_slot(m, "raw_counts", 'AnnData layers["counts"]')
        if m.identifiers(m.entry("raw_counts")["columns"]) != list(cells):
            raise ScrustError("the raw counts' cells differ from the output's")
        layers["counts"] = sp.csr_matrix(counts.T.tocsr())
        layers["counts"].sort_indices()
    else:
        var = None
    n_var = None
    if m.has_role("sct_counts"):
        tabs = sct_tables(m)
        ga = tabs["gene_attr"]
        if m.identifiers(m.entry("sct_counts")["rows"]) != ga["name"].tolist():
            raise ScrustError("the SCT features (sct_features) differ from the model's genes (gene_attr)")
        if var is None:
            var = pd.DataFrame({"gene_ids": pd.array(ga["id"].tolist(), dtype=object)}, index=_make_unique(ga["name"]))
            pos = np.arange(len(ga))
        else:
            pos = map_features(feat, ga["name"].tolist(), ga["id"].tolist(), "model genes (gene_attr)")
        n_var = len(var)
        cc = read_for_slot(m, "sct_counts", 'AnnData layers["sct_corrected_counts"]')
        dd = read_for_slot(m, "sct_data", "AnnData X")
        if not (np.array_equal(cc.indices, dd.indices) and np.array_equal(cc.indptr, dd.indptr)):
            raise ScrustError("the SCT data layer's entries differ from the counts layer's")
        expected = libm_log1p(cc.data)
        if not np.array_equal(dd.data, expected):
            raise ScrustError(f"the SCT data layer is not log1p of the counts layer in {int((dd.data != expected).sum())} values")
        X = _embed(dd, pos, n_var)
        layers["sct_corrected_counts"] = _embed(cc, pos, n_var)
        modelled = np.zeros(n_var, dtype=bool)
        modelled[pos] = True
        var["sct_modelled"] = modelled
        uns["log1p"] = {"base": None}
        vf = read_tsv(m.path(m.entry("sct_variable_features")))
        vf_pos = pos[pd.Index(ga["name"].tolist()).get_indexer(vf["name"].tolist())] if feat is None else \
            map_features(feat, vf["name"].tolist(), vf["id"].tolist(), "variable features")
        if (vf_pos < 0).any():
            raise ScrustError("variable features outside the model")
        vf_names = var.index[vf_pos]
        sd = read_for_slot(m, "sct_scale_data", 'AnnData obsm["sct_scale_data"]')
        if m.identifiers(m.entry("sct_scale_data")["rows"]) != vf["name"].tolist():
            raise ScrustError("scale.data's rows differ from the variable features")
        if m.identifiers(m.entry("sct_scale_data")["columns"]) != list(cells):
            raise ScrustError("scale.data's columns differ from the cells")
        obsm["sct_scale_data"] = pd.DataFrame(np.ascontiguousarray(sd.T), index=pd.Index(cells, dtype=object),
                                              columns=pd.Index(list(vf_names), dtype=object))
        variable = np.zeros(n_var, dtype=bool)
        variable[vf_pos] = True
        var["sct_variable"] = variable
        step1 = tabs["step1"]
        step1_pos = pos[pd.Index(ga["name"].tolist()).get_indexer(step1["name"].tolist())] if feat is None else \
            map_features(feat, step1["name"].tolist(), step1["id"].tolist(), "step-1 genes (model_pars)")
        parts = model_parts(m, tabs, list(var.index[pos]), list(var.index[step1_pos]), cells)
        parts["variable_features"] = np.array(list(vf_names), dtype=object)
        uns["sct_model"] = parts
    if m.has_role("pca_embeddings"):
        for s in ('AnnData obsm["X_pca"]', 'AnnData varm["PCs"]', 'AnnData uns["pca"]["stdev"]'):
            from .slots import SLOTS

            e = m.entry(SLOTS[s]["role"])
            check_slot(s, e["role"], e["scale"], e["path"], e["format"])
            m.check_shape(e)
        e_emb = m.entry("pca_embeddings")
        pcs = read_first_line(m.path(e_emb)).split("\t")[1:]
        emb = read_tsv(m.path(e_emb), numeric=pcs)
        if emb["barcode"].tolist() != list(cells):
            raise ScrustError("the PCA embeddings' cells differ")
        obsm["X_pca"] = np.ascontiguousarray(emb[pcs].to_numpy(dtype=np.float64))
        ld = read_tsv(m.path(m.entry("pca_loadings")), numeric=pcs)
        ld_pos = map_features(feat, ld["name"].tolist(), ld["id"].tolist(), "PCA features") if feat is not None else \
            pd.Index(ga["name"].tolist()).get_indexer(ld["name"].tolist())
        if var is None:
            raise ScrustError("a PCA without genes")
        pcs_mat = np.zeros((len(var), len(pcs)), dtype=np.float64)
        pcs_mat[ld_pos, :] = ld[pcs].to_numpy(dtype=np.float64)
        varm["PCs"] = pcs_mat
        sdev = read_tsv(m.path(m.entry("pca_stdev")), numeric=["stdev"])
        if sdev["pc"].tolist() != pcs:
            raise ScrustError("the standard deviations' components differ from the embeddings'")
        stdev = sdev["stdev"].to_numpy(dtype=np.float64)
        res = m.summary.get("result") or {}
        tv = (res.get("pca") or {}).get("total_variance", res.get("total_variance"))
        uns["pca"] = {"stdev": stdev, "variance": stdev * stdev,
                      "features": np.array(list(var.index[ld_pos]), dtype=object)}
        if tv is not None:
            uns["pca"]["total_variance"] = float(tv)
            uns["pca"]["variance_ratio"] = stdev * stdev / float(tv)
    for role, key in (("graph_nn", "sct_nn"), ("graph_snn", "sct_snn")):
        if not m.has_role(role):
            continue
        e = m.entry(role)
        if m.identifiers(e["rows"]) != list(cells) or m.identifiers(e["columns"]) != list(cells):
            raise ScrustError(f"graph {key}'s cells differ")
        g = read_for_slot(m, role, f'AnnData obsp["{key}"]')
        g = sp.csr_matrix(g.tocsr())
        g.sort_indices()
        obsp[key] = g
    if "sct_snn" in obsp:
        par = (m.provenance.get("parameters") or {}).get("neighbors") or {}
        params = {"method": "Seurat FindNeighbors: shared nearest neighbours (Jaccard), pruned",
                  "use_rep": "X_pca", "graph_name": "SCT_snn"}
        for k_src, k_dst in (("k_param", "n_neighbors"), ("metric", "metric"), ("prune_snn", "prune_snn")):
            if par.get(k_src) is not None:
                params[k_dst] = par[k_src]
        if par.get("dims") is not None:
            params["n_pcs"] = len(par["dims"])
        uns["sct_snn"] = {"connectivities_key": "sct_snn", "params": params}
    if var is None:
        var = pd.DataFrame(index=pd.Index([], dtype=object))
    uns["scrust"] = {"command": m.command, "scrust_version": str(m.provenance.get("scrust_version", "")),
                     "command_line": np.array([str(x) for x in m.provenance.get("command_line", [])], dtype=object)}
    adata = anndata.AnnData(X=X, obs=obs, var=var, layers=layers or None, obsm=obsm or None, varm=varm or None,
                            obsp=obsp or None, uns=uns)
    # Text columns with repeated values become categoricals, as anndata makes them when it writes
    # .h5ad, so that the object written and read back equals this one.
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata.strings_to_categoricals()
    return adata


def write_h5ad(path: str, out: str, check_sha256: bool = False):
    """Build the AnnData (load_anndata) and write it as .h5ad. As one line:

        python -m scrustio h5ad <scrust output> <object.h5ad>
    """
    adata = load_anndata(path, check_sha256=check_sha256)
    adata.write_h5ad(out)
    return adata
