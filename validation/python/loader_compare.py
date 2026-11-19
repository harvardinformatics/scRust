"""K11 for the Python loader (docs/PLAN.md C3g, C4 "Loader tests", C5; milestone M14).

Builds an AnnData (scrustio.load_anndata) from one scrust output and compares it with the saved
workflow object, through M8's export of that object (data/clustering/<object>/saved/: the SCT
layers, scale.data, the SCTModel, the PCA, the graphs and the metadata, in plain files; input/: its
RNA counts layer), component by component. Grading is M11-M13's: exact when the export's recorded
OpenBLAS kernel is the saved objects' (SapphireRapids); otherwise the BLAS-dependent components may
differ within [clustering.floor.kernel] of validation/tolerances.toml (PCA components aligned in
sign first), and the rest stay exact. scale.data and the log1p layer of the objects whose export
kept only their digest (M8's levels standard and light) are compared by that digest; in floor
mode a differing digest is reported, not graded (it has no magnitude).

Also: the loader loads exactly what the files hold (every Matrix Market file parsed again with
Python's own float(), every table's numbers likewise, .npy files as bytes), and the one-line
command's .h5ad reloads into an AnnData equal to the loader's (anndata's read_h5ad).

Run in the hermetic env (envs/loaders-py.yml), isolated from the user's site:
    PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I validation/python/loader_compare.py \\
        --object <name> --export <scrust output> --saved <data/clustering/<object>/saved> \\
        --input <data/clustering/<object>/input> --out <dir> [--kind cluster|doubletfinder] [--tmp <dir>]
Output in <dir>: k11_anndata.tsv, files_anndata.tsv, roundtrip_anndata.tsv, compare_anndata.json.
Exit status 1 if a component fails.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import resource
import sys
import time
import tomllib

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "loaders", "python", "scrustio"))
import scrustio  # noqa: E402
from scrustio.read import read_mm  # noqa: E402

CHAIN = ["nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "seurat_clusters",
         "silhouette_width", "neighborhood_purity"]
COL_FLOOR = {"silhouette_width": "K8_silhouette_width", "neighborhood_purity": "K9_purity"}
BLAS_FREE_COLS = {"nCount_RNA", "nFeature_RNA", "percent.mt"}


def peak_rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20


def num_stats(a, b, align=False) -> dict:
    a = np.array(a, dtype=np.float64, copy=True)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        return {"n": int(b.size), "comparable": False}
    flips = 0
    if align and a.ndim == 2:
        for k in range(a.shape[1]):
            if np.sum((a[:, k] + b[:, k]) ** 2) < np.sum((a[:, k] - b[:, k]) ** 2):
                a[:, k] = -a[:, k]
                flips += 1
    na_a, na_b = np.isnan(a), np.isnan(b)
    both = ~na_a & ~na_b
    d = np.abs(a[both] - b[both])
    differing = int((na_a != na_b).sum() + (a[both] != b[both]).sum())
    scale = float(np.max(np.abs(b[both]))) if both.any() else 0.0
    max_abs = float(d.max()) if d.size else 0.0
    return {"n": int(b.size), "comparable": True, "differing": differing, "max_abs": max_abs,
            "max_scaled": (max_abs / scale) if scale > 0 else (math.inf if max_abs > 0 else 0.0),
            "fraction_differing": differing / b.size if b.size else 0.0, "flipped": flips}


def same_array(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if a.dtype.kind == "f":
        return bool(np.array_equal(a, b, equal_nan=True)) and bool(np.array_equal(np.signbit(a), np.signbit(b)))
    return bool(np.array_equal(a, b))


def same_sparse(a, b) -> bool:
    a, b = sp.csr_matrix(a), sp.csr_matrix(b)
    return a.shape == b.shape and same_array(a.indptr, b.indptr) and same_array(a.indices, b.indices) and \
        np.array_equal(a.data.astype(np.float64), b.data.astype(np.float64))


class Grader:
    def __init__(self, obj: str, mode: str, floor: dict):
        self.obj, self.mode, self.floor, self.rows = obj, mode, floor, []

    def grade(self, ident, stats, floor_name, blas) -> str:
        if ident:
            return "identical"
        fe = self.floor.get(floor_name) if (floor_name and blas) else None
        if self.mode != "floor" or fe is None or stats is None or not stats.get("comparable", True):
            return "differing"
        ok = True
        for k, v in fe.items():
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                continue
            got = {"max_abs": stats.get("max_abs"), "max_scaled": stats.get("max_scaled"),
                   "max_fraction_differing": stats.get("fraction_differing"),
                   "values_differing": stats.get("differing"), "cells_differing": stats.get("differing")}.get(k)
            if got is None:
                continue
            if not math.isfinite(got) or got > v:
                ok = False
        return "within floor" if ok else "beyond floor"

    def add(self, component, ident, stats=None, floor_name="", blas=True, detail="", status=None):
        st = status if status is not None else self.grade(ident, stats, floor_name, blas)
        s = stats or {}
        self.rows.append({"object": self.obj, "loader": "anndata", "component": component,
                          "class": "blas_dependent" if blas else "blas_free", "floor": floor_name, "status": st,
                          "n": s.get("n"), "differing": s.get("differing"), "max_abs": s.get("max_abs"),
                          "max_scaled": s.get("max_scaled"), "detail": detail})


def read_saved_meta(saved: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(saved, "metadata.tsv"), sep="\t", dtype=str, keep_default_na=False, quoting=3)
    return df


def sha256_doubles(x: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(x, dtype="<f8").tobytes()).hexdigest()


def compare_cluster(g: Grader, a, saved: str, input_dir: str):
    obj_json = json.load(open(os.path.join(saved, "object.json")))
    # raw counts
    want, _, _ = read_mm(os.path.join(input_dir, "matrix.mtx.gz"))
    got = a.layers["counts"]
    g.add("layers counts", same_sparse(got, want.T), blas=False, detail="against the saved RNA counts (input/)")
    # obs: the chain's columns
    meta = read_saved_meta(saved)
    g.add("obs names (cells)", list(a.obs_names) == meta["barcode"].tolist(), blas=False)
    levels = {c: open(os.path.join(saved, f"{c}_levels.txt")).read().splitlines() for c in ("seurat_clusters", "SCT_snn_res.0.8")}
    # M15b: percent.mt, when the export computed it and the saved object's export holds it
    extra = ["percent.mt"] if "percent.mt" in a.obs and "percent.mt" in meta else []
    for col in CHAIN + extra:
        x = a.obs[col]
        if col in levels:
            ident = list(x.cat.categories) == levels[col] and x.astype(str).tolist() == meta[col].tolist()
            g.add(f"obs {col}", ident, {"n": len(x), "differing": int((x.astype(str).to_numpy() != meta[col].to_numpy()).sum())},
                  blas=col not in BLAS_FREE_COLS, detail="values and categories (levels) in order")
        else:
            y = pd.to_numeric(meta[col].replace("NA", np.nan)).to_numpy(dtype=np.float64)
            xv = x.to_numpy(dtype=np.float64)
            g.add(f"obs {col}", same_array(xv, y), num_stats(xv, y), floor_name=COL_FLOOR.get(col, ""),
                  blas=col not in BLAS_FREE_COLS)
    # SCT layers, in the model genes' columns
    modelled = a.var["sct_modelled"].to_numpy()
    sct_feat = [l for l in gzip.open(os.path.join(saved, "sct", "sct_features.tsv.gz"), "rt").read().splitlines()]
    g.add("var sct_modelled (the model genes)", list(a.var_names[modelled]) == sct_feat, blas=True)
    cc = sp.csc_matrix(a.layers["sct_corrected_counts"][:, modelled].T)
    cc.sort_indices()
    want_cc, _, _ = read_mm(os.path.join(saved, "sct", "corrected_counts.mtx.gz"))
    g.add("layers sct_corrected_counts", same_sparse(cc, want_cc), blas=True)
    g.add("layers sct_corrected_counts (R's digest)", sha256_doubles(cc.data) == obj_json["sct_counts_x_sha256"], blas=True)
    dd = sp.csc_matrix(a.X[:, modelled].T)
    dd.sort_indices()
    p = os.path.join(saved, "sct", "log1p_corrected_counts.mtx.gz")
    if os.path.exists(p):
        want_dd, _, _ = read_mm(p)
        g.add("X (log1p corrected counts)", same_sparse(dd, want_dd), blas=True)
    g.add("X (R's digest)", sha256_doubles(dd.data) == obj_json["sct_data_x_sha256"], blas=True)
    g.add("uns log1p", a.uns.get("log1p") == {"base": None}, blas=False, detail="scanpy's mark of a natural-log X")
    # scale.data: saved rows in the SCT assay's order (scale_data_features.tsv)
    rows = open(os.path.join(saved, "sct", "scale_data_features.tsv")).read().splitlines()
    sd = a.obsm["sct_scale_data"]
    vf_saved = open(os.path.join(saved, "sct", "variable_features.tsv")).read().splitlines()
    g.add("obsm sct_scale_data columns (VariableFeatures order)", list(sd.columns) == vf_saved, blas=False)
    g.add("var sct_variable", set(a.var_names[a.var["sct_variable"].to_numpy()]) == set(vf_saved), blas=False)
    sd_saved_order = sd[rows].to_numpy(dtype=np.float64)          # cells x features (SCT order)
    p = os.path.join(saved, "sct", "scale_data.npy")
    if os.path.exists(p):
        want = np.load(p).T                                         # features x cells -> cells x features
        g.add("obsm sct_scale_data", same_array(sd_saved_order, want), num_stats(sd_saved_order, want),
              floor_name="S2d_scale_data")
    else:
        ident = sha256_doubles(sd_saved_order) == obj_json["scale_data_sha256"]
        g.add("obsm sct_scale_data (R's digest)", ident, None, blas=True,
              status=("identical" if ident else ("digest differs (reported, not graded)" if g.mode == "floor" else "differing")),
              detail="M8 kept only the digest of this object's scale.data")
    # the model
    model = a.uns["sct_model"]
    fa = model["feature_attributes"]
    fa_saved = pd.read_csv(os.path.join(saved, "sct", "model", "feature_attributes.tsv"), sep="\t", dtype=str,
                           keep_default_na=False, quoting=3)
    g.add("sct_model feature_attributes rows and columns",
          list(fa.index) == fa_saved["name"].tolist() and list(fa.columns) == list(fa_saved.columns[1:]), blas=False)
    # the floor's S2d_residual_stats is one array of both columns, cbind(residual_mean, residual_variance)
    rcols = ["residual_mean", "residual_variance"]
    rs = num_stats(fa[rcols].to_numpy(dtype=np.float64),
                   np.column_stack([pd.to_numeric(fa_saved[c].replace("NA", np.nan)).to_numpy(dtype=np.float64) for c in rcols]))
    for col in fa_saved.columns[1:]:
        if col == "genes_log_gmean_step1":
            g.add(f"sct_model feature_attributes {col}", fa[col].to_numpy().tolist() == (fa_saved[col] == "TRUE").tolist())
            continue
        y = pd.to_numeric(fa_saved[col].replace("NA", np.nan)).to_numpy(dtype=np.float64)
        x = fa[col].to_numpy(dtype=np.float64)
        fl = "S2d_residual_stats" if col in rcols else ""
        g.add(f"sct_model feature_attributes {col}", same_array(x, y), rs if col in rcols else num_stats(x, y), floor_name=fl,
              blas=col not in ("detection_rate", "gmean", "variance"),
              detail="statistics over residual_mean and residual_variance together, as the floor's S2d_residual_stats" if col in rcols else "")
    ca = model["cell_attributes"]
    ca_saved = pd.read_csv(os.path.join(saved, "sct", "model", "cell_attributes.tsv"), sep="\t", dtype=str,
                           keep_default_na=False, quoting=3)
    ok = list(ca.index) == ca_saved["barcode"].tolist() and \
        same_array(ca["umi"].to_numpy(), pd.to_numeric(ca_saved["umi"]).to_numpy(dtype=np.float64)) and \
        same_array(ca["log_umi"].to_numpy(), pd.to_numeric(ca_saved["log_umi"]).to_numpy(dtype=np.float64))
    g.add("sct_model cell_attributes umi, log_umi", ok, blas=False)
    if "cells_step1" in ca.columns:
        g.add("sct_model cell_attributes cells_step1", ca["cells_step1"].tolist() == (ca_saved["cells_step1"] == "TRUE").tolist(), blas=False)
    else:
        g.add("sct_model cell_attributes cells_step1", None, status="not rebuilt", blas=False,
              detail="vst() drew the 5,000 step-1 cells with R's sample(); the Python loader does not reproduce R's generator")
    mj = json.load(open(os.path.join(saved, "sct", "model", "model.json")))
    clips_ok = all(same_array(np.asarray(model["clips"][k]), np.asarray(mj["clips"][k], dtype=np.float64)) for k in ("vst", "sct"))
    g.add("sct_model clips", clips_ok)
    g.add("sct_model umi_assay, model", model["umi_assay"] == mj["umi.assay"] and model["model"] == mj["model"], blas=False)
    g.add("sct_model median_umi", same_array(np.float64(model["median_umi"]), np.float64(mj["median_umi"])), blas=False)
    args_ok, bad = True, []
    if list(model["arguments"].keys()) != list(mj["arguments"].keys()):
        args_ok, bad = False, ["names or order"]
    for k, v in mj["arguments"].items():
        mine = model["arguments"].get(k)
        if v is None:
            same = mine is None
        elif isinstance(v, list):
            same = mine is not None and same_array(np.asarray(mine, dtype=np.float64), np.asarray(v, dtype=np.float64))
        elif isinstance(v, bool):
            same = mine is v
        elif isinstance(v, (int, float)):
            same = isinstance(mine, float) and mine == float(v)
        else:
            same = mine == v
        if not same:
            args_ok = False
            bad.append(k)
    g.add("sct_model arguments", args_ok, detail=", ".join(bad))
    # PCA
    emb = a.obsm["X_pca"]
    want = np.load(os.path.join(saved, "pca", "embeddings.npy"))
    k30 = slice(0, min(30, want.shape[1]))
    g.add("obsm X_pca, dims 1:30", same_array(emb[:, k30], want[:, k30]), num_stats(emb[:, k30], want[:, k30], align=True),
          floor_name="S3_embeddings_dims_1_30")
    if want.shape[1] > 30:
        g.add("obsm X_pca, dims 31:50", same_array(emb[:, 30:], want[:, 30:]), num_stats(emb[:, 30:], want[:, 30:], align=True),
              floor_name="S3_embeddings_dims_31_50")
    feats = list(a.uns["pca"]["features"])
    lf = open(os.path.join(saved, "pca", "loadings_features.tsv")).read().splitlines()
    g.add("uns pca features (the loadings' rows)", feats == lf, blas=False)
    pos = a.var_names.get_indexer(feats)
    ld = a.varm["PCs"][pos, :]
    want = np.load(os.path.join(saved, "pca", "loadings.npy"))
    g.add("varm PCs (the PCA's genes)", same_array(ld, want), num_stats(ld, want, align=True), floor_name="S3_loadings")
    others = np.delete(a.varm["PCs"], pos, axis=0)
    g.add("varm PCs (other genes: 0)", bool((others == 0).all()), blas=False)
    want = np.load(os.path.join(saved, "pca", "stdev.npy"))
    sd_ = a.uns["pca"]["stdev"]
    g.add("uns pca stdev", same_array(sd_, want), num_stats(sd_, want), floor_name="S3_sdev")
    pj = json.load(open(os.path.join(saved, "pca", "pca.json")))
    g.add("uns pca total_variance", a.uns["pca"]["total_variance"] == pj["total_variance"])
    # graphs
    for key, f in (("sct_nn", "nn"), ("sct_snn", "snn")):
        want, _, _ = read_mm(os.path.join(saved, "graphs", f"{f}.mtx.gz"))
        g.add(f"obsp {key}", same_sparse(a.obsp[key], want), blas=True)
    g.add("uns sct_snn (scanpy's neighbors_key)", a.uns["sct_snn"].get("connectivities_key") == "sct_snn", blas=False)
    # not rebuilt
    g.add("umap", None, status="not rebuilt", blas=False, detail="RunUMAP is not part of the chain scrust reproduces")
    g.add("other metadata columns", None, status="not rebuilt", blas=False, detail="columns of other workflow steps")


def compare_doubletfinder(g: Grader, a, saved: str, m):
    """A scrust doubletfinder output against M6b's/M7's export of the saved DoubletFinder object
    (data/doubletfinder/<input>/saved/: singlets/ and cells.tsv with the saved pANN and calls)."""
    want, _, _ = read_mm(os.path.join(saved, "singlets", "matrix.mtx.gz"))
    g.add("layers counts (the singlets)", same_sparse(a.layers["counts"], want.T), blas=False)
    cells = pd.read_csv(os.path.join(saved, "cells.tsv"), sep="\t", dtype=str, keep_default_na=False, quoting=3)
    res = m.summary.get("result") or {}
    pcol, ccol = res.get("pANN_column"), res.get("classification_column")
    cells = cells.set_index("barcode")
    g.add("obs names (the saved object's cells)", list(a.obs_names) == list(cells.index), blas=False)
    sub = cells.reindex(list(a.obs_names))
    x = a.obs[pcol].to_numpy(dtype=np.float64)
    y = pd.to_numeric(sub["pANN"]).to_numpy(dtype=np.float64)
    g.add(f"obs {pcol}", same_array(x, y), num_stats(x, y))
    g.add(f"obs {ccol}", set(a.obs[ccol].tolist()) == {"Singlet"}, blas=True,
          detail="the saved object keeps the singlets, so every call is Singlet")
    g.add("X", a.X is None, blas=False, detail="empty: only SCTransform's log1p data go into X")


def files_check(a, m, rows: list, obj: str):
    """The loader loaded exactly what the files hold: each value parsed again with Python's float()."""

    def add(path, component, n, differing):
        rows.append({"object": obj, "file": path, "component": component, "values": n, "differing": differing,
                     "status": "identical" if differing == 0 else "differing"})

    def mm_triplets(path):
        with gzip.open(path, "rt") as f:
            line = f.readline()
            while line.startswith("%"):
                line = f.readline()
            df = pd.read_csv(f, sep=" ", header=None, dtype={0: np.int64, 1: np.int64, 2: str}, quoting=3)
        x = np.fromiter((float(s) for s in df[2]), dtype=np.float64, count=len(df))
        return df[0].to_numpy() - 1, df[1].to_numpy() - 1, x

    def check_mm(role, mat_cells_by_genes_cols, label, col_pos=None):
        e = m.entry(role)
        i, j, x = mm_triplets(m.path(e))
        mat = mat_cells_by_genes_cols
        # the file is rows (genes or cells) x columns (cells); the loaded matrix is cells x rows
        rr = col_pos[i] if col_pos is not None else i
        got = np.asarray(mat[j, rr]).ravel()
        extra = abs(int(mat.nnz) - len(x))  # entries the loaded matrix holds that the file does not
        add(e["path"], label, len(x), int((got != x).sum() + (np.signbit(got) != np.signbit(x)).sum()) + extra)

    if m.command == "doubletfinder":
        check_mm("raw_counts", a.layers["counts"], "singlets' raw counts")
        return
    check_mm("raw_counts", a.layers["counts"], "layers counts")
    modelled = np.flatnonzero(a.var["sct_modelled"].to_numpy())
    check_mm("sct_counts", a.layers["sct_corrected_counts"], "layers sct_corrected_counts", modelled)
    check_mm("sct_data", a.X, "X", modelled)
    check_mm("graph_nn", a.obsp["sct_nn"].T, "obsp sct_nn")
    check_mm("graph_snn", a.obsp["sct_snn"].T, "obsp sct_snn")
    e = m.entry("sct_scale_data")
    raw = np.load(m.path(e))
    got = a.obsm["sct_scale_data"].to_numpy(dtype=np.float64).T
    add(e["path"], "obsm sct_scale_data", raw.size, int((raw.view(np.uint64) != np.ascontiguousarray(got).view(np.uint64)).sum()))

    def check_tsv(role, cols: dict, label):
        e = m.entry(role)
        df = pd.read_csv(m.path(e), sep="\t", dtype=str, keep_default_na=False, quoting=3)
        n = bad = 0
        for c, v in cols.items():
            x = np.fromiter((float(s) if s != "NA" else np.nan for s in df[c]), dtype=np.float64, count=len(df))
            v = np.asarray(v, dtype=np.float64)
            n += len(x)
            bad += int((~((x == v) | (np.isnan(x) & np.isnan(v)))).sum())
        add(e["path"], label, n, bad)

    fa = a.uns["sct_model"]["feature_attributes"]
    check_tsv("sct_feature_attributes", {c: fa[c] for c in ("detection_rate", "gmean", "variance", "residual_mean", "residual_variance")},
              "sct_model feature_attributes")
    check_tsv("sct_model_pars_fit", {c: fa[c] for c in ("theta", "(Intercept)", "log_umi")}, "sct_model theta, (Intercept), log_umi")
    ca = a.uns["sct_model"]["cell_attributes"]
    check_tsv("sct_cell_attributes", {"umi": ca["umi"], "log_umi": ca["log_umi"]}, "sct_model cell_attributes")
    num = [c for c in ("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "silhouette_width", "neighborhood_purity") if c in a.obs]
    check_tsv("cell_metadata", {c: a.obs[c] for c in num}, "obs (numeric columns)")
    emb = a.obsm["X_pca"]
    e = m.entry("pca_embeddings")
    pcs = pd.read_csv(m.path(e), sep="\t", nrows=0).columns[1:]
    check_tsv("pca_embeddings", {c: emb[:, k] for k, c in enumerate(pcs)}, "obsm X_pca")
    pos = a.var_names.get_indexer(list(a.uns["pca"]["features"]))
    ld = a.varm["PCs"][pos, :]
    check_tsv("pca_loadings", {c: ld[:, k] for k, c in enumerate(pcs)}, "varm PCs")
    check_tsv("pca_stdev", {"stdev": a.uns["pca"]["stdev"]}, "uns pca stdev")


def same_anndata(a, b) -> list:
    """Where two AnnData differ (empty if equal): X, layers, obs, var, obsm, varm, obsp, uns."""
    diffs = []

    def same_frame(x: pd.DataFrame, y: pd.DataFrame) -> bool:
        if list(x.index) != list(y.index) or list(x.columns) != list(y.columns):
            return False
        for c in x.columns:
            u, v = x[c], y[c]
            if isinstance(u.dtype, pd.CategoricalDtype) or isinstance(v.dtype, pd.CategoricalDtype):
                if not (isinstance(u.dtype, pd.CategoricalDtype) and isinstance(v.dtype, pd.CategoricalDtype)):
                    return False
                if list(u.cat.categories) != list(v.cat.categories) or not np.array_equal(u.cat.codes, v.cat.codes):
                    return False
            elif u.dtype.kind in "fiub" and v.dtype.kind in "fiub":
                if u.dtype != v.dtype or not same_array(u.to_numpy(), v.to_numpy()):
                    return False
            elif [str(t) for t in u.tolist()] != [str(t) for t in v.tolist()]:
                return False
        return True

    def same_value(x, y) -> bool:
        if isinstance(x, dict) or isinstance(y, dict):
            return isinstance(x, dict) and isinstance(y, dict) and set(x) == set(y) and all(same_value(x[k], y[k]) for k in x)
        if isinstance(x, pd.DataFrame) or isinstance(y, pd.DataFrame):
            return isinstance(x, pd.DataFrame) and isinstance(y, pd.DataFrame) and same_frame(x, y)
        if sp.issparse(x) or sp.issparse(y):
            return sp.issparse(x) and sp.issparse(y) and same_sparse(x, y)
        if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
            xa, ya = np.asarray(x), np.asarray(y)
            if xa.dtype.kind in "OUT" or ya.dtype.kind in "OUT":
                return xa.shape == ya.shape and [str(t) for t in xa.ravel()] == [str(t) for t in ya.ravel()]
            return same_array(xa, ya)
        if x is None or y is None:
            return x is None and y is None
        if isinstance(x, float) and isinstance(y, float) and math.isnan(x) and math.isnan(y):
            return True
        return type(x) is type(y) and x == y

    if not same_value(a.X, b.X):
        diffs.append("X")
    for name in ("layers", "obsm", "varm", "obsp"):
        x, y = getattr(a, name), getattr(b, name)
        if set(x.keys()) != set(y.keys()):
            diffs.append(f"{name} keys")
            continue
        for k in x.keys():
            if not same_value(x[k], y[k]):
                diffs.append(f"{name}[{k}]")
    if not same_frame(a.obs, b.obs):
        diffs.append("obs")
    if not same_frame(a.var, b.var):
        diffs.append("var")
    if not same_value(dict(a.uns), dict(b.uns)):
        diffs.append("uns")
    return diffs


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--object", required=True)
    p.add_argument("--export", required=True)
    p.add_argument("--saved", required=True)
    p.add_argument("--input", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--kind", default="cluster", choices=["cluster", "doubletfinder"])
    p.add_argument("--tmp", default=None)
    p.add_argument("--files", default="yes", choices=["yes", "no"])
    p.add_argument("--round-trip", default="yes", choices=["yes", "no"])
    p.add_argument("--tolerances", default=os.path.join(ROOT, "validation", "tolerances.toml"))
    p.add_argument("--saved-kernel", default="SapphireRapids")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    steps = {}
    m = scrustio.read_manifest(a.export)
    kernel = ((m.summary.get("run") or {}).get("fingerprint") or {}).get("openblas_corename")
    mode = "exact" if kernel == a.saved_kernel else "floor"
    with open(a.tolerances, "rb") as f:
        floor = tomllib.load(f)["clustering"]["floor"]["kernel"]
    g = Grader(a.object, mode, floor)
    ad = scrustio.load_anndata(m)
    steps["load"] = {"elapsed_s": time.time() - t0, "peak_rss_gib": peak_rss_gib()}
    frows = []
    if a.files == "yes":
        files_check(ad, m, frows, a.object)
        steps["files"] = {"elapsed_s": time.time() - t0, "peak_rss_gib": peak_rss_gib()}
    if a.kind == "cluster":
        compare_cluster(g, ad, a.saved, a.input)
    else:
        compare_doubletfinder(g, ad, a.saved, m)
    steps["compare"] = {"elapsed_s": time.time() - t0, "peak_rss_gib": peak_rss_gib()}
    rt = None
    if a.round_trip == "yes":
        import anndata

        tmp = a.tmp or a.out
        os.makedirs(tmp, exist_ok=True)
        f = os.path.join(tmp, f"{a.object}.h5ad")
        written = scrustio.write_h5ad(a.export, f)
        back = anndata.read_h5ad(f)
        d1, d2 = same_anndata(ad, back), same_anndata(written, back)
        rt = {"object": a.object, "loader": "anndata", "h5ad_bytes": os.path.getsize(f),
              "identical_to_loader": not d1, "identical_to_rebuilt": not d2, "differences": "; ".join(d1 + d2)}
        os.remove(f)
        steps["round_trip"] = {"elapsed_s": time.time() - t0, "peak_rss_gib": peak_rss_gib()}
    tab = pd.DataFrame(g.rows)
    tab.to_csv(os.path.join(a.out, "k11_anndata.tsv"), sep="\t", index=False, na_rep="")
    if frows:
        pd.DataFrame(frows).to_csv(os.path.join(a.out, "files_anndata.tsv"), sep="\t", index=False)
    if rt:
        pd.DataFrame([rt]).to_csv(os.path.join(a.out, "roundtrip_anndata.tsv"), sep="\t", index=False)
    failing = tab["status"].isin(["beyond floor", "differing"])
    files_ok = all(r["differing"] == 0 for r in frows) if frows else None
    rt_ok = (rt["identical_to_loader"] and rt["identical_to_rebuilt"]) if rt else None
    counts = tab["status"].value_counts().to_dict()
    summary = {"object": a.object, "kind": a.kind, "export": os.path.abspath(a.export), "saved": os.path.abspath(a.saved),
               "export_kernel": kernel, "saved_kernel": a.saved_kernel, "mode": mode,
               "pass": bool(not failing.any() and files_ok is not False and rt_ok is not False),
               "components": int(len(tab)), "status": counts,
               "failing": [r for r in tab.loc[failing, "component"]], "files_identical": files_ok,
               "round_trip_identical": rt_ok, "steps": steps, "elapsed_s": time.time() - t0,
               "peak_rss_gib": peak_rss_gib(), "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
               "versions": {p: __import__("importlib.metadata").metadata.version(p) for p in ("anndata", "numpy", "scipy", "pandas")}}
    with open(os.path.join(a.out, "compare_anndata.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"{a.object} ({a.kind}, {mode} mode): {len(tab)} components: {counts}; files "
          f"{'not checked' if files_ok is None else ('identical' if files_ok else 'DIFFERING')}; round trip "
          f"{'not run' if rt_ok is None else ('identical' if rt_ok else 'DIFFERING: ' + rt['differences'])}; "
          f"{time.time() - t0:.0f} s, peak RSS {peak_rss_gib():.2f} GiB")
    if failing.any():
        print(tab.loc[failing].to_string())
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
