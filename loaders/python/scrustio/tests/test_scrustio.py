"""scrustio's tests (docs/PLAN.md C3g; milestone M14): the slots, the manifest, the AnnData loader,
other commands' outputs and the one-line command's .h5ad."""

from __future__ import annotations

import gzip
import os

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import scrustio
from scrustio.cli import main as cli_main
from scrustio.manifest import ROLES, SCALES, ScrustError
from scrustio.slots import SLOTS, check_comments, check_slot

from conftest import FIXTURES, tamper

FORMATS = ("matrix_market", "tsv", "npy", "json", "text", "annoy")


# 1. The slots ------------------------------------------------------------------------------------

def test_every_wrong_combination_is_refused():
    wrong = refused = 0
    for slot, s in SLOTS.items():
        for role in ROLES:
            for scale in SCALES:
                want = role == s["role"] and scale == s["scale"]
                try:
                    check_slot(slot, role, scale, "f", s["format"])
                    ok = True
                except ScrustError:
                    ok = False
                assert ok == want, (slot, role, scale)
                if not want:
                    wrong += 1
                    refused += not ok
        for fmt in FORMATS:
            want = fmt == s["format"]
            try:
                check_slot(slot, s["role"], s["scale"], "f", fmt)
                ok = True
            except ScrustError:
                ok = False
            assert ok == want, (slot, fmt)
            if not want:
                wrong += 1
                refused += not ok
    assert wrong == refused == len(SLOTS) * (len(ROLES) * len(SCALES) - 1 + len(FORMATS) - 1)


def test_raw_count_slots_take_only_raw_counts():
    t = scrustio.slot_table()
    assert (t.loc[t.raw_counts, "scale"] == "raw_counts").all()
    assert (t.loc[~t.raw_counts, "scale"] != "raw_counts").all()


def test_the_comment_lines():
    with pytest.raises(ScrustError, match="NOT raw counts"):
        check_comments('AnnData layers["counts"]', ["% scrust: x. NOT raw counts: y"], "f")
    with pytest.raises(ScrustError, match="not so marked"):
        check_comments("AnnData X", ["%metadata_json: {}"], "f")


def test_tampered_manifests_are_refused(outputs, copy_output):
    d = copy_output(outputs["export"], "swap")
    tamper(d, "sct/corrected_counts.mtx.gz", "role", '"raw_counts"')
    tamper(d, "sct/corrected_counts.mtx.gz", "scale", '"raw_counts"')
    tamper(d, "counts/matrix.mtx.gz", "role", '"sct_counts"')
    tamper(d, "counts/matrix.mtx.gz", "scale", '"corrected_counts"')
    with pytest.raises(ScrustError, match="NOT raw counts"):
        scrustio.load_anndata(d)
    d = copy_output(outputs["export"], "scale_only")
    tamper(d, "sct/corrected_counts.mtx.gz", "scale", '"raw_counts"')
    with pytest.raises(ScrustError, match=r"role sct_counts, scale raw_counts"):
        scrustio.load_anndata(d)
    d = copy_output(outputs["export"], "residuals")
    tamper(d, "counts/matrix.mtx.gz", "role", '"stage_dump"')
    tamper(d, "sct/scale_data.npy", "role", '"raw_counts"')
    tamper(d, "sct/scale_data.npy", "scale", '"raw_counts"')
    with pytest.raises(ScrustError, match="format npy"):
        scrustio.load_anndata(d)
    d = copy_output(outputs["export"], "log1p_swap")
    tamper(d, "sct/corrected_counts.mtx.gz", "role", '"sct_data"')
    tamper(d, "sct/corrected_counts.mtx.gz", "scale", '"log1p_corrected_counts"')
    tamper(d, "sct/log1p_corrected_counts.mtx.gz", "role", '"sct_counts"')
    tamper(d, "sct/log1p_corrected_counts.mtx.gz", "scale", '"corrected_counts"')
    with pytest.raises(ScrustError, match="not log1p of the counts layer"):
        scrustio.load_anndata(d)


# 2. The manifest -----------------------------------------------------------------------------------

def _edit(path, old, new):
    s = open(path).read()
    assert old in s
    open(path, "w").write(s.replace(old, new, 1))


def test_manifest_rules(outputs, copy_output):
    d = copy_output(outputs["export"], "schema")
    _edit(os.path.join(d, "summary.json"), '"schema": "scrust-files"', '"schema": "other-files"')
    with pytest.raises(ScrustError, match="schema"):
        scrustio.read_manifest(d)
    d = copy_output(outputs["export"], "version")
    _edit(os.path.join(d, "summary.json"), '    "version": 1,', '    "version": 2,')
    with pytest.raises(ScrustError, match="version 2"):
        scrustio.read_manifest(d)
    d = copy_output(outputs["export"], "missing")
    os.remove(os.path.join(d, "sct", "corrected_counts.mtx.gz"))
    with pytest.raises(ScrustError, match="missing"):
        scrustio.read_manifest(d)
    d = copy_output(outputs["export"], "size")
    open(os.path.join(d, "pca", "stdev.tsv"), "a").write("x")
    with pytest.raises(ScrustError, match="bytes"):
        scrustio.read_manifest(d)
    d = copy_output(outputs["export"], "sha")
    p = os.path.join(d, "pca", "stdev.tsv")
    lines = open(p).read().split("\n")
    lines[1] = lines[1][:-1] + ("2" if lines[1][-1] != "2" else "3")
    open(p, "w").write("\n".join(lines))
    scrustio.read_manifest(d)
    with pytest.raises(ScrustError, match="SHA-256"):
        scrustio.read_manifest(d, check_sha256=True)
    d = copy_output(outputs["export"], "shape")
    tamper(d, "sct/corrected_counts.mtx.gz", "shape", "[820, 300]")
    with pytest.raises(ScrustError, match="shape gives 820 rows"):
        scrustio.load_anndata(d)
    d = copy_output(outputs["export"], "traversal")
    tamper(d, "pca/stdev.tsv", "path", '"../stdev.tsv"')
    with pytest.raises(ScrustError, match="not a relative path"):
        scrustio.read_manifest(d)


def test_files_are_found_by_role(outputs, copy_output):
    d = copy_output(outputs["export"], "renamed")
    os.rename(os.path.join(d, "sct", "corrected_counts.mtx.gz"), os.path.join(d, "sct", "layer_a.mtx.gz"))
    os.rename(os.path.join(d, "sct", "sct_features.tsv.gz"), os.path.join(d, "sct", "rows_a.tsv.gz"))
    p = os.path.join(d, "summary.json")
    s = open(p).read().replace("sct/corrected_counts.mtx.gz", "sct/layer_a.mtx.gz").replace("sct/sct_features.tsv.gz", "sct/rows_a.tsv.gz")
    open(p, "w").write(s)
    a, b = scrustio.load_anndata(outputs["export"]), scrustio.load_anndata(d)
    for key in ("counts", "sct_corrected_counts"):
        assert (a.layers[key] != b.layers[key]).nnz == 0
    assert (a.X != b.X).nnz == 0


# 3. AnnData ------------------------------------------------------------------------------------------

def test_anndata(outputs):
    a = scrustio.load_anndata(outputs["export"])
    m = scrustio.read_manifest(outputs["export"])
    assert a.shape == tuple(reversed(m.entry("raw_counts")["shape"]))
    levels = open(os.path.join(outputs["export"], "seurat_clusters_levels.txt")).read().split()
    assert list(a.obs["seurat_clusters"].cat.categories) == levels
    saved = pd.read_csv(os.path.join(FIXTURES, "clustering", "filtered_seurat_tenx_test", "saved", "metadata.tsv.gz"), sep="\t", dtype=str)
    assert a.obs["seurat_clusters"].astype(str).tolist() == saved["seurat_clusters"].tolist()
    assert a.obs["nFeature_RNA"].dtype == np.int32 and a.obs["nCount_RNA"].dtype == np.float64
    modelled = a.var["sct_modelled"].to_numpy()
    assert modelled.sum() == m.entry("sct_counts")["shape"][0]
    # genes outside the model are 0 in X and in the corrected counts
    assert a.X[:, ~modelled].nnz == 0 and a.layers["sct_corrected_counts"][:, ~modelled].nnz == 0
    assert np.array_equal(a.X.data, scrustio.anndata_loader.libm_log1p(a.layers["sct_corrected_counts"].data))
    assert a.uns["log1p"] == {"base": None}
    vf = pd.read_csv(os.path.join(outputs["export"], "sct", "variable_features.tsv"), sep="\t")["name"].tolist()
    assert list(a.obsm["sct_scale_data"].columns) == vf
    assert a.obsm["X_pca"].shape == (a.n_obs, 50)
    pos = a.var_names.get_indexer(vf)
    assert (a.varm["PCs"][np.setdiff1d(np.arange(a.n_vars), pos)] == 0).all()
    assert set(a.obsp) == {"sct_nn", "sct_snn"} and a.uns["sct_snn"]["connectivities_key"] == "sct_snn"
    assert set(a.uns["sct_model"]["feature_attributes"].columns) >= {"residual_variance", "theta", "step1_theta"}
    assert a.uns["sct_model"]["cell_attributes"]["cells_step1"].all()
    # the raw counts equal the input's
    raw, _, _ = scrustio.read.read_mm(os.path.join(outputs["export"], "counts", "matrix.mtx.gz"))
    assert (a.layers["counts"] != raw.T.tocsr()).nnz == 0


def test_numpy_log1p_is_not_r_s():
    """np.log1p may differ from C's log1p (R's) in the last bit: the loader reads Seurat's data
    layer and checks it with C's log1p (math.log1p)."""
    x = np.arange(1.0, 20000.0)
    import math

    assert np.array_equal(scrustio.anndata_loader.libm_log1p(x), np.array([math.log1p(v) for v in x]))


def test_h5ad_one_line_command(outputs, tmp_path):
    out = str(tmp_path / "object.h5ad")
    assert cli_main(["h5ad", outputs["export"], out]) == 0
    a = scrustio.load_anndata(outputs["export"])
    b = anndata.read_h5ad(out)
    assert list(a.obs.columns) == list(b.obs.columns)
    for c in a.obs.columns:
        if isinstance(a.obs[c].dtype, pd.CategoricalDtype):
            assert list(a.obs[c].cat.categories) == list(b.obs[c].cat.categories)
            assert np.array_equal(a.obs[c].cat.codes, b.obs[c].cat.codes)
        else:
            assert np.array_equal(a.obs[c].to_numpy(), b.obs[c].to_numpy(), equal_nan=True)
    for k in a.layers:
        assert (a.layers[k] != b.layers[k]).nnz == 0
    assert (a.X != b.X).nnz == 0
    assert np.array_equal(a.obsm["X_pca"], b.obsm["X_pca"])
    assert a.obsm["sct_scale_data"].equals(b.obsm["sct_scale_data"])
    assert np.array_equal(a.varm["PCs"], b.varm["PCs"])
    for k in a.obsp:
        assert (a.obsp[k] != b.obsp[k]).nnz == 0
    assert a.uns["sct_model"]["feature_attributes"].equals(b.uns["sct_model"]["feature_attributes"])


# 4. Other commands' outputs ---------------------------------------------------------------------------

def test_sctransform_output(outputs):
    s = scrustio.load_anndata(outputs["sctransform"])
    e = scrustio.load_anndata(outputs["export"])
    assert "counts" not in s.layers
    modelled = e.var["sct_modelled"].to_numpy()
    assert list(s.var_names) == list(e.var_names[modelled])
    assert (s.X != e.X[:, modelled]).nnz == 0
    assert {"nCount_RNA", "nFeature_RNA"} <= set(s.obs.columns)


def test_doubletfinder_output(outputs):
    a = scrustio.load_anndata(outputs["doubletfinder"])
    m = scrustio.read_manifest(outputs["doubletfinder"])
    res = m.summary["result"]
    singlets = gzip.open(os.path.join(outputs["doubletfinder"], "singlets", "barcodes.tsv.gz"), "rt").read().split()
    assert list(a.obs_names) == singlets
    assert a.X is None and "counts" in a.layers
    assert {res["pANN_column"], res["classification_column"]} <= set(a.obs.columns)
    assert set(a.obs[res["classification_column"]].astype(str)) == {"Singlet"}


def test_nested_outputs(outputs):
    """M15: scrust doubletfinder's clusters/ and recluster/ are nested outputs. The DoubletFinder
    output loads as before, by its own files' roles; each nested output loads on its own."""
    d = outputs["doubletfinder_nested"]
    m = scrustio.read_manifest(d)
    assert len(m.entries_by_role("nested_output")) == 2
    assert len(m.entries_by_role("raw_counts")) == 1
    a = scrustio.load_anndata(d)
    b = scrustio.load_anndata(outputs["doubletfinder"])
    assert list(a.obs_names) == list(b.obs_names)
    res = m.summary["result"]
    assert (a.obs[res["pANN_column"]].to_numpy() == b.obs[res["pANN_column"]].to_numpy()).all()
    r = scrustio.load_anndata(os.path.join(d, "recluster"))
    assert list(r.obs_names) == list(a.obs_names)
    assert "X_pca" in r.obsm and r.X is not None and "seurat_clusters" in r.obs.columns
    c = scrustio.read_manifest(os.path.join(d, "clusters"))
    assert c.command == "cluster" and c.has_role("cell_metadata")


def test_cluster_output_without_export_holds_no_matrices(outputs):
    m = scrustio.read_manifest(outputs["cluster"])
    assert not m.has_role("raw_counts") and not m.has_role("sct_counts")
    with pytest.raises(ScrustError, match="neither raw counts nor"):
        scrustio.load_anndata(outputs["cluster"])
