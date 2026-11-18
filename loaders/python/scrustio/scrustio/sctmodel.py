"""SCTransform's model, from the model files, as Seurat 5.5.1's SCTModel holds it.

The same tables as Seurat's PrepVSTResults() makes from sctransform 0.4.3's vst() (feature and cell
attributes), with the clip ranges, the model formula, vst()'s arguments as Seurat's SCTransform()
records them and the median UMI. AnnData has no model class: they go into
``uns["sct_model"]``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .manifest import Manifest, ScrustError
from .read import read_tsv

PARS = ["theta", "(Intercept)", "log_umi"]
FEATURE_COLS = ["detection_rate", "gmean", "variance", "residual_mean", "residual_variance"]


def sct_parameters(m: Manifest) -> dict:
    p = m.provenance.get("parameters") or {}
    return p.get("sctransform", p)


def sct_result(m: Manifest) -> dict:
    r = m.summary.get("result") or {}
    return r.get("sctransform", r)


def sct_tables(m: Manifest) -> dict:
    ga = read_tsv(m.path(m.entry("sct_feature_attributes")),
                  numeric=["detection_rate", "gmean", "amean", "variance", "residual_mean", "residual_variance"])
    fit = read_tsv(m.path(m.entry("sct_model_pars_fit")), numeric=PARS)
    step1 = read_tsv(m.path(m.entry("sct_model_pars")), numeric=PARS)
    ca = read_tsv(m.path(m.entry("sct_cell_attributes")), numeric=["umi", "log_umi"])
    if fit["name"].tolist() != ga["name"].tolist():
        raise ScrustError("model_pars_fit's genes differ from gene_attr's")
    return {"gene_attr": ga, "fit": fit, "step1": step1, "cell_attr": ca}


def cells_step1(m: Manifest, n_cells: int):
    """Which cells vst() fitted on in its first step, when the output determines it: every cell
    when there are at most n_cells. Otherwise vst() drew them with R's sample() after
    SCTransform's set.seed(), which this loader does not reproduce: None."""
    n1 = sct_result(m).get("n_cells_step1")
    if n1 is None:
        raise ScrustError("the summary does not give vst's number of step-1 cells")
    return np.ones(n_cells, dtype=bool) if n1 >= n_cells else None


def model_parts(m: Manifest, tabs: dict, gene_names, step1_names, cells) -> dict:
    """The SCTModel's slots: feature_attributes, cell_attributes, clips, umi_assay, model,
    arguments, median_umi."""
    ga = tabs["gene_attr"]
    fa = pd.DataFrame({c: ga[c].to_numpy(dtype=np.float64) for c in FEATURE_COLS}, index=pd.Index(gene_names))
    for c in PARS:
        fa[c] = tabs["fit"][c].to_numpy(dtype=np.float64)
    s1 = tabs["step1"]
    in_step1 = pd.Index(gene_names).isin(step1_names)
    fa["genes_log_gmean_step1"] = in_step1
    pos = pd.Index(step1_names)
    for c in PARS:
        col = np.full(len(fa), np.nan)
        k = pos.get_indexer(pd.Index(gene_names)[in_step1])
        col[in_step1] = s1[c].to_numpy(dtype=np.float64)[k]
        fa["step1_" + c] = col
    ca = pd.DataFrame({"umi": tabs["cell_attr"]["umi"].to_numpy(dtype=np.float64),
                       "log_umi": tabs["cell_attr"]["log_umi"].to_numpy(dtype=np.float64)}, index=pd.Index(cells))
    step1 = cells_step1(m, len(cells))
    if step1 is not None:
        ca["cells_step1"] = step1
    res = sct_result(m)
    par = sct_parameters(m)

    def num(x, what):
        if x is None:
            raise ScrustError(f"the summary does not give {what}")
        return x

    res_clip = np.asarray(num(res.get("res_clip_range"), "vst's res_clip_range"), dtype=np.float64)
    sct_clip = np.asarray(num(res.get("clip_range"), "SCTransform's clip.range"), dtype=np.float64)
    arguments = {
        "glmGamPoi_check": True,
        "latent_var": "log_umi",
        "batch_var": None,
        "latent_var_nonreg": None,
        "n_genes": float(num(par.get("n_genes"), "vst's n_genes")),
        "n_cells": float(min(float(num(par.get("ncells"), "SCTransform's ncells")), len(cells))),
        "method": "glmGamPoi_offset",
        "do_regularize": True,
        "theta_regularization": "od_factor",
        "res_clip_range": res_clip,
        "bin_size": float(num(par.get("bin_size"), "vst's bin_size")),
        "min_cells": float(num(par.get("min_cells"), "vst's min_cells")),
        "residual_type": "pearson",
        "return_cell_attr": True,
        "return_gene_attr": True,
        "return_corrected_umi": True,
        "min_variance": "umi_median",
        "bw_adjust": float(num(par.get("bw_adjust"), "vst's bw_adjust")),
        "gmean_eps": 1.0,
        "theta_estimation_fun": "theta.ml",
        "theta_given": None,
        "exclude_poisson": True,
        "use_geometric_mean": True,
        "use_geometric_mean_offset": False,
        "fix_intercept": False,
        "fix_slope": False,
        "scale_factor": None,
        "vst.flavor": "v2",
        "verbosity": 0.0,
        "set_min_var": float(num(res.get("min_variance"), "vst's set_min_var")),
        "sct.clip.range": sct_clip,
    }
    umi = ca["umi"].to_numpy()
    median_umi = float(np.median(umi))
    return {"feature_attributes": fa, "cell_attributes": ca, "clips": {"vst": res_clip, "sct": sct_clip},
            "umi_assay": "RNA", "model": "y ~ log_umi", "arguments": arguments, "median_umi": median_umi}
