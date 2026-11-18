"""The rule (docs/PLAN.md C3g): normalised data goes only into slots meant for normalised data.

The raw counts always come from the raw counts (counts/ of an export, singlets/ of scrust
doubletfinder), and the corrected counts never go into a raw-count slot (AnnData's
``layers["counts"]``). Every slot the loader fills from a matrix file takes exactly one role and one
scale. A file is put into a slot only through :func:`read_for_slot`, which refuses any other role
or scale, and also refuses, for a raw-count slot, a Matrix Market file whose comment lines say
"NOT raw counts" (scrust marks every normalised matrix so), and, for a slot of SCTransform's
layers, one whose comment lines do not say it.
"""

from __future__ import annotations

import pandas as pd

from .manifest import Manifest, ScrustError
from .read import read_mm, read_npy, read_tsv

SLOTS = {
    'AnnData layers["counts"]': {"role": "raw_counts", "format": "matrix_market", "scale": "raw_counts", "raw": True},
    "AnnData X": {"role": "sct_data", "format": "matrix_market", "scale": "log1p_corrected_counts", "normalised": True},
    'AnnData layers["sct_corrected_counts"]': {"role": "sct_counts", "format": "matrix_market", "scale": "corrected_counts", "normalised": True},
    'AnnData obsm["sct_scale_data"]': {"role": "sct_scale_data", "format": "npy", "scale": "residuals"},
    'AnnData obsm["X_pca"]': {"role": "pca_embeddings", "format": "tsv", "scale": "pca_embeddings"},
    'AnnData varm["PCs"]': {"role": "pca_loadings", "format": "tsv", "scale": "pca_loadings"},
    'AnnData uns["pca"]["stdev"]': {"role": "pca_stdev", "format": "tsv", "scale": "pca_stdev"},
    'AnnData obsp["sct_nn"]': {"role": "graph_nn", "format": "matrix_market", "scale": "graph_nn"},
    'AnnData obsp["sct_snn"]': {"role": "graph_snn", "format": "matrix_market", "scale": "graph_snn"},
}


def slot_table() -> pd.DataFrame:
    """The slots the loader fills from matrix files, with the one role, scale and format each takes."""
    return pd.DataFrame(
        [(k, v["role"], v["scale"], v["format"], bool(v.get("raw"))) for k, v in SLOTS.items()],
        columns=["slot", "role", "scale", "format", "raw_counts"],
    )


def check_slot(slot: str, role: str, scale: str, path: str = "a file", fmt: str | None = None) -> None:
    """Refuse a file of this role, scale (and format) for this slot unless the slot takes exactly them."""
    s = SLOTS.get(slot)
    if s is None:
        raise ScrustError(f"no slot named {slot}")
    if fmt is not None and fmt != s["format"]:
        raise ScrustError(f"refusing to put {path} (format {fmt}) into {slot}, which takes only {s['format']} files")
    if role != s["role"] or scale != s["scale"]:
        why = ": raw counts come only from the raw counts (counts/ or singlets/)" if s.get("raw") else ""
        raise ScrustError(
            f"refusing to put {path} (role {role}, scale {scale}) into {slot}, which takes only role "
            f"{s['role']} with scale {s['scale']}{why}"
        )


def check_comments(slot: str, comments: list[str], path: str) -> None:
    s = SLOTS[slot]
    marked = any("NOT raw counts" in c for c in comments)
    if s.get("raw") and marked:
        raise ScrustError(f"refusing to put {path} into {slot}: its comment lines say it is NOT raw counts")
    if s.get("normalised") and not marked:
        raise ScrustError(
            f'refusing to put {path} into {slot}: scrust marks SCTransform\'s layers "NOT raw counts" '
            "in their comment lines, and this file is not so marked"
        )


def read_for_slot(m: Manifest, entry, slot: str):
    """Read a listed file (an entry, or a role held by one file) for a slot, refusing it if the
    slot is not meant for it. Matrix Market files come back as CSC matrices (rows x columns as
    stored), .npy files as arrays, tables as data frames."""
    e = m.entry(entry) if isinstance(entry, str) else entry
    check_slot(slot, e["role"], e["scale"], e["path"], e["format"])
    m.check_shape(e)
    fmt = e["format"]
    if fmt == "matrix_market":
        mat, _field, comments = read_mm(m.path(e))
        check_comments(slot, comments, e["path"])
        shape = e.get("shape")
        if shape is not None and list(mat.shape) != list(shape[:2]):
            raise ScrustError(f"{e['path']} is {mat.shape}; the manifest's shape is {shape}")
        if e.get("nnz") is not None and mat.nnz != e["nnz"]:
            raise ScrustError(f"{e['path']} stores {mat.nnz} entries; the manifest lists {e['nnz']}")
        return mat
    if fmt == "npy":
        a = read_npy(m.path(e))
        shape = e.get("shape")
        if shape is not None and list(a.shape) != list(shape):
            raise ScrustError(f"{e['path']}: the array's shape {a.shape} differs from the manifest's {shape}")
        return a
    if fmt == "tsv":
        return read_tsv(m.path(e))
    raise ScrustError(f"{e['path']}: format {fmt} cannot go into a slot")
