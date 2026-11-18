"""Readers for the formats of docs/output_manifest.md.

Exactness: scrust writes doubles in text with 17 significant digits (C's %.17g), so a correctly
rounded reader recovers every double. Tables are read with pandas' ``float_precision="round_trip"``
(Python's own correctly rounded conversion), Matrix Market files with scipy.io.mmread (scipy's
fast_matrix_market, whose fast_float parser is correctly rounded), .npy files as bytes.
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .manifest import ScrustError


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def read_lines(path: str) -> list[str]:
    """The lines of a text file, gzipped or not, without their line ends."""
    with _open_text(path) as f:
        return f.read().splitlines()


def read_first_line(path: str) -> str:
    with _open_text(path) as f:
        return f.readline().rstrip("\n")


def read_id_column(path: str, header: bool, column) -> list[str]:
    """One column of a TSV file, as text: by header name (header true) or 1-based number."""
    if not header and column in (1, 1.0):
        return [line.split("\t", 1)[0] for line in read_lines(path)]
    df = pd.read_csv(path, sep="\t", header=0 if header else None, dtype=str, keep_default_na=False,
                     na_filter=False, quoting=3)
    if isinstance(column, str):
        if column not in df.columns:
            raise ScrustError(f"{path} has no column {column}")
        return df[column].tolist()
    k = int(column) - 1
    if k >= df.shape[1]:
        raise ScrustError(f"{path} has no column {column}")
    return df.iloc[:, k].tolist()


def read_tsv(path: str, numeric=(), integer=()) -> pd.DataFrame:
    """A TSV file with a header line: ``numeric`` columns as float64, ``integer`` columns as int32
    (R's integers), the others as text. "NA" is missing in numeric columns only."""
    header = read_first_line(path).split("\t")
    dtype = {c: str for c in header}
    na_values = {}
    for c in header:
        if c in numeric:
            dtype[c] = np.float64
            na_values[c] = ["NA"]
        elif c in integer:
            dtype[c] = np.int32
    df = pd.read_csv(path, sep="\t", header=0, dtype=dtype, keep_default_na=False, na_values=na_values,
                     float_precision="round_trip", quoting=3)
    if list(df.columns) != header:
        raise ScrustError(f"{path}: unexpected header")
    return df


def read_mm(path: str):
    """Matrix Market (coordinate, general): (CSC matrix, field, comment lines). Entries are kept as
    stored, explicit zeros included; a position stored twice is refused."""
    with _open_text(path) as f:
        banner = f.readline().strip()
        comments = []
        while True:
            line = f.readline()
            if not line:
                raise ScrustError(f"{path}: no size line")
            if not line.startswith("%"):
                break
            comments.append(line.rstrip("\n"))
    w = banner.split()
    if len(w) < 5 or w[0] != "%%MatrixMarket" or w[1].lower() != "matrix" or w[2].lower() != "coordinate" or w[4].lower() != "general":
        raise ScrustError(f"{path}: unsupported Matrix Market header: {banner}")
    field = w[3].lower()
    if field not in ("integer", "real", "pattern"):
        raise ScrustError(f"{path}: unsupported field {field}")
    dims = [int(x) for x in line.split()]
    if len(dims) != 3:
        raise ScrustError(f"{path}: bad size line: {line.strip()}")
    coo = scipy.io.mmread(path)
    if not sp.issparse(coo):
        coo = sp.coo_matrix(coo)
    coo = sp.coo_matrix(coo)
    if coo.shape != (dims[0], dims[1]) or coo.nnz != dims[2]:
        raise ScrustError(f"{path}: read {coo.shape} with {coo.nnz} entries; the size line says {dims}")
    data = coo.data.astype(np.float64, copy=False)
    m = sp.csc_matrix((data, (coo.row, coo.col)), shape=coo.shape)
    if m.nnz != dims[2]:
        raise ScrustError(f"{path}: a position is stored twice")
    m.has_sorted_indices = False
    m.sort_indices()
    return m, field, comments


def read_npy(path: str) -> np.ndarray:
    """A .npy file (little-endian <f8, <f4 or <i4; any order)."""
    a = np.load(path, allow_pickle=False)
    if a.dtype.byteorder == ">":
        raise ScrustError(f"{path}: big-endian arrays are not written by scrust")
    return a
