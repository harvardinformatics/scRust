"""The output manifest: summary.json's ``files`` section, schema ``scrust-files`` version 1.

The schema is docs/output_manifest.md in the scRust repository. Every file is found through the
manifest, by its role, never by its name.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Any

SCHEMA = "scrust-files"
SCHEMA_VERSION = 1

# Roles and scales of schema version 1 (docs/output_manifest.md, "Roles" and "Scales").
ROLES = (
    "summary", "readme", "raw_counts", "raw_counts_cells", "raw_counts_features", "rna_cell_metadata",
    "cell_metadata", "cluster_labels", "cluster_levels", "cluster_columns", "workflow_cluster_ids",
    "workflow_nclusters", "optimizer_ids", "sct_cell_attributes", "sct_feature_attributes",
    "sct_model_pars", "sct_model_pars_fit", "sct_variable_features", "sct_scale_data", "sct_counts",
    "sct_data", "sct_features", "sct_cells", "pca_embeddings", "pca_loadings", "pca_stdev", "rng_state",
    "knn_index", "knn_distance", "graph_nn", "graph_snn", "annoy_index", "doublet_calls",
    "doublet_sweep_stats", "doublet_bcmvn", "stage_dump", "stage_index", "stage_scalars",
    "start_partitions", "nested_output",
)
SCALES = (
    "raw_counts", "corrected_counts", "log1p_corrected_counts", "residuals", "pca_embeddings",
    "pca_loadings", "pca_stdev", "graph_nn", "graph_snn", "neighbour_indices", "neighbour_distances",
    "none",
)


class ScrustError(Exception):
    """A scrust output that cannot be read as asked: a refused manifest, a missing or changed
    file, a shape that does not match, or a file offered to a slot not meant for it."""


@dataclass
class Manifest:
    """A scrust output's manifest: its directory, the command that wrote it, the provenance,
    the entries (one per file) and the whole summary.json."""

    dir: str
    command: str
    provenance: dict
    entries: list
    summary: dict
    _cache: dict = field(default_factory=dict, repr=False)

    def entries_by_role(self, role: str) -> list:
        return [e for e in self.entries if e.get("role") == role]

    def has_role(self, role: str) -> bool:
        return bool(self.entries_by_role(role))

    def entry(self, role: str, required: bool = True) -> dict | None:
        """The one entry with this role; an error if there is none (when required) or several."""
        hits = self.entries_by_role(role)
        if len(hits) > 1:
            paths = ", ".join(e["path"] for e in hits)
            raise ScrustError(f"{len(hits)} files have role {role} ({paths}); expected one")
        if not hits:
            if required:
                raise ScrustError(f"this {self.command} output has no file with role {role}")
            return None
        return hits[0]

    def path(self, entry: dict) -> str:
        return os.path.join(self.dir, entry["path"])

    def entry_at(self, path: str) -> dict:
        for e in self.entries:
            if e["path"] == path:
                return e
        raise ScrustError(f"the manifest names {path} as an identifier file but does not list it")

    def identifiers(self, ids: dict | None) -> list | None:
        """The identifiers a ``rows`` or ``columns`` field points to (cached)."""
        if ids is None:
            return None
        key = (ids["file"], ids["column"])
        if key not in self._cache:
            from .read import read_id_column

            e = self.entry_at(ids["file"])
            self._cache[key] = read_id_column(self.path(e), bool(e.get("header")), ids["column"])
        return self._cache[key]

    def check_shape(self, entry: dict) -> None:
        """The number of row and column identifiers against the entry's shape."""
        shape = entry.get("shape")
        if shape is None:
            return
        for axis, field_name in ((0, "rows"), (1, "columns")):
            ids = entry.get(field_name)
            if ids is None:
                continue
            if len(shape) <= axis:
                raise ScrustError(f"{entry['path']}: {field_name} identifiers but a {len(shape)}-dimensional shape")
            n = len(self.identifiers(ids))
            if n != shape[axis]:
                raise ScrustError(
                    f"{entry['path']}: the manifest's shape gives {shape[axis]} {field_name}, but its "
                    f"identifiers ({ids['file']}, column {ids['column']}) number {n}"
                )

    def check(self, sha256: bool = False) -> None:
        """Every entry's identifiers against its shape, and with ``sha256`` every checksum."""
        for e in self.entries:
            self.check_shape(e)
        if sha256:
            verify_sha256(self)

    def __repr__(self) -> str:
        return f"Manifest({self.command!r}, {self.dir!r}, {len(self.entries)} files)"


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def verify_sha256(m: Manifest, entries: list | None = None) -> None:
    """Every listed file's SHA-256 against the manifest's."""
    for e in entries if entries is not None else m.entries:
        if e.get("sha256") is None:
            continue
        got = sha256_file(m.path(e))
        if got != e["sha256"]:
            raise ScrustError(f"{e['path']}: SHA-256 {got}, the manifest lists {e['sha256']}")


def read_manifest(path: str, check_sha256: bool = False) -> Manifest:
    """Read and check the manifest of a scrust output directory.

    Refuses a manifest whose ``schema`` is not "scrust-files" or whose ``version`` is not 1,
    and checks that every listed file exists with its listed size (and, with ``check_sha256``,
    its SHA-256).
    """
    if not os.path.isdir(path):
        raise ScrustError(f"no directory {path}")
    path = os.path.abspath(path)
    summary_path = os.path.join(path, "summary.json")
    if not os.path.isfile(summary_path):
        raise ScrustError(f"{path} has no summary.json: not a scrust output")
    with open(summary_path) as f:
        js: dict[str, Any] = json.load(f)
    files = js.get("files")
    if files is None:
        raise ScrustError(f"{summary_path} has no `files` section (the manifest)")
    if files.get("schema") != SCHEMA:
        raise ScrustError(f"{summary_path}: manifest schema {files.get('schema')!r}, not {SCHEMA!r}; refusing it")
    version = files.get("version")
    if isinstance(version, bool) or not isinstance(version, (int, float)) or version != SCHEMA_VERSION:
        raise ScrustError(
            f"{summary_path}: manifest schema version {version!r}; this loader reads version "
            f"{SCHEMA_VERSION} only. Refusing it"
        )
    entries = files.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ScrustError(f"{summary_path}: the manifest lists no files")
    seen = set()
    for k, e in enumerate(entries):
        for name in ("path", "role", "scale", "format"):
            if not isinstance(e.get(name), str) or not e[name]:
                raise ScrustError(f"{summary_path}: entry {k + 1} has no valid `{name}`")
        p = e["path"]
        parts = p.split("/")
        if p.startswith("/") or ".." in parts or "\\" in p:
            raise ScrustError(f"{summary_path}: entry path {p} is not a relative path inside the output")
        if p in seen:
            raise ScrustError(f"{summary_path}: {p} is listed twice")
        seen.add(p)
        full = os.path.join(path, p)
        if not os.path.isfile(full):
            raise ScrustError(f"{p} is listed in the manifest but missing from {path}")
        if e.get("bytes") is not None and os.path.getsize(full) != e["bytes"]:
            raise ScrustError(f"{p} has {os.path.getsize(full)} bytes; the manifest lists {e['bytes']}")
    m = Manifest(dir=path, command=files.get("command", ""), provenance=files.get("provenance") or {},
                 entries=entries, summary=js)
    if check_sha256:
        verify_sha256(m)
    return m
