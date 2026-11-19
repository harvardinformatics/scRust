"""Test outputs, made once per session by the scrust binary from the committed fixtures.

The binary is SCRUST_BIN, else <repository>/target/release/scrust (cargo build --release -p scrust).
"""

from __future__ import annotations

import os
import shutil
import subprocess

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
FIXTURES = os.path.join(ROOT, "tests", "fixtures")
COUNTS_IN = os.path.join(FIXTURES, "clustering", "filtered_seurat_tenx_test", "input")
DF_IN = os.path.join(FIXTURES, "doubletfinder", "soupx_tenx_test", "input")


def _bin() -> str:
    b = os.environ.get("SCRUST_BIN") or os.path.join(ROOT, "target", "release", "scrust")
    if not os.path.isfile(b):
        pytest.fail(f"no scrust binary at {b} (cargo build --release -p scrust, or set SCRUST_BIN)")
    return b


def _run(*args):
    r = subprocess.run([_bin(), *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"scrust {' '.join(args)} failed:\n{r.stdout}\n{r.stderr}")


@pytest.fixture(scope="session")
def outputs(tmp_path_factory):
    d = tmp_path_factory.mktemp("scrust")
    out = {k: str(d / k) for k in ("export", "sctransform", "cluster", "doubletfinder", "doubletfinder_nested")}
    _run("cluster", "--counts", COUNTS_IN, "--preset", "workflow", "--export", "--out", out["export"], "--threads", "2")
    _run("sctransform", "--counts", COUNTS_IN, "--seed", "12345", "--corrected-counts", "--out", out["sctransform"], "--threads", "2")
    _run("cluster", "--counts", COUNTS_IN, "--preset", "workflow", "--out", out["cluster"], "--threads", "2")
    _run("doubletfinder", "--counts", DF_IN, "--clusters", os.path.join(DF_IN, "seurat_clusters.tsv"),
         "--out", out["doubletfinder"], "--threads", "2")
    # M15: clustering on the fly (clusters/) and the singlets' re-clustering with the export
    # (recluster/), nested outputs
    _run("doubletfinder", "--counts", DF_IN, "--out", out["doubletfinder_nested"], "--cluster-stats", "silhouette,purity",
         "--recluster", "--recluster-export", "--threads", "2")
    return out


@pytest.fixture
def copy_output(tmp_path):
    """A writable copy of an output, to tamper with."""

    def make(src: str, name: str) -> str:
        dst = tmp_path / name
        shutil.copytree(src, dst)
        return str(dst)

    return make


def tamper(dir_: str, path: str, field: str, value: str) -> None:
    """Edit one field of the entry whose path is ``path`` in summary.json, as text (``value`` is
    JSON text), so that every other number stays as written."""
    f = os.path.join(dir_, "summary.json")
    lines = open(f).read().split("\n")
    k = [i for i, line in enumerate(lines) if f'"path": "{path}",' in line]
    assert len(k) == 1
    for j in range(k[0], k[0] + 31):
        s = lines[j].lstrip()
        if s.startswith(f'"{field}": '):
            indent = lines[j][: len(lines[j]) - len(s)]
            comma = "," if lines[j].rstrip().endswith(",") else ""
            lines[j] = f'{indent}"{field}": {value}{comma}'
            open(f, "w").write("\n".join(lines))
            return
    raise AssertionError(f"no field {field} after {path}")
