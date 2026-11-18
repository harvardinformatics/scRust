"""Command line: ``python -m scrustio h5ad <scrust output> <object.h5ad> [--check-sha256]``, and
``python -m scrustio manifest <scrust output>`` (the files by role)."""

from __future__ import annotations

import argparse
import sys
from collections import Counter

from .anndata_loader import write_h5ad
from .manifest import ScrustError, read_manifest


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m scrustio", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("h5ad", help="write a scrust output as an .h5ad file (AnnData)")
    h.add_argument("output", help="the scrust output directory (it holds summary.json)")
    h.add_argument("h5ad", help="the .h5ad file to write")
    h.add_argument("--check-sha256", action="store_true", help="check every listed file's SHA-256 first")
    mf = sub.add_parser("manifest", help="list a scrust output's files by role")
    mf.add_argument("output")
    mf.add_argument("--check-sha256", action="store_true")
    a = p.parse_args(argv)
    try:
        if a.cmd == "h5ad":
            adata = write_h5ad(a.output, a.h5ad, check_sha256=a.check_sha256)
            print(f"wrote {a.h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")
        else:
            m = read_manifest(a.output, check_sha256=a.check_sha256)
            m.check()
            roles = Counter(e["role"] for e in m.entries)
            print(f"scrust output ({m.command}), {m.dir}: {len(m.entries)} files")
            for role, n in sorted(roles.items()):
                print(f"  {role}: {n}")
    except ScrustError as e:
        print(f"scrustio: {e}", file=sys.stderr)
        return 1
    return 0
