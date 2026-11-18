"""The export's guard against being read as raw counts (docs/PLAN.md C3f, C5's K11; milestone M13):
scanpy's read_10x_mtx on every directory of a `scrust cluster --export` output and of a
`scrust sctransform --corrected-counts` output, without prefix= and with every prefix that could
name a file there. It must fail everywhere but on the export's counts/ (and on the top directory
with prefix="counts/", which also names counts/), which it must read as the input counts.

read_10x_mtx(path, prefix=p) reads <p>matrix.mtx.gz, <p>features.tsv.gz and <p>barcodes.tsv.gz
(or <p>matrix.mtx, <p>genes.tsv, <p>barcodes.tsv when <p>genes.tsv exists), so a prefix can name
a file only if it is a beginning of that file's path below the directory: every beginning of every
file's path is tried, which is exhaustive.

Run in the hermetic env (envs/loaders-py.yml), isolated from the user's site:
    PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I validation/python/guard_read10x.py \\
        <out.tsv> <export dir> <input counts dir> [<sct dir> ...]
Exit status 1 if any attempt behaves otherwise.
"""

import os
import sys
import warnings


def files_under(root):
    out = []
    for d, _, fs in os.walk(root):
        for f in fs:
            out.append(os.path.relpath(os.path.join(d, f), root))
    return sorted(out)


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    import numpy as np
    import scanpy as sc

    out, export, input_dir = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
    others = [os.path.abspath(a) for a in sys.argv[4:]]
    want = sc.read_10x_mtx(input_dir)
    rows = []
    ok_all = True
    warnings.simplefilter("ignore")
    for root in [export] + others:
        dirs = sorted({root} | {os.path.join(root, os.path.dirname(f)) for f in files_under(root) if os.path.dirname(f)})
        dirs = [d for d in dirs if "/stages" not in d[len(root):]]
        for d in dirs:
            rel = os.path.relpath(d, root)
            files = [f for f in files_under(d) if not f.startswith("stages/")]
            prefixes = {""}
            for f in files:
                prefixes.update(f[:k] for k in range(1, len(f) + 1))
            n_tried = n_read = 0
            for p in sorted(prefixes):
                n_tried += 1
                try:
                    a = sc.read_10x_mtx(d, prefix=p)
                    read, detail = True, f"read: {a.n_obs} cells x {a.n_vars} genes"
                except Exception as e:  # noqa: BLE001 - any failure is a refusal to read
                    read, detail = False, f"{type(e).__name__}: {str(e).splitlines()[0][:160] if str(e) else ''}"
                counts_dir = root == export and ((rel == "counts" and p == "") or (rel == "." and p == "counts/"))
                ok = read == counts_dir
                if read and counts_dir:
                    same = a.shape == want.shape and (a.obs_names == want.obs_names).all() and (a.var_names == want.var_names).all() and (a.X != want.X).nnz == 0
                    detail += "; identical to the input" if same else "; DIFFERS from the input"
                    ok = ok and bool(same)
                n_read += read
                if read or p == "" or not ok:
                    rows.append((os.path.basename(root), rel, repr(p), "read" if counts_dir else "fail", read, ok, detail))
                    print(f"{os.path.basename(root):28s} {rel:12s} prefix={p!r:40s} {'ok' if ok else 'WRONG'} {detail}")
                ok_all = ok_all and ok
            rows.append((os.path.basename(root), rel, f"all {n_tried} prefixes", "fail but counts", n_read > 0, True, f"{n_read} of {n_tried} read"))
            print(f"{os.path.basename(root):28s} {rel:12s} {n_tried} prefixes tried, {n_read} read")
    with open(out, "w") as f:
        f.write("output\tdir\tprefix\texpected\tread\tas_required\tdetail\n")
        for r in rows:
            f.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")
    from importlib.metadata import version
    print(f"scanpy {version('scanpy')} (anndata {version('anndata')}) read_10x_mtx: {'all as required' if ok_all else 'NOT as required'}")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
