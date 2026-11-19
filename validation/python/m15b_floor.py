"""Milestone M15b: the R-vs-R floor across OpenBLAS kernels with SCTransform's vars.to.regress
([clustering.floor.kernel_regressed], validation/tolerances.toml).

For every case of validation/m15b_cases.tsv that has both R reference runs
(validation/slurm/m15b_references.sbatch: the production profile, OpenBLAS's SapphireRapids kernel,
and the portable one, Haswell, both on Sapphire Rapids), compares the two runs' files with the
statistics kernel-gated grading uses (scio::grade::stats, ported here line for line: bitwise
differences, a column's sign flipped when its negation is closer in L1, max abs, max scaled = max
abs over the reference's largest absolute value, max rel, rows with a difference), on the arrays
`scrust validate cluster` grades against the floor in these runs:
  S2d_scale_data            sct/scale_data.npy, rows matched by name (features x cells)
  S3_embeddings_dims_1_30   pca/embeddings.npy, dims 1:30 (cells x 30), signs aligned
  S3_embeddings_dims_31_50  dims 31:50, signs aligned
  S3_loadings               pca/loadings.npy (features x 50), signs aligned
  S3_sdev                   pca/stdev.npy
  total_variance            pca/pca.json's total_variance
  K8_silhouette_width       metadata.tsv's silhouette_width
  K9_purity                 metadata.tsv's neighborhood_purity
each in both directions (either run as the reference), and the checks that stay exact: the labels
and levels, the graphs, percent.mt and the counts' columns, the variable features, the corrected
counts, the regression's rank and pivot. As [clustering.floor.kernel] was derived (M8, M11), each
bound is the maximum over the cases, rounded up (two significant digits), and a count
(cells_differing, values_differing) is the total over the cases, so it bounds one case too.

With --rust <label> <profile>, compares instead each case's production R run with scrust's run of
that label and profile (validation/slurm/m15b_cluster.sbatch's run/, scrust's layout): the
comparison kernel-gated grading makes.

  PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I validation/python/m15b_floor.py \\
      [--tsv <rows.tsv>] [--summary <summary.md>] [--rust <label> <profile>]
"""

import argparse
import csv
import glob
import gzip
import hashlib
import json
import math
import os

import numpy as np

ROOT = os.getcwd()
KEYS = ["S2d_scale_data", "S3_embeddings_dims_1_30", "S3_embeddings_dims_31_50", "S3_loadings", "S3_sdev",
        "total_variance", "K8_silhouette_width", "K9_purity"]
# the bounds each entry takes: those of [clustering.floor.kernel]'s entry of the same name
# (total_variance, identical across kernels without the regression, takes max_abs and max_scaled)
BOUNDS = {
    "S2d_scale_data": ["max_abs", "max_fraction_differing"],
    "S3_embeddings_dims_1_30": ["max_abs", "max_scaled"],
    "S3_embeddings_dims_31_50": ["max_abs"],
    "S3_loadings": ["max_abs", "max_scaled"],
    "S3_sdev": ["max_abs", "max_scaled"],
    "total_variance": ["max_abs", "max_scaled"],
    "K8_silhouette_width": ["max_abs"],
    "K9_purity": ["cells_differing", "max_abs"],
}
# each bound's statistic
STAT = {"max_abs": "max_abs", "max_scaled": "max_scaled", "max_fraction_differing": "fraction_differing",
        "cells_differing": "rows_differ", "values_differing": "n_differ"}


def stats(a, b, align):
    """scio::grade::stats on two arrays of shape (nrow, ncol), compared column by column."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    out = {"n": int(a.size), "n_differ": 0, "rows_differ": 0, "max_abs": 0.0, "max_scaled": 0.0, "max_rel": 0.0,
           "flipped": 0, "incomparable": False}
    if a.shape != b.shape:
        out["incomparable"] = True
        out["n_differ"] = max(a.size, b.size)
        return out
    if a.size == 0:
        return out
    absa = np.abs(a)
    max_a = float(np.max(np.where(np.isnan(absa), 0.0, absa)))
    if align:
        same = np.abs(a - b).sum(axis=0)
        opp = np.abs(a + b).sum(axis=0)
        flip = opp < same  # a NaN sum compares false: no flip, as in Rust
        b = np.where(flip[None, :], -b, b)
        out["flipped"] = int(flip.sum())
    both_nan = np.isnan(a) & np.isnan(b)
    differ = (a.view(np.uint64) != b.view(np.uint64)) & ~both_nan
    out["n_differ"] = int(differ.sum())
    out["rows_differ"] = int(differ.any(axis=1).sum())
    one_nan = differ & (np.isnan(a) != np.isnan(b))
    out["incomparable"] = bool(one_nan.any())
    num = differ & ~one_nan
    if num.any():
        x = a[num]
        d = np.abs(x - b[num])
        out["max_abs"] = float(d.max())
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(x == 0.0, np.where(d == 0.0, 0.0, np.inf), d / np.abs(x))
        out["max_rel"] = float(rel.max())
    m = out["max_abs"]
    out["max_scaled"] = m / max_a if max_a > 0 else (math.inf if m > 0 else 0.0)
    out["fraction_differing"] = out["n_differ"] / out["n"]
    return out


def table(path):
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        rows = list(r)
    return h, rows


def col(t, name):
    h, rows = t
    if name not in h:
        return None
    i = h.index(name)
    return [row[i] for row in rows]


def nums(v):
    return np.array([[np.nan if x in ("NA", "NaN", "") else float(x)] for x in v])


# readers for both layouts: R's runs (validation/R/lib/clustering.R) and scrust's
def scale_data(d):
    f = os.path.join(d, "sct", "scale_data.npy")
    if not os.path.isfile(f):
        return None, None
    fr = os.path.join(d, "sct", "scale_data_features.tsv")
    if os.path.isfile(fr):
        names = open(fr).read().split("\n")[:-1]
    else:
        names = [r[1] for r in table(os.path.join(d, "sct", "variable_features.tsv"))[1]]
    return names, np.load(f)


def pca_array(d, what):
    f = os.path.join(d, "pca", what + ".npy")
    if os.path.isfile(f):
        return np.load(f)
    h, rows = table(os.path.join(d, "pca", what + ".tsv"))
    skip = sum(1 for x in h if not x.startswith("PC_") and x != "stdev")
    m = np.array([[float(x) for x in row[skip:]] for row in rows])
    return m[:, 0] if what == "stdev" else m


def total_variance(d):
    f = os.path.join(d, "pca", "pca.json")
    if os.path.isfile(f):
        v = json.load(open(f)).get("total_variance")
        return v[0] if isinstance(v, list) else v
    return json.load(open(os.path.join(d, "summary.json")))["result"]["pca"]["total_variance"]


def same_file(fa, fb):
    if fa.endswith(".gz"):
        with gzip.open(fa) as x, gzip.open(fb) as y:
            return x.read() == y.read()
    h = [hashlib.sha256(open(f, "rb").read()).hexdigest() for f in (fa, fb)]
    return h[0] == h[1]


def exact_checks(a, b, rust):
    """The checks that stay exact: True where identical."""
    out = {}
    ma, mb = table(os.path.join(a, "metadata.tsv")), table(os.path.join(b, "metadata.tsv"))
    for c in ("seurat_clusters", "SCT_snn_res.0.8", "percent.mt", "nCount_RNA", "nFeature_RNA", "silhouette_other", "purity_maximum"):
        x, y = col(ma, c), col(mb, c)
        if x is None or y is None:
            continue
        if c in ("percent.mt", "nCount_RNA"):  # as doubles (R's %.17g and scrust's text)
            x = [v if v == "NA" else float(v).hex() for v in x]
            y = [v if v == "NA" else float(v).hex() for v in y]
        out["metadata.tsv " + c] = x == y
    for f in ("seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt"):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if os.path.isfile(fa) and os.path.isfile(fb):
            out[f] = open(fa, "rb").read() == open(fb, "rb").read()
    if not rust:
        for f in ("graphs/nn.mtx.gz", "graphs/snn.mtx.gz", "sct/corrected_counts.mtx.gz", "sct/variable_features.tsv"):
            fa, fb = os.path.join(a, f), os.path.join(b, f)
            if os.path.isfile(fa) and os.path.isfile(fb):
                out[f] = same_file(fa, fb)
        ra = json.load(open(os.path.join(a, "summary.json")))["result"].get("regression", {})
        rb = json.load(open(os.path.join(b, "summary.json")))["result"].get("regression", {})
        out["the regression's rank and pivot"] = (ra.get("rank"), ra.get("pivot")) == (rb.get("rank"), rb.get("pivot"))
    return out


def measure(a, b):
    """The statistics of each graded array, b against a (a the reference)."""
    out = {}
    na, sa = scale_data(a)
    nb, sb = scale_data(b)
    if sa is not None and sb is not None:
        pos = {x: i for i, x in enumerate(nb)}
        out["S2d_scale_data"] = stats(sa, sb[[pos[x] for x in na], :], False)  # features x cells
    ea, eb = pca_array(a, "embeddings"), pca_array(b, "embeddings")
    out["S3_embeddings_dims_1_30"] = stats(ea[:, :30], eb[:, :30], True)
    if ea.shape[1] > 30:
        out["S3_embeddings_dims_31_50"] = stats(ea[:, 30:], eb[:, 30:], True)
    out["S3_loadings"] = stats(pca_array(a, "loadings"), pca_array(b, "loadings"), True)
    out["S3_sdev"] = stats(pca_array(a, "stdev")[:, None], pca_array(b, "stdev")[:, None], False)
    out["total_variance"] = stats(np.array([[total_variance(a)]]), np.array([[total_variance(b)]]), False)
    ma, mb = table(os.path.join(a, "metadata.tsv")), table(os.path.join(b, "metadata.tsv"))
    for key, c in (("K8_silhouette_width", "silhouette_width"), ("K9_purity", "neighborhood_purity")):
        x, y = col(ma, c), col(mb, c)
        if x is not None and y is not None:
            out[key] = stats(nums(x), nums(y), False)
    return out


def runs(obj, case, rust=None):
    base = os.path.join(ROOT, "data", "clustering", "m15b", "runs", obj, case)
    prod = [d for d in sorted(glob.glob(os.path.join(base, "direct-production-*"))) if ".partial" not in d]
    if rust:
        d = os.path.join(ROOT, "data", "clustering", "m15b", "rust", rust[0], obj, case, rust[1], "run")
        port = [d] if os.path.isdir(d) else []
    else:
        port = [d for d in sorted(glob.glob(os.path.join(base, "direct-portable-*"))) if ".partial" not in d]
    return (prod[0] if prod else None), (port[0] if port else None)


def round_up(x):
    """x rounded up to two significant digits (never below x)."""
    if x == 0 or not math.isfinite(x):
        return x
    e = math.floor(math.log10(abs(x)))
    q = 10.0 ** (e - 1)
    v = float("%.2g" % (math.ceil(x / q - 1e-9) * q))
    while v < x:
        v = float("%.2g" % (v + q))
    return v


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv")
    p.add_argument("--summary")
    p.add_argument("--rust", nargs=2, metavar=("LABEL", "PROFILE"))
    args = p.parse_args()
    cases = list(csv.DictReader(open(os.path.join(ROOT, "validation", "m15b_cases.tsv")), delimiter="\t"))
    rows, exact, measured = [], {}, []
    for c in cases:
        prod, port = runs(c["object"], c["case"], args.rust)
        if not (prod and port):
            continue
        measured.append(c)
        exact[(c["object"], c["case"])] = exact_checks(prod, port, bool(args.rust))
        directions = [("production reference", prod, port)]
        if not args.rust:
            directions.append(("portable reference", port, prod))
        for direction, a, b in directions:
            for key, s in measure(a, b).items():
                rows.append({"object": c["object"], "tier": c["tier"], "case": c["case"], "vars_to_regress": c["vars_to_regress"],
                             "key": key, "direction": direction, **s})
    if args.tsv:
        cols = list(rows[0].keys())
        with open(args.tsv, "w") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    tiers = sorted({c["tier"] for c in measured})
    count = {t: sum(1 for c in measured if c["tier"] == t) for t in tiers}
    what = f"R's production runs against scrust {args.rust[0]} {args.rust[1]}" if args.rust else "R's production runs against R's portable runs, both directions"
    lines = [f"# M15b: {what}", "", f"{len(measured)} cases: " + ", ".join(f"{count[t]} tier {t}" for t in tiers) + ".", ""]
    lines.append("| key | statistic | maximum | set by | " + " | ".join(f"tier {t}" for t in tiers) + " | bound (rounded up) |")
    lines.append("|---|---|---|---|" + "---|" * len(tiers) + "---|")
    toml = []
    for key in KEYS:
        kr = [r for r in rows if r["key"] == key]
        if not kr:
            continue
        entry = []
        shown = BOUNDS[key] + (["values_differing"] if key == "K9_purity" else [])
        for bound in shown:
            st = STAT[bound]
            best = max(kr, key=lambda r: r[st])
            v = best[st]
            count = bound in ("cells_differing", "values_differing")
            if count:
                # the total over the cases (one direction: a count is the same either way)
                one = [r for r in kr if r["direction"] == "production reference"]
                b = sum(r[st] for r in one)
                per_tier = [str(sum(r[st] for r in one if r["tier"] == t)) for t in tiers]
                setby = f"total over {sum(1 for r in one if r[st])} case(s); at most {v} in one ({best['object']} {best['case']})"
                entry.append(f"{bound} = {b}") if bound in BOUNDS[key] else None
                lines.append(f"| {key} | {bound} | {b} | {setby} | " + " | ".join(per_tier) + f" | {b} |")
                continue
            per_tier = []
            for t in tiers:
                tv = [r[st] for r in kr if r["tier"] == t]
                per_tier.append("%.4g" % max(tv) if tv else "-")
            b = round_up(v)
            if bound in BOUNDS[key]:
                entry.append(f"{bound} = {b:.2g}")
            lines.append(f"| {key} | {bound} | {v:.4g} | {best['object']} {best['case']} ({best['direction']}) | "
                         + " | ".join(per_tier) + f" | {b:.2g} |")
        n_cases = len({(r["object"], r["case"]) for r in kr if r["n_differ"]})
        toml.append(f"{key} = {{ {', '.join(entry)}, cases_differing = {n_cases} }}")
    lines += ["", "Exact checks (identical in every case unless listed):", ""]
    for name in sorted({k for v in exact.values() for k in v}):
        bad = [f"{o} {c}" for (o, c), v in exact.items() if name in v and not v[name]]
        n = sum(1 for v in exact.values() if name in v)
        lines.append(f"- {name}: {n - len(bad)} of {n} identical" + (f"; differing: {', '.join(bad)}" if bad else ""))
    lines += ["", "The bounds, as TOML (counts as measured):", "", "```toml"] + toml + ["```", ""]
    text = "\n".join(lines)
    if args.summary:
        open(args.summary, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
