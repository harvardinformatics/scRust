"""Milestone M15b: scrust's percent.mt against the saved workflow objects'.

  one <metadata.tsv> <saved export dir> [--json <file>]
      compares one run's metadata.tsv (barcode, nCount_RNA, nFeature_RNA, percent.mt) with
      validation/R/export_percent_mt.R's percent_mt.tsv and percent_mt.json: the cells in order,
      nCount_RNA, nFeature_RNA and percent.mt as %.17g text (so bit for bit), the features the
      pattern matched (the run's summary.json record against R's grep on the saved object), and
      R's own facts (percent.mt computed again on the saved object equals its saved column; the
      saved counts equal the export scrust read). Exit status 1 if anything differs.
  report <saved_rust dir> [--tsv <file>]
      gathers every object's compare.json into one table and prints the totals.
"""

import argparse
import json
import os
import sys


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]
    return header, rows


def col(t, name):
    h, rows = t
    k = h.index(name)
    return [r[k] for r in rows]


def one(args):
    run = read_tsv(args.metadata)
    saved = read_tsv(os.path.join(args.saved, "percent_mt.tsv"))
    facts = json.load(open(os.path.join(args.saved, "percent_mt.json")))
    summary_dir = os.path.dirname(args.metadata)
    s = json.load(open(os.path.join(summary_dir, "summary.json")))
    rec = s["result"]["sctransform"]["percent_mt"]
    out = {"metadata": args.metadata, "saved": args.saved, "n_cells": len(run[1])}
    checks = {}
    checks["barcodes"] = col(run, "barcode") == col(saved, "barcode")
    for c in ("nCount_RNA", "nFeature_RNA", "percent.mt"):
        a, b = col(run, c), col(saved, c)
        checks[c] = a == b
        out[f"{c}_differing"] = sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else len(a)
    feats = facts.get("features") or []
    if isinstance(feats, str):
        feats = [feats]
    checks["features"] = rec["n_features"] == len(feats) and rec["features_first"] == feats[: len(rec["features_first"])]
    checks["r_recomputed_identical"] = bool(facts.get("recomputed_identical"))
    checks["saved_counts_are_the_export"] = bool(facts.get("counts_identical_to_export", True))
    out["n_features"] = rec["n_features"]
    out["features"] = feats
    pm = [float(v) for v in col(run, "percent.mt")]
    out["cells_nonzero"] = sum(v != 0.0 for v in pm)
    out["max"] = max(pm) if pm else None
    out["checks"] = checks
    out["identical"] = all(checks.values())
    if args.json:
        json.dump(out, open(args.json, "w"), indent=2)
    print(f"{args.metadata}: {out['n_cells']} cells, {out['n_features']} features match, {out['cells_nonzero']} cells non-zero; "
          + ("identical" if out["identical"] else "DIFFERS: " + ", ".join(k for k, v in checks.items() if not v)))
    return 0 if out["identical"] else 1


def report(args):
    rows = []
    for obj in sorted(os.listdir(args.dir)):
        p = os.path.join(args.dir, obj, "compare.json")
        if ".partial." in obj or not os.path.isfile(p):
            continue
        c = json.load(open(p))
        kind = "doubletfinder" if "singlets_clusters" in c["metadata"] else "clustering"
        rows.append((obj, kind, c["n_cells"], c["n_features"], c["cells_nonzero"], c["max"], c["identical"], ",".join(k for k, v in c["checks"].items() if not v)))
    if args.tsv:
        with open(args.tsv, "w") as f:
            f.write("object\tkind\tcells\tfeatures_matched\tcells_nonzero\tmax_percent_mt\tidentical\tdiffering\n")
            for r in rows:
                f.write("\t".join("" if v is None else (repr(v) if isinstance(v, float) else str(v)) for v in r) + "\n")
    for kind in ("clustering", "doubletfinder"):
        k = [r for r in rows if r[1] == kind]
        print(f"{kind}: {len(k)} objects, {sum(r[6] for r in k)} identical, {sum(r[2] for r in k)} cells, "
              f"{sum(1 for r in k if r[3] > 0)} with features matched, {sum(r[4] for r in k)} cells non-zero")
    return 0 if all(r[6] for r in rows) else 1


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("one")
    a.add_argument("metadata")
    a.add_argument("saved")
    a.add_argument("--json")
    b = sub.add_parser("report")
    b.add_argument("dir")
    b.add_argument("--tsv")
    args = p.parse_args()
    sys.exit(one(args) if args.cmd == "one" else report(args))


if __name__ == "__main__":
    main()
