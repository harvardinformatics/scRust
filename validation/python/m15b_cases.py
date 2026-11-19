"""Milestone M15b: the cases of the vars.to.regress validation and their inputs.

Writes, from each object's exported counts (no R, no .rds):

- data/clustering/m15b/inputs/<object>/mt_genes.txt: the mitochondrial features' ids, one per
  line (column 1 of the input's features.tsv.gz), built from the features' names (column 2):
    * mouse (the regex finds 13 features): the names the workflow's "(?i)^mt-" matches;
    * zebrafish (soupx_tenx_drerio): the two "mt-" rRNAs and the 13 protein-coding genes named
      ND1-ND6, ND4L, COX1-COX3, ATP6, ATP8, CYTB (the regex finds only the rRNAs);
    * bearded dragon (soupx_tenx_pvitticeps): the 13 protein-coding genes, named as above (the
      regex finds none);
    * tier T (the test data has no mitochondrial gene): a stand-in set, the ribosomal protein
      genes (names starting with Rpl or Rps), to exercise the list with non-zero values.
  Each list also has a comment line, a blank line, one id twice and two ids no feature has, so the
  parsing and the report of ids not found are exercised. R and Rust read the same file.
- zero_genes.txt: the ids of three features with no counts in any cell, for objects where the
  regex finds no feature: R's run regresses its all-zero percent.mt (the rank-deficient path),
  which scRust refuses with the regex; scRust's run gives the same all-zero column through this
  list, so the two compare file by file.
- metadata.tsv: barcode, score (a number per cell, written in mixed forms: 17 and 6 significant
  digits, exponents, integers, signs) and label (text), in a shuffled cell order, for the
  --metadata case.
- validation/m15b_cases.tsv: one row per (object, case): the R reference run's options and the
  Rust run's.

Run from the repository's top directory:
  PYTHONNOUSERSITE=1 python -I validation/python/m15b_cases.py
"""

import csv
import gzip
import os
import random
import sys

ROOT = os.getcwd()
OUT = os.path.join(ROOT, "data", "clustering", "m15b", "inputs")
MT_OTHER = ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"]
NOT_FOUND = ["NOT_A_FEATURE_1", "NOT_A_FEATURE_2"]


def read_features(inp):
    with gzip.open(os.path.join(inp, "features.tsv.gz"), "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f]


def read_barcodes(inp):
    with gzip.open(os.path.join(inp, "barcodes.tsv.gz"), "rt") as f:
        return [line.rstrip("\n") for line in f]


def rows_with_entries(inp, nrow):
    seen = [False] * nrow
    with gzip.open(os.path.join(inp, "matrix.mtx.gz"), "rt") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        for line in f:
            seen[int(line.split(" ", 1)[0]) - 1] = True
    return seen


def seurat_like(name):
    return name.replace("_", "-").replace("|", "-")


def matches_mt(name):
    n = seurat_like(name)
    return len(n) >= 3 and n[0] in "mM" and n[1] in "tT" and n[2] == "-"


def write_list(path, ids, what):
    dup = ids[0] if ids else None
    lines = [f"# {what}", *ids[: len(ids) // 2], "", *ids[len(ids) // 2 :]]
    if dup is not None:
        lines.append(dup)
    lines += NOT_FOUND
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def score_text(rng, i):
    x = rng.uniform(-3.0, 50.0)
    k = i % 6
    if k == 0:
        return "%.17g" % x
    if k == 1:
        return "%.6g" % x
    if k == 2:
        return "%.4e" % x
    if k == 3:
        return str(int(x))
    if k == 4:
        return "%.3f" % (x / 1000.0)
    return "+%.15g" % abs(x)


def main():
    datasets = list(csv.DictReader(open("validation/clustering_datasets.tsv"), delimiter="\t"))
    objects = []
    for d in datasets:
        if d["tier"] in ("T", "P"):
            objects.append((d["object"], d["tier"], f"data/clustering/{d['object']}/input", d["assay"], d["level"]))
    for o in ("soupx_tenx_drerio", "soupx_tenx_pvitticeps"):
        objects.append((o, "X", f"data/doubletfinder/{o}/input", "v5", "standard"))
    rows = []
    for obj, tier, inp, assay, level in objects:
        feats = read_features(inp)
        names = [f[1] if len(f) > 1 else f[0] for f in feats]
        ids = [f[0] for f in feats]
        regex = [i for i, n in enumerate(names) if matches_mt(n)]
        other = [i for i, n in enumerate(names) if n.upper() in MT_OTHER]
        d = os.path.join(OUT, obj)
        os.makedirs(d, exist_ok=True)
        if tier == "T":
            listed = [i for i, n in enumerate(names) if n.startswith("Rpl") or n.startswith("Rps")]
            what = "tier T has no mitochondrial gene: a stand-in set, the ribosomal protein genes (Rpl*, Rps*)"
        elif obj == "soupx_tenx_drerio":
            listed = sorted(set(regex) | set(other))
            what = "zebrafish: the two mt- rRNAs the regex finds, and ND1-6, ND4L, COX1-3, ATP6, ATP8, CYTB"
        elif obj == "soupx_tenx_pvitticeps":
            listed = other
            what = "bearded dragon: ND1-6, ND4L, COX1-3, ATP6, ATP8, CYTB (the regex finds none)"
        elif regex:
            listed = regex
            what = "mouse: the features the workflow's (?i)^mt- matches"
        else:
            listed = []
            what = ""
        mt_file = ""
        if listed:
            mt_file = os.path.join(d, "mt_genes.txt")
            write_list(mt_file, [ids[i] for i in listed], what)
        zero_file = ""
        if not regex:
            seen = rows_with_entries(inp, len(feats))
            zero = [ids[i] for i, s in enumerate(seen) if not s][:3]
            if len(zero) < 3:
                sys.exit(f"{obj}: fewer than three features without counts")
            zero_file = os.path.join(d, "zero_genes.txt")
            with open(zero_file, "w") as f:
                f.write("# three features with no counts in any cell: percent.mt is 0 in every cell\n" + "\n".join(zero) + "\n")
        bc = read_barcodes(inp)
        rng = random.Random(f"m15b-{obj}")
        order = list(range(len(bc)))
        rng.shuffle(order)
        md_file = os.path.join(d, "metadata.tsv")
        with open(md_file, "w") as f:
            f.write("barcode\tscore\tlabel\n")
            for i in order:
                f.write(f"{bc[i]}\t{score_text(rng, i)}\t{'ab'[i % 2]}\n")
        rel = lambda p: os.path.relpath(p, ROOT) if p else ""
        # the cases: R's options, and the Rust run's where they differ
        cases = []
        large_without_mt = tier == "P" and not regex
        cases.append(("pm", "percent.mt", "yes", "", "", rel(zero_file) if not regex else ""))
        if listed and not large_without_mt:
            cases.append(("ids", "percent.mt", "yes", rel(mt_file), "", rel(mt_file)))
            cases.append(("ids_ncount", "percent.mt,nCount_RNA", "yes", rel(mt_file), "", rel(mt_file)))
        if not large_without_mt or obj == "filtered_seurat_emptydrops_sscovelli":
            if tier != "X":
                cases.append(("meta", "score", "no", "", rel(md_file), ""))
        for case, vars_, pm, mt, md, rust_mt in cases:
            rows.append({
                "object": obj, "tier": tier, "input": inp, "assay": assay, "level": level, "case": case,
                "vars_to_regress": vars_, "percent_mt": pm, "mt_genes": mt, "metadata": md,
                "rust_mt_genes": rust_mt if rust_mt else mt,
                "regex_features": str(len(regex)), "listed_features": str(len(listed)),
            })
        print(f"{obj}: regex {len(regex)}, listed {len(listed)}, cases {', '.join(c[0] for c in cases)}")
    with open("validation/m15b_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} cases in validation/m15b_cases.tsv")


if __name__ == "__main__":
    main()
