#!/usr/bin/env python3
"""checks.py: each run's output checks, in one table (scRust data/benchmarks).

  PYTHONNOUSERSITE=1 python3 -I scripts/checks.py > checks.tsv

R runs: the number of clusters (and their ids) against the workflow's log (the clustering panel's
log_final_clusters; for DoubletFinder, Rust's clustering of the same singlets), and the size of
the output .rds against the size of the workflow's saved .rds.
Rust runs: same.txt (outputs against the reference runs: identical / identical data lines /
DIFFERENT) and the exit code.
"""
import csv, glob, json, os, re, sys

# the repository (validation/benchmark/ is two levels down), and the benchmark's data
S = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = os.path.join(S, "data", "benchmarks")
PIPE = "/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing"

panel = {r["object"]: r for r in csv.DictReader(open(os.path.join(S, "validation", "clustering_panel.tsv")), delimiter="\t")}
dfsets = {r["object"]: r for r in csv.DictReader(open(os.path.join(S, "validation", "doubletfinder_datasets.tsv")), delimiter="\t")}

print("\t".join(["side", "tool", "object", "rep", "variant", "check", "value", "reference", "result"]))
for d in sorted(glob.glob(os.path.join(B, "runs", "r", "*", "*", "rep*"))):
    if ".partial." in d:
        continue
    tool, obj, rep = d.split(os.sep)[-3:]
    outs = {}
    for line in open(os.path.join(d, "outputs.tsv")):
        size, rel, full = line.rstrip("\n").split("\t")
        outs[rel] = int(size)
    rds = [k for k in outs if k.endswith(".rds")]
    saved = panel[obj]["saved"] if tool == "cluster" else dfsets[obj]["saved"]
    saved_size = os.path.getsize(os.path.join(PIPE, saved))
    for k in rds:
        print("\t".join(["R", tool, obj, rep, "workflow script", "output .rds bytes", str(outs[k]), str(saved_size),
                         "same" if outs[k] == saved_size else "differs by %+d" % (outs[k] - saved_size)]))
    chk = open(os.path.join(d, "check.txt")).read()
    m = re.search(r"nclusters=(\d+)", chk)
    ncl = m.group(1) if m else "?"
    if tool == "cluster":
        want = panel[obj]["log_final_clusters"]
    else:
        # the workflow's DoubletFinder job log's final clustering: Rust's clustering of the same
        # singlets (identical singlets: cells.tsv identical to R's direct run) as the reference
        rr = sorted(glob.glob(os.path.join(B, "runs", "rust", "doubletfinder", obj, "rep*", "r1", "nclusters.txt")))
        want = open(rr[0]).read().strip() if rr else "?"
    print("\t".join(["R", tool, obj, rep, "workflow script", "clusters", ncl, want, "same" if ncl == want else "DIFFERENT"]))
for d in sorted(glob.glob(os.path.join(B, "runs", "rust", "*", "*", "rep*", "*"))):
    if ".partial." in d or not os.path.isdir(d):
        continue
    tool, obj, rep, v = d.split(os.sep)[-4:]
    pm = json.load(open(os.path.join(d, "procmon.json")))
    print("\t".join(["Rust", tool, obj, rep, v, "exit code", str(pm["exit_code"]), "0", "same" if pm["exit_code"] == 0 else "FAILED"]))
    sf = os.path.join(d, "same.txt")
    if os.path.exists(sf):
        lines = open(sf).read().splitlines()
        ident = sum(1 for l in lines if l.endswith(": identical"))
        data = sum(1 for l in lines if "identical data lines" in l)
        diff = sum(1 for l in lines if "DIFFERENT" in l)
        print("\t".join(["Rust", tool, obj, rep, v, "outputs vs reference", "%d identical, %d identical data lines" % (ident, data), "%d files" % len(lines), "same" if diff == 0 else "%d DIFFERENT" % diff]))
