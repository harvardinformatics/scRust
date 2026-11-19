"""Milestone M15's tables (docs/validation/clustering.md, "M15"), from the panel's runs
(validation/slurm/m15_panel.sbatch) and the T and P re-clustering checks
(validation/slurm/m15_recluster.sbatch). Reads plain files only.

  python -I validation/python/report_m15.py <label> [--tsv <per-object.tsv>] > tables.md

Per object of validation/doubletfinder_panel.tsv, from data/doubletfinder/<object>/rust/<label>/
production/: the cells and singlets; each check (c_step1, k12_labels, k12_files, k12_saved,
k12_logged, a_files, b_saved, manifests); the clusters found on the fly and in the re-clustering;
peak RSS and wall time of the four runs, and M7's; where the in-process run's and the hook's peaks
were reached (summary.json's phases); the node, CPU model and OpenBLAS kernel.
"""

import csv
import filecmp
import json
import os
import re
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CHECKS = ["c_step1", "k12_labels", "k12_files", "k12_saved", "k12_logged", "a_files", "b_saved", "manifests"]
GIB = 1024.0 * 1024.0


def jload(p):
    with open(p) as f:
        return json.load(f)


def lines(p):
    with open(p) as f:
        return f.read().splitlines()


def fmt_int(v):
    return f"{v:,}"


def gib(kib):
    return float(kib) / GIB


def check_files(j):
    """compare_runs.py's report: (identical, compared)."""
    rows = [r for r in j["rows"] if not r["status"].endswith("allowed")]
    return sum(1 for r in rows if r["status"].startswith("identical")), len(rows)


def obj_row(obj, label, logs):
    base = os.path.join(ROOT, "data", "doubletfinder", obj)
    d = os.path.join(base, "rust", label, "production")
    if not os.path.isdir(d):
        return None
    row = list(csv.DictReader(open(os.path.join(d, "row.tsv")), delimiter="\t"))[0]
    out = {"object": obj, "arm": re.sub(r"_[^_]+$", "", obj) if not obj.startswith("cellbender") else "cellbender_fromraw"}
    for k in ("soupx_tenx", "soupx_emptydrops", "decontx_raw_tenx", "decontx_raw_emptydrops", "cellbender_fromraw"):
        if obj.startswith(k + "_"):
            out["arm"] = k
    out["cells"] = int(row["cells"])
    out["singlets"] = int(row["singlets"])
    c = os.path.join(d, "checks")
    st = {}
    i, n = check_files(jload(os.path.join(c, "c_step1.json")))
    st["c_step1"] = (i == n and n > 0, f"{i} of {n} files")
    inp = os.path.join(base, "input")
    same = all(filecmp.cmp(os.path.join(d, "hook", "run", "clusters", f), os.path.join(inp, f), shallow=False)
               for f in ("seurat_clusters.tsv", "seurat_clusters_levels.txt"))
    st["k12_labels"] = (same, f"{out['cells']} cells")
    i, n = check_files(jload(os.path.join(c, "k12_files.json")))
    st["k12_files"] = (i == n and n > 0, f"{i} of {n} files")
    ks = jload(os.path.join(c, "k12_saved.json"))
    made = len(ks["checks"])
    ident = sum(1 for x in ks["checks"] if x["status"] == "identical")
    elements = sum(x.get("n", 0) for x in ks["checks"])
    st["k12_saved"] = (made == 9 and ident == 9, f"{ident} of {made} checks, {elements} values")
    out["k12_saved_values"] = elements
    kl = jload(os.path.join(c, "k12_logged.json"))
    st["k12_logged"] = (kl["pass"], "nodes, edges, modularity, communities, singletons, clusters")
    i, n = check_files(jload(os.path.join(c, "a_files.json")))
    st["a_files"] = (i == n and n > 0, f"{i} of {n} files")
    out["a_files_n"] = n
    bs = jload(os.path.join(c, "b_saved.json"))
    bvals = sum(x["compared"] for x in bs["checks"] if x["check"] in ("seurat_clusters", "SCT_snn_res.0.8", "silhouette_width", "neighborhood_purity"))
    st["b_saved"] = (bs["pass"], f"{sum(1 for x in bs['checks'] if x['status'] == 'identical')} of {len(bs['checks'])} checks")
    out["b_values"] = bvals
    log = logs.get(obj, "")
    st["manifests"] = ("manifests: identical" in log, "hook/, inprocess/")
    for k in CHECKS:
        out[k] = "identical" if st[k][0] else "DIFFERS"
        out[k + "_detail"] = st[k][1]
    out["pass"] = all(st[k][0] for k in CHECKS)
    out["clusters_on_the_fly"] = len(lines(os.path.join(d, "hook", "run", "clusters", "seurat_clusters_levels.txt")))
    out["clusters_recluster"] = len(lines(os.path.join(d, "inprocess", "run", "recluster", "seurat_clusters_levels.txt")))
    res = jload(os.path.join(d, "hook", "run", "summary.json"))["result"]
    out["pK"] = res["pK_label"]
    out["nExp"] = int(res["nExp"])
    for k in ("m7", "given", "hook", "recluster2", "inprocess"):
        out[k + "_peak_gib"] = gib(row[k + "_peak_kib"])
        out[k + "_wall_s"] = float(row[k + "_wall_s"])
    # where each in-process peak was reached: summary.json records the peak at the end of each phase
    ph = jload(os.path.join(d, "inprocess", "run", "summary.json"))["run"]["phases"]
    peaks = [v["peak_rss_kib"] for v in ph.values() if v.get("peak_rss_kib")]
    out["inprocess_df_peak_gib"] = gib(ph["singlets"]["peak_rss_kib"])
    out["inprocess_rose"] = max(peaks) > ph["singlets"]["peak_rss_kib"]
    ph = jload(os.path.join(d, "hook", "run", "summary.json"))["run"]["phases"]
    peaks = [v["peak_rss_kib"] for v in ph.values() if v.get("peak_rss_kib")]
    out["hook_cluster_peak_gib"] = gib(ph["clustering: write"]["peak_rss_kib"])
    out["hook_set_by_clustering"] = max(peaks) == ph["clustering: write"]["peak_rss_kib"]
    fp = jload(os.path.join(d, "fingerprint.json"))
    out["node"] = row["node"]
    out["cpu_model"] = fp.get("cpu_model", "")
    out["openblas"] = fp.get("openblas_corename", "")
    return out


TP = ["soupx_tenx_test", "soupx_emptydrops_test", "decontx_raw_tenx_test", "decontx_raw_emptydrops_test", "cellbender_fromraw_test",
      "soupx_tenx_SRR18256555", "soupx_tenx_10xmousebrain", "soupx_emptydrops_4500871lane", "soupx_emptydrops_sscovelli",
      "cellbender_fromraw_ghirsutum2"]


def tp_table(label):
    """The T and P re-clustering checks (m15_recluster.sbatch), both profiles."""
    out = ["| object | profile | (a) in process ÷ two commands | (a) with the export | (b) against the saved re-clustering | clustering on the fly against the saved input |",
           "|---|---|---|---|---|---|"]
    n = 0
    for obj in TP:
        for prof in ("production", "portable"):
            d = os.path.join(ROOT, "data", "doubletfinder", obj, "rust", label, prof + "_recluster", "checks")
            if not os.path.isdir(d):
                continue
            n += 1
            cells = []
            for f in ("a_files.json", "a_export.json"):
                i, k = check_files(jload(os.path.join(d, f)))
                cells.append(f"{'identical' if i == k else 'DIFFERS'}, {i} of {k} files")
            for f in ("b_saved.json", "hook_saved.json"):
                j = jload(os.path.join(d, f))
                diff = [c["check"] for c in j["checks"] if c["status"] == "differs"]
                made = [c for c in j["checks"] if c["status"] != "not made"]
                cells.append(f"identical, {len(made)} checks" if j["pass"] else "differs: " + ", ".join(diff))
            out.append(f"| {obj} | {prof} | " + " | ".join(cells) + " |")
    return "\n".join(out) if n else ""


def main(argv):
    args = argv[1:]
    tsv = None
    if "--tsv" in args:
        k = args.index("--tsv")
        tsv = args[k + 1]
        del args[k:k + 2]
    if len(args) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    label = args[0]
    panel = list(csv.DictReader(open(os.path.join(ROOT, "validation", "doubletfinder_panel.tsv")), delimiter="\t"))
    logs = {}
    logdir = os.path.join(ROOT, "data", "doubletfinder", "logs", "m15")
    for f in sorted(os.listdir(logdir)):
        if not f.startswith("scrust-m15-df-"):
            continue
        text = open(os.path.join(logdir, f)).read()
        m = re.search(r"== task \d+: (\S+) on", text)
        if m and f"scrust-{label} " in text:
            logs[m.group(1)] = text
    rows = [r for r in (obj_row(p["object"], label, logs) for p in panel) if r is not None]
    n = len(rows)
    print(f"Label `{label}`: {n} of {len(panel)} objects run; {sum(r['pass'] for r in rows)} passing every check.\n")
    print("| check | objects | identical | compared |")
    print("|---|---|---|---|")
    for k in CHECKS:
        ok = sum(1 for r in rows if r[k] == "identical")
        what = {
            "c_step1": "files (cells.tsv, the sweep tables, singlets/) and result, against M7's run",
            "k12_labels": "seurat_clusters.tsv and its levels, against the input's",
            "k12_files": "DoubletFinder's files and result, against M7's run",
            "k12_saved": f"M7's nine checks against the saved output: {fmt_int(sum(r['k12_saved_values'] for r in rows))} values",
            "k12_logged": "the input's job's logged numbers",
            "a_files": f"in process against two commands: {fmt_int(sum(r['a_files_n'] for r in rows))} files",
            "b_saved": f"the saved re-clustering: {fmt_int(sum(r['b_values'] for r in rows))} labels, silhouette widths and purities",
            "manifests": "scrust validate manifest, nested outputs",
        }[k]
        print(f"| {k} | {n} | {ok} | {what} |")
    print()
    hook_m7 = [(r["hook_peak_gib"] / r["m7_peak_gib"], r["object"]) for r in rows]
    hook_given = [(r["hook_peak_gib"] / r["given_peak_gib"], r["object"]) for r in rows]
    inproc = [(r["inprocess_peak_gib"] / max(r["given_peak_gib"], r["recluster2_peak_gib"]), r["object"]) for r in rows]
    inproc_sum = [(r["inprocess_peak_gib"] / (r["given_peak_gib"] + r["recluster2_peak_gib"]), r["object"]) for r in rows]
    inproc_df = [(r["inprocess_df_peak_gib"] / r["given_peak_gib"], r["object"]) for r in rows]
    t_ratio = [(r["inprocess_wall_s"] / (r["given_wall_s"] + r["recluster2_wall_s"]), r["object"]) for r in rows]
    h_ratio = [(r["hook_wall_s"] / r["given_wall_s"], r["object"]) for r in rows]

    def span(v):
        v = sorted(v)
        return f"{v[0][0]:.3f} ({v[0][1]}) to {v[-1][0]:.3f} ({v[-1][1]}), median {statistics.median(x for x, _ in v):.3f}"

    print("| measure (8 threads, the same node for the four runs) | ratio: smallest, largest, median |")
    print("|---|---|")
    print(f"| peak RSS, clustering on the fly ÷ M7's run | {span(hook_m7)} |")
    print(f"| peak RSS, clustering on the fly ÷ --clusters (same binary, same node) | {span(hook_given)} |")
    print(f"| peak RSS, whole rule in process ÷ the larger of the two commands' | {span(inproc)} |")
    print(f"| peak RSS, whole rule in process ÷ the sum of the two commands' | {span(inproc_sum)} |")
    print(f"| peak RSS, the in-process run's DoubletFinder phases (to the end of `singlets`) ÷ --clusters' run, the same code | {span(inproc_df)} |")
    print(f"| wall time, clustering on the fly ÷ --clusters | {span(h_ratio)} |")
    print(f"| wall time, whole rule in process ÷ the two commands' sum | {span(t_ratio)} |")
    print(f"| objects whose peak with clustering on the fly is more than 10% above M7's | {sum(1 for x, _ in hook_m7 if x > 1.10)} |")
    print(f"| objects whose peak with clustering on the fly was reached in the clustering's phases | {sum(1 for r in rows if r['hook_set_by_clustering'])} |")
    print(f"| objects whose in-process peak rose during the re-clustering | {sum(1 for r in rows if r['inprocess_rose'])} |")
    print()
    tot = lambda k: sum(r[k] for r in rows)
    print(f"Peak RSS, largest: M7 {max(r['m7_peak_gib'] for r in rows):.2f} GiB, --clusters {max(r['given_peak_gib'] for r in rows):.2f}, "
          f"on the fly {max(r['hook_peak_gib'] for r in rows):.2f}, the singlets' clustering alone {max(r['recluster2_peak_gib'] for r in rows):.2f}, "
          f"in process {max(r['inprocess_peak_gib'] for r in rows):.2f}. Wall time, total: M7 {tot('m7_wall_s'):.0f} s, --clusters {tot('given_wall_s'):.0f} s, "
          f"on the fly {tot('hook_wall_s'):.0f} s, two commands {tot('given_wall_s') + tot('recluster2_wall_s'):.0f} s, in process {tot('inprocess_wall_s'):.0f} s.\n")
    print("| arm | objects | passing | cells | singlets | clusters on the fly (input's) | clusters, re-clustering |")
    print("|---|---|---|---|---|---|---|")
    for arm in ("soupx_tenx", "soupx_emptydrops", "decontx_raw_tenx", "decontx_raw_emptydrops", "cellbender_fromraw"):
        a = [r for r in rows if r["arm"] == arm]
        if not a:
            continue
        print(f"| {arm} | {len(a)} | {sum(r['pass'] for r in a)} | {min(r['cells'] for r in a):,} to {max(r['cells'] for r in a):,} | "
              f"{min(r['singlets'] for r in a):,} to {max(r['singlets'] for r in a):,} | {min(r['clusters_on_the_fly'] for r in a)} to {max(r['clusters_on_the_fly'] for r in a)} | "
              f"{min(r['clusters_recluster'] for r in a)} to {max(r['clusters_recluster'] for r in a)} |")
    print()
    cpus = sorted(set(r["cpu_model"] for r in rows))
    cores = sorted(set(r["openblas"] for r in rows))
    print(f"Nodes: {len(set(r['node'] for r in rows))}; CPU models: {', '.join(cpus)}; OpenBLAS kernels: {', '.join(cores)}.\n")
    print("| object | cells | singlets | clusters (on the fly; re-clustering) | pK | nExp | checks | peak GiB: M7, --clusters, on the fly, cluster alone, in process | wall s: M7, --clusters, on the fly, cluster alone, in process |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['object']} | {r['cells']:,} | {r['singlets']:,} | {r['clusters_on_the_fly']}; {r['clusters_recluster']} | {r['pK']} | {r['nExp']:,} | "
              f"{'8 of 8' if r['pass'] else ', '.join(k for k in CHECKS if r[k] != 'identical')} | "
              f"{r['m7_peak_gib']:.2f}, {r['given_peak_gib']:.2f}, {r['hook_peak_gib']:.2f}, {r['recluster2_peak_gib']:.2f}, {r['inprocess_peak_gib']:.2f} | "
              f"{r['m7_wall_s']:.1f}, {r['given_wall_s']:.1f}, {r['hook_wall_s']:.1f}, {r['recluster2_wall_s']:.1f}, {r['inprocess_wall_s']:.1f} |")
    tp = tp_table(label)
    if tp:
        print()
        print(tp)
    if tsv:
        cols = ["object", "arm", "cells", "singlets", "clusters_on_the_fly", "clusters_recluster", "pK", "nExp"] + CHECKS + \
               [k + "_detail" for k in CHECKS] + ["pass"] + [f"{k}_{m}" for m in ("peak_gib", "wall_s") for k in ("m7", "given", "hook", "recluster2", "inprocess")] + \
               ["inprocess_df_peak_gib", "inprocess_rose", "hook_cluster_peak_gib", "hook_set_by_clustering", "node", "cpu_model", "openblas"]
        with open(tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
