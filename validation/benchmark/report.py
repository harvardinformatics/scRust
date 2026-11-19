#!/usr/bin/env python3
"""report.py: benchmark.md from analyze.py's tables (scRust data/benchmarks).

  PYTHONNOUSERSITE=1 python3 -I scripts/report.py

Reads benchmark_table.tsv, benchmark_kinds.tsv, benchmark_runs.tsv, sacct_workflow.tsv and
runs/calib/; writes benchmark.md (and prints the main table).
"""

import collections
import csv
import glob
import json
import os
import statistics

# the benchmark's data (validation/benchmark/ is two levels down from the repository)
B = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "benchmarks")
# where the tables are written (and report.py reads them): BENCH_OUT, else data/benchmarks/
OUT = os.environ.get("BENCH_OUT", B)


def rows(name):
    path = os.path.join(OUT, name)
    if not os.path.exists(path):
        path = os.path.join(B, name)  # inputs analyze.py does not write (sacct_workflow.tsv)
    return list(csv.DictReader(open(path), delimiter="\t"))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def f_s(v):
    """seconds, 3 significant digits"""
    if v != v:
        return "–"
    if v >= 100:
        return "%.0f" % v
    if v >= 10:
        return "%.1f" % v
    return "%.2f" % v


def f_gib(v):
    if v != v:
        return "–"
    if v >= 10:
        return "%.1f" % v
    if v >= 1:
        return "%.2f" % v
    return "%.3f" % v


def f_pct(v):
    if v != v:
        return "–"
    if v < 0.5:
        return "%.1f" % v if v >= 0.05 else "0"
    return "%.0f" % v


def f_x(v):
    if v != v:
        return "–"
    if v >= 10:
        return "%.0f×" % v
    return "%.1f×" % v


def f_int(v):
    return "{:,}".format(int(float(v)))


def f_big(v):
    v = float(v)
    if v >= 1e6:
        return "%.1f M" % (v / 1e6)
    return f_int(v)


def f_bytes(v):
    v = num(v)
    if v != v:
        return "–"
    if v >= 1e9:
        return "%.2f GB" % (v / 1e9)
    return "%.0f MB" % (v / 1e6) if v >= 1e7 else "%.1f MB" % (v / 1e6)


def side_cells(r, side):
    """'CPU / wall / peak' and 'I/O' cells for one side."""
    p = side + "_"
    if r.get(p + "reps", "0") in ("", "0"):
        return "–", "–"
    main = "%s / %s / %s" % (f_s(num(r[p + "cpu_total_s"])), f_s(num(r[p + "wall_total_s"])), f_gib(num(r[p + "peak_gib"])))
    no_write = side == "R" and r["rust_variant"] == "DoubletFinder alone"
    io = "%s + %s = %s%% · %s%%, %s" % (
        f_pct(num(r[p + "cpu_read_pct"])), "none" if no_write else f_pct(num(r[p + "cpu_write_pct"])),
        f_pct(num(r[p + "cpu_io_pct"])), f_pct(num(r[p + "mem_read_pct"])),
        "none" if no_write else f_pct(num(r[p + "mem_write_pct"])) + "%")
    return main, io


TOOL_LABEL = {
    ("cluster", "without --export"): "Clustering; Rust without `--export`",
    ("cluster", "with --export"): "Clustering; Rust with `--export`",
    # M15b: SCTransform with vars.to.regress = "percent.mt"
    ("cluster_regress", "without --export"): "M15b: clustering, vars.to.regress; Rust without `--export`",
    ("cluster_regress", "with --export"): "M15b: clustering, vars.to.regress; Rust with `--export`",
    ("cluster_regress", "without the regression"): "M15b: clustering, the same binary without vars.to.regress",
    ("doubletfinder", "DoubletFinder alone"): "DoubletFinder alone",
    ("doubletfinder", "whole rule"): "DoubletFinder rule, whole (no RunUMAP)",
    # M15, one binary for the four rows
    ("doubletfinder_m15", "DoubletFinder alone"): "M15: DoubletFinder alone, labels given",
    ("doubletfinder_m15", "DoubletFinder, clustering its own input"): "M15: DoubletFinder clustering its own input",
    ("doubletfinder_m15", "whole rule"): "M15: DoubletFinder rule, two commands",
    ("doubletfinder_m15", "whole rule, in process"): "M15: DoubletFinder rule, one process",
}


def m15_table(tab, t):
    """M15's rows at t threads: the hook's cost (clustering its own input against the labels
    given) and the in-process hand-off's saving (the whole rule in one process against the two
    commands), with R's whole rule."""
    by = {(r["object"], r["rust_variant"]): r for r in tab if r["tool"] == "doubletfinder_m15"}
    objs = []
    for r in tab:
        if r["tool"] == "doubletfinder_m15" and r["object"] not in objs:
            objs.append(r["object"])
    side = "Rust%s" % t

    def cell(r):
        if r is None or r.get(side + "_reps", "0") in ("", "0"):
            return "–"
        return "%s / %s / %s" % (f_s(num(r[side + "_cpu_total_s"])), f_s(num(r[side + "_wall_total_s"])), f_gib(num(r[side + "_peak_gib"])))

    def ratio(a, b, k):
        if a is None or b is None:
            return float("nan")
        x, y = num(a.get(side + "_" + k)), num(b.get(side + "_" + k))
        return x / y if x == x and y == y and y else float("nan")

    out = ["| object | DoubletFinder, labels given: CPU s / wall s / peak GiB | clustering its own input | its cost: CPU, wall, peak (÷ labels given) "
           "| rule, two commands | rule, one process | one process ÷ two commands: CPU, wall, peak | R's rule: CPU s / peak GiB | R ÷ one process: CPU, peak |",
           "|" + "---|" * 9]
    for obj in objs:
        a, h = by.get((obj, "DoubletFinder alone")), by.get((obj, "DoubletFinder, clustering its own input"))
        two, one = by.get((obj, "whole rule")), by.get((obj, "whole rule, in process"))
        rr = two if two is not None else one
        r_cpu, r_peak = num(rr.get("R_cpu_total_s")) if rr else float("nan"), num(rr.get("R_peak_gib")) if rr else float("nan")
        o_cpu, o_peak = (num(one.get(side + "_cpu_total_s")), num(one.get(side + "_peak_gib"))) if one else (float("nan"), float("nan"))
        out.append("| %s | %s | %s | %s, %s, %s | %s | %s | %s, %s, %s | %s / %s | %s, %s |" % (
            obj, cell(a), cell(h), f_x(ratio(h, a, "cpu_total_s")), f_x(ratio(h, a, "wall_total_s")), f_x(ratio(h, a, "peak_gib")),
            cell(two), cell(one), f_x(ratio(one, two, "cpu_total_s")), f_x(ratio(one, two, "wall_total_s")), f_x(ratio(one, two, "peak_gib")),
            f_s(r_cpu), f_gib(r_peak), f_x(r_cpu / o_cpu if o_cpu == o_cpu and o_cpu else float("nan")), f_x(r_peak / o_peak if o_peak == o_peak and o_peak else float("nan"))))
    return "\n".join(out)


def m15b_table(tab, t):
    """M15b's rows at t threads: the chain with SCTransform's vars.to.regress = "percent.mt"
    against R's workflow script as it was before the pipeline's commit 8d34de1, and what the
    regression costs each side: R's against the cluster rows' R runs (the current scripts, which
    differ in the argument only), Rust's against the same binary's chain without it (variant b)."""
    by = {(r["tool"], r["object"], r["rust_variant"]): r for r in tab}
    objs = []
    for r in tab:
        if r["tool"] == "cluster_regress" and r["object"] not in objs:
            objs.append(r["object"])
    side = "Rust%s" % t

    def get(r, s, k):
        return num(r.get("%s_%s" % (s, k))) if r is not None and r.get(s + "_reps", "0") not in ("", "0") else float("nan")

    def div(a, b):
        return a / b if a == a and b == b and b else float("nan")

    def f_r(v):
        # what the regression costs: a ratio near 1, to three decimals
        return "–" if v != v else "%.3f" % v

    out = ["| object | R with vars.to.regress: CPU s / peak GiB | R: ÷ without (CPU, peak) | Rust with: CPU s / wall s / peak GiB | Rust: ÷ without (CPU, peak) "
           "| with --export: CPU s / peak GiB | R ÷ Rust, with: CPU · peak | R ÷ Rust with --export: CPU · peak |",
           "|" + "---|" * 8]
    for obj in objs:
        rn, rx = by.get(("cluster_regress", obj, "without --export")), by.get(("cluster_regress", obj, "with --export"))
        bn = by.get(("cluster_regress", obj, "without the regression"))
        r_cpu, r_peak = get(rn, "R", "cpu_total_s"), get(rn, "R", "peak_gib")
        b_cpu, b_peak = get(bn, "R", "cpu_total_s"), get(bn, "R", "peak_gib")
        n_cpu, n_wall, n_peak = get(rn, side, "cpu_total_s"), get(rn, side, "wall_total_s"), get(rn, side, "peak_gib")
        x_cpu, x_peak = get(rx, side, "cpu_total_s"), get(rx, side, "peak_gib")
        bn_cpu, bn_peak = get(bn, side, "cpu_total_s"), get(bn, side, "peak_gib")
        out.append("| %s | %s / %s | %s, %s | %s / %s / %s | %s, %s | %s / %s | %s · %s | %s · %s |" % (
            obj, f_s(r_cpu), f_gib(r_peak), f_r(div(r_cpu, b_cpu)), f_r(div(r_peak, b_peak)),
            f_s(n_cpu), f_s(n_wall), f_gib(n_peak), f_r(div(n_cpu, bn_cpu)), f_r(div(n_peak, bn_peak)),
            f_s(x_cpu), f_gib(x_peak), f_x(div(r_cpu, n_cpu)), f_x(div(r_peak, n_peak)), f_x(div(r_cpu, x_cpu)), f_x(div(r_peak, x_peak))))
    return "\n".join(out)


def main_table(tab):
    head = ("| tool | object | input cells (genes; nonzeros) | input on disk: R; Rust | R: CPU s / wall s / peak GiB | R: I/O, CPU % (read + write) · memory % (reading, writing) "
            "| Rust 1 thread: CPU / wall / peak | Rust 1 thread: I/O | Rust 32 threads: CPU / wall / peak | Rust 32 threads: I/O "
            "| R ÷ Rust, CPU time (1; 32 threads) | R ÷ Rust, peak memory (1; 32 threads) |")
    out = [head, "|" + "---|" * 12]
    for r in tab:
        label = TOOL_LABEL[(r["tool"], r["rust_variant"])]
        rm, ri = side_cells(r, "R")
        m1, i1 = side_cells(r, "Rust1")
        m32, i32 = side_cells(r, "Rust32")
        out.append("| %s | %s | %s (%s; %s) | %s; %s | %s | %s | %s | %s | %s | %s | %s; %s | %s; %s |" % (
            label, r["object"], f_int(r["cells"]), f_int(r["genes"]), f_big(r["nonzeros"]),
            f_bytes(r.get("R_input_bytes")), f_bytes(r.get("Rust1_input_bytes")), rm, ri, m1, i1, m32, i32,
            f_x(num(r["cpu_ratio_R_over_Rust1"])), f_x(num(r["cpu_ratio_R_over_Rust32"])),
            f_x(num(r["peak_ratio_R_over_Rust1"])), f_x(num(r["peak_ratio_R_over_Rust32"]))))
    return "\n".join(out)


def ranges_table(tab):
    out = ["| tool | object | R: CPU s (range) | R: peak GiB (range) | Rust 1 thread: CPU s | Rust 1 thread: peak GiB | Rust 32 threads: CPU s | Rust 32 threads: peak GiB | repeats (R, Rust 1, Rust 32) |",
           "|" + "---|" * 9]
    for r in tab:
        label = TOOL_LABEL[(r["tool"], r["rust_variant"])]
        cells = []
        for side in ("R", "Rust1", "Rust32"):
            for k, f in (("cpu_total_s", f_s), ("peak_gib", f_gib)):
                p = "%s_%s" % (side, k)
                if r.get(side + "_reps", "0") in ("", "0"):
                    cells.append("–")
                else:
                    cells.append("%s (%s–%s)" % (f(num(r[p])), f(num(r[p + "_min"])), f(num(r[p + "_max"]))))
        out.append("| %s | %s | %s | %s, %s, %s |" % (label, r["object"], " | ".join(cells), r["R_reps"], r["Rust1_reps"], r["Rust32_reps"]))
    return "\n".join(out)


def soupx_inputs(runs):
    """R, rules soupx and soupx_emptydrops: CPU s of reading SoupX's own inputs (the raw matrix and
    the base object's .rds: phases read_raw and read_base_rds) and of the whole read, per object
    (medians over the runs that have these phases)."""
    import collections
    per = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in runs:
        if r["side"] != "R" or r["tool"] != "cluster" or r.get("included", "yes") != "yes":
            continue
        key = (r["object"], r["rep"])
        if r["kind"] == "read":
            per[key]["read"] += num(r["cpu_s"])
            if r["phase"] in ("read_raw", "read_base_rds"):
                per[key]["soupx"] += num(r["cpu_s"])
                per[key]["has"] = 1.0
    out = collections.defaultdict(list)
    for (obj, rep), d in per.items():
        if d.get("has"):
            out[obj].append((d["soupx"], d["read"]))
    return {o: (statistics.median([a for a, _ in v]), statistics.median([b for _, b in v]), len(v)) for o, v in out.items()}


def r_scripts_table(kinds, sacct, runs=None):
    """R: the whole script, what is left out (not ported), and the workflow job's sacct."""
    by = collections.defaultdict(list)
    for k in kinds:
        if k["side"] == "R" and k["variant"] == "workflow script":
            by[(k["tool"], k["object"])].append(k)
    sac = {}
    for s in sacct:
        if s["step"] == "step.0":
            sac[(s["tool"], s["object"])] = s
    sx = soupx_inputs(runs or [])
    out = ["| tool | object | rule | R's script, whole: CPU s / wall s / peak GiB | of which not ported, CPU s (shown apart) | of which R startup + packages, CPU s | of which reading SoupX's own inputs (raw matrix, base .rds), CPU s | like for like (main table): CPU s | workflow job (sacct, step .0): CPU s / elapsed s / MaxRSS GiB | job |",
           "|" + "---|" * 10]

    def med(lst, key):
        v = [num(x[key]) for x in lst]
        v = [x for x in v if x == x]
        return statistics.median(v) if v else float("nan")

    for (tool, obj), lst in sorted(by.items(), key=lambda kv: (kv[0][0], ORDER.index(kv[0][1]))):
        s = sac.get((tool, obj), {})
        rule = RULE.get((tool, obj), "")
        cpu = hms(s.get("total_cpu", ""))
        el = hms(s.get("elapsed", ""))
        rss = s.get("max_rss", "")
        rss_g = float(rss.rstrip("K")) / 1048576 if rss.endswith("K") else float("nan")
        sxv = sx.get(obj) if tool == "cluster" else None
        sxc = "%s of %s (%d runs)" % (f_s(sxv[0]), f_s(sxv[1]), sxv[2]) if sxv else "–"
        out.append("| %s | %s | %s | %s / %s / %s | %s | %s | %s | %s | %s / %s / %s | %s |" % (
            tool, obj, rule, f_s(med(lst, "cpu_whole_s")), f_s(med(lst, "wall_whole_s")), f_gib(med(lst, "peak_whole_gib")),
            f_s(med(lst, "cpu_unported_s")), f_s(med(lst, "cpu_startup_s")), sxc, f_s(med(lst, "cpu_total_s")),
            f_s(cpu), f_s(el), f_gib(rss_g), s.get("job", "")))
    return "\n".join(out)


def hms(t):
    """sacct's [D-]HH:MM:SS(.fff) or MM:SS.fff -> seconds"""
    if not t:
        return float("nan")
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)
    parts = [float(x) for x in t.split(":")]
    secs = 0.0
    for p in parts:
        secs = secs * 60 + p
    return days * 86400 + secs


ORDER = ["soupx_tenx_SRR18256555", "soupx_tenx_10xmousebrain", "filtered_seurat_tenx_10xmousebrain",
         "soupx_emptydrops_4500871lane", "soupx_emptydrops_sscovelli", "filtered_seurat_emptydrops_sscovelli",
         "cellbender_fromraw_ghirsutum2"]
RULE = {
    ("cluster", "soupx_tenx_SRR18256555"): "soupx (soupx.R)",
    ("cluster", "soupx_tenx_10xmousebrain"): "soupx (soupx.R)",
    ("cluster", "filtered_seurat_tenx_10xmousebrain"): "tenx2seuratrds",
    ("cluster", "soupx_emptydrops_4500871lane"): "soupx_emptydrops (soupx.R)",
    ("cluster", "soupx_emptydrops_sscovelli"): "soupx_emptydrops (soupx.R)",
    ("cluster", "filtered_seurat_emptydrops_sscovelli"): "emptydrops",
    ("cluster", "cellbender_fromraw_ghirsutum2"): "cellbender2seurat",
    ("doubletfinder", "cellbender_fromraw_ghirsutum2"): "doubletfinder_cellbender",
}
for o in ORDER:
    RULE.setdefault(("doubletfinder", o), "doubletfinder")
    # M15b: the same rules' scripts as they were before the pipeline's commit 8d34de1
    RULE[("cluster_regress", o)] = RULE[("cluster", o)] + ", at 8d34de1^ (vars.to.regress)"


def df_alone_table(runs, kinds):
    """DoubletFinder alone: R's DoubletFinder phases (paramSweep .. subset) against
    `scrust doubletfinder`; and the clustering of the singlets: R's against `scrust cluster`."""
    import collections
    seg = collections.defaultdict(list)
    for r in runs:
        if r["tool"] == "doubletfinder":
            seg[(r["side"], r["object"], r["variant"], r["rep"])].append(r)
    out = ["| object | R: DoubletFinder phases, CPU s / wall s / peak GiB | Rust `doubletfinder`, 1 thread: CPU / wall / peak | 32 threads | R ÷ Rust CPU (1; 32) | R: clustering of the singlets (no RunUMAP), CPU s | Rust `cluster --export` on singlets/, 1 thread: CPU s | 32 threads: CPU s |",
           "|" + "---|" * 8]

    def med(v):
        v = [x for x in v if x == x]
        return statistics.median(v) if v else float("nan")

    for obj in ORDER:
        rr = [v for k, v in seg.items() if k[0] == "R" and k[1] == obj]
        if not rr:
            continue
        r_cpu, r_wall, r_peak, r_rc = [], [], [], []
        for lst in rr:
            dfp = [x for x in lst if x["phase"].startswith("df_")]
            r_cpu.append(sum(num(x["cpu_s"]) for x in dfp))
            r_wall.append(sum(num(x["wall_s"]) for x in dfp))
            r_peak.append(max(num(x["rss_peak_kib"]) for x in dfp) / 1048576)
            r_rc.append(sum(num(x["cpu_s"]) for x in lst if x["kind"] == "recluster"))
        cells = []
        rust = {}
        for v in ("d1", "d32", "r1", "r32"):
            ks = [k for k in kinds if k["side"] == "Rust" and k["tool"] == "doubletfinder" and k["object"] == obj and k["variant"] == v]
            rust[v] = (med([num(k["cpu_whole_s"]) for k in ks]), med([num(k["wall_whole_s"]) for k in ks]), med([num(k["peak_whole_gib"]) for k in ks]))
        rc, rw, rp = med(r_cpu), med(r_wall), med(r_peak)
        out.append("| %s | %s / %s / %s | %s / %s / %s | %s / %s / %s | %s; %s | %s | %s | %s |" % (
            obj, f_s(rc), f_s(rw), f_gib(rp),
            f_s(rust["d1"][0]), f_s(rust["d1"][1]), f_gib(rust["d1"][2]),
            f_s(rust["d32"][0]), f_s(rust["d32"][1]), f_gib(rust["d32"][2]),
            f_x(rc / rust["d1"][0]) if rust["d1"][0] == rust["d1"][0] else "–",
            f_x(rc / rust["d32"][0]) if rust["d32"][0] == rust["d32"][0] else "–",
            f_s(med(r_rc)), f_s(rust["r1"][0]), f_s(rust["r32"][0])))
    return "\n".join(out)


def facts(kinds, other, runs):
    """Numbers the caveats quote, from the runs themselves."""
    import collections
    f = {}
    f["nodes_R"] = sorted(set(k["host"] for k in kinds if k["side"] == "R" and k["variant"] == "workflow script"))
    f["nodes_Rust"] = sorted(set(k["host"] for k in kinds if k["side"] == "Rust"))
    f["models"] = sorted(set(k["cpu_model"].replace("Intel(R)_Xeon(R)_Platinum_", "") for k in kinds))
    reps = collections.Counter()
    for k in kinds:
        if "+" in k["variant"]:
            continue
        reps[(k["side"], k["tool"], k["object"], k["variant"])] += 1
    f["reps_min"] = min(reps.values()) if reps else 0
    f["reps_max"] = max(reps.values()) if reps else 0
    # co-location of R runs with other runs of this benchmark
    rr = [k for k in kinds if k["side"] == "R" and k["variant"] == "workflow script"]
    f["r_runs"] = len(rr)
    f["r_colocated"] = sum(1 for k in rr if int(k.get("concurrent_bench_jobs", "0") or 0) > 0)
    f["r_coloc_max"] = max([int(k.get("concurrent_bench_jobs", "0") or 0) for k in rr] or [0])
    # spread of the repeats: the largest (max / min - 1) of the like-for-like CPU time, per side
    spread = {}
    for side in ("R", "Rust"):
        g = collections.defaultdict(list)
        for k in kinds:
            if k["side"] == side:
                g[(k["tool"], k["object"], k["variant"])].append(num(k["cpu_total_s"]))
        sp = [(max(v) / min(v) - 1.0, key) for key, v in g.items() if len(v) > 1 and min(v) > 0]
        spread[side] = max(sp) if sp else None
    f["spread"] = spread
    # the R runs against the workflow's own jobs (sacct step .0)
    sac = {(x["tool"], x["object"]): x for x in rows("sacct_workflow.tsv") if x["step"] == "step.0"}
    cmp_cpu, cmp_wall, cmp_peak = [], [], []
    g = collections.defaultdict(list)
    for k in rr:
        g[(k["tool"], k["object"])].append(k)
    for key, lst in g.items():
        x = sac.get(key)
        if not x:
            continue
        mc = statistics.median([num(k["cpu_whole_s"]) for k in lst])
        mw = statistics.median([num(k["wall_whole_s"]) for k in lst])
        mp = statistics.median([num(k["peak_whole_gib"]) for k in lst])
        wc, we = hms(x["total_cpu"]), hms(x["elapsed"])
        wr = float(x["max_rss"].rstrip("K")) / 1048576 if x["max_rss"].endswith("K") else float("nan")
        cmp_cpu.append((mc / wc, key))
        cmp_wall.append((mw / we, key))
        cmp_peak.append((mp - wr, key))
    f["vs_workflow"] = (cmp_cpu, cmp_wall, cmp_peak)
    # the same R runs on the other CPU model (8480+), against 8480CL
    by = collections.defaultdict(list)
    for k in kinds:
        if k["side"] == "R" and k["variant"] == "workflow script":
            by[(k["tool"], k["object"])].append(num(k["cpu_whole_s"]))
    ob = collections.defaultdict(list)
    for k in other:
        if k["side"] == "R" and k["variant"] == "workflow script":
            ob[(k["tool"], k["object"])].append(num(k["cpu_whole_s"]))
    ratios = []
    for key in ob:
        if key in by and by[key]:
            ratios.append(statistics.median(by[key]) / statistics.median(ob[key]))
    f["model_ratio"] = (min(ratios), max(ratios), len(ratios)) if ratios else None
    f["other_runs"] = len([k for k in other if k["variant"] == "workflow script"])
    # the pca writer's calibrated share
    moved = []
    for r in runs:
        if r["phase"] == "pca:write" and "calibrated" in r["peak_exact"]:
            try:
                moved.append(float(r["peak_exact"].split("+")[1].split(" s")[0]))
            except (IndexError, ValueError):
                pass
    f["pca_moved"] = (min(moved), max(moved)) if moved else None
    return f


def run_facts(f):
    """The section on nodes and repeats, from the runs."""
    lines = ["### Nodes and repeats", ""]
    lines.append("- **CPU model:** every run in the tables ran on an Intel Xeon Platinum %s (sapphire partition). "
                 "Node names map to models: holy8a24xxx and 26xxx are 8480+, the rest are 8480CL. "
                 "The jobs exclude the 8480+ nodes and check the model before running." % ", ".join(f["models"]))
    lines.append("- **Nodes:** R on %d (%s); Rust on %d (%s)." % (
        len(f["nodes_R"]), ", ".join(f["nodes_R"]), len(f["nodes_Rust"]), ", ".join(f["nodes_Rust"])))
    lines.append("- **Repeats:** %s per side and object (Rust: per variant and thread count), each in its own process, with its inputs evicted from the page cache first. The tables give medians; the ranges are in the ranges table and in `benchmark_table.tsv` (`_min`, `_max`)." % (str(f["reps_min"]) if f["reps_min"] == f["reps_max"] else "%d to %d" % (f["reps_min"], f["reps_max"])))
    if f["model_ratio"]:
        lo, hi, n = f["model_ratio"]
        lines.append("- **The 8480+ runs:** the first %d R runs landed on 8480+ nodes before the model was pinned. "
                     "They are kept in `benchmark_other_cpu.tsv` and left out of the tables. On %d objects, R's median CPU time on the 8480CL is %.2f to %.2f times its median on the 8480+." % (
                         f["other_runs"], n, lo, hi))
    lines.append("- **Shared nodes:** %d of %d R runs shared their node with other runs of this benchmark (at most %d at once), besides other users' jobs. "
                 "That is how the workflow's own 1-CPU jobs run. Rust's jobs had 33 CPUs each." % (
                     f["r_colocated"], f["r_runs"], f["r_coloc_max"]))
    for side in ("R", "Rust"):
        sp = f["spread"].get(side)
        if sp:
            lines.append("- **Repeat spread, %s:** the CPU time of one object's repeats differs by at most %.0f%% (%s, %s%s)." % (
                side, 100 * sp[0], sp[1][0], sp[1][1], "" if sp[1][2] in ("workflow script",) else ", " + sp[1][2]))
    cpu, wall, peak = f["vs_workflow"]
    if cpu:
        above = sum(1 for d, _ in peak if d > 0)
        lines.append("- **Against the workflow's own jobs (sacct, step .0):** the R scripts' whole-run CPU time is %.2f to %.2f times the workflow job's TotalCPU, and their wall time %.2f to %.2f times its Elapsed, on %d objects. "
                     "Their exact peak is above the job's MaxRSS on %d of %d objects, by up to %.1f GiB; on the others it is within %.1f GiB of it." % (
                         min(c for c, _ in cpu), max(c for c, _ in cpu), min(c for c, _ in wall), max(c for c, _ in wall), len(cpu),
                         above, len(peak), max(d for d, _ in peak), max(-d for d, _ in peak if d <= 0) if any(d <= 0 for d, _ in peak) else 0.0))
    if f["pca_moved"]:
        lines.append("- **The pca writer's calibration:** in Rust's clustering runs, %.2f to %.2f s of CPU were moved from compute to write (see the caveats)." % f["pca_moved"])
    return "\n".join(lines) + "\n"


def main():
    tab = rows("benchmark_table.tsv")
    kinds = rows("benchmark_kinds.tsv")
    sacct = rows("sacct_workflow.tsv")
    md = []
    md.append("# Benchmark: the workflow's R steps against the Rust ports\n")
    md.append(main_table(tab))
    if any(r["tool"] == "doubletfinder_m15" for r in tab):
        for t in ("1", "32"):
            md.append("\n## M15: clustering on the fly and the re-clustering in one process, %s thread%s\n" % (t, "" if t == "1" else "s"))
            md.append(m15_table(tab, t))
    if any(r["tool"] == "cluster_regress" for r in tab):
        for t in ("1", "32"):
            md.append("\n## M15b: the chain with SCTransform's vars.to.regress, %s thread%s\n" % (t, "" if t == "1" else "s"))
            md.append(m15b_table(tab, t))
    md.append("\n## Ranges\n")
    md.append(ranges_table(tab))
    md.append("\n## R's scripts, whole, and the workflow's jobs\n")
    md.append(r_scripts_table(kinds, sacct, rows("benchmark_runs.tsv")))
    text = "\n".join(md) + "\n"
    notes = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")
    if os.path.exists(notes):
        text += "\n" + open(notes).read()
    other = rows("benchmark_other_cpu.tsv") if os.path.exists(os.path.join(OUT, "benchmark_other_cpu.tsv")) else []
    text += "\n" + run_facts(facts(kinds, other, rows("benchmark_runs.tsv")))
    open(os.path.join(OUT, "benchmark.md"), "w").write(text)
    print(main_table(tab))


if __name__ == "__main__":
    main()
