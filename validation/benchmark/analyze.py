#!/usr/bin/env python3
"""analyze.py: the benchmark's tables from its runs (scRust data/benchmarks).

  PYTHONNOUSERSITE=1 python3 -I scripts/analyze.py [--include-pilot]

Reads runs/r/<tool>/<object>/rep<k>/ (run_r.sbatch) and runs/rust/<tool>/<object>/rep<k>/<variant>/
(run_rust.sbatch); rep0 is the pilot, left out unless --include-pilot. Writes:
  benchmark_runs.tsv   one row per run and segment (each phase, split into read / compute /
                       write where files were open), plus one "whole" row per run
  benchmark_kinds.tsv  one row per run: CPU, wall and memory by kind, the I/O shares
  benchmark_table.tsv  one row per tool x object (x Rust variant): medians and ranges
(benchmark.md is written by report.py from these.)

Measurement, per side:
- R: the phase markers in the R process (bench_markers.R): proc.time() (getrusage of the
  process and its waited-for children) for CPU; VmHWM, reset to the current RSS at each phase's
  start (5 > /proc/self/clear_refs), for the phase's own peak; VmRSS at its start and end;
  /proc/self/io. Startup = from the fork (the sampler's clock) to the first mark; exit = from the
  last mark to the process's end, with the rest of getrusage's CPU (wait4).
- Rust: the sampler (procmon.py, 10 ms): /proc/<pid>/stat (utime + stime of the whole thread
  group), status (VmRSS, VmHWM), io, and the open files (/proc/<pid>/fd, fdinfo), aligned with
  summary.json's run.phases (cumulative ends, from Phases::new at the start of the command:
  taken as the fork, +-10 ms). Each phase's default kind comes from its name; within a phase, the
  span from the first to the last sample with an output file open (any mode: written, or re-read
  for the manifest's digests) is write, and with an input file open (and no output) read.
  getrusage (wait4) gives the exact totals: CPU and peak RSS (= VmHWM; never reset).
"""

import argparse
import collections
import csv
import glob
import gzip
import json
import os
import statistics
import sys

# the repository (validation/benchmark/ is two levels down), and the benchmark's data
S = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = os.path.join(S, "data", "benchmarks")
# where the tables are written (and report.py reads them): BENCH_OUT, else data/benchmarks/
OUT = os.environ.get("BENCH_OUT", B)
CLK = 100.0

OBJECTS = [
    "soupx_tenx_SRR18256555",
    "soupx_tenx_10xmousebrain",
    "filtered_seurat_tenx_10xmousebrain",
    "soupx_emptydrops_4500871lane",
    "soupx_emptydrops_sscovelli",
    "filtered_seurat_emptydrops_sscovelli",
    "cellbender_fromraw_ghirsutum2",
]
DF_OBJECTS = [o for o in OBJECTS if o.startswith(("soupx_", "cellbender_"))]

# M15's DoubletFinder runs (run_rust.sbatch doubletfinder_m15): with the labels given (d), the
# two-command route's clustering of singlets/ (r), clustering on the fly (h), and the whole rule in
# one process, the singlets re-clustered from memory with the export (p)
DF_M15 = "doubletfinder_m15"
DF_TOOLS = ("doubletfinder", DF_M15)

# Rust phases' default kinds; M15's in-process clusterings record their stages as
# "clustering: <stage>" and "recluster: <stage>" (compute unless named here)
RUST_KIND = {
    "load": "read", "write": "write", "manifest": "write", "singlets": "write",
    "(head)": "startup", "(tail)": "write",
    "clustering: digest": "read", "clustering: write": "write",
    "recluster: load": "read", "recluster: write": "write", "recluster: manifest": "write", "recluster": "write",
}

# The binary whose chain wrote each stage's files right after the stage (sct/, pca/, neighbors/ in
# the stage's phase); later ones (M15 on: the chain's entry point returns its result, which is
# written in the write phase) write them all in the write phase.
STAGE_WRITER = "scrust-b6f6e4e"
KINDS = ["startup", "read", "compute", "recluster", "unported", "write", "exit"]


def fnum(x):
    return float(x.strip())


def parse_env(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        for tok in line.split():
            if "=" in tok:
                k, _, v = tok.partition("=")
                out.setdefault(k, v)
    return out


# ---------------------------------------------------------------------------------------------
# R runs


def parse_r(run_dir):
    env = parse_env(os.path.join(run_dir, "env.txt"))
    pm = json.load(open(os.path.join(run_dir, "procmon.json")))
    rows = [{k: (v or "").strip() for k, v in r.items()}
            for r in csv.DictReader(open(os.path.join(run_dir, "markers.tsv")), delimiter="\t")]
    t0 = pm["start_epoch_s"]
    cpu_keys = ("user_s", "sys_s", "child_user_s", "child_sys_s")
    io_keys = ("rchar", "wchar", "read_bytes", "write_bytes")

    def cpu(r):
        return sum(fnum(r[k]) for k in cpu_keys)

    segs = []
    prev_end = None  # the previous phase's end row: CPU, time and io continue from it
    start_row = None
    for r in rows:
        if r["event"] == "start":
            start_row = r
            continue
        # an end row closes the phase
        if prev_end is None:
            a_t, a_cpu, a_io, rss0, a_sys = 0.0, 0.0, {k: 0.0 for k in io_keys}, 0.0, 0.0
        else:
            a_t = fnum(prev_end["epoch_s"]) - t0
            a_cpu = cpu(prev_end)
            a_sys = fnum(prev_end["sys_s"]) + fnum(prev_end["child_sys_s"])
            a_io = {k: fnum(prev_end[k]) for k in io_keys}
            rss0 = fnum(start_row["vmrss_kib"]) if start_row is not None else fnum(prev_end["vmrss_kib"])
        b_t = fnum(r["epoch_s"]) - t0
        seg = {
            "phase": r["phase"], "kind": r["kind"], "t_start_s": a_t, "t_end_s": b_t,
            "wall_s": b_t - a_t, "cpu_s": cpu(r) - a_cpu,
            "sys_s": fnum(r["sys_s"]) + fnum(r["child_sys_s"]) - a_sys,
            "rss_start_kib": rss0, "rss_peak_kib": fnum(r["vmhwm_kib"]), "rss_end_kib": fnum(r["vmrss_kib"]),
            "peak_exact": "yes",
        }
        for k in io_keys:
            seg[k] = fnum(r[k]) - a_io[k]
        if r["phase"] == "startup":
            seg["rss_start_kib"] = 0.0
        segs.append(seg)
        prev_end = r
        start_row = None
    # exit: the last mark to the end of the process (the rest of getrusage's CPU)
    end_t = pm["wall_s"]
    last = segs[-1]
    segs.append({
        "phase": "exit", "kind": "exit", "t_start_s": last["t_end_s"], "t_end_s": end_t,
        "wall_s": end_t - last["t_end_s"], "cpu_s": pm["rusage"]["cpu_s"] - cpu(prev_end),
        "sys_s": pm["rusage"]["stime_s"] - fnum(prev_end["sys_s"]) - fnum(prev_end["child_sys_s"]),
        "rss_start_kib": last["rss_end_kib"], "rss_peak_kib": float("nan"), "rss_end_kib": float("nan"),
        "peak_exact": "no", "rchar": float("nan"), "wchar": float("nan"),
        "read_bytes": float("nan"), "write_bytes": float("nan"),
    })
    peak = max(s["rss_peak_kib"] for s in segs if s["rss_peak_kib"] == s["rss_peak_kib"])
    whole = {"wall_s": pm["wall_s"], "cpu_s": pm["rusage"]["cpu_s"], "peak_kib": peak,
             "rusage_maxrss_kib": pm["rusage"]["maxrss_kib"], "sampler_cpu_s": pm["sampler_cpu_s"],
             "exit_code": pm["exit_code"]}
    inputs = sum(e["bytes"] for e in pm.get("evicted", []))
    return env, pm, segs, whole, inputs


# ---------------------------------------------------------------------------------------------
# Rust runs


def load_samples(path):
    T, C, R, H, IO, OPEN, SY = [], [], [], [], [], [], []
    with gzip.open(path, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        ix = {k: i for i, k in enumerate(head)}
        for row in rd:
            T.append(float(row[ix["t_s"]]))
            C.append((int(row[ix["utime"]]) + int(row[ix["stime"]]) + int(row[ix["cutime"]]) + int(row[ix["cstime"]])) / CLK)
            SY.append((int(row[ix["stime"]]) + int(row[ix["cstime"]])) / CLK)
            R.append(float(row[ix["vmrss_kib"]]))
            H.append(float(row[ix["vmhwm_kib"]]))
            IO.append(tuple(float(row[ix[k]]) for k in ("rchar", "wchar", "read_bytes", "write_bytes")))
            OPEN.append((row[ix["open_r"]], row[ix["open_w"]]))
    return T, C, R, H, IO, OPEN, SY


def interp(T, V, t):
    """V at time t, linear between samples (V given at the times T, ascending)."""
    if t <= T[0]:
        return V[0]
    if t >= T[-1]:
        return V[-1]
    lo, hi = 0, len(T) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if T[mid] <= t:
            lo = mid
        else:
            hi = mid
    if T[hi] == T[lo]:
        return V[hi]
    f = (t - T[lo]) / (T[hi] - T[lo])
    return V[lo] + f * (V[hi] - V[lo])


def roots_of(open_field):
    """r<i>:path@pos;... -> set of root indices."""
    if open_field == "-":
        return set()
    out = set()
    for item in open_field.split(";"):
        head = item.split(":", 1)[0]
        if head.startswith("r") and head[1:].isdigit():
            out.add(int(head[1:]))
    return out


# which phase writes (w) or re-reads (r) each output file, by its path in --out
def output_phase(tool, rel, mode, late=False):
    """The phase that writes (or re-reads) an output file. late: the binary writes the chain's
    stage files in the write phase (STAGE_WRITER's did it after each stage)."""
    if rel == "summary.json":
        return "(tail)"  # written after the last phase mark
    if tool == "cluster":
        if rel.startswith("counts/"):
            return "load"
        if "w" not in mode:
            return "manifest"
        if rel == "README.txt":
            return "manifest"
        if not late:
            for prefix, phase in (("sct/", "sctransform"), ("pca/", "pca"), ("neighbors/", "neighbors")):
                if rel.startswith(prefix):
                    return phase
        return "write"
    # doubletfinder; M15: clusters/ (clustering on the fly) and recluster/ (--recluster)
    if rel.startswith("clusters/"):
        return "clustering: write"
    if rel.startswith("recluster/"):
        sub = rel[len("recluster/"):]
        if sub == "summary.json":
            return "recluster"
        if sub.startswith("counts/"):
            return "recluster: load"
        if "w" not in mode or sub == "README.txt":
            return "recluster: manifest"
        return "recluster: write"
    if "w" not in mode:
        return "(tail)"
    if rel.startswith("singlets/"):
        return "singlets"
    return "write"


WCHAR_JUMP = 64 * 1024  # bytes written between two samples that mark output being written
GUARD = 0.03  # s: write signals this close to a stage's start are the previous stage's


def parse_rust(var_dir, tool):
    env = parse_env(os.path.join(var_dir, "env.txt"))
    pm = json.load(open(os.path.join(var_dir, "procmon.json")))
    summ = json.load(open(os.path.join(var_dir, "summary.json")), object_pairs_hook=collections.OrderedDict)
    T, C, R, H, IO, OPEN, SY = load_samples(os.path.join(var_dir, "procmon.samples.tsv.gz"))
    wall = pm["wall_s"]
    cpu_total = pm["rusage"]["cpu_s"]
    peak_total = float(pm["rusage"]["maxrss_kib"])
    cmd = pm["cmd"]
    out_dir = cmd[cmd.index("--out") + 1]
    out_paths = {os.path.abspath(out_dir).rstrip("/"), os.path.realpath(out_dir).rstrip("/")}
    roots = pm["watch"]
    out_ix = {i for i, r in enumerate(roots) if r in out_paths}
    sub = "cluster" if cmd[1] == "cluster" else "doubletfinder"
    late = STAGE_WRITER not in env.get("binary", "")

    # the timeline: the fork (t = 0), the samples, the end (getrusage's CPU and peak; io and RSS
    # as last sampled)
    TT = [0.0] + T + [wall]
    CC = [0.0] + C + [cpu_total]
    SS = [0.0] + SY + [pm["rusage"]["stime_s"]]
    IOs = [(0.0, 0.0, 0.0, 0.0)] + IO + [IO[-1] if IO else (0.0, 0.0, 0.0, 0.0)]
    RR = [0.0] + R + [R[-1] if R else 0.0]
    HH = [0.0] + H + [peak_total]
    n = len(TT)

    def entries(field):
        if field == "-":
            return []
        out = []
        for item in field.split(";"):
            head, _, rest = item.partition(":")
            if not (head.startswith("r") and head[1:].isdigit()):
                continue
            out.append((int(head[1:]), rest.rsplit("@", 1)[0]))
        return out

    out_open = [False] * n
    in_open = [False] * n
    out_seen = []  # (k, rel, mode)
    for k, (fr, fw) in enumerate(OPEN, start=1):
        for mode, field in (("r", fr), ("w", fw)):
            for i, rel in entries(field):
                if i in out_ix:
                    out_open[k] = True
                    out_seen.append((k, rel, mode))
                else:
                    in_open[k] = True
    wjump = [False] * n
    for k in range(1, n - 1):
        if IOs[k][1] - IOs[k - 1][1] >= WCHAR_JUMP:
            wjump[k] = True

    # phases: cumulative ends from summary.json, on the clock of Phases::new (at the start of the
    # command), which starts delta after the fork: delta from the output files, each of which can
    # only be open within the phase that writes (or re-reads) it
    names, ends = [], []
    for name, rec in summ["run"]["phases"].items():
        names.append(name)
        ends.append(float(rec["elapsed_s"]))
    lo_b, hi_b = 0.0, 0.05
    first_in = [k for k in range(n) if in_open[k]]
    if first_in:
        hi_b = min(hi_b, TT[first_in[0]])  # Phases::new comes before the first input is opened
    for k, rel, mode in out_seen:
        if "w" not in mode:
            continue  # re-reads (verification, digests) can fall in more than one phase
        ph = output_phase(sub, rel, mode, late)
        if ph not in names:
            continue
        i = names.index(ph)
        lo_b = max(lo_b, TT[k] - ends[i])
        start_i = ends[i - 1] if i > 0 else 0.0
        hi_b = min(hi_b, TT[k] - start_i)
    delta = min(lo_b + 0.002, hi_b) if lo_b <= hi_b else lo_b

    phases = [("(head)", 0.0, delta)]
    prev = delta
    for name, e in zip(names, ends):
        phases.append((name, prev, e + delta))
        prev = e + delta
    phases.append(("(tail)", prev, wall))

    def kind_of(name):
        return RUST_KIND.get(name, "compute")

    def span(flags, a, b):
        idx = [k for k in range(1, n - 1) if a < TT[k] < b and flags[k]]
        if not idx:
            return None
        k1, k2 = idx[0], idx[-1]
        lo = max(a, (TT[k1 - 1] + TT[k1]) / 2.0)
        hi = min(b, (TT[k2] + TT[k2 + 1]) / 2.0)
        return (lo, hi)

    def measure(a, b):
        cpu = interp(TT, CC, b) - interp(TT, CC, a)
        sys_ = interp(TT, SS, b) - interp(TT, SS, a)
        io = [interp(TT, [x[i] for x in IOs], b) - interp(TT, [x[i] for x in IOs], a) for i in range(4)]
        rss_start = interp(TT, RR, a)
        rss_end = interp(TT, RR, b)
        inside = [RR[k] for k in range(n) if a <= TT[k] <= b]
        peak = max(inside + [rss_start, rss_end])
        # VmHWM: if it rose between the samples around the segment, the new high-water mark was
        # reached in the segment or within one sample of it
        before = [k for k in range(n) if TT[k] <= a]
        after = [k for k in range(n) if TT[k] >= b]
        h_a = HH[before[-1]] if before else 0.0
        h_b = HH[after[0]] if after else HH[-1]
        exact = "sampled"
        if h_b > h_a and h_b > peak:
            peak = h_b
            exact = "VmHWM"
        return cpu, io, rss_start, peak, rss_end, exact, sys_

    # Within a phase (the code's structure, docs in the module string):
    # - sctransform, pca, neighbors (cluster): the stage's files are written after its
    #   computation, until the phase ends: write = from the first write signal (an output file
    #   open, or >= 64 KiB written since the previous sample) at least GUARD after the phase's
    #   start (signals closer to the start belong to the previous phase's last writes) to the end;
    # - load (export: the counts/ copy and its verification come first): write = from the first
    #   to the last write signal;
    # - singlets (doubletfinder: the counts read again, then singlets/ written): read = from the
    #   phase's start to the first write signal (GUARD after the start or later);
    # - the tail (doubletfinder: the inputs read again for their digests, then the outputs'
    #   digests and summary.json): read = from its start to the last sample with an input open;
    # - every other phase: its default kind throughout.
    write_flags = [o or j for o, j in zip(out_open, wjump)]
    only_in = [i and not o for i, o in zip(in_open, out_open)]
    segs = []
    for name, a, b in phases:
        if b <= a:
            continue
        default = kind_of(name)
        cuts = []
        if sub == "cluster" and name in ("sctransform", "pca", "neighbors"):
            w = span(write_flags, a + GUARD, b)
            if w is not None:
                cuts.append((w[0], b, "write", name + ":write"))
        elif name in ("load", "recluster: load"):
            w = span(write_flags, a, b)
            if w is not None:
                cuts.append((w[0], w[1], "write", name + ":write"))
        elif sub == "doubletfinder" and name == "singlets":
            # write_singlets reads the counts again (whole files, decompressed and parsed in
            # memory, so the files are open only briefly), then writes singlets/: read = up to
            # the first write signal
            w = span(write_flags, a + GUARD, b)
            if w is not None and w[0] > a:
                cuts.append((a, w[0], "read", name + ":read"))
        elif sub == "doubletfinder" and name == "(tail)":
            # after the singlets mark: the inputs' digests (read again), then the outputs' shapes
            # and digests and summary.json: read = up to the last sample with an input open
            r_span = span(only_in, a, b)
            if r_span is not None:
                cuts.append((a, r_span[1], "read", name + ":read"))
        cuts.sort()
        t = a
        pieces = []
        for lo, hi, k, label in cuts:
            if lo > t:
                pieces.append((t, lo, default, name))
            pieces.append((max(lo, t), hi, k, label))
            t = max(t, hi)
        if t < b:
            pieces.append((t, b, default, name))
        for lo, hi, k, label in pieces:
            if hi <= lo:
                continue
            cpu, io, rs, pk, re_, exact, sys_ = measure(lo, hi)
            segs.append({
                "phase": label, "kind": k, "t_start_s": lo, "t_end_s": hi, "wall_s": hi - lo,
                "cpu_s": cpu, "sys_s": sys_, "rss_start_kib": rs, "rss_peak_kib": pk, "rss_end_kib": re_,
                "peak_exact": exact, "rchar": io[0], "wchar": io[1], "read_bytes": io[2], "write_bytes": io[3],
            })
    ph = summ["run"]["phases"]
    load_peak = float(ph["load"]["peak_rss_kib"]) if "load" in ph and ph["load"].get("peak_rss_kib") is not None else None
    whole = {"wall_s": wall, "cpu_s": cpu_total, "peak_kib": peak_total,
             "rusage_maxrss_kib": peak_total, "sampler_cpu_s": pm["sampler_cpu_s"],
             "exit_code": pm["exit_code"], "load_peak_kib": load_peak,
             "samples": pm["samples"], "delta_s": delta, "late_writer": late}
    inputs = sum(e["bytes"] for e in pm.get("evicted", []))
    return env, pm, segs, whole, inputs


# ---------------------------------------------------------------------------------------------
# the pca writer's calibration (pca_calib.sbatch): `scrust pca`'s write phase, the same function
# as the chain's pca/ writer, which formats each file before opening it


def load_pca_calibration():
    """object -> (cells, median seconds of the write phase); and a linear fit seconds ~ cells."""
    out = {}
    for d in sorted(glob.glob(os.path.join(B, "runs", "calib", "pca", "*", "rep*"))):
        obj = d.split(os.sep)[-2]
        sj = os.path.join(d, "summary.json")
        if not os.path.exists(sj):
            continue
        summ = json.load(open(sj), object_pairs_hook=collections.OrderedDict)
        ph = summ["run"]["phases"]
        names = list(ph.keys())
        i = names.index("write")
        w = float(ph["write"]["elapsed_s"]) - float(ph[names[i - 1]]["elapsed_s"])
        cells = summ["result"]["n_cells"]
        out.setdefault(obj, [cells, []])[1].append(w)
    cal = {o: (c, statistics.median(v)) for o, (c, v) in out.items()}
    fit = None
    if len(cal) >= 2:
        xs = [c for c, _ in cal.values()]
        ys = [w for _, w in cal.values()]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        sxx = sum((x - mx) ** 2 for x in xs)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx if sxx else 0.0
        fit = (my - b * mx, b)
    return cal, fit


PCA_CAL, PCA_FIT = None, None


def pca_writer_seconds(obj, cells):
    global PCA_CAL, PCA_FIT
    if PCA_CAL is None:
        PCA_CAL, PCA_FIT = load_pca_calibration()
    if obj in PCA_CAL and PCA_CAL[obj][0] == cells:
        return PCA_CAL[obj][1]
    if PCA_FIT is not None:
        return PCA_FIT[0] + PCA_FIT[1] * cells
    return None


def adjust_pca_write(segs, obj, cells):
    """Moves the pca writer's unseen formatting time (calibrated, minus what the sampler saw as
    writing) from the pca stage's compute segment to its write segment. Returns the seconds moved."""
    target = pca_writer_seconds(obj, cells)
    if target is None:
        return 0.0
    comp = [x for x in segs if x["phase"] == "pca" and x["kind"] == "compute"]
    wr = [x for x in segs if x["phase"] == "pca:write"]
    if not comp:
        return 0.0
    seen = sum(x["cpu_s"] for x in wr)
    move = max(0.0, min(target - seen, comp[-1]["cpu_s"]))
    if move <= 0:
        return 0.0
    c = comp[-1]
    c["cpu_s"] -= move
    c["wall_s"] -= move
    c["t_end_s"] -= move
    if wr:
        w = wr[0]
        w["cpu_s"] += move
        w["wall_s"] += move
        w["t_start_s"] -= move
        w["peak_exact"] += "; +%.3f s calibrated" % move
    else:
        i = segs.index(c)
        segs.insert(i + 1, {"phase": "pca:write", "kind": "write", "t_start_s": c["t_end_s"],
                            "t_end_s": c["t_end_s"] + move, "wall_s": move, "cpu_s": move, "sys_s": 0.0,
                            "rss_start_kib": c["rss_end_kib"], "rss_peak_kib": c["rss_end_kib"],
                            "rss_end_kib": c["rss_end_kib"], "peak_exact": "calibrated",
                            "rchar": 0.0, "wchar": 0.0, "read_bytes": 0.0, "write_bytes": 0.0})
    return move


# ---------------------------------------------------------------------------------------------
# per-run summaries by kind


def kind_summary(segs, whole, side, included):
    """CPU, wall by kind; the memory shares. included: kinds counted in the like-for-like total."""
    cpu = collections.OrderedDict((k, 0.0) for k in KINDS)
    wall = collections.OrderedDict((k, 0.0) for k in KINDS)
    for s in segs:
        k = s["kind"]
        if k not in cpu:
            cpu[k] = wall[k] = 0.0
        cpu[k] += s["cpu_s"]
        wall[k] += s["wall_s"]
    tot_cpu = sum(v for k, v in cpu.items() if k in included)
    tot_wall = sum(v for k, v in wall.items() if k in included)
    peaks = [s["rss_peak_kib"] for s in segs if s["kind"] in included and s["rss_peak_kib"] == s["rss_peak_kib"]]
    peak_incl = max(peaks) if peaks else float("nan")
    # memory due to reading: the peak reached by the end of the read phase (R: the read phases,
    # read, read_raw, read_base_rds, which follow each other; Rust: taken from its own record below)
    mem_read = 0.0
    seen_read = False
    for s in segs:
        if seen_read and s["kind"] != "read":
            break
        if s["kind"] == "read":
            seen_read = True
        if s["rss_peak_kib"] == s["rss_peak_kib"]:
            mem_read = max(mem_read, s["rss_peak_kib"])
    # memory added by writing: per contiguous write stretch, its peak minus the RSS when it began
    mem_write = 0.0
    run_start = None
    run_peak = 0.0
    for s in segs + [{"kind": "(end)"}]:
        if s["kind"] == "write":
            if run_start is None:
                run_start = s["rss_start_kib"]
                run_peak = 0.0
            if s["rss_peak_kib"] == s["rss_peak_kib"]:
                run_peak = max(run_peak, s["rss_peak_kib"])
        else:
            if run_start is not None:
                mem_write = max(mem_write, run_peak - run_start)
            run_start = None
    if whole.get("load_peak_kib") is not None:
        # Rust: VmHWM at the end of load, from the run's own record (exact)
        mem_read = whole["load_peak_kib"]
    out = {
        "cpu": cpu, "wall": wall, "cpu_total": tot_cpu, "wall_total": tot_wall,
        "peak_kib": peak_incl, "peak_whole_kib": whole["peak_kib"],
        "mem_read_kib": mem_read, "mem_write_kib": mem_write,
    }
    for k in ("rchar", "wchar", "read_bytes", "write_bytes"):
        out[k] = sum(s[k] for s in segs if s[k] == s[k])
    out["sys_write"] = sum(s.get("sys_s", 0.0) for s in segs if s["kind"] == "write")
    out["sys_total"] = sum(s.get("sys_s", 0.0) for s in segs if s["kind"] in included)
    return out


# ---------------------------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-pilot", action="store_true")
    ap.add_argument("--cpu-model", default="8480CL",
                    help="runs on other CPU models go to benchmark_other_cpu.tsv, not into the tables")
    a = ap.parse_args()
    min_rep = 0 if a.include_pilot else 1
    other = collections.defaultdict(list)

    runs_rows = []
    kinds_rows = []
    per = collections.defaultdict(list)  # (tool, object, side, variant) -> [summary]
    parsed = {}  # Rust: (tool, object, rep, variant) -> (meta, segs, whole)

    def add_rows(meta, segs, whole, ks):
        for s in segs:
            row = collections.OrderedDict(meta)
            row.update(s)
            runs_rows.append(row)
        row = collections.OrderedDict(meta)
        row.update({"phase": "(whole run)", "kind": "whole", "t_start_s": 0.0, "t_end_s": whole["wall_s"],
                    "wall_s": whole["wall_s"], "cpu_s": whole["cpu_s"], "sys_s": ks["sys_total"], "rss_start_kib": 0.0,
                    "rss_peak_kib": whole["peak_kib"], "rss_end_kib": float("nan"), "peak_exact": "yes",
                    "rchar": ks["rchar"], "wchar": ks["wchar"], "read_bytes": ks["read_bytes"],
                    "write_bytes": ks["write_bytes"]})
        runs_rows.append(row)

    # R
    for run_dir in sorted(glob.glob(os.path.join(B, "runs", "r", "*", "*", "rep*"))):
        if ".partial." in run_dir:
            continue
        tool, obj, rep = run_dir.split(os.sep)[-3:]
        repn = int(rep[3:])
        if repn < min_rep:
            continue
        env, pm, segs, whole, inputs = parse_r(run_dir)
        included = {"startup", "read", "compute", "write", "exit"}
        if tool in DF_TOOLS:
            included = included | {"recluster"}
        ks = kind_summary(segs, whole, "R", included)
        meta = collections.OrderedDict([
            ("side", "R"), ("tool", tool), ("object", obj), ("variant", "workflow script"),
            ("rep", repn), ("threads", 1), ("job", env.get("job", "")), ("host", env.get("host", "")),
            ("cpu_model", env.get("cpu", "")), ("input_bytes", inputs),
            ("start_epoch_s", pm["start_epoch_s"]), ("end_epoch_s", pm["start_epoch_s"] + pm["wall_s"])])
        meta["included"] = "yes" if a.cpu_model in meta["cpu_model"] else "no (CPU model)"
        add_rows(meta, segs, whole, ks)
        target = per if meta["included"] == "yes" else other
        target[(tool, obj, "R", "R")].append((meta, ks, whole))
        if tool == "doubletfinder":
            # DoubletFinder alone: R's startup, readRDS and its DoubletFinder phases (paramSweep ..
            # the subset to the singlets); R writes nothing in this part (its outputs are saved
            # after the re-clustering)
            segs_a = [x for x in segs if x["kind"] in ("startup", "read", "compute")]
            whole_a = {"wall_s": sum(x["wall_s"] for x in segs_a), "cpu_s": sum(x["cpu_s"] for x in segs_a),
                       "peak_kib": max(x["rss_peak_kib"] for x in segs_a)}
            ks_a = kind_summary(segs_a, whole_a, "R", {"startup", "read", "compute"})
            m_a = collections.OrderedDict(meta)
            m_a["variant"] = "DoubletFinder alone"
            target[(tool, obj, "R", "R-alone")].append((m_a, ks_a, whole_a))

    # Rust
    for var_dir in sorted(glob.glob(os.path.join(B, "runs", "rust", "*", "*", "rep*", "*"))):
        if ".partial." in var_dir or not os.path.isdir(var_dir):
            continue
        tool, obj, rep, variant = var_dir.split(os.sep)[-4:]
        repn = int(rep[3:])
        if repn < min_rep:
            continue
        if not os.path.exists(os.path.join(var_dir, "summary.json")):
            print("skipped (no summary.json): %s" % var_dir, file=sys.stderr)
            continue
        env, pm, segs, whole, inputs = parse_rust(var_dir, tool)
        threads = int(variant[1:])
        summ = json.load(open(os.path.join(var_dir, "summary.json")))
        whole["pca_write_moved_s"] = 0.0
        # the pca writer's calibration is for the binary that wrote pca/ inside the pca phase
        if (summ.get("tool") == "cluster" or pm["cmd"][1] == "cluster") and not whole["late_writer"]:
            whole["pca_write_moved_s"] = adjust_pca_write(segs, obj, summ["result"].get("n_cells") if isinstance(summ.get("result"), dict) else None)
        ks = kind_summary(segs, whole, "Rust", {"startup", "read", "compute", "write", "exit"})
        meta = collections.OrderedDict([
            ("side", "Rust"), ("tool", tool), ("object", obj), ("variant", variant),
            ("rep", repn), ("threads", threads), ("job", env.get("job", "")), ("host", env.get("host", "")),
            ("cpu_model", env.get("cpu", "")), ("input_bytes", inputs),
            ("start_epoch_s", pm["start_epoch_s"]), ("end_epoch_s", pm["start_epoch_s"] + pm["wall_s"])])
        meta["included"] = "yes" if a.cpu_model in meta["cpu_model"] else "no (CPU model)"
        add_rows(meta, segs, whole, ks)
        if meta["included"] != "yes":
            other[(tool, obj, "Rust", variant)].append((meta, ks, whole))
            continue
        per[(tool, obj, "Rust", variant)].append((meta, ks, whole))
        parsed[(tool, obj, repn, variant)] = (meta, segs, whole)

    # doubletfinder in Rust: `scrust doubletfinder` then `scrust cluster` on its singlets/ (the
    # workflow's doubletfinder.R clusters the singlets after DoubletFinder), as one pipeline
    for (tool, obj, repn, variant), (meta, segs, whole) in sorted(parsed.items()):
        if tool not in DF_TOOLS or not variant.startswith("d"):
            continue
        t = variant[1:]
        rrun = parsed.get((tool, obj, repn, "r" + t))
        if rrun is None:
            continue
        meta2, segs2, whole2 = rrun
        shift = whole["wall_s"]
        segs_c = [dict(x) for x in segs] + [dict(x, t_start_s=x["t_start_s"] + shift, t_end_s=x["t_end_s"] + shift,
                                                 phase="cluster:" + x["phase"]) for x in segs2]
        whole_c = {"wall_s": whole["wall_s"] + whole2["wall_s"], "cpu_s": whole["cpu_s"] + whole2["cpu_s"],
                   "peak_kib": max(whole["peak_kib"], whole2["peak_kib"]),
                   "load_peak_kib": max(whole["load_peak_kib"] or 0.0, whole2["load_peak_kib"] or 0.0)}
        ks = kind_summary(segs_c, whole_c, "Rust", {"startup", "read", "compute", "write", "exit"})
        m = collections.OrderedDict(meta)
        m["variant"] = "d%s+r%s" % (t, t)
        m["job"] = meta["job"]
        per[(tool, obj, "Rust", m["variant"])].append((m, ks, whole_c))

    cols = ["side", "tool", "object", "variant", "rep", "threads", "job", "host", "cpu_model", "included", "input_bytes",
            "phase", "kind", "t_start_s", "t_end_s", "wall_s", "cpu_s", "sys_s", "rss_start_kib", "rss_peak_kib",
            "rss_end_kib", "peak_exact", "rchar", "wchar", "read_bytes", "write_bytes"]
    with open(os.path.join(OUT, "benchmark_runs.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in runs_rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")

    # co-location: the other measured runs of this benchmark on the same node, overlapping in time
    allruns = [(m, ks, w) for tab in (per, other) for lst in tab.values() for (m, ks, w) in lst
               if "+" not in m["variant"] and m["variant"] != "DoubletFinder alone"]
    for lst in list(per.values()) + list(other.values()):
        for m, ks, w in lst:
            ids = set()
            for m2, _, _ in allruns:
                if m2 is m or m2["host"] != m["host"] or m2["job"] == m["job"]:
                    continue
                if m2["start_epoch_s"] < m["end_epoch_s"] and m2["end_epoch_s"] > m["start_epoch_s"]:
                    ids.add(m2["job"])
            m["concurrent_bench_jobs"] = len(ids)

    # per-run kinds
    kcols = ["side", "tool", "object", "variant", "rep", "threads", "job", "host", "cpu_model", "input_bytes",
             "concurrent_bench_jobs", "wall_whole_s", "cpu_whole_s", "peak_whole_gib",
             "cpu_total_s", "wall_total_s", "peak_gib"] + \
            ["cpu_%s_s" % k for k in KINDS] + ["wall_%s_s" % k for k in KINDS] + \
            ["cpu_read_pct", "cpu_write_pct", "cpu_io_pct", "mem_read_gib", "mem_write_gib",
             "mem_read_pct", "mem_write_pct", "sys_total_s", "sys_write_s", "rchar", "wchar", "read_bytes", "write_bytes"]
    for name, table in (("benchmark_kinds.tsv", per), ("benchmark_other_cpu.tsv", other)):
        with open(os.path.join(OUT, name), "w") as fh:
            fh.write("\t".join(kcols) + "\n")
            for key in sorted(table):
                for meta, ks, whole in sorted(table[key], key=lambda x: x[0]["rep"]):
                    row = dict(meta)
                    row.update(derived(ks, whole))
                    fh.write("\t".join(fmt(row.get(c, "")) for c in kcols) + "\n")
    print("wrote benchmark_runs.tsv (%d rows), benchmark_kinds.tsv (%d runs)" % (
        len(runs_rows), sum(len(v) for v in per.values())))
    write_table(per)


METRICS = ["cpu_total_s", "wall_total_s", "peak_gib", "cpu_read_pct", "cpu_write_pct", "cpu_io_pct",
           "mem_read_pct", "mem_write_pct", "cpu_read_s", "cpu_write_s", "mem_read_gib", "mem_write_gib"]
R_EXTRA = ["cpu_startup_s", "cpu_unported_s", "cpu_whole_s", "wall_whole_s", "peak_whole_gib", "cpu_compute_s",
           "cpu_recluster_s"]


def input_info(obj):
    d = json.load(open(os.path.join(S, "data", "clustering", obj, "input", "export.json")))
    return d["n_cells"], d["n_features"], d["nnz"]


def stats(vals):
    vals = [v for v in vals if v == v]
    if not vals:
        return (float("nan"),) * 3 + (0,)
    return statistics.median(vals), min(vals), max(vals), len(vals)


def write_table(per):
    rows = []
    for tool in ("cluster", "doubletfinder", DF_M15):
        objs = OBJECTS if tool == "cluster" else DF_OBJECTS
        variants = {
            "cluster": [("without --export", "n"), ("with --export", "x")],
            "doubletfinder": [("DoubletFinder alone", "alone"), ("whole rule", "rule")],
            # M15 (one binary for all four): the hook's cost (h against alone) and the in-process
            # hand-off's (inprocess against rule)
            DF_M15: [("DoubletFinder alone", "alone"), ("DoubletFinder, clustering its own input", "hook"),
                     ("whole rule", "rule"), ("whole rule, in process", "inprocess")],
        }[tool]
        for obj in objs:
            if tool == DF_M15 and not any(k[0] == DF_M15 and k[1] == obj for k in per):
                continue
            for label, v in variants:
                # R's DoubletFinder runs serve both DoubletFinder tools (no new R runs for M15)
                r_tool = "doubletfinder" if tool == DF_M15 else tool
                r_runs = per.get((r_tool, obj, "R", "R-alone" if v in ("alone", "hook") else "R"), [])
                v1, v32 = {"n": ("n1", "n32"), "x": ("x1", "x32"), "alone": ("d1", "d32"), "rule": ("d1+r1", "d32+r32"),
                           "hook": ("h1", "h32"), "inprocess": ("p1", "p32")}[v]
                units = {"R": r_runs, "Rust1": per.get((tool, obj, "Rust", v1), []),
                         "Rust32": per.get((tool, obj, "Rust", v32), [])}
                cells, genes, nnz = input_info(obj)
                row = collections.OrderedDict([("tool", tool), ("object", obj), ("rust_variant", label),
                                               ("cells", cells), ("genes", genes), ("nonzeros", nnz)])
                for side, lst in units.items():
                    ds = [derived(ks, whole) for _, ks, whole in lst]
                    row["%s_reps" % side] = len(ds)
                    row["%s_input_bytes" % side] = lst[0][0]["input_bytes"] if lst else ""
                    row["%s_nodes" % side] = ",".join(sorted(set(m["host"] for m, _, _ in lst)))
                    row["%s_cpus" % side] = ",".join(sorted(set(m["cpu_model"].replace("Intel(R)_Xeon(R)_Platinum_", "") for m, _, _ in lst)))
                    for k in METRICS + (R_EXTRA if side == "R" else []):
                        med, lo, hi, nn = stats([d.get(k, float("nan")) for d in ds])
                        row["%s_%s" % (side, k)] = med
                        row["%s_%s_min" % (side, k)] = lo
                        row["%s_%s_max" % (side, k)] = hi
                for t in ("1", "32"):
                    rc, uc = row.get("R_cpu_total_s"), row.get("Rust%s_cpu_total_s" % t)
                    rp, up = row.get("R_peak_gib"), row.get("Rust%s_peak_gib" % t)
                    rw, uw = row.get("R_wall_total_s"), row.get("Rust%s_wall_total_s" % t)
                    row["cpu_ratio_R_over_Rust%s" % t] = rc / uc if rc == rc and uc and uc == uc else float("nan")
                    row["wall_ratio_R_over_Rust%s" % t] = rw / uw if rw == rw and uw and uw == uw else float("nan")
                    row["peak_ratio_R_over_Rust%s" % t] = rp / up if rp == rp and up and up == up else float("nan")
                rows.append(row)
    cols = list(rows[0].keys())
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(os.path.join(OUT, "benchmark_table.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")
    print("wrote benchmark_table.tsv (%d rows)" % len(rows))


def derived(ks, whole):
    g = 1024.0 * 1024.0
    d = {
        "wall_whole_s": whole["wall_s"], "cpu_whole_s": whole["cpu_s"], "peak_whole_gib": whole["peak_kib"] / g,
        "cpu_total_s": ks["cpu_total"], "wall_total_s": ks["wall_total"], "peak_gib": ks["peak_kib"] / g,
        "mem_read_gib": ks["mem_read_kib"] / g, "mem_write_gib": ks["mem_write_kib"] / g,
        "rchar": ks["rchar"], "wchar": ks["wchar"], "read_bytes": ks["read_bytes"], "write_bytes": ks["write_bytes"],
    }
    for k in KINDS:
        d["cpu_%s_s" % k] = ks["cpu"].get(k, 0.0)
        d["wall_%s_s" % k] = ks["wall"].get(k, 0.0)
    d["cpu_whole_s"] = whole["cpu_s"]
    d["sys_write_s"] = ks["sys_write"]
    d["sys_total_s"] = ks["sys_total"]
    tot = ks["cpu_total"]
    d["cpu_read_pct"] = 100.0 * ks["cpu"]["read"] / tot if tot else float("nan")
    d["cpu_write_pct"] = 100.0 * ks["cpu"]["write"] / tot if tot else float("nan")
    d["cpu_io_pct"] = d["cpu_read_pct"] + d["cpu_write_pct"]
    d["mem_read_pct"] = 100.0 * ks["mem_read_kib"] / ks["peak_kib"] if ks["peak_kib"] else float("nan")
    d["mem_write_pct"] = 100.0 * ks["mem_write_kib"] / ks["peak_kib"] if ks["peak_kib"] else float("nan")
    return d


def fmt(v):
    if isinstance(v, float):
        if v != v:
            return "NA"
        if abs(v) >= 1e6:
            return "%.0f" % v
        return "%.4f" % v
    return str(v)


if __name__ == "__main__":
    main()
