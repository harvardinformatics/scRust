#!/bin/bash
# Milestone M14: the tables of docs/validation/clustering.md ("M14") from the loader checks' plain
# files (no R, no Rust). From the repository's top directory:
#   bash validation/slurm/m14_report.sh <label> > docs/validation/clustering_m14.tsv          (K11)
#   bash validation/slurm/m14_report.sh <label> k13 > docs/validation/clustering_m14_k13.tsv
#   bash validation/slurm/m14_report.sh <label> m17 > docs/validation/clustering_m14_m17.tsv
#   bash validation/slurm/m14_report.sh <label> components    (every component's status, summed)
# K11 (m14_loaders.sbatch): one row per object, profile and loader (seurat, sce, anndata): the
# grading mode, pass, the components by status, the files checked (values, differing), the round
# trip, the step's wall time and peak RSS, and for the R loaders the 10x readers on the export.
# K13 (m14_markers.sbatch): one row per object and profile: clusters, rows, identical, time, RSS.
# M17 (m14_m17.sbatch): one row per object: the comparison's time and peak RSS, by step.
set -euo pipefail
label=${1:?label}
what=${2:-k11}
python3 -I - "$label" "$what" <<'PY'
import csv, glob, json, os, re, sys
label, what = sys.argv[1], sys.argv[2]
tiers = {}
for f in ("validation/clustering_datasets.tsv", "validation/doubletfinder_datasets.tsv"):
    for row in csv.DictReader(open(f), delimiter="\t"):
        tiers[row["object"]] = row["tier"]

def time_v(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        if "Elapsed (wall clock)" in line:
            t = line.rsplit(" ", 1)[1].strip().split(":")
            s = 0.0
            for x in t:
                s = s * 60 + float(x)
            out["wall_s"] = s
        elif "Maximum resident set size" in line:
            out["peak_gib"] = int(line.rsplit(" ", 1)[1]) / 2**20
    return out

def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []

def guard(d):
    sr = tsv(os.path.join(d, "read10x_seurat.tsv"))
    sc = tsv(os.path.join(d, "read10x_scanpy.tsv"))
    attempts = 0
    for r in sc:
        m = re.match(r"all (\d+) prefixes", r["prefix"])
        if m:
            attempts += int(m.group(1))
    return {"read10x_dirs": len(sr), "read10x_as_required": sum(r["as_required"] == "TRUE" for r in sr),
            "read_10x_mtx_attempts": attempts, "read_10x_mtx_as_required": all(r["as_required"] == "True" for r in sc) if sc else ""}

STATUSES = ["identical", "within floor", "beyond floor", "differing", "not rebuilt", "added"]
if what == "k11":
    cols = ["object", "tier", "kind", "profile", "loader", "mode", "pass", "components"] + [s.replace(" ", "_") for s in STATUSES] + \
           ["files_values", "files_differing", "round_trip", "wall_s", "peak_gib", "read10x_dirs", "read10x_as_required",
            "read_10x_mtx_attempts", "read_10x_mtx_as_required"]
    print("\t".join(cols))
    dirs = sorted(glob.glob(f"data/clustering/*/loaders/{label}/*/r")) + sorted(glob.glob(f"data/doubletfinder/*/loaders/{label}/*/r"))
    for rdir in dirs:
        parts = rdir.split("/")
        obj, profile = parts[2], parts[5]
        base = os.path.dirname(rdir)
        kind = "doubletfinder" if parts[1] == "doubletfinder" else "cluster"
        name = obj if kind == "cluster" else f"doubletfinder_{obj}"
        g = guard(os.path.join(base, "guard")) if kind == "cluster" else {}
        cj = os.path.join(rdir, "compare.json")
        k11 = tsv(os.path.join(rdir, "k11.tsv"))
        rt = {r["loader"]: r for r in tsv(os.path.join(rdir, "roundtrip.tsv"))}
        files = tsv(os.path.join(rdir, "files.tsv"))
        tv = time_v(os.path.join(rdir, "time_v.txt"))
        c = json.load(open(cj)) if os.path.exists(cj) else {}
        for loader in ("seurat", "sce"):
            rows = [r for r in k11 if r["loader"] == loader]
            if not rows:
                continue
            counts = {s: sum(r["status"] == s for r in rows) for s in STATUSES}
            fv = sum(int(r["values"]) for r in files) if loader == "seurat" and files else ""
            fd = sum(int(r["differing"]) for r in files) if loader == "seurat" and files else ""
            r = rt.get(loader)
            rtx = "" if r is None else ("identical" if r["identical_to_loader"] == "TRUE" and r["identical_to_rebuilt"] == "TRUE" else "DIFFERS")
            ok = counts["beyond floor"] == 0 and counts["differing"] == 0 and (fd in ("", 0)) and rtx != "DIFFERS"
            vals = [name, tiers.get(obj, ""), kind, profile, loader, c.get("mode", ""), str(ok).upper(), str(len(rows))] + \
                   [str(counts[s]) for s in STATUSES] + [str(fv), str(fd), rtx, f"{tv.get('wall_s', 0):.0f}", f"{tv.get('peak_gib', 0):.2f}"] + \
                   [str(g.get(k, "")) for k in ("read10x_dirs", "read10x_as_required", "read_10x_mtx_attempts", "read_10x_mtx_as_required")]
            print("\t".join(vals))
        pdir = os.path.join(base, "python")
        pj = os.path.join(pdir, "compare_anndata.json")
        if os.path.exists(pj):
            p = json.load(open(pj))
            rows = tsv(os.path.join(pdir, "k11_anndata.tsv"))
            counts = {s: sum(r["status"] == s for r in rows) for s in STATUSES}
            files = tsv(os.path.join(pdir, "files_anndata.tsv"))
            fv = sum(int(r["values"]) for r in files)
            fd = sum(int(r["differing"]) for r in files)
            r = tsv(os.path.join(pdir, "roundtrip_anndata.tsv"))
            rtx = "" if not r else ("identical" if r[0]["identical_to_loader"] == "True" and r[0]["identical_to_rebuilt"] == "True" else "DIFFERS")
            tv = time_v(os.path.join(pdir, "time_v.txt"))
            other = len(rows) - sum(counts.values())
            vals = [name, tiers.get(obj, ""), kind, profile, "anndata", p.get("mode", ""), str(bool(p.get("pass"))).upper(), str(len(rows))] + \
                   [str(counts[s]) for s in STATUSES] + [str(fv), str(fd), rtx, f"{tv.get('wall_s', 0):.0f}", f"{tv.get('peak_gib', 0):.2f}"] + \
                   ["", "", "", ""]
            if other:
                vals[STATUSES.index("not rebuilt") + 8] += f" (+{other} reported)"
            print("\t".join(vals))
elif what == "components":
    # every component's status over all objects, profiles and loaders
    from collections import Counter
    agg = Counter()
    for f in glob.glob(f"data/*/*/loaders/{label}/*/r/k11.tsv"):
        for r in tsv(f):
            agg[(r["loader"], r["component"], r["status"])] += 1
    for f in glob.glob(f"data/*/*/loaders/{label}/*/python/k11_anndata.tsv"):
        for r in tsv(f):
            agg[("anndata", r["component"], r["status"])] += 1
    print("loader\tcomponent\tstatus\tcount")
    for (l, comp, st), n in sorted(agg.items()):
        print(f"{l}\t{comp}\t{st}\t{n}")
elif what == "floor":
    # the components graded within the floor (portable profile), per loader: objects, the largest
    # statistics, and the floor's bounds
    import tomllib
    floor = tomllib.load(open("validation/tolerances.toml", "rb"))["clustering"]["floor"]["kernel"]
    agg = {}
    files = [(f, None) for f in glob.glob(f"data/*/*/loaders/{label}/*/r/k11.tsv")] + \
            [(f, "anndata") for f in glob.glob(f"data/*/*/loaders/{label}/*/python/k11_anndata.tsv")]
    for f, loader in files:
        tier = tiers.get(f.split("/")[2], "")
        for r in tsv(f):
            if r["status"] not in ("within floor", "beyond floor"):
                continue
            key = (loader or r["loader"], r["component"], r["floor"], tier)
            a = agg.setdefault(key, {"objects": 0, "beyond": 0, "differing": 0, "values": 0, "max_abs": 0.0, "max_scaled": 0.0})
            a["objects"] += 1
            a["beyond"] += r["status"] == "beyond floor"
            for k in ("differing", "n"):
                if r.get(k) not in (None, ""):
                    a["values" if k == "n" else k] += int(float(r[k]))
            for k in ("max_abs", "max_scaled"):
                if r.get(k) not in (None, ""):
                    a[k] = max(a[k], float(r[k]))
    print("loader\tcomponent\ttier\tfloor\tobjects\tbeyond\tvalues_differing\tvalues\tmax_abs\tmax_scaled\tfloor_bounds")
    for (l, comp, fl, tier), a in sorted(agg.items()):
        bounds = ", ".join(f"{k} {v}" for k, v in floor.get(fl, {}).items() if isinstance(v, (int, float)) and not isinstance(v, bool))
        print(f"{l}\t{comp}\t{tier}\t{fl}\t{a['objects']}\t{a['beyond']}\t{a['differing']}\t{a['values']}\t{a['max_abs']:.3g}\t{a['max_scaled']:.3g}\t{bounds}")
elif what == "notrebuilt":
    from collections import defaultdict
    det = defaultdict(set)
    for f in glob.glob(f"data/*/*/loaders/{label}/*/r/k11.tsv") + glob.glob(f"data/*/*/loaders/{label}/*/python/k11_anndata.tsv"):
        obj = f.split("/")[2]
        for r in tsv(f):
            if r["status"] in ("not rebuilt", "added") or r["status"].startswith("digest"):
                det[(r["loader"], r["component"], r["status"], r["detail"])].add(obj)
    print("loader\tcomponent\tstatus\tobjects\tdetail")
    for (l, comp, st, d), objs in sorted(det.items()):
        print(f"{l}\t{comp}\t{st}\t{len(objs)}\t{d}")
elif what == "k13":
    print("\t".join(["object", "tier", "profile", "clusters", "rows", "identical", "wall_s", "load_s", "peak_gib"]))
    for f in sorted(glob.glob(f"data/clustering/*/loaders/{label}/*/markers/markers.json")):
        parts = f.split("/")
        obj, profile = parts[2], parts[5]
        j = json.load(open(f))
        tv = time_v(os.path.join(os.path.dirname(f), "time_v.txt"))
        print("\t".join([obj, tiers.get(obj, ""), profile, str(len(j["cluster_order_used"])), str(j["rows"]),
                         str(j["comparison"]["identical"]).upper(), f"{tv.get('wall_s', 0):.0f}",
                         f"{j['timing_s']['load_seurat']:.0f}", f"{tv.get('peak_gib', 0):.2f}"]))
elif what == "m17":
    print("\t".join(["object", "tier", "profile", "cells", "features", "pass", "components", "identical", "not_rebuilt",
                     "wall_s", "peak_gib", "manifest_s", "load_seurat_s", "read_rds_s", "compare_s"]))
    for f in sorted(glob.glob(f"data/clustering/*/loaders/{label}/*/m17/compare.json")):
        parts = f.split("/")
        obj, profile = parts[2], parts[5]
        j = json.load(open(f))
        tv = time_v(os.path.join(os.path.dirname(f), "time_v.txt"))
        st = j["steps"]
        prev = 0.0
        times = {}
        for k in ("manifest", "load_seurat", "read_rds", "compare"):
            times[k] = st[k]["elapsed_s"] - prev
            prev = st[k]["elapsed_s"]
        summ = json.load(open(f"data/clustering/{obj}/rust/{label}/{profile}/run/summary.json"))
        cells = summ["result"]["n_cells"]
        feats = summ["result"]["sctransform"]["n_features"]
        print("\t".join([obj, tiers.get(obj, ""), profile, str(cells), str(feats), str(bool(j["pass"])).upper(), str(j["components"]),
                         str(j["status"].get("identical", 0)), str(j["status"].get("not rebuilt", 0)),
                         f"{tv.get('wall_s', 0):.0f}", f"{tv.get('peak_gib', 0):.2f}"] + [f"{times[k]:.0f}" for k in times]))
PY
