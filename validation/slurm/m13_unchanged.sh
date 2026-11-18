#!/bin/bash
# Milestone M13: everything M11 and M12 checked keeps its status. Compares every `scrust validate
# cluster` report of an earlier label (M12's, m12-4231635) with the same comparison's report of
# M13's label, check by check (stage and name): the status, the values compared and the values
# differing must be the same. The earlier reports' placeholders for M13 (stage "later": the K10
# arrays, nCount_SCT and nFeature_SCT, "not compared") are listed apart: M13 compares them (stage
# K10, and the saved object's). Reads plain files only. From the repository's top directory:
#   bash validation/slurm/m13_unchanged.sh <old label> <new label> [regrade] > data/clustering/compare/m13_unchanged_<new label>.tsv
# The new reports are the new label's own runs' (m13_cluster.sbatch), or with "regrade" the old
# runs re-graded by the new label's binary (m13_regrade.sbatch: <old run dir>/<new label>/).
# One row per comparison: object, profile, reference, the old report's checks, those unchanged,
# changed, gone (absent from the new report: only the placeholders may be), the new report's new
# checks by stage. Exit status 1 if any check changed or went missing.
set -euo pipefail
old=${1:?old label}
new=${2:?new label}
mode=${3:-runs}
python3 -I - "$old" "$new" "$mode" <<'EOF'
import glob, json, os, sys
old, new, mode = sys.argv[1], sys.argv[2], sys.argv[3]
print("object\tprofile\treference\told_checks\tunchanged\tchanged\tgone_placeholders\tgone_other\tnew_checks_by_stage")
bad = 0
tot = {"old": 0, "same": 0, "changed": 0, "ph": 0, "gone": 0}
for f in sorted(glob.glob(f"data/clustering/*/rust/{old}/*/vs_*.json")):
    parts = f.split("/")
    obj, profile, name = parts[2], parts[5], os.path.basename(f)
    g = f"data/clustering/{obj}/rust/{old}/{profile}/{new}/{name}" if mode == "regrade" else f"data/clustering/{obj}/rust/{new}/{profile}/{name}"
    if not os.path.exists(g):
        continue
    a, b = json.load(open(f)), json.load(open(g))
    key = lambda c: (c["stage"], c["name"])
    bk = {key(c): c for c in b["checks"]}
    same = changed = ph = gone = 0
    for c in a["checks"]:
        if c["stage"] == "environment":
            # the environment line names the kernels and the mode: compared apart (the mode)
            same += 1
            continue
        d = bk.get(key(c))
        if d is None:
            if c["stage"] == "later":
                ph += 1
            else:
                gone += 1
                print(f"# {obj} {profile} {name}: {key(c)} gone", file=sys.stderr)
            continue
        if (c["status"], c["n"], c["n_differ"]) == (d["status"], d["n"], d["n_differ"]):
            same += 1
        else:
            changed += 1
            print(f"# {obj} {profile} {name}: {key(c)}: {c['status']} {c['n']} {c['n_differ']} -> {d['status']} {d['n']} {d['n_differ']}", file=sys.stderr)
    ak = {key(c) for c in a["checks"]}
    news = {}
    for c in b["checks"]:
        if key(c) not in ak:
            news[c["stage"]] = news.get(c["stage"], 0) + 1
    if a["grading"]["mode"] != b["grading"]["mode"]:
        changed += 1
        print(f"# {obj} {profile} {name}: grading mode {a['grading']['mode']} -> {b['grading']['mode']}", file=sys.stderr)
    bad += changed + gone
    for k, v in (("old", len(a["checks"])), ("same", same), ("changed", changed), ("ph", ph), ("gone", gone)):
        tot[k] += v
    ref = name[3:-5]
    print(f"{obj}\t{profile}\t{ref}\t{len(a['checks'])}\t{same}\t{changed}\t{ph}\t{gone}\t" + ",".join(f"{k} {v}" for k, v in sorted(news.items())))
print(f"# total: {tot['old']} checks of {old}'s reports; {tot['same']} unchanged, {tot['changed']} changed, {tot['ph']} placeholders compared now, {tot['gone']} gone", file=sys.stderr)
sys.exit(1 if bad else 0)
EOF
