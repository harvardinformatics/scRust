"""Milestone M15b: scrust's validation runs against R's, case by case, in one table.

Reads validation/slurm/m15b_revalidate.sbatch's reports (data/clustering/m15b/revalidate/<label>/
<object>/<case>/<profile>/validate_<R profile>.json: scrust's run of each case in each profile
against R's runs of the case in both profiles) and writes one row per case, scrust's profile and
R's profile: the grading (exact when the two runs' OpenBLAS kernels match, else tolerance), the
floor, whether it passes, the checks by status, the values compared (the sum of every check's n)
and those differing, and the CPU models of the runs. Prints the totals: the main criterion first
(the same profile: exact), then across kernels.

  PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I validation/python/report_m15b.py <label> [--tsv <file>]
"""

import argparse
import csv
import json
import os
import re

ROOT = os.getcwd()


def cpu(path):
    try:
        text = open(path).read()
    except OSError:
        return ""
    m = re.search(r'"cpu_model": *"([^"]*)"', text)
    return m.group(1).replace("Intel(R) Xeon(R) Platinum ", "") if m else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("label")
    p.add_argument("--runs", default="m15b-25fe846", help="the validation runs' label")
    p.add_argument("--tsv")
    a = p.parse_args()
    cases = list(csv.DictReader(open(os.path.join(ROOT, "validation", "m15b_cases.tsv")), delimiter="\t"))
    rows = []
    for c in cases:
        for prof in ("production", "portable"):
            d = os.path.join(ROOT, "data", "clustering", "m15b", "revalidate", a.label, c["object"], c["case"], prof)
            run_fp = os.path.join(ROOT, "data", "clustering", "m15b", "rust", a.runs, c["object"], c["case"], prof, "fingerprint.json")
            for rp in ("production", "portable"):
                f = os.path.join(d, f"validate_{rp}.json")
                if not os.path.isfile(f):
                    continue
                j = json.load(open(f))
                st = {}
                for ch in j["checks"]:
                    st[ch["status"]] = st.get(ch["status"], 0) + 1
                # values compared and differing: the graded checks (not the n/a ones, nor the
                # purity certificate's pairs, which are reported, not compared)
                graded = [ch for ch in j["checks"] if ch["status"] in ("identical", "within floor", "differ", "missing")]
                n = sum(ch.get("n", 0) for ch in graded)
                nd = sum(ch.get("n_differ", 0) for ch in graded)
                bad = [f"{ch['stage']}: {ch['name']}" for ch in j["checks"] if ch["status"] in ("differ", "missing")]
                rdir = os.path.dirname(j["reference"]) if j["reference"].endswith("/") else j["reference"]
                rows.append({
                    "object": c["object"], "tier": c["tier"], "case": c["case"], "vars_to_regress": c["vars_to_regress"],
                    "scrust_profile": prof, "r_profile": rp, "grading": j["grading"]["mode"], "floor": j.get("floor", ""),
                    "pass": "yes" if j["pass"] else "NO", "identical": st.get("identical", 0), "within_floor": st.get("within floor", 0),
                    "differing": st.get("differ", 0), "missing": st.get("missing", 0), "not_applicable": st.get("n/a", 0),
                    "reported": st.get("flagged, reported", 0), "values_compared": n, "values_differing": nd,
                    "scrust_cpu": cpu(run_fp), "r_cpu": cpu(os.path.join(rdir, "fingerprint_wrapper.json")),
                    "failing": "; ".join(bad),
                })
    if a.tsv:
        with open(a.tsv, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    def total(sel, title):
        rs = [r for r in rows if sel(r)]
        if not rs:
            return
        ok = sum(r["pass"] == "yes" for r in rs)
        tiers = {}
        for r in rs:
            tiers.setdefault(r["tier"], [0, 0])
            tiers[r["tier"]][0] += 1
            tiers[r["tier"]][1] += r["pass"] == "yes"
        print(f"{title}: {ok} of {len(rs)} pass ({', '.join(f'{t} {v[1]}/{v[0]}' for t, v in sorted(tiers.items()))}); "
              f"checks identical {sum(r['identical'] for r in rs)}, within floor {sum(r['within_floor'] for r in rs)}, "
              f"differing {sum(r['differing'] for r in rs)}, missing {sum(r['missing'] for r in rs)}, reported {sum(r['reported'] for r in rs)}; "
              f"values compared {sum(r['values_compared'] for r in rs):,}, differing {sum(r['values_differing'] for r in rs):,}; "
              f"gradings {sorted({r['grading'] for r in rs})}, floors {sorted({r['floor'] for r in rs})}")
        for r in rs:
            if r["pass"] != "yes":
                print(f"   {r['object']} {r['case']} scrust {r['scrust_profile']} vs R {r['r_profile']}: {r['failing']}")
    total(lambda r: r["scrust_profile"] == r["r_profile"] == "production", "the same profile, production (exact)")
    total(lambda r: r["scrust_profile"] == r["r_profile"] == "portable", "the same profile, portable (exact)")
    total(lambda r: r["scrust_profile"] == "portable" and r["r_profile"] == "production", "across kernels, scrust portable against R production")
    total(lambda r: r["scrust_profile"] == "production" and r["r_profile"] == "portable", "across kernels, scrust production against R portable")
    pairs = sorted({(r["scrust_cpu"], r["r_cpu"]) for r in rows if r["scrust_profile"] == r["r_profile"] == "portable"})
    print("CPU models of the portable pairs (scrust, R):", pairs)


if __name__ == "__main__":
    main()
