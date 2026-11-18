"""Compare a scrust clustering with a saved workflow object's (milestone M15).

The candidate is a `scrust cluster` output directory (here the singlets' re-clustering,
`scrust doubletfinder --recluster`'s `recluster/`, or `clusters/` of clustering on the fly); the
reference is an export's `saved/` directory (validation/R/export_clustering.R, any level), which
holds the saved object's metadata.tsv and its factors' levels. Checks, each identical or not:

- the cells: the barcodes, in order;
- every metadata column both have (seurat_clusters, SCT_snn_res.0.8, silhouette_width,
  neighborhood_purity, nCount_RNA, nFeature_RNA, and nCount_SCT and nFeature_SCT when the candidate
  exported them), as text: both sides write doubles with 17 significant digits, so equal text is
  equal bits;
- the levels of seurat_clusters and SCT_snn_res.0.8, byte for byte, and idents_levels.txt
  against seurat_clusters' levels;
- with --panel and --object, the numbers the workflow's job logged (validation/clustering_panel.tsv):
  nodes, edges, the maximum modularity as printed (%.4f), communities, singletons and the final
  clusters, against the candidate's summary.json result.

With "-" for the saved directory, only the logged numbers are compared (an object whose
metadata was not exported, such as a DoubletFinder input, which clustering on the fly reclusters).

Usage:
  python -I compare_saved_clustering.py <candidate dir> <saved dir>|- [--panel <tsv> --object <name>]
                                         [--columns a,b,...] [--json <report.json>]

Exit status 0 when every check is made and identical.
"""

import csv
import json
import os
import sys

COLUMNS = ["seurat_clusters", "SCT_snn_res.0.8", "silhouette_width", "neighborhood_purity", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT"]


def read_tsv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    return rows[0], rows[1:]


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def compare_metadata(cand, saved, columns, check):
    """The metadata columns and the levels, against the saved object's export."""
    if saved is None:
        return
    ch, crows = read_tsv(os.path.join(cand, "metadata.tsv"))
    sh, srows = read_tsv(os.path.join(saved, "metadata.tsv"))
    n = len(srows)
    cb, sb = ch.index("barcode"), sh.index("barcode")
    same_n = len(crows) == n
    differ = max(n, len(crows)) if not same_n else sum(1 for a, b in zip(crows, srows) if a[cb] != b[sb])
    check("barcodes", n, differ, f"{len(crows)} cells in the candidate, {n} saved")
    for col in columns or COLUMNS:
        if col not in ch or col not in sh:
            if columns:
                check(col, 0, 0, "absent from " + ("the candidate" if col not in ch else "the saved object"))
            continue
        a, b = ch.index(col), sh.index(col)
        if not same_n:
            check(col, n, n, "different cells")
            continue
        bad = [k for k in range(n) if crows[k][a] != srows[k][b]]
        check(col, n, len(bad), "" if not bad else f"first at cell {bad[0] + 1}: candidate {crows[bad[0]][a]}, saved {srows[bad[0]][b]}")
    for f in ["seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt"]:
        p, q = os.path.join(cand, f), os.path.join(saved, f)
        if os.path.exists(p) and os.path.exists(q):
            check(f, 1, int(read_bytes(p) != read_bytes(q)))
    p, q = os.path.join(cand, "seurat_clusters_levels.txt"), os.path.join(saved, "idents_levels.txt")
    if os.path.exists(p) and os.path.exists(q):
        check("idents_levels.txt (the identities' levels)", 1, int(read_bytes(p) != read_bytes(q)))


def main(argv):
    args = argv[1:]
    report = None
    panel = obj = None
    columns = None
    pos = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("--json", "--panel", "--object", "--columns"):
            v = args[i + 1]
            if a == "--json":
                report = v
            elif a == "--panel":
                panel = v
            elif a == "--object":
                obj = v
            else:
                columns = v.split(",")
            i += 2
        else:
            pos.append(a)
            i += 1
    if len(pos) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    cand, saved = pos
    checks = []

    def check(name, compared, differing, detail=""):
        checks.append({"check": name, "compared": compared, "differing": differing, "status": "identical" if differing == 0 and compared > 0 else ("differs" if differing else "not made"), "detail": detail})

    if saved == "-":
        compare_metadata(cand, None, columns, check)
    else:
        compare_metadata(cand, saved, columns, check)
    if panel and obj:

        with open(panel, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
        row = next((r for r in rows if r["object"] == obj), None)
        with open(os.path.join(cand, "summary.json")) as f:
            res = json.load(f)["result"]
        r0 = res["resolutions"][-1]
        got = {
            "log_nodes": str(res["nodes"]),
            "log_edges": str(res["edges"]),
            "log_max_modularity": r0["max_modularity_printed"],
            "log_communities": str(r0["communities"]),
            "log_singletons": str(r0["n_singletons"]),
            "log_final_clusters": str(r0["final_clusters"]),
        }
        if row is None:
            check("logged numbers", 0, 0, f"{obj} is not in {panel}")
        else:
            bad = [f"{k}: logged {row[k]}, candidate {v}" for k, v in got.items() if row[k] != v]
            check("logged numbers (nodes, edges, maximum modularity, communities, singletons, clusters)", len(got), len(bad), "; ".join(bad))
    made = [c for c in checks if c["status"] != "not made"]
    ok = all(c["status"] == "identical" for c in checks) and len(made) > 0
    out = {"candidate": cand, "saved": saved, "object": obj, "pass": ok, "checks": checks}
    if report:
        with open(report, "w") as f:
            json.dump(out, f, indent=1)
    for c in checks:
        print(f"{c['status']:10} {c['check']} ({c['compared']}) {c['detail']}")
    print("pass" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
