"""Compare two scrust output directories file by file (milestone M15).

Every file must be byte-identical, except `summary.json` (at any depth), which is compared as
JSON with the members that record how the run went set aside: its timings, its start time and
peak memory, and the provenance's command line (``files.provenance.command_line``). Optionally
more members can be set aside (``--ignore <file>:<dotted path>``, repeatable, where ``<file>`` is
the summary's relative path, e.g. ``summary.json:result.clustering``, or ``*`` for every
summary), and files present on one side only can be allowed under given directories
(``--only-left``, ``--only-right``): so are their entries in that side's manifests and lists of
files written (``files.entries`` and ``run.files`` of each ``summary.json`` above them), which the
other side cannot have.

``--nested-summaries`` (M15b) also sets aside, in each manifest (``files.entries`` of a
``summary.json``), the ``sha256`` and ``bytes`` of the entries that list a nested output's
``summary.json``: that file holds timings, so its digest differs between any two runs, while the file
itself is compared as JSON like every summary. Without it, those members are compared, as M15 did.

Used for:
- K12: ``scrust doubletfinder`` clustering on the fly against M7's runs (``--only-right
  clusters`` and ``--ignore`` for the clustering's own members of summary.json), and its
  ``clusters/`` against ``scrust cluster --counts <input>``'s output (``--only-left pca
  --only-left neighbors``: only the labels are kept);
- the re-clustering's check (a): ``scrust doubletfinder --recluster`` against the two-command
  route (``scrust doubletfinder``, then ``scrust cluster --counts <out>/singlets``);
- check (c): ``scrust doubletfinder`` without the new options against M7's runs.

Usage:
  python -I compare_runs.py <left dir> <right dir> [--json <report.json>] [--ignore <file>:a.b]...
                            [--only-left <dir>]... [--only-right <dir>]... [--result-only]
                            [--nested-summaries]

``--result-only`` compares only ``result`` of the top-level ``summary.json`` (M7's runs record
another command line, inputs and paths). Exit status 0 when nothing differs.
"""

import hashlib
import json
import os
import sys

# members of summary.json that record how a run went, not what it computed (any depth)
VOLATILE = {
    "start_time", "elapsed_s", "phases", "peak_rss_kib", "printed_elapsed",
    "sctransform_s", "pca_s", "neighbors_s", "cluster_s", "silhouette_s", "purity_s",
    "network_s", "resolutions_s", "optimizer_s", "groupsingletons_s",
    "export_counts_copy_s", "manifest_s",
}


def files_of(root):
    out = {}
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            out[os.path.relpath(p, root).replace(os.sep, "/")] = p
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def strip(j, ignore, path=""):
    """j without the volatile members and the ignored dotted paths."""
    if isinstance(j, dict):
        out = {}
        for k, v in j.items():
            p = f"{path}.{k}" if path else k
            if k in VOLATILE or p in ignore:
                continue
            out[k] = strip(v, ignore, p)
        return out
    if isinstance(j, list):
        return [strip(v, ignore, path) for v in j]
    return j


def diff_json(a, b, path=""):
    """The dotted paths where a and b differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            p = f"{path}.{k}" if path else k
            if k not in a or k not in b:
                out.append(p + (" (right only)" if k not in a else " (left only)"))
            else:
                out.extend(diff_json(a[k], b[k], p))
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path} (lengths {len(a)}, {len(b)})"]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out.extend(diff_json(x, y, f"{path}[{i}]"))
        return out
    return [] if a == b else [path]


def drop_entries(j, rel, dirs):
    """The summary.json ``j`` (at relative path ``rel``) without its manifest's entries
    (``files.entries``) and its list of files written (``run.files``) under ``dirs``; and how many
    were dropped."""
    if not dirs or not isinstance(j, dict):
        return j, 0
    base = os.path.dirname(rel)
    base = base + "/" if base else ""
    under = lambda p: any((base + str(p)).startswith(d) for d in dirs)
    out, n = dict(j), 0
    f = j.get("files")
    if isinstance(f, dict) and isinstance(f.get("entries"), list):
        keep = [e for e in f["entries"] if not (isinstance(e, dict) and under(e.get("path", "")))]
        n += len(f["entries"]) - len(keep)
        out["files"] = {**f, "entries": keep}
    r = j.get("run")
    if isinstance(r, dict) and isinstance(r.get("files"), list):
        keep = [p for p in r["files"] if not (isinstance(p, str) and under(p))]
        n += len(r["files"]) - len(keep)
        out["run"] = {**r, "files": keep}
    return out, n


def drop_nested_digests(j):
    """The summary.json ``j`` with the ``sha256`` and ``bytes`` of its manifest's entries for nested
    ``summary.json`` files set aside (``--nested-summaries``)."""
    f = j.get("files") if isinstance(j, dict) else None
    if not (isinstance(f, dict) and isinstance(f.get("entries"), list)):
        return j
    nested = lambda e: isinstance(e, dict) and str(e.get("path", "")).endswith("/summary.json")
    entries = [{k: v for k, v in e.items() if k not in ("sha256", "bytes")} if nested(e) else e for e in f["entries"]]
    return {**j, "files": {**f, "entries": entries}}


def main(argv):
    args = argv[1:]
    report = None
    ignore = set(["*:files.provenance.command_line"])
    only_left, only_right = [], []
    result_only = False
    nested_summaries = False
    pos = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--json":
            report = args[i + 1]
            i += 2
        elif a == "--ignore":
            if ":" not in args[i + 1]:
                print("--ignore takes <file>:<dotted path>", file=sys.stderr)
                return 2
            ignore.add(args[i + 1])
            i += 2
        elif a == "--only-left":
            only_left.append(args[i + 1].rstrip("/") + "/")
            i += 2
        elif a == "--only-right":
            only_right.append(args[i + 1].rstrip("/") + "/")
            i += 2
        elif a == "--result-only":
            result_only = True
            i += 1
        elif a == "--nested-summaries":
            nested_summaries = True
            i += 1
        else:
            pos.append(a)
            i += 1
    if len(pos) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    left, right = pos
    fl, fr = files_of(left), files_of(right)
    rows = []
    differ = 0
    for rel in sorted(set(fl) | set(fr)):
        if rel not in fr:
            ok = any(rel.startswith(d) for d in only_left)
            rows.append({"file": rel, "status": "left only" + (", allowed" if ok else "")})
            differ += not ok
            continue
        if rel not in fl:
            ok = any(rel.startswith(d) for d in only_right)
            rows.append({"file": rel, "status": "right only" + (", allowed" if ok else "")})
            differ += not ok
            continue
        if os.path.basename(rel) == "summary.json":
            with open(fl[rel]) as f:
                a = json.load(f)
            with open(fr[rel]) as f:
                b = json.load(f)
            if result_only and rel == "summary.json":
                a, b = {"result": a.get("result")}, {"result": b.get("result")}
            a, na = drop_entries(a, rel, only_left)
            b, nb = drop_entries(b, rel, only_right)
            if nested_summaries:
                a, b = drop_nested_digests(a), drop_nested_digests(b)
            ign = {p.split(":", 1)[1] for p in ignore if p.split(":", 1)[0] in ("*", rel)}
            d = diff_json(strip(a, ign), strip(b, ign))
            row = {"file": rel, "status": "identical (as JSON, timings and command line aside)" if not d else "differs", "members": d}
            if na or nb:
                row["entries_set_aside"] = {"left": na, "right": nb}
            rows.append(row)
            differ += bool(d)
            continue
        same = os.path.getsize(fl[rel]) == os.path.getsize(fr[rel]) and sha256(fl[rel]) == sha256(fr[rel])
        rows.append({"file": rel, "status": "identical" if same else "differs", "bytes": os.path.getsize(fl[rel])})
        differ += not same
    out = {"left": left, "right": right, "files": len(rows), "differing": differ, "rows": rows}
    if report:
        with open(report, "w") as f:
            json.dump(out, f, indent=1)
    for r in rows:
        if not r["status"].startswith("identical"):
            print(f"{r['status']}: {r['file']} {r.get('members', '')}")
    print(f"{len(rows) - differ} of {len(rows)} files identical; {differ} differing")
    return 0 if differ == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
