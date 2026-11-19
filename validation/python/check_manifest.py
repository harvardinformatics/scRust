"""The whole manifest (validation/manifest.toml) against the files on disk, every milestone's
entries alike (AGENTS.md: every reference dataset and validation output with its size and SHA-256).

Checks that the manifest parses (tomllib), that no path is listed twice, that every [[file]] entry
has a path, a sha256 and a size (bytes), that every listed file exists, and that every file's size
is the one its entry records. SHA-256 is checked for every committed file the manifest lists (git
ls-files: the fixtures, small, which a later milestone may change in place) and, with --hash, for
every entry of the milestones named ('all': every entry, which reads every listed file in full).

From the repository's top directory, on a compute node (it stats every listed file):

  PYTHONNOUSERSITE=1 .envs/loaders-py/bin/python -I validation/python/check_manifest.py [--hash <milestone>]...

Prints the counts and every failing entry; exits 1 if any entry fails. Milestone M15b
(validation/python/m15b_manifest.py writes M15b's section; this checks the whole manifest after it).
"""

import argparse
import collections
import hashlib
import os
import subprocess
import sys
import tomllib


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="validation/manifest.toml")
    p.add_argument("--hash", action="append", default=[], metavar="MILESTONE",
                   help="also check the SHA-256 of every entry of this milestone ('all': every entry)")
    a = p.parse_args()
    with open(a.manifest, "rb") as f:
        manifest = tomllib.load(f)
    files = manifest.get("file", [])
    ls = subprocess.run(["git", "--no-optional-locks", "ls-files", "-z"], check=True, capture_output=True)
    committed = set(ls.stdout.decode().split("\0")) - {""}
    failing = []
    for path, n in collections.Counter(e.get("path") for e in files).items():
        if n > 1:
            failing.append(f"{path}: listed {n} times")
    n_sized = n_hashed = 0
    for e in files:
        path, ms = e.get("path"), e.get("milestone", "no milestone")
        absent = [k for k in ("path", "sha256", "bytes") if k not in e]
        if absent:
            failing.append(f"{path} ({ms}): no {', '.join(absent)}")
            continue
        if not os.path.isfile(path):
            failing.append(f"{path} ({ms}): missing")
            continue
        n_sized += 1
        size = os.path.getsize(path)
        if size != e["bytes"]:
            failing.append(f"{path} ({ms}): {size} bytes, the entry records {e['bytes']}")
        if path in committed or ms in a.hash or "all" in a.hash:
            n_hashed += 1
            digest = sha256(path)
            if digest != e["sha256"]:
                failing.append(f"{path} ({ms}): sha256 {digest}, the entry records {e['sha256']}")
    hashed = "every entry" if "all" in a.hash else ", ".join(["the committed files"] + [f"{m}'s entries" for m in a.hash])
    print(f"{a.manifest}: {len(manifest.get('run', []))} runs, {len(files)} file entries "
          f"({sum(e.get('path') in committed for e in files)} committed files); sizes checked {n_sized}, "
          f"SHA-256 checked {n_hashed} ({hashed}); failing {len(failing)}")
    for x in failing:
        print("  " + x)
    sys.exit(1 if failing else 0)


if __name__ == "__main__":
    main()
