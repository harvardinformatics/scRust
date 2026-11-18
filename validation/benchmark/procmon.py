#!/usr/bin/env python3
"""procmon: run one command and sample its /proc counters at a fixed interval.

Benchmark tool for data/benchmarks (scRust): R's workflow scripts against the Rust ports.
It changes nothing in the measured program. Standard library only; runs under the system
python3 (3.6) with `PYTHONNOUSERSITE=1 python3 -I`.

usage:
  procmon.py --out <prefix> [--interval 0.01] [--watch <dir-or-file> ...] [--evict <path> ...]
             [--cpu <list>] [--sampler-cpu <list>] [--child-ncpus N] -- <command> [<args> ...]

Writes:
  <prefix>.samples.tsv  one row per sample of the measured process (and its descendants):
      t_s (monotonic seconds since just before the fork), utime/stime (clock ticks, the whole
      thread group, live and exited threads), cutime/cstime (waited-for children), threads,
      vmrss_kib, vmhwm_kib, rssanon_kib, rssfile_kib, rchar, wchar, syscr, syscw, read_bytes,
      write_bytes, desc_cpu_ticks / desc_rss_kib (live descendants other than the process
      itself), open_r / open_w (watched files open for reading / writing, with their offsets,
      as r<i>:<path relative to the i-th watched root (the .json's "watch" list)>)
  <prefix>.fds.tsv      every watched path seen open, per access mode (r, w, rw): the first and
                        last sample it was seen open in, and the largest offset seen
  <prefix>.json         the run: command, start (UTC and monotonic), exit status, wall time,
                        and getrusage of the child from wait4 (user, sys, maxrss: exact totals),
                        the sampler's own CPU time, the eviction report
--evict: before the fork, posix_fadvise(POSIX_FADV_DONTNEED) on each file (directories: each file
in them, recursively), so the run reads them cold from Lustre rather than from this node's page
cache; the report says how many bytes were resident before and after (mincore via mmap).
--cpu / --sampler-cpu: pin the measured process / this sampler to CPUs (sched_setaffinity);
--child-ncpus N: the measured process gets the first N CPUs the job may use, the sampler the rest.
"""

import argparse
import ctypes
import ctypes.util
import json
import mmap
import os
import resource
import signal
import stat as statmod
import sys
import time

CLK_TCK = os.sysconf("SC_CLK_TCK")
PAGE = os.sysconf("SC_PAGE_SIZE")

O_ACCMODE = 0o3


def parse_cpus(text):
    out = set()
    for part in text.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        elif part:
            out.add(int(part))
    return out


def read_file(path):
    with open(path, "rb") as fh:
        return fh.read()


def read_stat(pid):
    raw = read_file("/proc/%d/stat" % pid).decode()
    # the command name is in parentheses and may contain spaces
    rest = raw[raw.rindex(")") + 2:].split()
    # rest[0] is field 3 (state); field n is rest[n - 3]
    return {
        "minflt": int(rest[7]),
        "majflt": int(rest[9]),
        "utime": int(rest[11]),
        "stime": int(rest[12]),
        "cutime": int(rest[13]),
        "cstime": int(rest[14]),
        "threads": int(rest[17]),
    }


def read_status(pid):
    out = {}
    for line in read_file("/proc/%d/status" % pid).decode().splitlines():
        key, _, val = line.partition(":")
        if key in ("VmRSS", "VmHWM", "RssAnon", "RssFile", "VmSwap"):
            out[key] = int(val.split()[0])
    return out


def read_io(pid):
    out = {}
    for line in read_file("/proc/%d/io" % pid).decode().splitlines():
        key, _, val = line.partition(":")
        out[key.strip()] = int(val)
    return out


def children(pid):
    """Live descendants of pid (not pid itself), via /proc/<pid>/task/<tid>/children."""
    found = []
    stack = [pid]
    while stack:
        p = stack.pop()
        try:
            tids = os.listdir("/proc/%d/task" % p)
        except OSError:
            continue
        for tid in tids:
            try:
                kids = read_file("/proc/%d/task/%s/children" % (p, tid)).split()
            except OSError:
                continue
            for k in kids:
                k = int(k)
                found.append(k)
                stack.append(k)
    return found


def open_watched(pid, roots):
    """Watched paths the process has open: {path: (mode 'r'|'w'|'rw', pos)}."""
    out = {}
    fd_dir = "/proc/%d/fd" % pid
    try:
        fds = os.listdir(fd_dir)
    except OSError:
        return out
    for fd in fds:
        try:
            target = os.readlink(fd_dir + "/" + fd)
        except OSError:
            continue
        if not target.startswith("/"):
            continue
        if not any(target == r or target.startswith(r + "/") for r in roots):
            continue
        try:
            info = read_file("/proc/%d/fdinfo/%s" % (pid, fd)).decode()
        except OSError:
            continue
        pos = flags = 0
        for line in info.splitlines():
            if line.startswith("pos:"):
                pos = int(line.split()[1])
            elif line.startswith("flags:"):
                flags = int(line.split()[1], 8)
        acc = flags & O_ACCMODE
        mode = {os.O_RDONLY: "r", os.O_WRONLY: "w", os.O_RDWR: "rw"}.get(acc, "?")
        prev = out.get(target)
        if prev is None or pos > prev[1]:
            out[target] = (mode if prev is None else (prev[0] if prev[0] == mode else "rw"), pos)
    return out


def short(path, roots):
    """A watched path relative to its watched root, prefixed with the root's index (r0:, r1:)."""
    for i, r in enumerate(roots):
        if path == r:
            return "r%d:." % i
        if path.startswith(r + "/"):
            return "r%d:%s" % (i, path[len(r) + 1:])
    return path


_libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
_libc.mmap.restype = ctypes.c_void_p
_libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
_libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
_libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_char_p]
MAP_FAILED = ctypes.c_void_p(-1).value


def resident_bytes(path):
    """Bytes of the file resident in this node's page cache (mincore on a read-only mapping),
    or None when it cannot be told."""
    size = os.path.getsize(path)
    if size == 0:
        return 0
    fd = os.open(path, os.O_RDONLY)
    try:
        addr = _libc.mmap(None, size, mmap.PROT_READ, mmap.MAP_SHARED, fd, 0)
    finally:
        os.close(fd)
    if addr is None or addr == MAP_FAILED:
        return None
    try:
        npages = (size + PAGE - 1) // PAGE
        vec = ctypes.create_string_buffer(npages)
        if _libc.mincore(ctypes.c_void_p(addr), ctypes.c_size_t(size), vec) != 0:
            return None
        return sum(1 for b in vec.raw[:npages] if b & 1) * PAGE
    finally:
        _libc.munmap(ctypes.c_void_p(addr), ctypes.c_size_t(size))


def files_under(path):
    if os.path.isdir(path):
        for root, _dirs, names in os.walk(path, followlinks=True):
            for n in sorted(names):
                yield os.path.join(root, n)
    else:
        yield path


def evict(paths):
    report = []
    for p in paths:
        for f in files_under(p):
            real = os.path.realpath(f)
            if not statmod.S_ISREG(os.stat(real).st_mode):
                continue
            before = resident_bytes(real)
            fd = os.open(real, os.O_RDONLY)
            try:
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(fd)
            after = resident_bytes(real)
            report.append({"path": f, "real": real, "bytes": os.path.getsize(real),
                           "resident_before": before, "resident_after": after})
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--interval", type=float, default=0.01)
    ap.add_argument("--watch", action="append", default=[])
    ap.add_argument("--evict", action="append", default=[])
    ap.add_argument("--cpu")
    ap.add_argument("--sampler-cpu")
    ap.add_argument("--child-ncpus", type=int,
                    help="pin the command to the first N CPUs this job may use and the sampler to the rest")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        ap.error("no command")
    # the watched paths, resolved (/proc/<pid>/fd shows resolved paths), plus the targets of
    # symlinks inside watched directories (the inputs' 10x files can be symlinks)
    roots = set()
    for w in a.watch:
        roots.add(os.path.abspath(w))
        roots.add(os.path.realpath(w))
        if os.path.isdir(w):
            for root, _dirs, names in os.walk(w):
                for name in names:
                    q = os.path.join(root, name)
                    if os.path.islink(q):
                        roots.add(os.path.realpath(q))
    roots = sorted(set(r.rstrip("/") for r in roots))

    child_cpus = parse_cpus(a.cpu) if a.cpu else None
    sampler_cpus = parse_cpus(a.sampler_cpu) if a.sampler_cpu else None
    if a.child_ncpus:
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < a.child_ncpus:
            ap.error("--child-ncpus %d but only %d CPUs allowed" % (a.child_ncpus, len(allowed)))
        child_cpus = set(allowed[:a.child_ncpus])
        if len(allowed) > a.child_ncpus:
            sampler_cpus = set(allowed[a.child_ncpus:])
    if sampler_cpus:
        os.sched_setaffinity(0, sampler_cpus)
    evicted = evict(a.evict) if a.evict else []

    t_wall0 = time.time()
    t0 = time.monotonic()
    pid = os.fork()
    if pid == 0:
        try:
            if child_cpus is not None:
                os.sched_setaffinity(0, child_cpus)
            os.execvp(cmd[0], cmd)
        except BaseException as e:  # noqa: B902
            sys.stderr.write("procmon: exec failed: %s\n" % e)
        os._exit(127)

    cols = ["t_s", "utime", "stime", "cutime", "cstime", "threads", "minflt", "majflt",
            "vmrss_kib", "vmhwm_kib", "rssanon_kib", "rssfile_kib", "vmswap_kib",
            "rchar", "wchar", "syscr", "syscw", "read_bytes", "write_bytes",
            "desc_n", "desc_cpu_ticks", "desc_rss_kib", "open_r", "open_w"]
    samples = open(a.out + ".samples.tsv", "w")
    samples.write("\t".join(cols) + "\n")
    seen = {}  # path -> [first_t, last_t, mode, max_pos]
    n = 0
    status = None
    rusage = None
    next_t = t0
    while True:
        wpid, st, ru = os.wait4(pid, os.WNOHANG)
        if wpid == pid:
            status, rusage = st, ru
            break
        now = time.monotonic()
        try:
            s = read_stat(pid)
            v = read_status(pid)
            io = read_io(pid)
            kids = children(pid)
            kcpu = krss = 0
            for k in kids:
                try:
                    ks = read_stat(k)
                    kcpu += ks["utime"] + ks["stime"]
                    krss += read_status(k).get("VmRSS", 0)
                except OSError:
                    pass
            ow = open_watched(pid, roots) if roots else {}
        except (OSError, ValueError, IndexError):
            # the process is exiting: retry wait4 at once
            time.sleep(0.001)
            continue
        t = now - t0
        r_list = []
        w_list = []
        for path, (mode, pos) in sorted(ow.items()):
            key = (path, mode)
            rec = seen.get(key)
            if rec is None:
                seen[key] = [t, t, pos, n]
            else:
                rec[1] = t
                rec[2] = max(rec[2], pos)
            (w_list if "w" in mode else r_list).append("%s@%d" % (short(path, roots), pos))
        samples.write("\t".join(str(x) for x in [
            "%.4f" % t, s["utime"], s["stime"], s["cutime"], s["cstime"], s["threads"],
            s["minflt"], s["majflt"],
            v.get("VmRSS", 0), v.get("VmHWM", 0), v.get("RssAnon", 0), v.get("RssFile", 0),
            v.get("VmSwap", 0),
            io.get("rchar", 0), io.get("wchar", 0), io.get("syscr", 0), io.get("syscw", 0),
            io.get("read_bytes", 0), io.get("write_bytes", 0),
            len(kids), kcpu, krss, ";".join(r_list) or "-", ";".join(w_list) or "-"]) + "\n")
        n += 1
        next_t += a.interval
        delay = next_t - time.monotonic()
        if delay > 0:
            time.sleep(delay)
        else:
            next_t = time.monotonic()
    t_end = time.monotonic()
    samples.close()

    with open(a.out + ".fds.tsv", "w") as fh:
        fh.write("path\tmode\tfirst_t_s\tlast_t_s\tmax_pos\tfirst_sample\n")
        for (path, mode), (f, l, pos, first_n) in sorted(seen.items(), key=lambda kv: (kv[1][0], kv[0])):
            fh.write("%s\t%s\t%.4f\t%.4f\t%d\t%d\n" % (path, mode, f, l, pos, first_n))

    me = resource.getrusage(resource.RUSAGE_SELF)
    if os.WIFEXITED(status):
        exit_code, sig = os.WEXITSTATUS(status), None
    else:
        exit_code, sig = None, os.WTERMSIG(status)
    rec = {
        "cmd": cmd,
        "start_utc": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t_wall0)) + ("%.6f" % (t_wall0 % 1))[1:] + "Z",
        "start_epoch_s": t_wall0,
        "start_monotonic_s": t0,
        "wall_s": t_end - t0,
        "exit_code": exit_code,
        "signal": sig,
        "rusage": {"utime_s": rusage.ru_utime, "stime_s": rusage.ru_stime,
                   "cpu_s": rusage.ru_utime + rusage.ru_stime, "maxrss_kib": rusage.ru_maxrss,
                   "minflt": rusage.ru_minflt, "majflt": rusage.ru_majflt,
                   "inblock": rusage.ru_inblock, "oublock": rusage.ru_oublock,
                   "nvcsw": rusage.ru_nvcsw, "nivcsw": rusage.ru_nivcsw},
        "samples": n,
        "interval_s": a.interval,
        "clk_tck": CLK_TCK,
        "sampler_cpu_s": me.ru_utime + me.ru_stime,
        "child_cpus": sorted(child_cpus) if child_cpus is not None else None,
        "sampler_cpus": sorted(os.sched_getaffinity(0)),
        "watch": roots,
        "evicted": evicted,
    }
    with open(a.out + ".json", "w") as fh:
        json.dump(rec, fh, indent=1)
        fh.write("\n")
    if sig is not None:
        sys.exit(128 + sig)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
