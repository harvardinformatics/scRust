#!/bin/bash
# Benchmark (scRust data/benchmarks): checks, on a Sapphire Rapids node, the measurement features
# the benchmark relies on:
#  1. writing 5 to /proc/self/clear_refs resets VmHWM (kernel 4.18), from Python and from R;
#  2. /proc/<pid>/io on Lustre reads: rchar and read_bytes;
#  3. posix_fadvise(DONTNEED) evicts a Lustre file from this node's page cache (procmon --evict);
#  4. procmon.py on a toy command (samples, fds, rusage) and the sampler's own cost at 10 ms.
#SBATCH --job-name=bench-test
#SBATCH -A informatics
#SBATCH --partition=sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-00:30:00
#SBATCH --output=data/benchmarks/logs/%x-%j.log

set -uo pipefail
source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh
S=$(realpath "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")  # the repository (submit from its top directory)
B=$S/data/benchmarks
REnv=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
T=$B/tests/$SLURM_JOB_ID
mkdir -p "$T"
cd "$T"
PY="env PYTHONNOUSERSITE=1 /usr/bin/python3 -I"
echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' kernel=$(uname -r) cpus=$(taskset -cp $$)"

echo "== 1a. clear_refs from python"
$PY - <<'EOF'
import os
def hwm():
    for l in open("/proc/self/status"):
        if l.startswith(("VmHWM", "VmRSS")):
            print("   ", l.strip())
x = bytearray(800 * 1024 * 1024)
for i in range(0, len(x), 4096): x[i] = 1
hwm()
del x
print("  after del:"); hwm()
with open("/proc/self/clear_refs", "w") as fh: fh.write("5\n")
print("  after clear_refs 5:"); hwm()
EOF

echo "== 1b. clear_refs from R (production R env)"
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 BENCH_MARKERS="$T/r_markers.tsv" \
  "$REnv/bin/Rscript" -e '
source("'"$S"'/validation/benchmark/workflow/scripts/bench_markers.R")
bench_phase("alloc", "compute")
x <- numeric(1e8); x[] <- 1
bench_phase("free", "compute")
rm(x); invisible(gc())
bench_phase("small", "compute")
y <- numeric(1e6); y[] <- 2
bench_end()
'
column -t "$T/r_markers.tsv"

echo "== 2/3. Lustre read counters and eviction"
f=$B/../doubletfinder/soupx_emptydrops_sscovelli/input/matrix.mtx.gz
ls -la "$(readlink -f "$f")"
command -v fincore && fincore "$(readlink -f "$f")"
PROCMON="$S/validation/benchmark/procmon.py" $PY - "$f" <<'EOF'
import os, sys, time
sys.path.insert(0, os.environ.get("B_SCRIPTS", ""))
f = os.path.realpath(sys.argv[1])
def io():
    return {k: int(v) for k, v in (l.split(":") for l in open("/proc/self/io"))}
import importlib.util
spec = importlib.util.spec_from_file_location("procmon", os.environ["PROCMON"])
pm = importlib.util.module_from_spec(spec); spec.loader.exec_module(pm)
print("  resident before:", pm.resident_bytes(f), "of", os.path.getsize(f))
for attempt in ("first read", "second read (cached?)", "after evict"):
    if attempt == "after evict":
        print("  evict:", pm.evict([f]))
    a = io(); t = time.monotonic()
    with open(f, "rb") as fh:
        while fh.read(1 << 20): pass
    b = io(); dt = time.monotonic() - t
    print("  %-22s %.2f s  rchar +%d  read_bytes +%d  resident now %s" % (attempt, dt, b["rchar"] - a["rchar"], b["read_bytes"] - a["read_bytes"], pm.resident_bytes(f)))
EOF
command -v fincore && fincore "$(readlink -f "$f")"

echo "== 4. procmon on a toy command (python allocating, writing and reading a file)"
mkdir -p "$T/toy"
$PY "$S/validation/benchmark/procmon.py" --out "$T/toy/run" --interval 0.01 --watch "$T/toy" --child-ncpus 1 -- \
  /usr/bin/python3 -I -c '
import time
x = bytearray(300 * 1024 * 1024)
for i in range(0, len(x), 4096): x[i] = 1
time.sleep(0.3)
with open("'"$T/toy"'/out.bin", "wb") as fh:
    for i in range(50): fh.write(b"x" * (1 << 20)); time.sleep(0.005)
del x
time.sleep(0.2)
with open("'"$T/toy"'/out.bin", "rb") as fh:
    while fh.read(1 << 20): time.sleep(0.002)
t = time.process_time()
while time.process_time() - t < 0.5: pass
'
echo "exit=$?"
cat "$T/toy/run.json"
cat "$T/toy/run.fds.tsv"
wc -l "$T/toy/run.samples.tsv"
awk -F'\t' 'NR==1 || NR%10==2' "$T/toy/run.samples.tsv" | cut -f1-3,9,10,14,15,23,24 | column -t | head -40
