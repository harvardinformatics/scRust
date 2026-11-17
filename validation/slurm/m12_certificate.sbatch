#!/bin/bash
# Milestone M12 (docs/PLAN.md C3e, C5: K9; C6: the decision point): neighborPurity teacher-forced
# on R's own embeddings and labels, with the certificate of its search, for every object and
# both OpenBLAS profiles' instrumented R runs. BLAS-free (purity reads only the embeddings and the
# labels), so it runs on any node. From the repository's top directory, after
# `cargo build --release -p scrust` and copying the binary to a file of its own:
#   SCRUST_BIN=<copy> sbatch validation/slurm/m12_certificate.sbatch <label> <object> [<object> ...]
#
# For each object (data/clustering/<object>/) and each instrumented R run on Sapphire Rapids
# (portable and production), `scrust cluster --snn <R run>/graphs --pca <R run>/pca --preset
# workflow` clusters R's SNN graph (M10: R's labels exactly) and scores R's embeddings, then
# `scrust validate cluster` compares it with the R run (its K9 arrays, radius and purities, and the
# teacher-forced K9 on R's embeddings, with that search's certificate). One row per run in
# data/clustering/m12_certificate/<label>/certificate.<job>.tsv: the cells, the radius, the
# flagged pairs (cells, step, kind where the binary labels it, how far below the radius squared
# in units in the last place), the width of the flagged band and how close the nearest unflagged
# pair came, the purity search's time, and the report's other checks (a staged start has no stage
# dumps; anything else not identical is named). This produced the decision point's evidence
# (binary scrust-m12-survey, before flags had kinds). Output,
# data/clustering/m12_certificate/<label>/<object>/<profile>/: run/, vs.json, vs.md, stderr.txt.
#
#SBATCH --job-name=scrust-m12-cert
#SBATCH -A informatics
#SBATCH --partition=shared
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh

cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -ge 2 ] || { echo "usage: m12_certificate.sbatch <label> <object> [<object> ...]" >&2; exit 2; }
label=$1
shift
source ./env.sh

echo "host=$(hostname) cpu='$(lscpu | sed -n 's/^Model name: *//p')' $(ldd --version | head -n 1)"
echo "OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE:-auto} $(rustc --version) commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
bin=${SCRUST_BIN:-target/release/scrust}
[ -x "$bin" ] || { echo "no $bin: run cargo build --release -p scrust first" >&2; exit 1; }
echo "binary=$bin sha256=$(sha256sum "$bin" | cut -c1-16) $("$bin" --version)"
threads=${SLURM_CPUS_PER_TASK:-1}
job=${SLURM_JOB_ID:-local}
top=data/clustering/m12_certificate/$label
mkdir -p "$top"
tsv=$top/certificate.$job.tsv
printf 'object\tprofile\tcells\tparity\tradius\tflagged\tflagged_kth\tflagged_radius\tpairs\tradius_band_ulps\tradius_gap_below_ulps\tkth_gap_ratio\tpurity_s\treport\n' > "$tsv"
status=0

for o in "$@"; do
    d=data/clustering/$o
    for profile in portable production; do
        ref=$d/runs/instrumented-$profile-sapphirerapids
        [ -d "$ref/graphs" ] && [ -d "$ref/pca" ] || { echo "$o: no $ref/graphs or pca" >&2; status=1; continue; }
        out=$top/$o/$profile
        if [ -e "$out" ]; then echo "$o $profile: $out exists, skipped"; continue; fi
        work=$out.partial.$job
        rm -rf "$work"
        mkdir -p "$work"
        if ! "$bin" cluster --snn "$ref/graphs" --pca "$ref/pca" --preset workflow --out "$work/run" --threads "$threads" 2> "$work/stderr.txt"; then
            echo "$o $profile: scrust cluster failed" >&2
            tail -5 "$work/stderr.txt" >&2
            status=1
            continue
        fi
        if "$bin" validate cluster --reference "$ref" --candidate "$work/run" --optimizer "$d/optimizer" \
                --json "$work/vs.json" --markdown "$work/vs.md" 2> "$work/vs.log"; then
            verdict=pass
        else
            verdict=FAIL
        fi
        mv "$work" "$out"
        python3 -I - "$o" "$profile" "$out" "$verdict" >> "$tsv" <<'EOF'
import json, sys
o, profile, out, verdict = sys.argv[1:5]
s = json.load(open(out + "/run/summary.json"))
p = s["result"]["purity"]
c = p["certificate"]
n = s["result"]["n_cells"]
pairs = ";".join(f'{q.get("kind", q["step"])}:{q["cell"]}-{q["other"]}@{q["ulps_from_threshold"]:.0f}' for q in c["pairs"])
r = json.load(open(out + "/vs.json"))
# the checks other than the certificates' (the decision point's) and the stage dumps a staged
# start cannot have
bad = [ch["name"] for ch in r["checks"] if ch["status"] in ("differ", "missing") and "certificate" not in ch["name"] and ch["name"] != "stage dumps"]
rep = "identical" if not bad else "differ: " + ",".join(bad)
fmt = lambda v: "-" if v is None else f"{v:.6g}"
print("\t".join([o, profile, str(n), "odd" if n % 2 else "even", f'{p["radius"]!r}', str(c["flagged_pairs"]), str(c["flagged_kth_distance"]),
                 str(c["flagged_radius"]), pairs or "-", fmt(c["radius_band_ulps"]), fmt(c["radius_gap_below_ulps"]), fmt(c["kth_gap_ratio"]),
                 f'{s["run"]["purity_s"]:.2f}', rep]))
EOF
        tail -1 "$tsv"
    done
done
exit $status
