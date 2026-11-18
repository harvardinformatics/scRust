#!/bin/bash
# Milestone M13's committed fixtures (docs/PLAN.md A1, C5's K10): what cargo test needs to check
# the corrected counts against M8's portable runs, without committing the matrices themselves.
# Reads plain files only, so a login node is enough. From the repository's top directory:
#   bash validation/slurm/fixtures_m13.sh
#
# For each of the seven tier-T objects cargo test covers (tests/fixtures/clustering/<object>/):
#   reference/k10.tsv   key and value: R's corrected counts as the portable instrumented run wrote
#                       them (data/clustering/<object>/runs/instrumented-portable-sapphirerapids/):
#                       the SHA-256 of sct/corrected_counts.mtx.gz's and
#                       sct/log1p_corrected_counts.mtx.gz's data lines (decompressed; the banner,
#                       the size line and the entries, each ending in a newline; R's comment lines
#                       left out), of sct_features.tsv.gz and sct_cells.tsv.gz (decompressed), R's
#                       own digests of the values (summary.json: sct_counts_x_sha256,
#                       sct_data_x_sha256), the number of entries, and correct()'s median log_umi
#                       and residuals beyond the clip range (stages/stages.json). The counts and
#                       nCount_SCT, nFeature_SCT are already committed (input/, reference/metadata.tsv.gz).
# For the two base objects also R's model for K10 teacher-forced (gzip -n -9):
#   reference/stages/vst/model_pars_fit.tsv.gz, reference/stages/vst/cell_attr.tsv.gz,
#   reference/stages/K10_genes.npy.gz
# The kernel of every file is the run's, Haswell (tests/fixtures/kernels.tsv: <object>/reference).
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/../.."
src_root=data/clustering
dst_root=tests/fixtures/clustering
primary=runs/instrumented-portable-sapphirerapids
all="filtered_seurat_tenx_test filtered_seurat_emptydrops_test soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test"
base_objects="filtered_seurat_tenx_test filtered_seurat_emptydrops_test"

for obj in $all; do
    r=$src_root/$obj/$primary
    [ -f "$r/sct/log1p_corrected_counts.mtx.gz" ] || { echo "$obj: no full-level portable run" >&2; exit 2; }
    python3 -I - "$r" "$dst_root/$obj/reference/k10.tsv" <<'EOF'
import gzip, hashlib, json, sys
run, out = sys.argv[1], sys.argv[2]
def data_lines(path):
    h, n = hashlib.sha256(), 0
    with gzip.open(path, "rb") as f:
        for i, line in enumerate(f):
            if i > 0 and line.startswith(b"%"):
                continue
            h.update(line if line.endswith(b"\n") else line + b"\n")
            n += 1
    return h.hexdigest(), n - 2
def plain(path):
    with gzip.open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()
res = json.load(open(run + "/summary.json"))["result"]
sc = json.load(open(run + "/stages/stages.json"))["scalars"]["chain"]
cnt, n1 = data_lines(run + "/sct/corrected_counts.mtx.gz")
dat, n2 = data_lines(run + "/sct/log1p_corrected_counts.mtx.gz")
assert n1 == n2 == res["sct_nnz"], (n1, n2, res["sct_nnz"])
rows = [
    ("source", run),
    ("kernel", json.load(open(run + "/summary.json"))["run"]["fingerprint"]["openblas_corename"]),
    ("corrected_counts_data_lines_sha256", cnt),
    ("log1p_data_lines_sha256", dat),
    ("sct_features_sha256", plain(run + "/sct/sct_features.tsv.gz")),
    ("sct_cells_sha256", plain(run + "/sct/sct_cells.tsv.gz")),
    ("sct_counts_x_sha256", res["sct_counts_x_sha256"]),
    ("sct_data_x_sha256", res["sct_data_x_sha256"]),
    ("sct_n_features", res["sct_n_features"]),
    ("sct_nnz", res["sct_nnz"]),
    ("correct_median_log_umi", repr(float(sc["correct_median_log_umi"]))),
    ("correct_residuals_beyond_clip", sc["correct_residuals_beyond_clip"]),
]
with open(out, "w") as f:
    f.write("key\tvalue\n")
    for k, v in rows:
        f.write(f"{k}\t{v}\n")
EOF
done
for obj in $base_objects; do
    r=$src_root/$obj/$primary
    d=$dst_root/$obj/reference/stages
    mkdir -p "$d/vst"
    gzip -n -9 -c "$r/stages/vst/model_pars_fit.tsv" > "$d/vst/model_pars_fit.tsv.gz"
    gzip -n -9 -c "$r/stages/vst/cell_attr.tsv" > "$d/vst/cell_attr.tsv.gz"
    gzip -n -9 -c "$r/stages/K10_genes.npy" > "$d/K10_genes.npy.gz"
done
du -ch "$dst_root"/*/reference/k10.tsv "$dst_root"/*/reference/stages/vst "$dst_root"/*/reference/stages/K10_genes.npy.gz | tail -1
