#!/bin/bash
# Milestone M12 (docs/PLAN.md C5: K9; validation/tolerances.toml [clustering.floor.kernel]
# K9_neighbors): the R-vs-R floor of neighborPurity's neighbour lists, cell by cell. M11's floor
# measurement (m11_floor.sbatch) compares arrays block by block, so where a list's length changes it
# can only call K9_neighbors_idx incomparable; scrust validate cluster compares the lists cell by
# cell, and this measures the same on R's own runs: the portable (Haswell) and production
# (SapphireRapids) instrumented runs on Sapphire Rapids of every object whose runs kept the lists
# (level full, and the small tier-P object). Reads plain .npy files only (no R, no Rust; numpy, in
# seconds). From the repository's top directory:
#   bash validation/slurm/m12_lists_floor.sh > data/clustering/floor_kernel_m12_lists.tsv
# One row per object: its cells, the cells whose list differs between the kernels, those that hold
# the same cells in another order, and the largest change in a list's length; a last line sums them.
set -euo pipefail
python3 -I - <<'EOF'
import os
import numpy as np

root = "data/clustering"
print("object\tcells\tcells_differing\treordered_only\tmax_length_change")
tot = [0, 0, 0, 0]
for o in sorted(os.listdir(root)):
    a = f"{root}/{o}/runs/instrumented-portable-sapphirerapids/stages"
    b = f"{root}/{o}/runs/instrumented-production-sapphirerapids/stages"
    if not (os.path.exists(f"{a}/K9_neighbors_ptr.npy") and os.path.exists(f"{b}/K9_neighbors_ptr.npy")):
        continue
    pa, ia = np.load(f"{a}/K9_neighbors_ptr.npy"), np.load(f"{a}/K9_neighbors_idx.npy")
    pb, ib = np.load(f"{b}/K9_neighbors_ptr.npy"), np.load(f"{b}/K9_neighbors_idx.npy")
    n = len(pa) - 1
    differ = reordered = maxlen = 0
    for i in range(n):
        x, y = ia[pa[i]:pa[i + 1]], ib[pb[i]:pb[i + 1]]
        if len(x) == len(y) and (x == y).all():
            continue
        differ += 1
        maxlen = max(maxlen, abs(len(x) - len(y)))
        if sorted(x) == sorted(y):
            reordered += 1
    print(f"{o}\t{n}\t{differ}\t{reordered}\t{maxlen}")
    tot = [tot[0] + 1, tot[1] + differ, tot[2] + reordered, max(tot[3], maxlen)]
print(f"total: {tot[0]} objects\t-\t{tot[1]}\t{tot[2]}\t{tot[3]}")
EOF
