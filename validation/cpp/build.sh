#!/usr/bin/env bash
# Builds milestone M8's C++ reference programs (docs/PLAN.md C4) into data/clustering/bin/
# (gitignored), with the production R env's own compiler and the flags of its Makeconf (`R CMD
# config`), as R builds a package's C++17 code:
#   ModularityOptimizer   Seurat 5.5.1's src/ModularityOptimizer.cpp in its standalone mode
#                         (-DSTANDALONE, the file's own main()), with ModularityOptimizer.patch,
#                         which only prints each random start's iterations, modularity (%.17g) and
#                         partition on a line starting with SCRUST_START
#   javarandom_vectors    javarandom_vectors.cpp linked with the unmodified ModularityOptimizer.cpp:
#                         the K0 vectors of the optimizer's JavaRandom
#   kiss64_vectors        kiss64_vectors.cpp with RcppAnnoy 0.0.23's kissrandom.h, as installed in
#                         the R env: the K0 vectors of Annoy's Kiss64Random
# and records the compiler, flags and source checksums in data/clustering/bin/build.json.
#
# Usage, from the repository's top directory (a login node is enough; it reads no .rds):
#   bash validation/cpp/build.sh
# Needs Seurat's sources from third_party/fetch_sources.sh (third_party/src/Seurat-5.5.1/).
# The R env is read-only: nothing is written there.
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
seurat_src=$root/third_party/src/Seurat-5.5.1/src
annoy_inc=$R_ENV/lib/R/library/RcppAnnoy/include
bin=$root/data/clustering/bin
build=$bin/build
here=$root/validation/cpp

[ -f "$seurat_src/ModularityOptimizer.cpp" ] || { echo "run third_party/fetch_sources.sh first" >&2; exit 2; }
want=$(awk '$2 == "Seurat_5.5.1.tar.gz" { print $1 }' "$root/third_party/SHA256SUMS")
[ "$(cat "$root/third_party/src/.Seurat-5.5.1.sha256")" = "$want" ] || { echo "third_party/src/Seurat-5.5.1 is not the verified tarball" >&2; exit 2; }

rconfig() { env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB "$R_ENV/bin/R" CMD config "$1"; }
CXX="$R_ENV/bin/$(rconfig CXX17)"
CXXSTD=$(rconfig CXX17STD)
CPPFLAGS=$(rconfig CPPFLAGS)
CXXFLAGS=$(rconfig CXX17FLAGS)
LDFLAGS=$(rconfig LDFLAGS)
read -r -a cppflags <<< "$CPPFLAGS"
read -r -a cxxflags <<< "$CXXFLAGS"
read -r -a ldflags <<< "$LDFLAGS"

rm -rf "$build"
mkdir -p "$build/orig" "$build/patched"
cp "$seurat_src/ModularityOptimizer.cpp" "$seurat_src/ModularityOptimizer.h" "$build/orig/"
cp "$seurat_src/ModularityOptimizer.cpp" "$seurat_src/ModularityOptimizer.h" "$build/patched/"
(cd "$build/patched" && patch -p1 --quiet --no-backup-if-mismatch -d . ModularityOptimizer.cpp < "$here/ModularityOptimizer.patch")

compile() {  # output, then sources and extra flags
    local out=$1
    shift
    echo "building $(basename "$out")"
    "$CXX" "$CXXSTD" "${cppflags[@]}" "${cxxflags[@]}" "$@" -o "$out" "${ldflags[@]}"
}
compile "$bin/ModularityOptimizer" -DSTANDALONE -I "$build/patched" "$build/patched/ModularityOptimizer.cpp"
compile "$bin/javarandom_vectors" -I "$build/orig" "$here/javarandom_vectors.cpp" "$build/orig/ModularityOptimizer.cpp"
compile "$bin/kiss64_vectors" -I "$annoy_inc" "$here/kiss64_vectors.cpp"

sha() { sha256sum "$1" | cut -d' ' -f1; }
cat > "$bin/build.json" <<EOF
{
  "compiler": "$CXX",
  "compiler_version": "$("$CXX" --version | head -n 1)",
  "cxxstd": "$CXXSTD",
  "cppflags": "$CPPFLAGS",
  "cxxflags": "$CXXFLAGS",
  "ldflags": "$LDFLAGS",
  "seurat_tarball_sha256": "$want",
  "ModularityOptimizer.cpp_sha256": "$(sha "$build/orig/ModularityOptimizer.cpp")",
  "ModularityOptimizer.h_sha256": "$(sha "$build/orig/ModularityOptimizer.h")",
  "patch_sha256": "$(sha "$here/ModularityOptimizer.patch")",
  "patched_cpp_sha256": "$(sha "$build/patched/ModularityOptimizer.cpp")",
  "kissrandom.h_sha256": "$(sha "$annoy_inc/kissrandom.h")",
  "javarandom_vectors.cpp_sha256": "$(sha "$here/javarandom_vectors.cpp")",
  "kiss64_vectors.cpp_sha256": "$(sha "$here/kiss64_vectors.cpp")",
  "ModularityOptimizer_sha256": "$(sha "$bin/ModularityOptimizer")",
  "javarandom_vectors_sha256": "$(sha "$bin/javarandom_vectors")",
  "kiss64_vectors_sha256": "$(sha "$bin/kiss64_vectors")",
  "built": "$(date -Iseconds)",
  "host": "$(hostname)"
}
EOF
echo "built into ${bin#"$root"/}/"
