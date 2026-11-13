// JavaRandom reference vectors (docs/PLAN.md C2 and C5, stage K0; milestone M8) from Seurat
// 5.5.1's own C++ (src/ModularityOptimizer.cpp, compiled unmodified and linked in): the
// java.util.Random generator that FindClusters' modularity optimizer seeds once per resolution,
// nextInt() with its power-of-two branch and its rejection loop, and
// Arrays2::generateRandomPermutation(), which each local-moving pass draws.
//
// Built by validation/cpp/build.sh with the R env's compiler and Makeconf flags; run as
//   data/clustering/bin/javarandom_vectors <out dir>
// Writes, tab-separated, lines starting with "#" naming the columns:
//   javarandom_nextint.tsv      seed, n, count, the first `count` values of nextInt(n) from a fresh
//                               generator seeded with `seed`
//   javarandom_permutation.tsv  seed, the sizes of successive permutations drawn from ONE generator
//                               (the stream continues across them, as in the optimizer), and each
//                               permutation, ';'-separated
// Seeds are given as the int the R wrapper passes (RunModularityClusteringCpp's `int randomSeed`),
// converted to uint64_t as the call does; negative ints therefore sign-extend.

#include "ModularityOptimizer.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace ModularityOptimizer;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: javarandom_vectors <out dir>\n");
    return 2;
  }
  std::string dir = argv[1];
  const std::vector<int> seeds = {0, 42, 12345, -1, 2147483647, -2147483647 - 1};
  // powers of two (multiply-and-shift branch), small and large n, n near 2^31 (rejection loop
  // with signed wrap-around), and sizes of real networks
  const std::vector<int> ns = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 17, 20, 31, 32, 33, 64, 100,
                               127, 128, 255, 256, 267, 276, 300, 1000, 1023, 1024, 1025, 3108,
                               11831, 26598, 59323, 65536, 79409, 1048576, 1073741823, 1073741824,
                               1073741825, 1500000000, 2147483646, 2147483647};
  const int count = 100;

  FILE* f = std::fopen((dir + "/javarandom_nextint.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tn\tcount\tnextInt(n) values from a fresh JavaRandom(seed), comma-separated\n");
  for (int seed : seeds) {
    for (int n : ns) {
      JavaRandom random(seed);
      std::fprintf(f, "%d\t%d\t%d\t", seed, n, count);
      for (int k = 0; k < count; k++) std::fprintf(f, k ? ",%d" : "%d", random.nextInt(n));
      std::fprintf(f, "\n");
    }
  }
  std::fclose(f);

  f = std::fopen((dir + "/javarandom_permutation.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tsizes\tpermutations (0-based), one per size in order, drawn from one JavaRandom(seed); ';' between permutations, ',' within\n");
  const std::vector<std::vector<int>> plans = {
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
      {300, 25, 7, 4, 300, 23, 6},
      {16, 16, 8, 4, 2, 1},
      {267, 40, 12, 267, 38, 11, 4, 267},
      {1000, 128, 64, 1000, 127, 63}};
  for (int seed : {0, 12345, 42}) {
    for (const auto& sizes : plans) {
      JavaRandom random(seed);
      std::fprintf(f, "%d\t", seed);
      for (size_t s = 0; s < sizes.size(); s++) std::fprintf(f, s ? ",%d" : "%d", sizes[s]);
      std::fprintf(f, "\t");
      for (size_t s = 0; s < sizes.size(); s++) {
        IVector p = Arrays2::generateRandomPermutation(sizes[s], random);
        if (s) std::fprintf(f, ";");
        for (size_t k = 0; k < p.size(); k++) std::fprintf(f, k ? ",%d" : "%d", p[k]);
      }
      std::fprintf(f, "\n");
    }
  }
  std::fclose(f);
  return 0;
}
