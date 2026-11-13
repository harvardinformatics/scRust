// Kiss64Random reference vectors (docs/PLAN.md C2 and C5, stage K0; milestone M8) from RcppAnnoy
// 0.0.23's own header, kissrandom.h (Annoy 1.17.3), included unmodified from the R env's installed
// RcppAnnoy: the generator Annoy's index build draws from. Seurat never calls setSeed, and
// RcppAnnoy builds single-threaded, so each build uses one Kiss64Random(default_seed + 0) for all
// of its trees (annoylib.h: thread_build's `Random _random(_seed + thread_idx)`).
//
// Built by validation/cpp/build.sh with the R env's compiler and Makeconf flags; run as
//   data/clustering/bin/kiss64_vectors <out dir>
// Writes, tab-separated, lines starting with "#" naming the columns:
//   kiss64_kiss.tsv   seed, count, the first `count` kiss() values (uint64, decimal)
//   kiss64_index.tsv  seed, n, count, the first `count` index(n) values (kiss() % n)
//   kiss64_flip.tsv   seed, count, the first `count` flip() values (kiss() & 1)
//   kiss64_mixed.tsv  seed, a sequence of calls as Annoy's two_means and split make them:
//                     "i<n>" for index(n), "f" for flip(), and the values returned
// Seeds: the default seed (1234567890987654321) and a few others, including 0's neighbours.

// kissrandom.h uses size_t without including <cstddef>; annoylib.h includes it after other headers.
#include <cinttypes>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include <kissrandom.h>

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: kiss64_vectors <out dir>\n");
    return 2;
  }
  std::string dir = argv[1];
  const std::vector<uint64_t> seeds = {Kiss64Random::default_seed, Kiss64Random::default_seed + 1,
                                       1, 42, 18446744073709551615ULL};
  // index(n): the default seed (the one every build uses) and one other
  const std::vector<uint64_t> index_seeds = {Kiss64Random::default_seed, 42};
  const int count = 200;

  FILE* f = std::fopen((dir + "/kiss64_kiss.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tcount\tkiss() values, comma-separated\n");
  for (uint64_t seed : seeds) {
    Kiss64Random r(seed);
    std::fprintf(f, "%" PRIu64 "\t%d\t", seed, count);
    for (int k = 0; k < count; k++) std::fprintf(f, k ? ",%" PRIu64 : "%" PRIu64, r.kiss());
    std::fprintf(f, "\n");
  }
  std::fclose(f);

  const std::vector<size_t> ns = {1, 2, 3, 7, 10, 20, 30, 32, 33, 100, 267, 300, 1000, 3108, 11831,
                                  26598, 59323, 79409, 1000000, 4294967296ULL, 4294967297ULL};
  f = std::fopen((dir + "/kiss64_index.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tn\tcount\tindex(n) values from a fresh Kiss64Random(seed), comma-separated\n");
  for (uint64_t seed : index_seeds) {
    for (size_t n : ns) {
      Kiss64Random r(seed);
      std::fprintf(f, "%" PRIu64 "\t%zu\t%d\t", seed, n, count);
      for (int k = 0; k < count; k++) std::fprintf(f, k ? ",%zu" : "%zu", r.index(n));
      std::fprintf(f, "\n");
    }
  }
  std::fclose(f);

  f = std::fopen((dir + "/kiss64_flip.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tcount\tflip() values, comma-separated\n");
  for (uint64_t seed : seeds) {
    Kiss64Random r(seed);
    std::fprintf(f, "%" PRIu64 "\t%d\t", seed, count);
    for (int k = 0; k < count; k++) std::fprintf(f, k ? ",%d" : "%d", r.flip());
    std::fprintf(f, "\n");
  }
  std::fclose(f);

  // One stream through calls of mixed kinds, as a build makes them: per split attempt two
  // index(n) for the starting pair, then 200 index(n) for the two-means iterations, then flips.
  f = std::fopen((dir + "/kiss64_mixed.tsv").c_str(), "w");
  if (!f) return 1;
  std::fprintf(f, "# seed\tcalls (i<n> = index(n), f = flip()), comma-separated\tvalues, comma-separated\n");
  for (uint64_t seed : {Kiss64Random::default_seed, (uint64_t)42}) {
    Kiss64Random r(seed);
    std::vector<std::string> calls;
    std::vector<uint64_t> vals;
    const std::vector<size_t> sizes = {300, 299, 160, 140, 79409, 37, 2};
    for (size_t n : sizes) {
      for (int k = 0; k < 2; k++) { calls.push_back("i" + std::to_string(n)); vals.push_back(r.index(n)); }
      for (int k = 0; k < 200; k++) { calls.push_back("i" + std::to_string(n)); vals.push_back(r.index(n)); }
      for (int k = 0; k < 5; k++) { calls.push_back("f"); vals.push_back(r.flip()); }
    }
    std::fprintf(f, "%" PRIu64 "\t", seed);
    for (size_t k = 0; k < calls.size(); k++) std::fprintf(f, k ? ",%s" : "%s", calls[k].c_str());
    std::fprintf(f, "\t");
    for (size_t k = 0; k < vals.size(); k++) std::fprintf(f, k ? ",%" PRIu64 : "%" PRIu64, vals[k]);
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  return 0;
}
