# tests/fixtures

Small fixtures committed to the repository, which `cargo test` checks on every change (docs/PLAN.md A1, A4 and A5b): tier-T data exported from the workflow's test data, and R-generated vectors for `rcompat` (stage S0).

- **Size:** about 20 MB in total at most. Larger fixtures go to `data/` (gitignored) and are listed in `validation/manifest.toml`.
- **Layout:** one directory per crate or tool, for example `rcompat/` (R vectors, M1), `doubletfinder/<object>/` (inputs and stage references for the five test objects, M2) `clustering/` (Step 2's K0 vectors and tier-T references, M8; [its README](clustering/README.md)) and `neighbors/` (the neighbour stage's constructed edge cases, M9; [its README](neighbors/README.md)).
- **Provenance:** every file is listed in `validation/manifest.toml` with its sha256, size, the script that made it, and the env and profile it was made in. Regenerate fixtures from those recipes; never edit them by hand.
- **Profile:** fixtures that depend on BLAS are made in the portable profile (`OPENBLAS_CORETYPE=Haswell`, one thread), so tests can check them exactly on any x86-64 CPU with AVX2 (A5b).
- **From a crate's tests:** `concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures")`.

## rcompat/

R 4.5.3's results for the `rcompat` crate (stage S0), made by [`validation/R/rcompat_vectors.R`](../../validation/R/rcompat_vectors.R) and checked by `crates/rcompat/tests/r_vectors.rs`. One tab-separated file per family (M1: `unif`, `random_seed`, `random_seed_after`, `sample`, `rnorm`, `qnorm`, `round`, `pow`, `seq_by`, `seq_len`, `order`, `summary`, `format`; M3: `strtod` (`as.numeric()` of strings, R's own parser, including strings it does not round correctly), `fft`, `dnorm`, `approx`); lines starting with `#` name the columns. Doubles are written as their IEEE-754 bit patterns (16 hex digits), so the comparison is bitwise; lists are comma-separated, `-` for an empty one.

## doubletfinder/

The five tier-T objects (the workflow's test data, plan B4), one directory each: `soupx_tenx_test`, `soupx_emptydrops_test`, `decontx_raw_tenx_test`, `decontx_raw_emptydrops_test`, `cellbender_fromraw_test`. Each is a copy of files under `data/doubletfinder/<object>/` (milestone M2; recipes and checks in [`docs/validation/doubletfinder_m2.md`](../../docs/validation/doubletfinder_m2.md)), made in the portable profile:

- `input/`: what DoubletFinder reads, exported by [`export_doubletfinder.R`](../../validation/R/export_doubletfinder.R): `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` (10x Matrix Market, entries exactly as stored in the object's RNA `counts` layer, **explicit zeros included**: the SoupX objects have them, and sctransform's `min_variance` is `median(umi@x)`, so a reader must keep them), `seurat_clusters.tsv` and `seurat_clusters_levels.txt` (the factor's level order), `export.json`.
- `saved/`: the saved workflow output's compared fields: `cells.tsv` (kept barcodes and their pANN), `doublets.txt` (the doublet set), `saved.json` (column names, pK, nExp).
- `reference/`: the direct run's `cells.tsv` (barcode, pANN, classification for every cell), `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json` (compare its `result` section; the `run` section holds timings, node and job).
- `stages/`: the instrumented run's stage arrays S1–S7 (`.npy`), `index.tsv` (which rows of which file belong to which context and stage), `rng.npy`/`rng.tsv` (`.Random.seed` checkpoints) and `stages.json`. Conventions: [`validation/R/lib/doubletfinder_hooks.R`](../../validation/R/lib/doubletfinder_hooks.R). Not copied, to stay within the size budget (they remain in `data/` and in `index.tsv`): the per-context matrices `S2d_y.<ctx>.npy`, `S2d_scale_data.<ctx>.npy` and `S4_positions.<ctx>.npy`, and `S6_bkde.npy`.

The whole tool is checked on them end to end (milestone M6, `crates/doubletfinder/tests/m6_tier_t.rs` and `crates/scrust/tests/cli.rs`): `scrust doubletfinder` from `input/` must reproduce `reference/` byte for byte (the digest of the final `.Random.seed` included), every array, RNG checkpoint and scalar in `stages/`, and `saved/`. Its `singlets/` (milestone M6b) must be `input/`'s counts restricted to `saved/cells.tsv`'s kept cells, every stored entry as it is, with the input's feature lines. No fixture was added for this. The saved objects' own counts are compared on `data/`, where `export_doubletfinder_singlets.R` exported them.

Doubles in the TSV and JSON files have 17 significant digits, so Rust's `str::parse::<f64>` gives back the exact value; `crates/doubletfinder/tests/m2_fixtures.rs` checks this against the `.npy` copies.

## clustering/

Step 2's committed reference data (milestone M8; [`docs/validation/clustering_m8.md`](../../docs/validation/clustering_m8.md)): the K0 vectors of the clustering chain's generators (`k0/`: Seurat's `JavaRandom`, RcppAnnoy's `Kiss64Random`, R's draw in GroupSingletons), and, for seven tier-T objects (the two base objects and the five DoubletFinder inputs, whose counts are the `doubletfinder/` fixtures above), the expected labels and scores of the portable profile and of the saved object and the instrumented optimizer's per-start results; for the two base objects also their counts and the stage references from R's own embeddings (the Annoy index, the neighbours, both graphs, silhouette's and purity's intermediates). Some files are gzipped (`gzip -n -9`) to stay near the size budget. Details in [clustering/README.md](clustering/README.md); the other 20 tier-T objects and the tier-P objects are in `data/clustering/`, listed in `validation/manifest.toml`.

## neighbors/

The neighbour stage's constructed edge cases (milestone M9; [`docs/validation/clustering_m9.md`](../../docs/validation/clustering_m9.md)): 25 inputs that reach Annoy's and FindNeighbors' rare branches, and the size and SHA-256 of every file R wrote for them with each of the three metrics; R's files themselves are in `data/neighbors/edge/`. Details in [neighbors/README.md](neighbors/README.md). The real objects' neighbour references are M8's `clustering/` fixtures.
