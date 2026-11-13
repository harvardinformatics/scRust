# DoubletFinder in Rust (Step 1; milestones M6 and M6b)

M6 puts the pieces of milestones M1 to M5 together into one command, `scrust doubletfinder`: DoubletFinder 2.0.6 as the workflow's `doubletfinder.R` calls it (docs/PLAN.md B1), from a 10x counts directory and the cells' cluster labels to the doublet calls, every stage computed by Rust. This is the report the plan names for Step 1. M6b adds `singlets/`, the counts without the doublets as a new 10x directory, on top of the outputs M6 reports on ([M6b: the singlets' counts](#m6b-the-singlets-counts)). The full panel of 150 production objects (M7) comes next; substituting `scrust` into the Snakemake workflow is a later, separate change.

## Outcome

Every criterion of plan B6's row M6 and of B5 is met.

| Done when (plan B6 M6, B5, and the milestone's brief) | Result |
|---|---|
| The command runs all of B1 itself | **Yes.** `scrust doubletfinder --counts <10x dir> --clusters <tsv> --out <dir>` runs the whole tool from the counts, the cluster labels and the seed:<br>- the sweep subsample and the doublet draws, from Rust's generator in R's stream;<br>- SCTransform, RunPCA, the neighbours and pANN in all six sweep iterations and the final run;<br>- `summarizeSweep`, `find.pK`, nExp with the homotypic adjustment, and the calls.<br>It writes `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json`, and takes the plan's options. `--from-stages` stays as M3's teacher-forcing mode. |
| Instrumented runs, all 24, each under its own OpenBLAS profile: S1-S7 and every RNG checkpoint identical; the only inputs the counts, the labels and the seed | **Yes, all 24.** Every block of every stage array R dumped, S1 to S7, is identical: `scale.data` and the residual matrices wherever R dumped them, and the S4 order positions on tier T. So are every `.Random.seed` checkpoint (2,646 of 2,646) and their sequence, every scalar, the final files and the `result` section: 13,396 checks, 1.58 billion values, none differing. |
| Direct reference runs, all 23: `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` byte-identical; `summary.json` identical but for fields that cannot match | **Yes, all 23.** The three files are byte-identical, and all 20 fields of the `result` section are identical, the digest of R's final `.Random.seed` included. What cannot match is the `run` section, a record of the run (listed [below](#summaryjson-what-is-compared-and-what-cannot-match)); its OpenBLAS kernel is compared and matches. |
| Saved workflow outputs (T and P): the doublet set, singlet pANN, pK and nExp as R writes them | **Yes**, for all 10 objects in both profiles (20 Rust runs): the same doublet sets, every kept cell's pANN identical as R wrote it, the same pK labels, nExp and column names. |
| Settled and tested | `tolerances.toml` settles the end-to-end, saved-output and performance rows. `cargo test` runs the tier-T end-to-end checks (`m6_tier_t.rs`: all five objects, every stage, RNG checkpoint and scalar the fixtures hold, the reference files, the saved outputs, thread counts 1 to 7; `cli.rs`: the command and `scrust validate`). No new fixtures were needed: the tests use M2's, already in the manifest. `cargo test --workspace` passes in a fresh clone with its own target directory. |
| Thread count doesn't matter | **Yes.** At 1 and at 32 threads the outputs of the two largest objects are byte-identical to each other and to R's; the validation ran 8 threads; tier T is checked at 1, 2, 3 and 7. |
| Benchmark: time and peak RSS at 1 and 32 threads against the direct R run; less time and memory than R; sscovelli at most 10 GiB | **Yes** ([Benchmark](#benchmark)).<br>- soupx_emptydrops_sscovelli: R 36:15 and 126.9 GiB; Rust 6:21 and 4.41 GiB on one thread (5.7 times faster, a 29th of the memory), 43.5 s and 4.51 GiB on 32 (50 times faster).<br>- cellbender_fromraw_ghirsutum2: R 33:53 and 150.4 GiB; Rust 8:58 and 4.02 GiB on one thread (3.8 times faster, a 37th of the memory), 60.9 s and 4.08 GiB on 32 (33 times faster).<br>Every Rust output is identical to R's. |
| Report and README | This report; the README's status and usage. |

R's peaks are `/usr/bin/time`'s 133,113,204 and 157,716,268 KiB, i.e. 126.9 and 150.4 GiB (136.3 and 161.5 GB); the brief's "133.1 GB" and "157.7 GB" read those KiB as kB. This report uses GiB throughout.

## Step 1 at a glance

| Milestone | What it established | Report |
|---|---|---|
| M1 | `rcompat`: R 4.5.3's random numbers (Mersenne-Twister, `sample`, `rnorm`), `round`, `seq`, `order`, long-double `sum`/`mean`, number formatting, bit for bit | `crates/rcompat`, `tests/fixtures/rcompat/` |
| M2 | the R reference data: 10 objects (tier T: the five workflow test objects; tier P: five production objects), 23 direct and 24 instrumented runs in the portable and production OpenBLAS profiles, each reproducing its saved workflow output; stage dumps S1-S7; the R-vs-R floor | [doubletfinder_m2.md](doubletfinder_m2.md) |
| M3 | the core after PCA (neighbours, pANN, `summarizeSweep`, `find.pK`, nExp, calls), exact on R's own embeddings, without the N x N distance matrix | [doubletfinder_m3.md](doubletfinder_m3.md) |
| M4 | SCTransform v2 (sctransform 0.4.3, glmGamPoi 1.22.0's fits, R's `nlminb`), bit-identical in all 168 SCTransform calls, teacher-forced and chained | [doubletfinder_m4.md](doubletfinder_m4.md) |
| M5 | RunPCA (irlba 2.4.1), bit-identical in all 168 PCA calls, signs and the generator's state included | [doubletfinder_m5.md](doubletfinder_m5.md) |
| M6 | the whole tool from the counts, the cluster labels and the seed: every stage, every RNG checkpoint and every output identical to R's; the benchmark | this report |

## What was built

| Crate | M6 contents |
|---|---|
| `doubletfinder` | `chain::RustSource`: the embeddings of every context computed from the counts. It draws the sweep subsample and the doublet pairs from R's random-number stream, builds DoubletFinder's input (real cells, then doublets), and runs SCTransform and `RunPCA(npcs = 10)` on it, in R's order with R's `set.seed` calls. `run::run` (M3) does the rest, as before. `record::Recorder`: a run's stage dumps S1 to S7, `.Random.seed` checkpoints and per-context scalars, in the instrumented R runs' layout and levels. `input`: the counts directory, the cluster labels and the factor's levels. `validate`: every block of every stage array, every RNG checkpoint and their sequence, every scalar, the OpenBLAS kernel; `compare_saved` against the saved workflow output. `sct::sct_input` builds the doublet-augmented input in one allocation of the exact size. |
| `sctransform` | `start` and `finish`: `sctransform()` split where it stops reading the counts (after `set.seed` and S2a), so the caller can free a large input before the fits. The residual stage (S2d) keeps the rows of the running best 3,000 genes as slots of one buffer, which becomes `scale.data` in place: reordered into R's order (`gather_rows`), clipped and centred. Before, the rows were separate allocations copied into a new matrix, which held `scale.data` twice for a moment (below). The arithmetic is unchanged in both. |
| `scio` | `sha256` (standard library only), for the digest of the final `.Random.seed` that R's runs record; `npy::write_f64_rows`, a row-major matrix streamed to R's column-major `.npy` layout. |
| `scrust` | `scrust doubletfinder --counts <dir> --clusters <tsv> --out <dir>` runs the whole tool, with the plan's options (`--seed`, `--pcs`, `--pn`, `--nexp-rate`, `--pk auto\|<value>`, `--threads`) plus `--cluster-levels` and `--dump-stages none\|light\|standard\|full`. `--from-stages` stays as M3's teacher-forcing mode. `scrust validate doubletfinder` gains `--saved`. |

No external crate was added: the workspace uses the Rust standard library and, through FFI, the OpenBLAS (with its LAPACK) that R uses.

**The command.**

```bash
scrust doubletfinder --counts <10x dir> --clusters <seurat_clusters.tsv> --out <dir> [--threads N]
```

- **Inputs.** `--counts` is a 10x Matrix Market directory (`matrix.mtx(.gz)`, `barcodes.tsv(.gz)`, `features.tsv(.gz)`; stored zeros are kept, as R keeps them). `--clusters` is a TSV with a header and the barcode and cluster label of every cell, in the matrix's column order. The factor's level order, which the homotypic adjustment tabulates, comes from `--cluster-levels`, else from `seurat_clusters_levels.txt` next to the clusters file, else from the order of first appearance.
- **Outputs.** `cells.tsv` (barcode, pANN and classification of every cell), `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json`. Its `result` section holds pK and its label as R writes it, nExp_poi, the homotypic proportion, nExp, the doublets made, k, the column names R writes, and the digest of the final `.Random.seed`. Its `run` section holds parameters, versions, per-phase times and peak memory, and the environment fingerprint.
- **Exactness.** SCTransform's and RunPCA's products are the OpenBLAS calls R makes, so the outputs are R's bit for bit under the kernel R ran with: `OPENBLAS_CORETYPE=Haswell` (env.sh's default, the portable profile), or the variable unset on Sapphire Rapids (the production profile, the kernel of the workflow's runs).

## How the run follows R's random-number stream

One generator, `RRng`, passes through the whole run in R's order (plan B2):

1. `set.seed(12345)`.
2. `paramSweep`: above 10,000 cells, `sample(1:n, 10000)`, a partial Fisher-Yates in draw order.
3. For each pN in turn (0.05 to 0.3):
   1. two `sample(real.cells, n_doublets, replace = TRUE)` draws;
   2. `SCTransform()`: three `set.seed(1448145)`, the cell sample (above 5,000 cells), the density-weighted gene sample (above 2,000 fitting genes), and two draws per fitting bin (`future.seed`);
   3. `RunPCA()`: `set.seed(42)`, then irlba's start vector, `rnorm(nfeatures)`.

   The next pN's doublets continue from the state this PCA leaves.
4. `set.seed(12345)`.
5. `doubletFinder`: the same three steps on all cells.

Nothing is read from an R run: the inputs are the counts, the labels and the seed. The instrumented R runs record the state at 102 to 120 points per run (`rng.tsv`); the Rust run records the same points under the same labels, and every one is compared.

**A consequence worth knowing.** Every RunPCA resets the generator to `set.seed(42)` and draws `rnorm(nfeatures)`. Where the previous context's PCA had the same number of features, the doublets of pN 0.1 to 0.3 are therefore drawn from the same state. That holds for every P object, where every sweep context has 3,000 features (R's own checkpoints: one distinct `before_doublets` state for pN 0.1 to 0.3 in each of the five P objects). The first draws of `real.cells1` then coincide across those five iterations. This is DoubletFinder's behaviour, reproduced, not changed.

## How the comparison works

- **The Rust runs.** For each object, [`validation/slurm/m6_doubletfinder.sbatch`](../../validation/slurm/m6_doubletfinder.sbatch) makes one Rust run per OpenBLAS profile its reference runs use (job 60007108, 8 threads, binary `scrust-m6-bc1b48e`). The portable profile sets `OPENBLAS_CORETYPE=Haswell`; the production profile leaves it unset, so that OpenBLAS detects SapphireRapids on the `sapphire` node. Each run uses `--dump-stages` at the level of the profile's instrumented R run: full for tier T, standard for the small and medium objects, light for the large ones. The core is checked with `scrust fingerprint` before each run.
- **Comparisons** (`scrust validate doubletfinder`):
  - against each instrumented R run of the profile: the final files, the `result` section, every block of every stage array R dumped, every `.Random.seed` checkpoint and their sequence, every scalar, S6 and S7 recomputed from R's S5, and S1 redrawn from R's own states;
  - against each direct R run of the profile: the final files and the `result` section;
  - against the saved workflow output: the doublet set, the kept cells' pANN and the column names.
- **One Rust run per profile, not per CPU type.** The portable R runs on Cascade Lake are compared with the portable Rust run made on Sapphire Rapids: the kernel, not the CPU, sets the bits. M2 showed R's portable runs identical on both CPU types.
- **Bit for bit.** Every comparison is bitwise: doubles by their bit patterns, indices and counts by value, files byte for byte, scalars from their `%.17g` text. A check is one array of one context (S6's 198 summary contexts gather into one check per array), one RNG checkpoint, one scalar of one context, one file or one `result` field.

The tier-T end-to-end checks run in `cargo test` ([`crates/doubletfinder/tests/m6_tier_t.rs`](../../crates/doubletfinder/tests/m6_tier_t.rs); the command itself in [`crates/scrust/tests/cli.rs`](../../crates/scrust/tests/cli.rs)), on the committed fixtures, under `OPENBLAS_CORETYPE=Haswell`.


## Results

Validation label `m6-bc1b48e` (`data/doubletfinder/<object>/rust/m6-bc1b48e/<profile>/`), job 60007108, binary `data/doubletfinder/logs/bin/scrust-m6-bc1b48e` (sha256 `ef3612dc7d08035f`), built from commit bc1b48e; Intel Xeon Platinum 8480+/8480CL (`sapphire`), 8 threads. Two earlier full runs gave the same results: `m6-7a8d7c5` (job 60002845, before the validator compared the OpenBLAS kernels) and `m6-f85b69a` (job 60005537, before `scale.data` was assembled in place). In all three, every Rust output is byte-identical: 1,496 files (the tables, `rng.npy`, every stage array but the per-context matrices) compared between runs.

### Instrumented runs

One row per instrumented R run. Each stage cell gives the checks identical and the checks made in that stage: one per array and context, plus that stage's scalars and, for S1, S6 and S7, the checks recomputed from R's own draws or S5. *RNG checkpoints* are R's `.Random.seed` states (102 to 120 per run), each against Rust's of the same context and label. *n/a* are R's structural checks of its own containers (five per context, four in light runs), which do not apply to Rust. *missing* counts anything R dumped that the Rust run did not: none. Not in the table but in the totals: one check per comparison of the OpenBLAS kernel each run recorded. Every run's report (`vs_<R run>.md`) lists every check.

| object | R run | dumps | S1 | S2a | S2b | S2c | S2d | S3 | S4 | S5 | S6 | S7 | RNG checkpoints | files | result | values compared | differing | n/a | missing | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_test | instrumented-portable-cascadelake | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 5,287,444 | 0 | 35 | 0 | yes |
| soupx_tenx_test | instrumented-portable-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 5,287,444 | 0 | 35 | 0 | yes |
| soupx_tenx_test | instrumented-production-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 5,287,444 | 0 | 35 | 0 | yes |
| soupx_emptydrops_test | instrumented-portable-cascadelake | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,700,445 | 0 | 35 | 0 | yes |
| soupx_emptydrops_test | instrumented-portable-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,700,445 | 0 | 35 | 0 | yes |
| soupx_emptydrops_test | instrumented-production-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,700,445 | 0 | 35 | 0 | yes |
| decontx_raw_tenx_test | instrumented-portable-cascadelake | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 5,057,213 | 0 | 35 | 0 | yes |
| decontx_raw_tenx_test | instrumented-portable-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 5,057,213 | 0 | 35 | 0 | yes |
| decontx_raw_tenx_test | instrumented-production-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 5,057,213 | 0 | 35 | 0 | yes |
| decontx_raw_emptydrops_test | instrumented-portable-cascadelake | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 4,570,495 | 0 | 35 | 0 | yes |
| decontx_raw_emptydrops_test | instrumented-portable-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 4,570,495 | 0 | 35 | 0 | yes |
| decontx_raw_emptydrops_test | instrumented-production-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (102) | = (3) | = (20) | 4,570,495 | 0 | 35 | 0 | yes |
| cellbender_fromraw_test | instrumented-portable-cascadelake | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,939,730 | 0 | 35 | 0 | yes |
| cellbender_fromraw_test | instrumented-portable-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,939,730 | 0 | 35 | 0 | yes |
| cellbender_fromraw_test | instrumented-production-sapphirerapids | full | = (34) | = (98) | = (42) | = (63) | = (63) | = (56) | = (35) | = (20) | = (4) | = (3) | = (106) | = (3) | = (20) | 4,939,730 | 0 | 35 | 0 | yes |
| soupx_tenx_SRR18256555 | instrumented-portable-cascadelake | standard | = (34) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 85,242,608 | 0 | 35 | 0 | yes |
| soupx_tenx_SRR18256555 | instrumented-portable-sapphirerapids | standard | = (34) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 85,242,608 | 0 | 35 | 0 | yes |
| soupx_tenx_SRR18256555 | instrumented-production-sapphirerapids | standard | = (34) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 85,242,608 | 0 | 35 | 0 | yes |
| soupx_tenx_10xmousebrain | instrumented-portable-sapphirerapids | standard | = (35) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 277,640,316 | 0 | 35 | 0 | yes |
| soupx_tenx_10xmousebrain | instrumented-production-sapphirerapids | standard | = (35) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 277,640,316 | 0 | 35 | 0 | yes |
| soupx_emptydrops_4500871lane | instrumented-portable-sapphirerapids | standard | = (35) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 337,266,987 | 0 | 35 | 0 | yes |
| soupx_emptydrops_4500871lane | instrumented-production-sapphirerapids | standard | = (35) | = (126) | = (42) | = (63) | = (56) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 337,266,987 | 0 | 35 | 0 | yes |
| soupx_emptydrops_sscovelli | instrumented-portable-sapphirerapids | light | = (35) | = (126) | = (42) | = (63) | = (49) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 11,481,122 | 0 | 28 | 0 | yes |
| cellbender_fromraw_ghirsutum2 | instrumented-portable-sapphirerapids | light | = (35) | = (126) | = (42) | = (63) | = (49) | = (56) | = (28) | = (20) | = (4) | = (3) | = (120) | = (3) | = (20) | 13,658,206 | 0 | 28 | 0 | yes |

### Direct runs

One row per direct R run, against the Rust run of its profile: the three files byte for byte and the 20 fields of the `result` section (the OpenBLAS kernel is compared too).

| object | R run | Rust run | cells.tsv | sweep_stats.tsv | bcmvn.tsv | result fields | pass |
|---|---|---|---|---|---|---|---|
| soupx_tenx_test | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| soupx_tenx_test | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| soupx_tenx_test | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| soupx_emptydrops_test | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| soupx_emptydrops_test | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| soupx_emptydrops_test | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| decontx_raw_tenx_test | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| decontx_raw_tenx_test | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| decontx_raw_tenx_test | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| decontx_raw_emptydrops_test | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| decontx_raw_emptydrops_test | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| decontx_raw_emptydrops_test | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| cellbender_fromraw_test | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| cellbender_fromraw_test | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| cellbender_fromraw_test | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| soupx_tenx_SRR18256555 | direct-portable-cascadelake | portable | = | = | = | = (20) | yes |
| soupx_tenx_SRR18256555 | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| soupx_tenx_SRR18256555 | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| soupx_tenx_10xmousebrain | direct-portable-sapphirerapids | portable | = | = | = | = (20) | yes |
| soupx_tenx_10xmousebrain | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| soupx_emptydrops_4500871lane | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| soupx_emptydrops_sscovelli | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |
| cellbender_fromraw_ghirsutum2 | direct-production-sapphirerapids | production | = | = | = | = (20) | yes |

### Saved workflow outputs

Each Rust run against the saved workflow output of its object: the doublet set (the input barcodes minus those the workflow kept), the pANN of every kept cell as R wrote it (17 significant digits, in order), pK and nExp, and the two column names they come from.

| object | Rust run | doublet set | singlet pANN | pK | nExp | columns | pass |
|---|---|---|---|---|---|---|---|
| soupx_tenx_test | portable | = (24) | = (276) | 0.29 | 24 | `pANN_0.25_0.29_24` | yes |
| soupx_tenx_test | production | = (24) | = (276) | 0.29 | 24 | `pANN_0.25_0.29_24` | yes |
| soupx_emptydrops_test | portable | = (30) | = (237) | 0.28 | 30 | `pANN_0.25_0.28_30` | yes |
| soupx_emptydrops_test | production | = (30) | = (237) | 0.28 | 30 | `pANN_0.25_0.28_30` | yes |
| decontx_raw_tenx_test | portable | = (27) | = (273) | 0.09 | 27 | `pANN_0.25_0.09_27` | yes |
| decontx_raw_tenx_test | production | = (27) | = (273) | 0.09 | 27 | `pANN_0.25_0.09_27` | yes |
| decontx_raw_emptydrops_test | portable | = (30) | = (237) | 0.22 | 30 | `pANN_0.25_0.22_30` | yes |
| decontx_raw_emptydrops_test | production | = (30) | = (237) | 0.22 | 30 | `pANN_0.25_0.22_30` | yes |
| cellbender_fromraw_test | portable | = (31) | = (247) | 0.23 | 31 | `pANN_0.25_0.23_31` | yes |
| cellbender_fromraw_test | production | = (31) | = (247) | 0.23 | 31 | `pANN_0.25_0.23_31` | yes |
| soupx_tenx_SRR18256555 | portable | = (420) | = (2,688) | 0.005 | 420 | `pANN_0.25_0.005_420` | yes |
| soupx_tenx_SRR18256555 | production | = (420) | = (2,688) | 0.005 | 420 | `pANN_0.25_0.005_420` | yes |
| soupx_tenx_10xmousebrain | portable | = (1,646) | = (10,185) | 0.07 | 1,646 | `pANN_0.25_0.07_1646` | yes |
| soupx_tenx_10xmousebrain | production | = (1,646) | = (10,185) | 0.07 | 1,646 | `pANN_0.25_0.07_1646` | yes |
| soupx_emptydrops_4500871lane | portable | = (3,359) | = (23,239) | 0.02 | 3,359 | `pANN_0.25_0.02_3359` | yes |
| soupx_emptydrops_4500871lane | production | = (3,359) | = (23,239) | 0.02 | 3,359 | `pANN_0.25_0.02_3359` | yes |
| soupx_emptydrops_sscovelli | portable | = (9,909) | = (69,500) | 0.19 | 9,909 | `pANN_0.25_0.19_9909` | yes |
| soupx_emptydrops_sscovelli | production | = (9,909) | = (69,500) | 0.19 | 9,909 | `pANN_0.25_0.19_9909` | yes |
| cellbender_fromraw_ghirsutum2 | portable | = (6,769) | = (52,554) | 0.03 | 6,769 | `pANN_0.25_0.03_6769` | yes |
| cellbender_fromraw_ghirsutum2 | production | = (6,769) | = (52,554) | 0.03 | 6,769 | `pANN_0.25_0.03_6769` | yes |

### Totals

| comparisons | runs | passed | checks identical | values compared | values differing | RNG checkpoints identical | not applicable | missing |
|---|---|---|---|---|---|---|---|---|
| instrumented | 24 | 24 | 13,396 of 13,396 | 1,584,347,739 | 0 | 2,646 of 2,646 | 826 | 0 |
| direct | 23 | 23 | 552 of 552 | 207,662 | 0 | - | 23 | 0 |
| saved | 20 | 20 | 120 of 120 | 363,442 | 0 | - | 0 | 0 |

### `summary.json`: what is compared and what cannot match

`summary.json` has two sections, in R's layout and Rust's alike.

**`result` is compared field by field, as text.** All 20 fields are identical in every comparison:

- `n_cells`, `n_genes`, `seed`, `PCs`, `pN`;
- the sweep's `sweep_pK_labels` and `sweep_pN_labels`, `bcmvn_max_row`;
- `pK` and `pK_label`, `nExp_poi`, `homotypic_prop`, `nExp`;
- `n_doublets`, `n_cells_with_doublets`, `k`;
- `pANN_column` and `classification_column`, `n_doublet_calls`;
- `rng_end_sha256`, the SHA-256 of R's final `.Random.seed` as text.

**`run` records the run, so it cannot match.**

- R's direct runs write:
  - `kind`, `script`, `input`, `input_source`, `outdir`, `start_time`, `slurm_job_id`, `profile`;
  - `load_input_s`, `peak_rss_after_load_kib`, `phases` (elapsed time and peak RSS after `paramSweep`, `summarizeSweep`, `find.pK`, the cleanup and `doubletFinder`), `peak_rss_kib`;
  - the environment `fingerprint`, which adds R's own view (R and package versions, BLAS and LAPACK library paths, the BLAS md5, RNG kind, locale, `sizeof(long double)`).
- Rust writes:
  - `kind`, `mode`, `scrust_version`, `reproduces` (the versions of DoubletFinder, sctransform and glmGamPoi, and irlba it reproduces), `counts`, `clusters`, `cluster_levels`, `from_stages`, `outdir`;
  - `parameters` (seed, PCs, pN, nExp rate, pK), `threads`, `dump_stages`, `start_time`, `slurm_job_id`, `elapsed_s`, `phases` (time and peak RSS after loading, after each context's SCTransform and RunPCA, each pN, each sweep step, the final run and the writing), `peak_rss_kib`;
  - `fingerprint` (`scrust fingerprint`: node, kernel, CPU model and flags, glibc, OpenBLAS core, configuration, threads and library, the OpenBLAS and SLURM environment variables).

Paths, times, memory, job numbers and versions differ by nature.

The fingerprint's BLAS-relevant fields agree: the same OpenBLAS build (md5 `d77a762c...`, the R env's and the Rust env's files), the same core (compared: `environment / OpenBLAS kernel`, identical in all 47 comparisons with R runs), one OpenBLAS thread, and the same CPU model and glibc wherever the CPU type is the same.

### Threads

The results do not depend on the number of threads:

- the benchmark's runs of the two largest objects at 1 and 32 threads are byte-identical to each other and to R's direct runs;
- the validation ran 8 threads, with the same results;
- `cargo test` checks tier T at 1, 2, 3 and 7 threads (the final `.Random.seed` included).

Where the work is split:

- SCTransform: the fitting bins' genes and the residuals' genes;
- RunPCA: irlba's two large products, only where OpenBLAS's kernels keep R's bits (M5);
- the neighbour search: the cells;
- the sweep summary: the (pN, pK) pairs.

Every result is placed by index, never by finishing order.

### Validation runs

The Rust runs of the validation, 8 threads, with stage dumps: standard and full dumps also write `scale.data` per context (full also the residual matrices), and the dumps' small arrays are held until the end, so these peaks are a little above the benchmark's.

| object | profile | cells | with doublets | k | dumps | threads | wall | peak RSS GiB |
|---|---|---|---|---|---|---|---|---|
| soupx_tenx_test | portable | 300 | 400 | 116 | full | 8 | 0.7 s | 0.04 |
| soupx_tenx_test | production | 300 | 400 | 116 | full | 8 | 0.6 s | 0.04 |
| soupx_emptydrops_test | portable | 267 | 356 | 100 | full | 8 | 0.7 s | 0.04 |
| soupx_emptydrops_test | production | 267 | 356 | 100 | full | 8 | 0.5 s | 0.03 |
| decontx_raw_tenx_test | portable | 300 | 400 | 36 | full | 8 | 0.6 s | 0.04 |
| decontx_raw_tenx_test | production | 300 | 400 | 36 | full | 8 | 0.5 s | 0.04 |
| decontx_raw_emptydrops_test | portable | 267 | 356 | 78 | full | 8 | 0.6 s | 0.04 |
| decontx_raw_emptydrops_test | production | 267 | 356 | 78 | full | 8 | 0.5 s | 0.03 |
| cellbender_fromraw_test | portable | 278 | 371 | 85 | full | 8 | 0.6 s | 0.04 |
| cellbender_fromraw_test | production | 278 | 371 | 85 | full | 8 | 0.5 s | 0.04 |
| soupx_tenx_SRR18256555 | portable | 3,108 | 4,144 | 21 | standard | 8 | 14.2 s | 0.42 |
| soupx_tenx_SRR18256555 | production | 3,108 | 4,144 | 21 | standard | 8 | 14.0 s | 0.42 |
| soupx_tenx_10xmousebrain | portable | 11,831 | 15,775 | 1,104 | standard | 8 | 55.4 s | 2.57 |
| soupx_tenx_10xmousebrain | production | 11,831 | 15,775 | 1,104 | standard | 8 | 56.3 s | 2.57 |
| soupx_emptydrops_4500871lane | portable | 26,598 | 35,464 | 709 | standard | 8 | 55.7 s | 2.45 |
| soupx_emptydrops_4500871lane | production | 26,598 | 35,464 | 709 | standard | 8 | 52.7 s | 2.45 |
| soupx_emptydrops_sscovelli | portable | 79,409 | 105,879 | 20,117 | light | 8 | 86.2 s | 4.48 |
| soupx_emptydrops_sscovelli | production | 79,409 | 105,879 | 20,117 | light | 8 | 79.4 s | 4.47 |
| cellbender_fromraw_ghirsutum2 | portable | 59,323 | 79,097 | 2,373 | light | 8 | 105.0 s | 4.09 |
| cellbender_fromraw_ghirsutum2 | production | 59,323 | 79,097 | 2,373 | light | 8 | 99.3 s | 4.09 |

## Benchmark

**What was run.**

- **Rust:** `scrust doubletfinder` without stage dumps, on each object's exported input: the counts as `matrix.mtx.gz` and the cluster labels, the same files R's direct run read. It ran at 1 and at 32 threads, in the production profile (`OPENBLAS_CORETYPE` unset, SapphireRapids detected) on Intel Xeon Platinum 8480CL nodes (`sapphire`, 32 CPUs allocated, shared nodes). The jobs were 60007089 (sscovelli) and 60007102 (ghirsutum2), [`validation/slurm/m6_benchmark.sbatch`](../../validation/slurm/m6_benchmark.sbatch), binary `scrust-m6-bc1b48e`.
- **R:** the direct production run of M2 on the same input and CPU type, with the same OpenBLAS kernel (one CPU, jobs 59639609 and 59639633).
- **Measures.** Wall time, CPU time and peak RSS come from `/usr/bin/time -v`. The phase times are R's `summary.json` marks and Rust's: *load* reads the input; *paramSweep* covers the six sweep iterations; *doubletFinder* is the final run (doublets, SCTransform, PCA, neighbours, pANN, calls).
- **Outputs.** Every Rust run's `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv` are byte-identical to R's, and its `result` section is identical.

| object | side | threads | wall | CPU time | peak RSS GiB | load | paramSweep | summary, pK | doubletFinder | vs R |
|---|---|---|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | R (direct run, M2) | 1 | 36:15.3 | 36:07.6 | 126.95 | 45.5 s | 10:41.0 | 0.4 s | 24:36.8 |  |
| soupx_emptydrops_sscovelli | Rust | 1 | 6:20.6 | 6:19.3 | 4.41 | 8.1 s | 2:21.7 | 1.0 s | 3:49.8 | identical: files 3/3, result 20/20 |
| soupx_emptydrops_sscovelli | Rust | 32 | 43.5 s | 7:12.2 | 4.51 | 8.1 s | 16.1 s | 0.1 s | 19.1 s | identical: files 3/3, result 20/20 |
| cellbender_fromraw_ghirsutum2 | R (direct run, M2) | 1 | 33:53.0 | 33:42.3 | 150.41 | 33.0 s | 12:12.7 | 0.3 s | 20:57.5 |  |
| cellbender_fromraw_ghirsutum2 | Rust | 1 | 8:58.4 | 8:55.4 | 4.02 | 9.1 s | 3:47.2 | 1.0 s | 5:01.0 | identical: files 3/3, result 20/20 |
| cellbender_fromraw_ghirsutum2 | Rust | 32 | 60.9 s | 10:57.1 | 4.08 | 9.0 s | 26.9 s | 0.1 s | 24.8 s | identical: files 3/3, result 20/20 |

Phase by phase in Rust: the time of each step of the final run, with the high-water mark of RSS after it (GiB).

| object | threads | after load | after the sweep | final SCTransform | final RunPCA | final kNN, pANN, calls | write |
|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | 1 | 1.84 | 1.84 | 2:18.5 (4.41 GiB) | 21.7 s (4.41 GiB) | 69.5 s (4.41 GiB) | 0.0 s |
| soupx_emptydrops_sscovelli | 32 | 1.84 | 1.84 | 13.6 s (4.51 GiB) | 3.1 s (4.51 GiB) | 2.4 s (4.51 GiB) | 0.1 s |
| cellbender_fromraw_ghirsutum2 | 1 | 2.04 | 2.04 | 4:05.6 (4.02 GiB) | 19.2 s (4.02 GiB) | 36.1 s (4.02 GiB) | 0.0 s |
| cellbender_fromraw_ghirsutum2 | 32 | 2.04 | 2.04 | 20.8 s (4.08 GiB) | 2.7 s (4.08 GiB) | 1.3 s (4.08 GiB) | 0.1 s |

**Time.**

| | sscovelli | ghirsutum2 |
|---|---|---|
| one thread | 6:21, 5.7 times faster than R | 8:58, 3.8 times faster |
| 32 threads | 43.5 s, 50 times faster | 60.9 s, 33 times faster |

On one thread, most of Rust's time goes where R's does, into the SCTransforms. The final run's takes 2:19 (sscovelli) and 4:06 (ghirsutum2): about 2.6 times faster than R's SCTransform on the same input (354 and 642 s, measured in M4). Also:

- PCA takes 19 to 22 s on one thread (M5).
- The neighbour search replaces R's distance matrices (89.7 GB for sscovelli, 50.1 GB for ghirsutum2): 70 and 36 s on one thread, 2.4 and 1.3 s on 32.
- Loading the counts takes 8 to 9 s, where R took 33 to 46 s.
- ghirsutum2 takes longer than sscovelli despite fewer cells: its 99,452 features leave 48,000 to 58,000 genes per SCTransform.

On 32 threads the final SCTransform runs 10 to 12 times faster than on one, PCA 7 times and the neighbour search 28 to 29 times.

**Memory.**

| | sscovelli | ghirsutum2 | R |
|---|---|---|---|
| peak RSS | 4.41 GiB (one thread), 4.51 GiB (32) | 4.02 and 4.08 GiB | 126.9 and 150.4 GiB |
| fraction of R's | 1/28 to 1/29 | 1/37 | |

sscovelli stays well within the plan's 10 GiB; 32 threads add 0.1 GiB of per-thread rows and buffers.

The high-water marks show the phases:

- **Loading:** 1.8 to 2.0 GiB. The decompressed `matrix.mtx` text and its triplets are read whole.
- **The sweep:** no higher. Each iteration's 10,000 real cells and their doublets are freed with the iteration.
- **The final run's SCTransform:** the peak. The sizes, measured on the final run's input:

  | | sscovelli | ghirsutum2 |
  |---|---|---|
  | input: real cells, then doublets | 106.1 million entries, 1.19 GiB | 118.1 million entries, 1.32 GiB |
  | the kept genes' rows | 22,247 genes, 1.19 GiB | 57,899 genes, 1.32 GiB |
  | `scale.data`, 3,000 x cells | 2.37 GiB | 1.77 GiB |
  | two bins' products | 0.79 GiB | 0.59 GiB |

  On sscovelli the peak comes at the residual stage: the kept genes' rows, the `scale.data` buffer and two products, 4.35 GiB plus a residual and a counts row per thread. On ghirsutum2 it comes at the setup: the input, its row-major copy and the kept genes' rows held at once, 3.96 GiB. The counts are freed before the final run's input is built, and the input once the setup has kept its rows.
- **RunPCA:** reads `scale.data` in place.
- **The neighbour search:** needs the embeddings and about 4 MB per thread.

**The change that lowered the peak.** The first full benchmark, with binary `scrust-m6-7a8d7c5` (jobs 60003079 and 60003081), peaked higher at the end of the residual stage:

| | first benchmark | now |
|---|---|---|
| sscovelli, 1 and 32 threads | 5.97 and 5.78 GiB | 4.41 and 4.51 GiB |
| ghirsutum2, 1 and 32 threads | 4.89 and 4.91 GiB | 4.02 and 4.08 GiB |

Times were within 4% of today's, and the outputs identical. At that peak, two copies of `scale.data` were held at once:

- the 3,000 kept rows, each a separate allocation;
- the matrix being assembled from them.

The rows were freed as they were copied, but glibc kept their pages. The kept genes' rows plus twice `scale.data` is 5.93 and 4.86 GiB, within 0.04 GiB of the peaks measured. The residual stage now keeps the rows as slots of one buffer that becomes `scale.data` itself (commit bc1b48e). `scrust sctransform` (M4) runs the same residual stage, so its assembly copy is gone too (not re-measured).

What remains could shrink further. Building the kept genes' rows straight from the input, without the full row-major copy, would take 1.2 to 1.3 GiB off the setup. Reading `matrix.mtx` in a stream would lower the loading. Neither is needed for the criterion.

## Findings

- **pN 0.1 to 0.3 share their starting state.** Each RunPCA ends with the generator at `set.seed(42)` plus `rnorm(nfeatures)`. So whenever two sweep contexts' PCAs have the same number of features, the next iterations' doublet draws start from the same state. On every P object all sweep PCAs have 3,000 features; R's own checkpoints show a single `before_doublets` state for pN 0.1 to 0.3 on each of them. The port reproduces this, as it must: the plan's B2 says the draws "continue from the previous pN's PCA", which is right, and in practice that state is nearly always the same one.
- **Three seeds, one stream.** The workflow's `set.seed(12345)` matters only for the sweep subsample, the first pN's doublets and the final run's doublets. Everything else is governed by SCTransform's `set.seed(1448145)` (three times per call) and RunPCA's `set.seed(42)`. Reproducing the stream needs every draw in between counted exactly: the sweep subsample (10,000 draws, partial Fisher-Yates), 2 x n_doublets per context, the 5,000-cell and 2,000-gene samples, two `future.seed` draws per fitting bin, and `rnorm(nfeatures)` per PCA.
- **The final `.Random.seed`.** R's reference runs record `digest::digest(paste(.Random.seed, collapse = ","), algo = "sha256", serialize = FALSE)` after `doubletFinder`. Rust's run reproduces the digest (all 47 comparisons), so the stream ends in R's state. Under teacher forcing (M3) this field was not applicable.
- **Freed memory is not always returned.** glibc's allocator raises its mmap threshold after a large block is freed, so later row-sized allocations come from the heap. A freed heap row stays resident unless it lies at the top of the heap. That made `scale.data`'s assembly from 3,000 separately allocated rows cost a second copy of the matrix ([Benchmark](#benchmark)). Keeping the rows in one buffer that becomes the matrix avoids it, whatever the allocator does.
- **M4's account of `scrust sctransform`'s peak** attributed it to the `.npy` write (a column-major copy of `scale.data`). The assembly from the kept rows doubled `scale.data` before that write too, and no longer does.

## Plan corrections

- **B2, "inside SCTransform ... then 1 draw per 500-gene bin"**: two draws per bin (`make_rng_seeds()` and `future_xapply()` each call `next_random_seed()`), as M4 found. The plan's table still says one. The end-to-end run confirms two in all 24 runs: every bin's checkpoint is identical.
- **B2, the consequence ("the doublets for pN 2-6 ... come from `set.seed(42)` plus whatever the previous pN's irlba consumed")** is right, and more uniform than it sounds. irlba consumes `rnorm(nfeatures)` and nothing more in every context of the reference data (M5), so pN 0.1 to 0.3 start from one and the same state whenever their predecessors' PCAs have the same number of features: on every production object.
- **B3, the inputs**: besides the counts and the labels, the order of the cluster factor's levels is an input (`seurat_clusters_levels.txt`, or `--cluster-levels`). `modelHomotypic` sums the squared cluster frequencies in level order, in long double. The command reads it as the M3 core did.
- **B3, `summary.json`**: "versions, seed, timings, peak RSS" sit in its `run` section, which records the run and cannot be compared. The comparable part is the `result` section, which also carries the digest of the final `.Random.seed` (SHA-256 of the state's decimal text joined by commas, as the reference scripts compute it).

## Deviations

- **One Rust run per profile.** The portable Rust run of each object is compared with the portable R runs on both CPU types (Sapphire Rapids and Cascade Lake). M2 showed that the kernel, not the CPU, sets R's bits, and Rust's own arithmetic does not depend on the CPU (plan A5b).
- **The large objects' production runs** have no instrumented R counterpart. Their Rust runs dump at the light level and are compared with the direct R run and the saved output. The portable runs are compared stage by stage, as in M3 to M5.
- **Validation and benchmark binaries.** The final validation (job 60007108) and benchmark (jobs 60007089, 60007102) ran binary `scrust-m6-bc1b48e`, built from commit bc1b48e. The commits after it change the documentation only. Two earlier validations and one earlier benchmark are kept as records (labels `m6-7a8d7c5`, `m6-f85b69a`, `m6-bench-7a8d7c5`).
- **The benchmark's R side is M2's direct run** (`runs/direct-production-sapphirerapids`, one CPU), on the same input files and CPU type, with the same OpenBLAS kernel. R was not rerun.
- **No new fixtures.** The tier-T end-to-end tests use the M2 fixtures already in the manifest (`tests/fixtures/doubletfinder/<object>/{input,reference,stages,saved}`). The Rust runs' outputs under `data/` are validation outputs, not references.

## M6b: the singlets' counts

M6b adds one output to `scrust doubletfinder`: `singlets/`, the counts without the doublets, as a new 10x directory (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) that the next tool can read. It is written in addition to `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json`, which are unchanged. This follows two rules of plan A3:
- every port keeps all of the original tool's outputs;
- every tool that processes counts also writes its result as new 10x count data, so the tools chain on plain files.

### Outcome

| Done when (plan B6 M6b, B5) | Result |
|---|---|
| `singlets/` identical to the saved objects' RNA counts, on T and P, in both profiles | **Yes, all 20 runs.** Identical: the barcodes in order, the feature lines, the header field, and every stored entry and value. That is 409,874,942 entries, explicit zeros included. |
| The other outputs unchanged | **Yes.** All 1,020 files of M6's validated runs (label `m6-bc1b48e`) are byte-identical in the new runs. The totals are M6's: 24 of 24 instrumented runs (13,396 checks, 2,646 RNG checkpoints, 1.58 billion values) and 23 of 23 direct runs. The saved-output comparisons, 20 of 20, now make 180 checks over 411,027,378 values, none differing. |
| Independent of the thread count | **Yes.** The benchmark's runs at 1 and 32 threads are byte-identical, `singlets/` included. `cargo test` writes it at 1, 2 and 7 threads. |
| The memory peak unchanged | **Yes, on the two largest objects:** 4.41 and 4.49 GiB on sscovelli, 4.02 and 4.08 GiB on ghirsutum2, as in M6. The counts are read again once the calls are made. |
| R reads it as the saved counts | **Yes, all 10 objects** (production runs, [`read10x_singlets.R`](../../validation/R/read10x_singlets.R), job 60139710). Both `Seurat::Read10X()` and the exact reader of `lib/io.R` give the saved object's RNA counts slot for slot, every stored entry and explicit zero included. |
| Tested and settled | `cargo test --workspace`: 160 tests (11 new), no warnings, in a fresh clone. The fresh clone's binary is byte-identical to the validated one, `scrust-m6b-13592a6` (sha256 `dfd3b9c697dcf263`). `tolerances.toml` gains `[doubletfinder.singlets]`. The manifest lists the 10 exports (60 new entries; no existing entry changed). |

Validation label `m6b-13592a6` (jobs 60114423 and 60117195, 8 threads), benchmark label `m6b-bench-13592a6` (jobs 60114424 and 60114425), binary `scrust-m6b-13592a6`, built from commit 13592a6.

### What the workflow's step leaves in the counts

`doubletfinder.R` keeps the singlets with `subset(seurat, cells = ...)`. It then runs SCTransform and the clustering on them, which add assays and leave the RNA counts as they are, and saves the object. So the saved object's RNA `counts` layer is the step's result at the level of the counts.

[`export_doubletfinder_singlets.R`](../../validation/R/export_doubletfinder_singlets.R) exports that layer into `data/doubletfinder/<object>/saved/singlets/`, next to M2's export. It first checks two things:
- the `.rds` is the one M2 exported (sha256);
- its cells are the kept cells, in order.

It reads the directory back `identical()`. It also compares the layer with the input export's counts restricted to the kept cells, slot for slot. They are identical for all 10 objects, explicit zeros included: R's `subset()` keeps the stored entries as they are. So "the input without the doublets, every entry as stored" is exactly R's result.

| object | singlets | features | entries | explicit zeros | equals the input restricted to the kept cells |
|---|---|---|---|---|---|
| soupx_tenx_test | 276 | 2,000 | 59,103 | 2,453 | yes |
| soupx_emptydrops_test | 237 | 2,000 | 54,893 | 684 | yes |
| decontx_raw_tenx_test | 273 | 2,000 | 52,332 | 0 | yes |
| decontx_raw_emptydrops_test | 237 | 2,000 | 50,896 | 0 | yes |
| cellbender_fromraw_test | 247 | 2,000 | 54,936 | 0 | yes |
| soupx_tenx_SRR18256555 | 2,688 | 32,285 | 2,461,479 | 64,274 | yes |
| soupx_tenx_10xmousebrain | 10,185 | 32,285 | 36,923,827 | 386,679 | yes |
| soupx_emptydrops_4500871lane | 23,239 | 32,285 | 35,907,669 | 1,223,791 | yes |
| soupx_emptydrops_sscovelli | 69,500 | 28,754 | 62,191,630 | 3,443,385 | yes |
| cellbender_fromraw_ghirsutum2 | 52,554 | 99,452 | 67,180,706 | 0 | yes |

The exports (jobs 60108830, 60108897, 60108900, 60108904, 60108908 and 60108910) took 1 s to 5.6 minutes each, with peaks of at most 5.1 GiB.

### How `singlets/` is written

- **The counts, read again.** The run frees the counts when it builds the final run's input, which keeps the peak down (M6). Once the calls are made, [`doubletfinder::singlets`](../../crates/doubletfinder/src/singlets.rs) reads them again from `--counts`. It keeps the singlets' columns in input order, each exactly as stored, and the input's feature lines as they are.
- **The directory.** [`scio::tenx::write_columns`](../../crates/scio/src/tenx.rs) writes a Cell Ranger v3 directory, which Seurat's `Read10X`, `DropletUtils` and scanpy read. Values are written as `%.17g` writes them, integers as integers, so a correctly rounded reader recovers every double. The header says `integer` or `real` according to the values, and a `%metadata_json` line names the writer, as in Cell Ranger's files. The matrix text is cut into pieces of 131,072 entries, formatted and compressed on the run's threads. The bytes do not depend on the number of threads.
- **gzip.** [`scio::deflate`](../../crates/scio/src/deflate.rs) is a gzip writer in the standard library only, the counterpart of the reader `scio::gz`:
  - LZ77 with hash chains, then dynamic Huffman codes per block, limited to 15 bits by package-merge;
  - each piece coded on its own and ended by an empty stored block, as zlib's `Z_SYNC_FLUSH` does, with the pieces' CRC-32s combined.

  The files are 0.04% to 1% larger than R's zlib output for the same matrices: 0.04% to 0.23% for the four files above 100 MB, up to 1% for the test objects' 150 KB. Chains of 16 candidates were chosen on sscovelli's 66 million entries:

  | candidates per position | one thread | size against `gzip -6` |
  |---|---|---|
  | 32 | 39.2 s | +0.05% |
  | **16** | **27.8 s** | **+0.24%** |
  | 8 | 17.0 s | +6.3% |
  | 4 | 13.2 s | +8.3% |

### Cost

Production profile, Intel Xeon Platinum 8480CL, the same inputs as M6's benchmark. *Until the calls* is the run up to the calls, as in M6 (whose benchmark is the comparison). *Writing `singlets/`* re-reads the counts and writes the directory.

| object | threads | M6 total | until the calls | writing `singlets/` | total | peak RSS GiB | vs R's direct run |
|---|---|---|---|---|---|---|---|
| soupx_emptydrops_sscovelli | 1 | 6:20.6 | 6:37.9 | 34.2 s | 7:12.1 | 4.41 | 5.0 times faster, 1/29 of the memory |
| soupx_emptydrops_sscovelli | 32 | 43.5 s | 44.1 s | 9.9 s | 54.0 s | 4.49 | 40 times faster |
| cellbender_fromraw_ghirsutum2 | 1 | 8:58.4 | 8:55.5 | 37.9 s | 9:33.4 | 4.02 | 3.5 times faster, 1/37 of the memory |
| cellbender_fromraw_ghirsutum2 | 32 | 60.9 s | 58.0 s | 11.0 s | 69.0 s | 4.08 | 29 times faster |

- **The new phase.** It reads the counts again, 8.4 to 9.4 s as at the start, and writes the directory, about 26 to 29 s on one thread and 1 to 2 s on 32.
- **The run before it** is within 5% of M6's on other nodes.
- **The memory peak** stays where M6 put it, in the final run's SCTransform. The re-read needs what loading needed, 1.8 to 2.0 GiB.
- **Small objects.** In the validation runs of soupx_tenx_SRR18256555 (standard dumps), the re-read is the high-water mark in one of the two runs: 0.45 GiB, against 0.42 GiB in M6.
- **The comparison is conservative.** R's direct run writes no count matrix.
- **Keeping the counts in memory through the final run** instead would save the re-read at the cost of 0.74 GiB (sscovelli) to 0.82 GiB (ghirsutum2) on the peak.

## Reproducing

From the repository's top directory, after `source env.sh`:

```bash
cargo build --release -p scrust
c=$(git rev-parse --short HEAD)
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m6-$c
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6-$c sbatch --export=ALL validation/slurm/m6_doubletfinder.sbatch m6-$c soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6-$c sbatch --export=ALL validation/slurm/m6_benchmark.sbatch m6-bench-$c soupx_emptydrops_sscovelli 1 32
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6-$c sbatch --export=ALL validation/slurm/m6_benchmark.sbatch m6-bench-$c cellbender_fromraw_ghirsutum2 1 32
```

The tables come from [`validation/R/report_doubletfinder_m6.R`](../../validation/R/report_doubletfinder_m6.R), which reads plain files only. Run it with the production R env and the hermetic flags:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" validation/R/report_doubletfinder_m6.R m6-<commit> m6-bench-<commit>
```

SLURM logs of the M6 jobs are in `data/doubletfinder/logs/`.

M6b uses the same wrappers. First comes a one-time export of the saved outputs' RNA counts, which reads `.rds` files and needs no particular CPU type. Then come a validation and a benchmark with the new binary. Last, R's own readers read `singlets/` back. The logs are in `data/doubletfinder/logs/m6b/`.

```bash
O="soupx_tenx_test soupx_emptydrops_test decontx_raw_tenx_test decontx_raw_emptydrops_test cellbender_fromraw_test soupx_tenx_SRR18256555 soupx_tenx_10xmousebrain soupx_emptydrops_4500871lane soupx_emptydrops_sscovelli cellbender_fromraw_ghirsutum2"
sbatch -p shared,sapphire validation/slurm/export_doubletfinder_singlets.sbatch $O
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m6b-$c
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6b-$c sbatch --export=ALL validation/slurm/m6_doubletfinder.sbatch m6b-$c $O
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6b-$c sbatch --export=ALL --cpus-per-task=32 --mem=32G validation/slurm/m6_benchmark.sbatch m6b-bench-$c soupx_emptydrops_sscovelli 1 32
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m6b-$c sbatch --export=ALL --cpus-per-task=32 --mem=32G validation/slurm/m6_benchmark.sbatch m6b-bench-$c cellbender_fromraw_ghirsutum2 1 32
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" validation/R/read10x_singlets.R m6b-$c production $O   # in a batch job: 5 minutes, 5.2 GiB
```

One object by hand, on a tier-T fixture (seconds, any x86-64 machine with AVX2, under env.sh's portable profile):

```bash
f=tests/fixtures/doubletfinder/soupx_tenx_test
target/release/scrust doubletfinder --counts $f/input --clusters $f/input/seurat_clusters.tsv --out /path/to/out --dump-stages light
target/release/scrust validate doubletfinder --reference $f/reference --reference-stages $f/stages --saved $f/saved --candidate /path/to/out
```
