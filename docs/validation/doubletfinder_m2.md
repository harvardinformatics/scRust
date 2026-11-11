# DoubletFinder reference data (milestone M2)

M2 builds the R reference data that the Rust port of DoubletFinder (M3–M7) is checked against (docs/PLAN.md B6). It ports nothing to Rust. This report gives the checks behind each "done when" criterion. Recipes: [`validation/R/README.md`](../../validation/R/README.md) and [`validation/slurm/README.md`](../../validation/slurm/README.md); every file is listed in [`validation/manifest.toml`](../../validation/manifest.toml).

## Outcome

All of plan B6's M2 criteria are met.

| Criterion | Result |
|---|---|
| Each reference run reproduces its saved workflow output | **Yes, all 47 runs of all 10 objects**, in both profiles: identical doublet sets, singlet pANN identical bit for bit, identical column names (so pK and nExp). |
| Portable fixtures have identical checksums on at least 2 CPU types | **Yes.** Portable runs on Sapphire Rapids (Xeon Platinum 8480CL) and Cascade Lake (Xeon Platinum 8268), tier T and soupx_tenx_SRR18256555, direct and instrumented: every deterministic file has the same sha256 (`identical_on` in the manifest), every stage array and RNG checkpoint is bitwise identical. |
| Production core name recorded | **`SapphireRapids`** (OpenBLAS 0.3.34 detection, `OPENBLAS_CORETYPE` unset), in every production run's `summary.json` and in the manifest. |
| R-vs-R floor measured | **Yes**, per stage (below; settled as measurements in `validation/tolerances.toml`). Same kernel on another CPU type: no difference at all. Haswell against SapphireRapids kernels: only S2d and S3 differ, in the last bits; S4–S7 and every final output are identical. |
| DoubletFinder's own memory and time measured | **Yes**, `/usr/bin/time -v` and sacct for every run (below and in the manifest's `[[run]]` tables). |

The instrumented runs change nothing: their final outputs equal the direct runs' byte for byte in all 21 pairs (tier T in three configurations, soupx_tenx_SRR18256555 in three, soupx_tenx_10xmousebrain in two, soupx_emptydrops_4500871lane in one). `cargo test` passes, including a new check of the tier-T fixtures.

## What was run

**Reference.** DoubletFinder 2.0.6 (commit `1b244d8`) as the workflow's `workflow/scripts/doubletfinder.R` calls it, seed 12345, in the production R env (R 4.5.3, Seurat 5.5.1, sctransform 0.4.3, glmGamPoi 1.22.0, irlba 2.4.1; OpenBLAS 0.3.34, one thread). The call sequence lives in one place, [`validation/R/lib/doubletfinder.R`](../../validation/R/lib/doubletfinder.R), shared by the direct and instrumented runs: the workflow's library calls and option, `set.seed(12345)`, `paramSweep(PCs = 1:10, sct = TRUE)`, `summarizeSweep(GT = FALSE)`, `find.pK` (its plot goes to `find_pK.pdf`), the pK through `as.numeric(as.character(...))`, the homotypic nExp, the workflow's `rm()`/`gc()`, `set.seed(12345)`, `doubletFinder(pN = 0.25, ...)`. What the workflow does afterwards is not run.

**Inputs.** The reference runs read the exported plain files, not the `.rds`: the counts are rebuilt as the identical dgCMatrix and wrapped with `CreateSeuratObject()`. With `sct = TRUE`, DoubletFinder 2.0.6 reads only the RNA assay's `counts` layer (`LayerData(seu, assay = "RNA", layer = "counts")`, since the objects are SeuratObject 5.4.0), the meta.data row names, and `seurat_clusters` for the homotypic adjustment. The export verified, for every object, that re-reading the Matrix Market files gives back the object's own matrix bit for bit (`identical(..., num.eq = FALSE)`: slots `i`, `p`, `x` and the dimnames). Passing the cross-check from these files also shows the exports are complete.

**Objects.**

| object | tier | features | cells | entries | explicit zeros | clusters | saved |
|---|---|---|---|---|---|---|---|
| soupx_tenx_test | T | 2,000 | 300 | 64,934 | 2,686 | 3 | `pANN_0.25_0.29_24` (0.29, 24) |
| soupx_emptydrops_test | T | 2,000 | 267 | 62,981 | 809 | 4 | `pANN_0.25_0.28_30` (0.28, 30) |
| decontx_raw_tenx_test | T | 2,000 | 300 | 57,708 | 0 | 3 | `pANN_0.25_0.09_27` (0.09, 27) |
| decontx_raw_emptydrops_test | T | 2,000 | 267 | 58,001 | 0 | 4 | `pANN_0.25_0.22_30` (0.22, 30) |
| cellbender_fromraw_test | T | 2,000 | 278 | 63,782 | 0 | 4 | `pANN_0.25_0.23_31` (0.23, 31) |
| soupx_tenx_SRR18256555 | P | 32,285 | 3,108 | 2,817,955 | 74,901 | 15 | `pANN_0.25_0.005_420` (0.005, 420) |
| soupx_tenx_10xmousebrain | P | 32,285 | 11,831 | 42,440,187 | 445,249 | 21 | `pANN_0.25_0.07_1646` (0.07, 1646) |
| soupx_emptydrops_4500871lane | P | 32,285 | 26,598 | 43,674,465 | 1,467,924 | 31 | `pANN_0.25_0.02_3359` (0.02, 3359) |
| soupx_emptydrops_sscovelli | P | 28,754 | 79,409 | 66,036,718 | 3,602,702 | 36 | `pANN_0.25_0.19_9909` (0.19, 9909) |
| cellbender_fromraw_ghirsutum2 | P | 99,452 | 59,323 | 73,602,177 | 0 | 22 | `pANN_0.25_0.03_6769` (0.03, 6769) |

**Profiles (A5b).** *Portable*: `OPENBLAS_CORETYPE=Haswell`. *Production*: `OPENBLAS_CORETYPE` unset, so OpenBLAS detects the CPU, on `sapphire` nodes; every production run reported core `SapphireRapids` (OpenBLAS 0.3.34 DYNAMIC_ARCH, one thread, on Intel Xeon Platinum 8480CL). Each job checked the core with `scrust fingerprint` before running and recorded the fingerprint in `summary.json`.

Setting the variable is not the same as detection: on these Sapphire Rapids CPUs, `OPENBLAS_CORETYPE=SapphireRapids` and a literal `OPENBLAS_CORETYPE=auto` ("Core not found: auto") both select `Cooperlake`. The production profile must leave the variable unset, as `env.sh` does for `auto`.

**Runs.** Run directories are `data/doubletfinder/<object>/runs/<mode>-<profile>-<cpu arch>/`.

- Tier T, all five objects: direct and instrumented runs in both profiles on Sapphire Rapids, and in the portable profile on Cascade Lake (`shared`, Xeon Platinum 8268). Runs on AMD Genoa (`serial_requeue --constraint=genoa`) were submitted but never scheduled and were cancelled.
- soupx_tenx_SRR18256555 (small): as tier T, including Cascade Lake.
- soupx_tenx_10xmousebrain (medium): direct and instrumented, both profiles.
- soupx_emptydrops_4500871lane (medium 2): direct production; instrumented in both profiles.
- soupx_emptydrops_sscovelli and cellbender_fromraw_ghirsutum2 (large): direct production (the cross-check and a clean memory and time baseline) and instrumented portable at the light dump level, which stands in for the portable direct run once instrumented runs are proven identical to direct ones.

## Where the test references were made

The tier-T references (`tests/reference_outputs/testdata/results/`, files dated 2026-11-06 14:40–14:48) come from the test-workflow run started at 14:36 (Snakemake log `.snakemake/log/2026-11-06T143629.726589.snakemake.log`, runner on holy8a24106); `tests/update_reference_outputs.py` copies them with `shutil.copy2`, which keeps those times. A second run at 15:01 rewrote `testdata/results/`, not the references.

| Object | Job | Partition | Node | Elapsed | sacct MaxRSS |
|---|---|---|---|---|---|
| cellbender_fromraw_test | 58374065 | sapphire | holy8a32608 | 1:20 | 0.63 GiB |
| soupx_tenx_test | 58374068 | sapphire | holy8a32608 | 1:20 | 0.62 GiB |
| soupx_emptydrops_test | 58374079 | sapphire | holy8a32608 | 1:18 | 0.56 GiB |
| decontx_raw_emptydrops_test | 58374282 | sapphire | holy8a32606 | 1:21 | 0.80 GiB |
| decontx_raw_tenx_test | 58374301 | sapphire | holy8a32606 | 1:19 | 0.82 GiB |

Both nodes are Sapphire Rapids (`sapphirerapids` feature; holy8a32608 is one of the two nodes now excluded for the lost OST). The run used the conda env `.snakemake/conda/b11ef12478f5805db716efa9cb089f6b_`, which differs from the production env only in xgboost (3.4.2 against 3.2.0; DoubletFinder does not use it): every other package build is the same, the lazy-load databases of Seurat, SeuratObject, sctransform, glmGamPoi, irlba, fields, KernSmooth, Matrix and matrixStats are byte-identical, and so is `libopenblasp-r0.3.34.so` (md5 `d77a762c…`, also in `.envs/rust-dev`). Its DoubletFinder was installed from the same commit `1b244d8` (its lazy-load database differs, as locally built packages do). Nothing set `OPENBLAS_CORETYPE`, so OpenBLAS detected the CPU: **the test references were made in the production profile**. Both profiles reproduce them exactly (below).

## Production jobs (tier P)

The saved outputs come from the production run of 2026-11-08/09 (runner 59050095, relaunch 59391964). The input `.rds` of sscovelli and ghirsutum2 were rewritten by the relaunch seconds before their DoubletFinder jobs started; the earlier jobs (58831227, 58828690, from the 2026-11-07 run) are superseded.

| Object | Job | Node | CPU | Elapsed | sacct MaxRSS |
|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | 59054994 | holy8a30604 | Sapphire Rapids | 4:36 | 5.2 GiB |
| soupx_tenx_10xmousebrain | 59123277 | holy8a32604 | Sapphire Rapids | 14:18 | 19.7 GiB |
| soupx_emptydrops_4500871lane | 59160818 | holy8a32606 | Sapphire Rapids | 15:31 | 38.6 GiB |
| soupx_emptydrops_sscovelli | 59409415 | holy8a24105 | Sapphire Rapids | 41:19 | 125.4 GiB |
| cellbender_fromraw_ghirsutum2 | 59406621 | holy8a32603 | Sapphire Rapids | 44:08 | 137.1 GiB |

All on the `sapphire` partition (features `intel,holyndr,sapphirerapids,avx,avx2,avx512`), one CPU each. These peaks include the wrapper's Seurat work after DoubletFinder, and sacct samples every 30 s.

## Cross-check against the saved workflow outputs

Each run against the saved output's fields (`data/doubletfinder/<object>/saved/`, exported by `export_doubletfinder.R`): the doublet set (input barcodes minus output barcodes), the pANN of every kept cell (compared as 17-digit text, i.e. bitwise), and the column names. Per-run JSON: `data/doubletfinder/compare/<object>/crosscheck_<run>.json`.

| object | run | pK | nExp | columns | doublet set | singlet pANN identical | result |
|---|---|---|---|---|---|---|---|
| soupx_tenx_test | direct-portable-cascadelake | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_tenx_test | direct-portable-sapphirerapids | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_tenx_test | direct-production-sapphirerapids | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_tenx_test | instrumented-portable-cascadelake | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_tenx_test | instrumented-portable-sapphirerapids | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_tenx_test | instrumented-production-sapphirerapids | 0.29 | 24 | identical | identical (24) | 276 / 276 | pass |
| soupx_emptydrops_test | direct-portable-cascadelake | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| soupx_emptydrops_test | direct-portable-sapphirerapids | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| soupx_emptydrops_test | direct-production-sapphirerapids | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| soupx_emptydrops_test | instrumented-portable-cascadelake | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| soupx_emptydrops_test | instrumented-portable-sapphirerapids | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| soupx_emptydrops_test | instrumented-production-sapphirerapids | 0.28 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_tenx_test | direct-portable-cascadelake | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_tenx_test | direct-portable-sapphirerapids | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_tenx_test | direct-production-sapphirerapids | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_tenx_test | instrumented-portable-cascadelake | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_tenx_test | instrumented-portable-sapphirerapids | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_tenx_test | instrumented-production-sapphirerapids | 0.09 | 27 | identical | identical (27) | 273 / 273 | pass |
| decontx_raw_emptydrops_test | direct-portable-cascadelake | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_emptydrops_test | direct-portable-sapphirerapids | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_emptydrops_test | direct-production-sapphirerapids | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_emptydrops_test | instrumented-portable-cascadelake | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_emptydrops_test | instrumented-portable-sapphirerapids | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| decontx_raw_emptydrops_test | instrumented-production-sapphirerapids | 0.22 | 30 | identical | identical (30) | 237 / 237 | pass |
| cellbender_fromraw_test | direct-portable-cascadelake | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| cellbender_fromraw_test | direct-portable-sapphirerapids | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| cellbender_fromraw_test | direct-production-sapphirerapids | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| cellbender_fromraw_test | instrumented-portable-cascadelake | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| cellbender_fromraw_test | instrumented-portable-sapphirerapids | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| cellbender_fromraw_test | instrumented-production-sapphirerapids | 0.23 | 31 | identical | identical (31) | 247 / 247 | pass |
| soupx_tenx_SRR18256555 | direct-portable-cascadelake | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_SRR18256555 | direct-portable-sapphirerapids | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_SRR18256555 | direct-production-sapphirerapids | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_SRR18256555 | instrumented-portable-cascadelake | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_SRR18256555 | instrumented-portable-sapphirerapids | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_SRR18256555 | instrumented-production-sapphirerapids | 0.005 | 420 | identical | identical (420) | 2688 / 2688 | pass |
| soupx_tenx_10xmousebrain | direct-portable-sapphirerapids | 0.07 | 1646 | identical | identical (1646) | 10185 / 10185 | pass |
| soupx_tenx_10xmousebrain | direct-production-sapphirerapids | 0.07 | 1646 | identical | identical (1646) | 10185 / 10185 | pass |
| soupx_tenx_10xmousebrain | instrumented-portable-sapphirerapids | 0.07 | 1646 | identical | identical (1646) | 10185 / 10185 | pass |
| soupx_tenx_10xmousebrain | instrumented-production-sapphirerapids | 0.07 | 1646 | identical | identical (1646) | 10185 / 10185 | pass |
| soupx_emptydrops_4500871lane | direct-production-sapphirerapids | 0.02 | 3359 | identical | identical (3359) | 23239 / 23239 | pass |
| soupx_emptydrops_4500871lane | instrumented-portable-sapphirerapids | 0.02 | 3359 | identical | identical (3359) | 23239 / 23239 | pass |
| soupx_emptydrops_4500871lane | instrumented-production-sapphirerapids | 0.02 | 3359 | identical | identical (3359) | 23239 / 23239 | pass |
| soupx_emptydrops_sscovelli | direct-production-sapphirerapids | 0.19 | 9909 | identical | identical (9909) | 69500 / 69500 | pass |
| soupx_emptydrops_sscovelli | instrumented-portable-sapphirerapids | 0.19 | 9909 | identical | identical (9909) | 69500 / 69500 | pass |
| cellbender_fromraw_ghirsutum2 | direct-production-sapphirerapids | 0.03 | 6769 | identical | identical (6769) | 52554 / 52554 | pass |
| cellbender_fromraw_ghirsutum2 | instrumented-portable-sapphirerapids | 0.03 | 6769 | identical | identical (6769) | 52554 / 52554 | pass |

The test references were made in the production profile, but the portable profile reproduces them too: on every object, the two profiles give identical final outputs.

## Run against run

Final files: `cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv` (byte comparison). Result: the `result` section of `summary.json` (pK, nExp, homotypic proportion, column names, digest of the final `.Random.seed`, ...). Stage arrays and RNG checkpoints: every block of every array, bitwise (`data/doubletfinder/compare/<object>/*_stages.tsv`). `instr_vs_direct_*`: the instrumented run against the direct run. `*_cascadelake_vs_sapphirerapids`: the portable profile on two CPU types. `portable_vs_production`: Haswell against SapphireRapids kernels on Sapphire Rapids (for the large objects, the portable instrumented run against the production direct run, final outputs only).

| object | comparison | final files identical | result identical | stage arrays identical | RNG checkpoints identical | calls differ | pANN differ |
|---|---|---|---|---|---|---|---|
| soupx_tenx_test | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_test | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 521 / 521 | TRUE | 0 | 0 |
| soupx_tenx_test | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_test | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_test | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_test | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 493 / 521 | TRUE | 0 | 0 |
| soupx_emptydrops_test | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_emptydrops_test | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 521 / 521 | TRUE | 0 | 0 |
| soupx_emptydrops_test | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_emptydrops_test | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_emptydrops_test | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_emptydrops_test | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 493 / 521 | TRUE | 0 | 0 |
| decontx_raw_tenx_test | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_tenx_test | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 521 / 521 | TRUE | 0 | 0 |
| decontx_raw_tenx_test | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_tenx_test | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_tenx_test | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_tenx_test | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 493 / 521 | TRUE | 0 | 0 |
| decontx_raw_emptydrops_test | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_emptydrops_test | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 521 / 521 | TRUE | 0 | 0 |
| decontx_raw_emptydrops_test | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_emptydrops_test | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_emptydrops_test | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| decontx_raw_emptydrops_test | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 493 / 521 | TRUE | 0 | 0 |
| cellbender_fromraw_test | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| cellbender_fromraw_test | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 521 / 521 | TRUE | 0 | 0 |
| cellbender_fromraw_test | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| cellbender_fromraw_test | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| cellbender_fromraw_test | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| cellbender_fromraw_test | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 493 / 521 | TRUE | 0 | 0 |
| soupx_tenx_SRR18256555 | direct_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_SRR18256555 | instrumented_cascadelake_vs_sapphirerapids | TRUE/TRUE/TRUE | TRUE | 569 / 569 | TRUE | 0 | 0 |
| soupx_tenx_SRR18256555 | instr_vs_direct_portable-cascadelake | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_SRR18256555 | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_SRR18256555 | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_SRR18256555 | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 527 / 569 | TRUE | 0 | 0 |
| soupx_tenx_10xmousebrain | direct_portable_vs_production | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_10xmousebrain | instr_vs_direct_portable-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_10xmousebrain | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_tenx_10xmousebrain | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 540 / 582 | TRUE | 0 | 0 |
| soupx_emptydrops_4500871lane | instr_vs_direct_production-sapphirerapids | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| soupx_emptydrops_4500871lane | portable_vs_production | TRUE/TRUE/TRUE | TRUE | 552 / 594 | TRUE | 0 | 0 |
| soupx_emptydrops_sscovelli | portable_vs_production | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |
| cellbender_fromraw_ghirsutum2 | portable_vs_production | TRUE/TRUE/TRUE | TRUE |  |  | 0 | 0 |

## CPU types (portable profile)

Tier T and soupx_tenx_SRR18256555 ran in the portable profile on Sapphire Rapids (`sapphire`) and Cascade Lake (`shared`), direct and instrumented. Every deterministic file has the same sha256 on both (all files except `summary.json`, `find_pK.pdf`, `time_v.txt` and `fingerprint_wrapper.json`, which record the run; the `result` sections of the `summary.json` files are identical), and every stage array and RNG checkpoint is bitwise identical. The manifest lists the Cascade Lake runs as `[[run]]` tables and under `identical_on` in the Sapphire Rapids files' entries. So the portable references can be regenerated on another CPU type and checked by checksum, as plan A5b intends. An AMD CPU (Genoa) would have been a stronger test of that claim, and of GitHub's runners, but the two jobs never got a node and were cancelled; that check remains open.

## R-vs-R floor

Same CPU type (Sapphire Rapids), OpenBLAS kernel Haswell (portable) against SapphireRapids (production), on the eight objects with instrumented runs in both profiles. Per array: how many objects differ at all, how many elements differ, the largest absolute difference, and the largest difference relative to the block's largest absolute value. S3 embeddings and loadings are compared after aligning each component's sign. The same kernel on another CPU type gives no difference at all (previous section).

| stage | quantity | objects | objects with differences | elements differing | max abs diff | max diff / max abs value | note |
|---|---|---|---|---|---|---|---|
| S1 | S1_doublets | 8 | 0 | 0 / 94944 | 0 | 0 |  |
| S1 | S1_subsample | 2 | 0 | 0 / 20000 | 0 | 0 |  |
| S2a | S2a_cells | 8 | 0 | 0 / 475082 | 0 | 0 |  |
| S2a | S2a_cells_step1 | 8 | 0 | 0 / 109260 | 0 | 0 |  |
| S2a | S2a_density | 3 | 0 | 0 / 21504 | 0 | 0 |  |
| S2a | S2a_gene_stats | 8 | 0 | 0 / 2219605 | 0 | 0 |  |
| S2a | S2a_genes | 8 | 0 | 0 / 443921 | 0 | 0 |  |
| S2a | S2a_genes_step1 | 8 | 0 | 0 / 62072 | 0 | 0 |  |
| S2a | S2a_genes_step1_pool | 8 | 0 | 0 / 290280 | 0 | 0 |  |
| S2a | S2a_log_gmean_step1 | 8 | 0 | 0 / 62072 | 0 | 0 |  |
| S2a | S2a_sampling_prob | 3 | 0 | 0 / 270208 | 0 | 0 |  |
| S2b | S2b_glm | 8 | 0 | 0 / 558648 | 0 | 0 |  |
| S2b | S2b_model_pars | 8 | 0 | 0 / 186216 | 0 | 0 |  |
| S2b | S2b_model_pars_raw | 8 | 0 | 0 / 186216 | 0 | 0 |  |
| S2c | S2c_model_pars_fit | 8 | 0 | 0 / 1331763 | 0 | 0 |  |
| S2c | S2c_outliers | 8 | 0 | 0 / 62072 | 0 | 0 |  |
| S2c | S2c_poisson_genes | 8 | 0 | 0 / 153160 | 0 | 0 |  |
| S2c | S2c_smoothing_genes | 8 | 0 | 0 / 60278 | 0 | 0 |  |
| S2d | S2d_pca_features | 8 | 0 | 0 / 92682 | 0 | 0 |  |
| S2d | S2d_residual_stats | 8 | 3 | 4065 / 887842 | 5.68e-14 | 1.36e-16 |  |
| S2d | S2d_scale_data | 8 | 3 | 3538894 / 686252033 | 2.13e-14 | 1.18e-15 |  |
| S2d | S2d_top_features | 8 | 0 | 0 / 92682 | 0 | 0 |  |
| S2d | S2d_variable_features | 8 | 0 | 0 / 92682 | 0 | 0 |  |
| S2d | S2d_y | 5 | 0 | 0 / 10427033 | 0 | 0 |  |
| S3 | S3_embeddings | 8 | 8 | 2344210 / 2375410 | 3.55e-11 | 2.41e-13 | after aligning signs; 53 of 56 blocks had 1 to 8 components flipped |
| S3 | S3_irlba_d | 8 | 8 | 534 / 560 | 7.5e-12 | 3.8e-15 |  |
| S3 | S3_loadings | 8 | 8 | 915739 / 926820 | 8.03e-14 | 2.97e-13 | after aligning signs; 53 of 56 blocks had 1 to 8 components flipped |
| S3 | S3_sdev | 8 | 8 | 529 / 560 | 6.75e-14 | 3.84e-15 |  |
| S4 | S4_digest | 8 | 0 | 0 / 380138 | 0 | 0 |  |
| S4 | S4_positions | 5 | 0 | 0 / 1042134 | 0 | 0 |  |
| S5 | S5_pann_final | 8 | 0 | 0 / 42949 | 0 | 0 |  |
| S5 | S5_pann_sweep | 8 | 0 | 0 / 4706832 | 0 | 0 |  |
| S6 | S6_bkde | 8 | 0 | 0 / 1111572 | 0 | 0 |  |
| S6 | S6_gkde | 8 | 0 | 0 / 4706832 | 0 | 0 |  |
| S7 | S7_doublet | 8 | 0 | 0 / 42949 | 0 | 0 |  |

- **S1, S2a, S2b, S2c: identical.** No BLAS in the sampling, filters, gene statistics, glmGamPoi fits (including the Cox–Reid term, a scalar for the intercept-only design) or regularisation.
- **S2d: last bits, larger objects only.** The residuals (`exp(tcrossprod(model_pars, regressors))`, OpenBLAS `dgemm` with an inner dimension of 2) differ in at most 0.8% of a block by at most 2.1e-14 on soupx_tenx_SRR18256555 and both medium objects. On tier T, `y`, the residual statistics and `scale.data` are bit-identical across kernels; OpenBLAS's kernels plausibly share code for such small matrices. The variable features and their order are identical everywhere.
- **S3: last bits, plus signs.** sdev differ by at most 7.3e-15 relative; embeddings by at most 3.6e-11 absolute, 2.4e-13 of the largest value; most blocks had some components with flipped sign, which leaves distances unchanged.
- **S4–S7: identical** on every object, including the large ones end to end, and so are all RNG checkpoints: irlba consumed the same number of random numbers under both kernels.

## Memory and time

DoubletFinder's own cost: the direct run in the production profile, with the input loaded from the exported files (`load min`), against the workflow job, which also ran Seurat after DoubletFinder and whose sacct MaxRSS is sampled every 30 s. Peak RSS is `/usr/bin/time -v`'s maximum resident set size (the same as the R process's VmHWM); `peak after sweep` is VmHWM after `paramSweep`, `summarizeSweep` and the workflow's `gc()`; `doubletFinder min` is the final run alone. Tier-T sacct values are per job (five objects per job). The last column is the instrumented run, for the overhead of the dumps.

| object | job | peak RSS GiB | sacct MaxRSS GiB | wall min | load min | paramSweep min | doubletFinder min | peak after sweep GiB | workflow job | workflow MaxRSS GiB | workflow min | instrumented peak GiB / min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_test | 59637971 | 1.0 | 1.0 | 0.5 | 0.0 | 0.4 | 0.0 | 1.0 | 58374068 | 0.6 | 1.3 | 1.1 / 0.5 (production) |
| soupx_emptydrops_test | 59637971 | 1.0 | 1.0 | 0.5 | 0.0 | 0.4 | 0.0 | 1.0 | 58374079 | 0.6 | 1.3 | 1.1 / 0.5 (production) |
| decontx_raw_tenx_test | 59637971 | 0.9 | 1.0 | 0.5 | 0.0 | 0.3 | 0.0 | 0.9 | 58374301 | 0.8 | 1.3 | 1.1 / 0.5 (production) |
| decontx_raw_emptydrops_test | 59637971 | 1.0 | 1.0 | 0.5 | 0.0 | 0.3 | 0.0 | 1.0 | 58374282 | 0.8 | 1.4 | 1.1 / 0.5 (production) |
| cellbender_fromraw_test | 59637971 | 1.0 | 1.0 | 0.5 | 0.0 | 0.4 | 0.0 | 1.0 | 58374065 | 0.6 | 1.3 | 1.1 / 0.5 (production) |
| soupx_tenx_SRR18256555 | 59639111 | 5.1 | 3.9 | 3.3 | 0.0 | 2.7 | 0.5 | 5.1 | 59054994 | 5.2 | 4.6 | 5.4 / 3.5 (production) |
| soupx_tenx_10xmousebrain | 59639138 | 22.8 | 19.7 | 12.4 | 0.4 | 9.9 | 2.0 | 22.8 | 59123277 | 19.7 | 14.3 | 22.5 / 17.9 (production) |
| soupx_emptydrops_4500871lane | 59639148 | 38.1 | 34.9 | 18.0 | 0.5 | 11.8 | 5.5 | 19.0 | 59160818 | 38.6 | 15.5 | 39.0 / 20.1 (production) |
| soupx_emptydrops_sscovelli | 59639609 | 126.9 | 127.2 | 36.3 | 0.8 | 10.7 | 24.6 | 18.2 | 59409415 | 125.4 | 41.3 | 127.5 / 37.0 (portable) |
| cellbender_fromraw_ghirsutum2 | 59639633 | 150.4 | 140.6 | 33.9 | 0.6 | 12.2 | 21.0 | 30.9 | 59406621 | 137.1 | 44.1 | 153.1 / 38.8 (portable) |

- The final run sets the peak on the large objects. For soupx_emptydrops_sscovelli it is the full `fields::rdist` matrix, 89.7 GB for 105,879 cells (with doublets). For cellbender_fromraw_ghirsutum2 it is the final SCTransform: its dense residual matrix over 57,899 genes × 79,097 cells is 36.6 GB, vst holds more than one such matrix at a time, and sstat already showed 106 GiB during it, above the 50.1 GB distance matrix. With the 10,000-cell subsample, the sweep peaks at 18–31 GiB.
- The direct runs' exact peaks exceed the workflow jobs' sampled MaxRSS for ghirsutum2 (150.4 GiB against 137.1) and sscovelli (126.9 against 125.4), so the workflow's true peaks were at least as high. Wall times are within minutes of the workflow jobs' (which also ran Seurat afterwards): from 2.5 minutes more (4500871lane) to 10 minutes less (ghirsutum2).
- The instrumented runs cost at most 2.7 GiB more memory and up to 5.5 minutes more time.

## Instrumented runs

[`instrumented_doubletfinder.R`](../../validation/R/instrumented_doubletfinder.R) runs the same code as the direct run with passive `trace()` hooks ([`lib/doubletfinder_hooks.R`](../../validation/R/lib/doubletfinder_hooks.R)) in the installed functions: DoubletFinder's `paramSweep`, `parallel_paramSweep`, `summarizeSweep`, `bimodality_coefficient` and `doubletFinder`; Seurat's `SCTransform.Seurat`, `SCTransform.default` and `RunPCA.default`; sctransform's `vst` (at entry, before its gene sampling and at exit), `get_model_pars` and `reg_model_pars`; future.apply's `future_lapply` (one RNG checkpoint per fitting bin); glmGamPoi's `glm_gp` and its group-wise estimators; `irlba`; KernSmooth's `bkde`. `trace()` also updates the copies other packages import and the S3 method tables, so the call path is the installed code. The hooks read values and write files; they draw no random numbers and change nothing. In a namespace, `trace()` does not reach the copy exported to `package:DoubletFinder`, so the instrumented run calls `DoubletFinder::paramSweep` etc.; the direct run calls the attached functions, like the workflow.

**Proof that the hooks change nothing.** For every object with both runs (tier T and soupx_tenx_SRR18256555 in three configurations, soupx_tenx_10xmousebrain in two, soupx_emptydrops_4500871lane in one: 21 pairs), the instrumented run's `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv` are byte-identical to the direct run's, and so is the `result` section of `summary.json`, which includes a digest of the final `.Random.seed`. Every instrumented run also checks that each hook fired as often as expected (`stages.json`: `hook_counts`, `problems`).

**What is dumped** (`stages/index.tsv` maps every block of every array to its context and stage; contexts `sweep_pN<pN>`, `summary:<sweep item>`, `final`):

| Stage | Arrays |
|---|---|
| S1 | sweep subsample (n > 10,000), doublet index pairs per pN and final; `.Random.seed` before and after the doublet draws |
| S2a | cell attributes (umi, log_umi), filtered genes, log geometric means, means, variances, detection rates, fitting cells and genes in draw order, the gene-sampling pool, probabilities and density; RNG at vst entry, after cell sampling, at the fits' entry |
| S2b | model_pars before and after the theta floor; per bin from glmGamPoi: initial dispersion and beta, first Newton fit, overdispersion MLE, dispersion trend, second fit, iterations, deviance; RNG at every bin |
| S2c | model_pars_fit, outliers, Poisson genes, smoothing genes; bandwidth and mean cell sum (`stages.json`) |
| S2d | residual mean and variance (all genes), top features and VariableFeatures in order, the PCA input (`scale.data` rows in VariableFeatures order; tier T and the small and medium objects), the clipped residual matrix y (tier T) |
| S3 | embeddings, loadings, sdev, irlba's d, iterations; RNG at irlba entry and at RunPCA exit |
| S4 | per real cell, a digest of its first `ind` (sweep) or k+1 (final) `order()` positions; the positions themselves for tier T |
| S5 | pANN per pK for each sweep pN; final pANN |
| S6 | bkde grid and density and the density passed to `bimodality_coefficient`, for every (pN, pK) |
| S7 | final Doublet calls |

S4 digest, per cell with positions `pos[1..L]` (R's 1-based `order()` result): `h1 = Σ r·pos[r]` and `h2 = (Σ (pos[r]·w[r] mod P)) mod P`, with `P = 2^31 − 1` and `w[r] = 16807^r mod P`; every intermediate is an integer below 2^53, so the digest is exact in doubles and in u64.

Structural facts the dumps confirmed on every run (`stages.json`): vst sees all input features in input order; the final SCT `scale.data` equals `SCTransform.default`'s clipped and centred `y`; `VariableFeatures` equal the top 3,000 by residual variance (Seurat 5.5.1 re-derives them in `SCTransform.StdAssay`, with the same order); the PCA input's cells are vst's cells.

## Findings for the port

- **Explicit zeros.** The SoupX objects store zeros as entries of the counts matrix (tier T: 2,686 and 809; tier P: 74,901 to 3,602,702). sctransform sets `min_variance = (median(umi@x)/5)^2`, and `umi@x` includes them, so the reader in `scio` must keep stored zeros, and the doublet construction must keep whatever structure Matrix's `(x1 + x2)/2` and `cbind()` produce. In these data the median stayed 1 (`min_variance` 0.04 in every context), but the port must reproduce the rule, not the value.
- **Counts are integer-valued** in every workflow input (SoupX `roundToInt = TRUE`, decontX `round()`, CellBender), so the Matrix Market files are exact as integers; only the doublets are half-integers.
- **SCTransform in Seurat 5.5.1.** `SCTransform.StdAssay` runs `SCTransform.Assay` → `SCTransform.default` on a v3 assay built from the counts layer, then recomputes `scale.data` through `FetchResiduals` and re-derives `VariableFeatures` through `VariableFeatures.SCTAssay`. The dumps show both are no-ops here: `scale.data` is `SCTransform.default`'s clipped and centred `y`, and `VariableFeatures` are the top 3,000 genes by decreasing residual variance with ties in gene order. The PCA input is `scale.data` restricted to `VariableFeatures`, in that order (rows with zero variance would be dropped; none were).
- **glmGamPoi takes its group-wise path.** With sctransform's design `~ 1` and offset `log(10^log_umi)`, glmGamPoi 1.22.0 uses `estimate_betas_roughly_group_wise` and `fitBeta_one_group` (Newton) twice, not Fisher scoring; `S2b_glm.npy` holds every step per gene. The Cox–Reid term is a scalar for an intercept-only design, and no S2b array differed between OpenBLAS kernels.
- **Where BLAS matters.** Between the Haswell and SapphireRapids kernels, only S2d (the residuals, `exp(tcrossprod(model_pars, regressors))`) and S3 (irlba) differ, in the last bits. irlba's components can change sign between kernels; distances, and so S4–S7, do not depend on the sign.
- **RNG.** All `.Random.seed` checkpoints (up to 120 per run: before and after each doublet draw, at vst entry, after cell sampling, at every fitting bin, around irlba) are identical across kernels and CPU types, so irlba's draw count did not depend on the kernel in any run.
- **Position 1 is not always the cell itself.** In cellbender_fromraw_ghirsutum2, 11 real cells in 5 groups have bitwise-identical PCA embeddings in the final run (cells with 1–4 UMIs whose counts fall outside the variable features, so their `scale.data` columns coincide; some have different counts). For the 6 of them with a lower-index twin, `order()`'s stable tie-break puts the twin at position 1 and the cell itself at position 2, where it counts as a neighbour (`first_position_is_self` = 59,317 of 59,323; the same happened to one cell in four sweep iterations). Everywhere else, position 1 was the cell. The port's selection by (distance, index) reproduces this; ghirsutum2's S4 digests check it.
- **R's decimal parser** (`R_strtod`) is not correctly rounded (long-double accumulation, then rounding to double). The text outputs use 17 significant digits, which Rust's parser reads back exactly (`crates/doubletfinder/tests/m2_fixtures.rs` checks every final pANN against its `.npy` copy); R-side checks compare text, not R-parsed non-integers.
- **OpenBLAS kernel selection.** `OPENBLAS_CORETYPE=SapphireRapids` and a literal `auto` select `Cooperlake` on these CPUs; detection (variable unset) selects `SapphireRapids`. Rust runs in the production profile must leave the variable unset as well.

## Deviations from the plan

- **Reference inputs are the exports, not the `.rds`.** The reference runs rebuild the Seurat object from the exported Matrix Market and TSV files, so references can be regenerated without the workflow's objects (plan A5b, "fixtures are recipes"). This is equivalent: the export checks the counts bit for bit, DoubletFinder with `sct = TRUE` reads nothing else from the object but the cell names and `seurat_clusters`, and every run reproduces its saved workflow output. `reference_doubletfinder.R --input <file>.rds` runs on the workflow's object instead, for diagnosis.
- **Instrumented runs stand in for some direct runs.** Having shown that instrumented runs equal direct runs (21 pairs, both profiles, four object sizes), the portable profile of soupx_emptydrops_4500871lane, soupx_emptydrops_sscovelli and cellbender_fromraw_ghirsutum2 comes from instrumented runs only; their production profile has a direct run (the cross-check and the memory and time baseline). The large objects have no production-profile instrumented run: their floor is measured end to end (portable instrumented against production direct), and per stage on the smaller objects.
- **Dump levels.** The clipped residual matrix `y` (all genes × all cells) is dumped for tier T only; the PCA input (`scale.data`, 3,000 genes × all cells) for tier T and the small and both medium objects (including medium 2, 4500871lane); for the large objects, residual statistics per gene but no residual-sized matrix. S4 order positions are dumped in full for tier T; elsewhere a per-cell digest (the plan's "compact digest").
- **Fixtures.** `S6_bkde.npy` and the per-context matrices stay in `data/` to keep `tests/fixtures/` within 20 MB (17.4 MB in total now).
- **Rust code.** M2 ports nothing, but adds two pieces of infrastructure: `scrust fingerprint` (the A5b fingerprint, also for future Rust runs) and a fixture test, `crates/doubletfinder/tests/m2_fixtures.rs`.
- **Plan details corrected by the sources and dumps.** `SCTransform.StdAssay` does not call `set.seed` itself (`.Seurat`, `.Assay` and `.default` do: three calls, as B2 says); glmGamPoi fits the intercept-only model group-wise, not by Fisher scoring; Seurat 5.5.1 re-derives `VariableFeatures` and `scale.data` after `SCTransform.default` (with no effect, see above). cellbender_fromraw_ghirsutum2's counts matrix has 99,452 features; plan B4's "57,899 genes" is the number SCTransform keeps (detected in at least 5 cells) in the final run.
- **CPU types.** Two (Sapphire Rapids, Cascade Lake), as the criterion asks; the AMD Genoa runs never got a node (cancelled after about 1.5 hours pending).

## Reproducing

From the repository's top directory, after `source env.sh` and `cargo build --release -p scrust`:

```bash
sbatch validation/slurm/export_doubletfinder.sbatch soupx_tenx_test
sbatch validation/slurm/reference_doubletfinder.sbatch portable direct soupx_tenx_test
sbatch validation/slurm/reference_doubletfinder.sbatch portable instrumented soupx_tenx_test
bash validation/slurm/compare_doubletfinder.sbatch soupx_tenx_test
```

Memory and time limits per object are in `validation/slurm/README.md`; production-profile runs keep the default `--partition=sapphire`.

The tables in this report come from `Rscript validation/R/report_doubletfinder.R`, and the manifest's DoubletFinder part from `Rscript validation/R/manifest_doubletfinder.R` (both read plain files only). SLURM logs of every M2 job are in `data/doubletfinder/logs/`.
