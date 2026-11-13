# validation/R

R scripts that make scRust's reference data (docs/PLAN.md A4). They run in the R reference env ([`envs/r-reference/`](../../envs/r-reference/README.md)). Scripts that read `.rds` files run only through sbatch wrappers in [`validation/slurm/`](../slurm/README.md), never on a login node.

[`rcompat_vectors.R`](rcompat_vectors.R) makes the R vectors for the `rcompat` crate (stage S0): R's random-number stream and `.Random.seed` checkpoints, `sample`, `rnorm`/`qnorm`, `round`, `^`, `seq`, `order`, `sum`/`mean`/`var`/`sd` and number formatting (M1); `as.numeric()` of strings (R's own parser, `R_strtod`, including strings it does not round correctly), `fft`, `dnorm` and `approxfun` (M3); over many seeds, sizes and edge cases. It reads no `.rds`, so it may run on a login node:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
  "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
```

Three kinds of script per tool:

- `export_<tool>.R`: one-time conversion of workflow artifacts into plain files: the tool's inputs, and the fields of the workflow's saved outputs that are compared (for DoubletFinder: the doublet set, singlet pANN, pK and nExp). Counts become Matrix Market, tables TSV, and dense arrays and `.Random.seed` checkpoints `.npy` (written by a small base-R helper).
- `reference_<tool>.R`: a direct run of the original tool through a thin script with the workflow's parameters and seeds, writing complete outputs as plain files. Each run must reproduce the corresponding saved workflow output (or test reference output), which shows the script calls the tool exactly as the workflow did.
- `instrumented_<tool>.R`: the same run, also dumping the tool's internal stages (for DoubletFinder: doublet indices, SCTransform, PCA, neighbours and RNG checkpoints), so a difference can be traced to the stage that caused it.

## DoubletFinder (milestone M2)

Datasets (tier T: the five workflow test objects; tier P: five production objects) and their pipeline paths are in [`../doubletfinder_datasets.tsv`](../doubletfinder_datasets.tsv). Results: [`docs/validation/doubletfinder_m2.md`](../../docs/validation/doubletfinder_m2.md).

- [`export_doubletfinder.R`](export_doubletfinder.R): from each input object, the counts DoubletFinder 2.0.6 reads with `sct = TRUE` (the RNA assay's `counts` layer) as 10x Matrix Market, entries exactly as stored, explicit zeros included, cells and features in object order; `seurat_clusters` and its level order. From each saved workflow output: the kept barcodes with their pANN, the doublet set (input minus output barcodes), and pK and nExp parsed from the column names. It re-reads the matrix with the reference runs' reader and checks `identical(..., num.eq = FALSE)` against the object's own matrix.
- [`reference_doubletfinder.R`](reference_doubletfinder.R): the direct run, with the workflow's calls and seeds (`set.seed(12345)`, `paramSweep(PCs = 1:10, sct = TRUE)`, `summarizeSweep(GT = FALSE)`, `find.pK`, the pK via `as.numeric(as.character(...))`, the homotypic nExp, `set.seed(12345)`, `doubletFinder(pN = 0.25, ...)`). Writes plan B3's `cells.tsv` (every input cell), `sweep_stats.tsv`, `bcmvn.tsv` and `summary.json` (a deterministic `result` section and a `run` section with timings, peak RSS and the environment fingerprint), plus `find_pK.pdf`, the plot `find.pK` draws. Its input is an export directory (the Seurat object is rebuilt with `CreateSeuratObject()` from the identical counts) or, for diagnosis, the workflow's `.rds`.
- [`instrumented_doubletfinder.R`](instrumented_doubletfinder.R): the same run, also dumping stages S1–S7 to `stages/` (`.npy`, with `index.tsv`, `rng.npy`/`rng.tsv` and `stages.json`). Mechanism: `trace()` hooks in the installed functions that only read values and write files ([`lib/doubletfinder_hooks.R`](lib/doubletfinder_hooks.R) documents the hooks, the array layout and the index conventions). Its final outputs are byte-identical to the direct run's, checked for every object that has both.
- [`compare_doubletfinder.R`](compare_doubletfinder.R): the M2 checks on plain files: a run against the saved workflow output (doublet set, singlet pANN bitwise, column names), and two runs against each other (final files, the `result` section, every stage array block by block, RNG checkpoints), driven by [`../slurm/compare_doubletfinder.sbatch`](../slurm/compare_doubletfinder.sbatch).
- [`manifest_doubletfinder.R`](manifest_doubletfinder.R): regenerates the DoubletFinder part of `validation/manifest.toml` from the files on disk.
- [`lib/io.R`](lib/io.R): the `.npy` writer and reader, exact text output (`%.17g`), JSON, the 10x Matrix Market writer and exact reader, peak RSS and the environment fingerprint (which calls `scrust fingerprint` for the CPU, OpenBLAS core and threads, and glibc). [`lib/doubletfinder.R`](lib/doubletfinder.R): the workflow's call sequence, shared by the direct and instrumented runs.

Exactness: doubles are written as text with 17 significant digits (glibc's correctly rounded printf; a correctly rounded reader, such as Rust's, gets the exact double back). R's own decimal parser is not correctly rounded, so checks compare the written text, never R-parsed non-integers. The workflow's counts are integer-valued (SoupX `roundToInt = TRUE`, decontX `round()`, CellBender), which R parses exactly.

Rules:

- R runs with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB` and `PYTHONNOUSERSITE=1`; sourcing `env.sh`, as every wrapper does, sets this up.
- OpenBLAS's kernel and thread count come from `env.sh` (portable profile by default, see the profiles in [`validation/slurm/README.md`](../slurm/README.md)).
- Seeds and parameters are the workflow's, as written in its scripts (for DoubletFinder, `workflow/scripts/doubletfinder.R` with seed 12345).
- Outputs go to `data/` (large, gitignored) or `tests/fixtures/` (small, committed), and every file is listed in `validation/manifest.toml` with its sha256.

## DoubletFinder core (milestone M3)

The Rust core runs on the instrumented runs' own embeddings and is compared with them by `scrust validate doubletfinder` ([`../slurm/m3_doubletfinder.sbatch`](../slurm/m3_doubletfinder.sbatch)); no new R run was needed beyond the `rcompat` vectors. Results: [`docs/validation/doubletfinder_m3.md`](../../docs/validation/doubletfinder_m3.md).

- [`report_doubletfinder_m3.R`](report_doubletfinder_m3.R): the report's tables from the M3 outputs (the `validate.json` of every Rust run, the direct R runs, the saved workflow outputs and the benchmark's `/usr/bin/time -v`). Reads only plain files.

## The whole tool (milestone M6)

`scrust doubletfinder` runs from the counts, the cluster labels and the seed; it is compared with every direct and instrumented reference run of M2 and with the saved workflow outputs by `scrust validate doubletfinder` ([`../slurm/m6_doubletfinder.sbatch`](../slurm/m6_doubletfinder.sbatch)), and benchmarked against M2's direct production runs ([`../slurm/m6_benchmark.sbatch`](../slurm/m6_benchmark.sbatch)). No new R run was needed. Results: [`docs/validation/doubletfinder.md`](../../docs/validation/doubletfinder.md).

- [`report_doubletfinder_m6.R`](report_doubletfinder_m6.R) `<label> <bench label>`: the report's tables from the M6 outputs (each Rust run's `vs_<R run>.json` and `vs_saved.json`, its `summary.json` and `time_v.txt`; the benchmark's runs and the direct R runs' `time_v.txt` and `summary.json`). Reads only plain files.

## The singlets' counts (milestone M6b)

`scrust doubletfinder` also writes `singlets/`: the counts without the doublets, as a new 10x directory (docs/PLAN.md A3 and B3). The workflow's equivalent is the RNA `counts` layer of its saved DoubletFinder output, since `doubletfinder.R` keeps the singlets with `subset()`.

- [`export_doubletfinder_singlets.R`](export_doubletfinder_singlets.R): exports that layer from each saved output into `data/doubletfinder/<object>/saved/singlets/`, next to the M2 export: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` (`lib/io.R`'s writer, with the input export's feature lines) and `export.json`. It checks three things:
  - the `.rds` is the one M2 exported (sha256);
  - its cells are the kept cells, in order;
  - the directory reads back `identical()`.

  It also records whether the layer equals the input export's counts restricted to the kept cells, slot for slot. Run it through [`../slurm/export_doubletfinder_singlets.sbatch`](../slurm/export_doubletfinder_singlets.sbatch).
- `scrust validate doubletfinder --saved` compares `singlets/` with that export: the barcodes in order, the feature lines, and every stored entry and value of the matrix. [`../slurm/m6_doubletfinder.sbatch`](../slurm/m6_doubletfinder.sbatch) runs this check with the others.
- [`read10x_singlets.R`](read10x_singlets.R) `<label> <profile> <object>...` reads the Rust runs' `singlets/` back with R's readers, as a downstream R step would: `lib/io.R`'s exact reader and `Seurat::Read10X()`. It compares each with the export, slot for slot, and records whether `Read10X` keeps every stored entry. It reads only plain files.
