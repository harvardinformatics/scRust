# validation/R

R scripts that make scRust's reference data (docs/PLAN.md A4). They run in the R reference env ([`envs/r-reference/`](../../envs/r-reference/README.md)). Scripts that read `.rds` files run only through sbatch wrappers in [`validation/slurm/`](../slurm/README.md), never on a login node.

[`rcompat_vectors.R`](rcompat_vectors.R) makes the R vectors for the `rcompat` crate (stage S0, milestone M1): R's random-number stream and `.Random.seed` checkpoints, `sample`, `rnorm`/`qnorm`, `round`, `^`, `seq`, `order`, `sum`/`mean`/`var`/`sd` and number formatting, over many seeds, sizes and edge cases. It reads no `.rds`, so it may run on a login node:

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
  "$R_ENV/bin/Rscript" validation/R/rcompat_vectors.R tests/fixtures/rcompat
```

Three kinds of script per tool:

- `export_<tool>.R`: one-time conversion of workflow artifacts into plain files: the tool's inputs, and the fields of the workflow's saved outputs that are compared (for DoubletFinder: the doublet set, singlet pANN, pK and nExp). Counts become Matrix Market, tables TSV, and dense arrays and `.Random.seed` checkpoints `.npy` (written by a small base-R helper).
- `reference_<tool>.R`: a direct run of the original tool through a thin script with the workflow's parameters and seeds, writing complete outputs as plain files. Each run must reproduce the corresponding saved workflow output (or test reference output), which shows the script calls the tool exactly as the workflow did.
- `instrumented_<tool>.R`: the same run, also dumping the tool's internal stages (for DoubletFinder: doublet indices, SCTransform, PCA, neighbours and RNG checkpoints), so a difference can be traced to the stage that caused it.

Rules:

- R runs with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB` and `PYTHONNOUSERSITE=1`; sourcing `env.sh`, as every wrapper does, sets this up.
- OpenBLAS's kernel and thread count come from `env.sh` (portable profile by default, see the profiles in [`validation/slurm/README.md`](../slurm/README.md)).
- Seeds and parameters are the workflow's, as written in its scripts (for DoubletFinder, `workflow/scripts/doubletfinder.R` with seed 12345).
- Outputs go to `data/` (large, gitignored) or `tests/fixtures/` (small, committed), and every file is listed in `validation/manifest.toml` with its sha256.
