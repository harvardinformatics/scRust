# AGENTS.md

Instructions for AI coding agents working in scRust; a summary for human contributors too. [docs/PLAN.md](docs/PLAN.md) is the source of truth. Read its Part A and the part for the tool at hand before changing anything. This page lists the rules every change must keep, each with the section of the plan that explains it.

## What scRust is

- Standalone Rust reimplementations of the tools in the scRNA-seq preprocessing workflow. They are validated against the original tools, not against a rebuilt workflow.
- One binary, `scrust <tool>`. GPL-3.0 ([LICENSE](LICENSE); upstream licenses in [NOTICE](NOTICE)).
- Status and the order of the tools: [README.md](README.md) and PLAN A7.

## Reading order for a task

1. This page.
2. PLAN Part A.
3. The part for the tool (B: DoubletFinder; C: the clustering crate).
4. The latest report in [docs/validation/](docs/validation/).
5. The brief or issue you were given.

## Standing rules for every tool and milestone

1. **Every output of the original, plus count data out** (A3). A port keeps all of the original tool's outputs. Every tool that processes counts also writes its result as new 10x count data. Never remove or rename an output.
2. **Normalised data is never written as 10x count data** (A3, C3f). It goes into the export, marked as what it is. The loaders put it only into slots meant for its scale (C3g).
3. **In-memory hand-offs** (A3):
   - each tool is one library function, from typed inputs in memory to typed outputs in memory, with no file I/O inside;
   - what passes between tools uses `scio`'s shared types;
   - the output value holds everything the files hold;
   - the command is a thin read–call–write layer;
   - a chain can free each stage;
   - the in-process route gives files byte-identical to the file route.
4. **Exact against the original** (A4, A5, A5b).
   - Under the reference's OpenBLAS kernel, outputs must equal the original tool's bit for bit.
   - Across kernels, checks that pass through BLAS are graded against the R-vs-R floor in `validation/tolerances.toml`. Everything else stays exact.
5. **The cycle for each tool** (A7):
   1. Read the exact code path in the pinned env.
   2. Write the reference scripts.
   3. Port it as a library.
   4. Validate stage by stage.
   5. Validate on T, P and the full panel.
   6. Benchmark, adding the rows to [docs/benchmarks.md](docs/benchmarks.md).
   7. Update the user documentation (M14b).
6. **Records go with the work:**
   - every reference dataset and validation output in `validation/manifest.toml` (sizes, SHA-256);
   - every settled criterion in `validation/tolerances.toml`;
   - every result in `docs/validation/`;
   - every adapted upstream source in `NOTICE` and `third_party/SHA256SUMS`.

## Decisions that are the user's

Ask before any of these:
- loosening a tolerance, or accepting a difference you cannot trace;
- changing the output manifest's schema (`scrust-files`; [docs/output_manifest.md](docs/output_manifest.md));
- making new R reference runs;
- removing or renaming an output;
- starting work beyond the task you were given.

## Environment and HPC (A2, A5b)

- **Envs:** hermetic, project-local conda envs made from the yaml files in `envs/`. Run `source env.sh` for the Rust toolchain. Nothing is installed into the home directory or `~/.local`.
- **Python:** run with `PYTHONNOUSERSITE=1` and `python -I`.
- **R:**
  - Always run with `env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB`.
  - The production R env is read-only. R packages under test install into `.rlib/`.
- **`.rds` files:** anything that reads one runs through sbatch or srun, never on a login node. The sbatch wrappers in `validation/slurm/` source the node guard and exclude the two nodes known to lose a Lustre OST.
- **BLAS kernel:** `env.sh` pins `OPENBLAS_CORETYPE` (Haswell, the portable profile) and `OPENBLAS_NUM_THREADS=1`. The production profile uses OpenBLAS's auto-detection (`OPENBLAS_CORETYPE=auto`) on Sapphire Rapids nodes. Don't force `SapphireRapids` by name: in a test it reported Cooperlake.
- **The workflow's repository** (scRNAseq-preprocessing) is read-only: never write there or run Snakemake there.

## Tests and git

- **Tests:** `cargo test --workspace` must pass in debug and in release, with no warnings, before every push.
  - Tests whose fixtures depend on BLAS check the kernel and stop under another one.
  - The loaders' tests run through `validation/slurm/m14_tests.sbatch`.
- **Data:** large data lives in `data/` (gitignored) and is listed in `validation/manifest.toml`. Committed fixtures in `tests/fixtures/` stay small (about 20 MB in all).
- **Git:**
  - One branch per milestone (`mN_dev`) from `main`, merged by pull request.
  - Commit only the files your task changed.
  - Never rewrite pushed history.
