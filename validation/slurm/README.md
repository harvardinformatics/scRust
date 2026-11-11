# validation/slurm

sbatch wrappers for everything that runs on compute nodes (docs/PLAN.md A2 and A4):

- the export scripts and the direct and instrumented reference runs in `validation/R/`; anything that reads `.rds` files runs this way, never on a login node;
- tier-P and full-panel validation (`scrust validate ...`, from M2/M3 on);
- memory and time benchmarks.

Start each wrapper from [`template.sbatch`](template.sbatch), which:

- sources the node guard, `/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh`. It reads one byte of a file stored on holylfs05 OST 51 and exits with code 75 if that hasn't returned within 30 s, instead of letting the job hang unkillably on a node that lost that storage target;
- excludes the nodes known to have lost it (`--exclude=holy8a32608,holy8a24102`) and runs under the `informatics` account (`-A informatics`);
- sources [`env.sh`](../../env.sh): the project env, OpenBLAS pinned to one kernel and one thread for R and Rust alike, and no user-level R, Perl or Python libraries;
- prints an environment fingerprint (node, CPU model, glibc, OpenBLAS kernel and threads, rustc, commit) to the log.

## Reference profiles (A5b)

| Profile | `OPENBLAS_CORETYPE` | Nodes | Use |
|---|---|---|---|
| portable | `Haswell` (env.sh's default) | any x86-64 node with AVX2 | development, CI, fixtures that can be regenerated on any cluster |
| production | `auto` (OpenBLAS detects `SapphireRapids`) | Sapphire Rapids: `--partition=sapphire`, or `--constraint=sapphirerapids` elsewhere | comparisons with the saved production outputs; all 150 production DoubletFinder jobs ran on Sapphire Rapids |

Use one CPU per job for anything that calls BLAS (`--cpus-per-task=1`), as the production jobs did.

Logs go to `<job name>-<job id>.sbatch.log` in the submit directory, which `.gitignore` excludes.
