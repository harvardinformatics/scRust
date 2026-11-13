# DoubletFinder in Rust: the full panel (milestone M7)

M7 is the plan's full-panel confirmation of Step 1 (docs/PLAN.md B4 to B6). `scrust doubletfinder` ran on all 150 production DoubletFinder objects of the preprocessing workflow, and each run was compared with the output the workflow saved. R was not run: the saved outputs are the references, after a one-time export. The comparison is M6's against the saved outputs (the doublet set, the kept cells' pANN, pK and nExp) plus M6b's singlets' counts (`singlets/` against the saved object's RNA counts). The five tier-P objects of M2 to M6b are among the 150; the other 145 are new to the validation.

## Outcome

Every criterion of plan B6's row M7 and of B5 is met, on all 150 objects.

| Done when (plan B6 M7, B5, and the milestone's brief) | Result |
|---|---|
| All 150 objects compared | **Yes.** 150 exports (145 new, 5 from M2 and M6b) and 150 Rust runs, each compared with its saved workflow output by `scrust validate doubletfinder --saved`. |
| The doublet set, the kept cells' pANN, pK and nExp exact | **Yes, all 150.** The same 353,762 doublets. The pANN of all 2,217,328 kept cells identical as R wrote them (17 significant digits, in order). The same pK labels, nExp and column names. |
| The singlets' counts identical to the saved objects' | **Yes, all 150.** The barcodes in order, the feature lines, the header field, and every stored entry and value of the saved RNA counts: 3,285,047,618 entries, 120,288,482 of them explicit zeros. |
| Each run on the CPU type of its workflow job, with the kernel OpenBLAS detects there | **Yes.** All 150 workflow jobs ran on Sapphire Rapids nodes (38 nodes; sacct and the nodes' SLURM features). Every Rust run ran on a Sapphire Rapids node, checked against its job's CPU type, with `OPENBLAS_CORETYPE` unset: OpenBLAS reported `SapphireRapids` in all 150. |
| A summary table | [Results](#results): per object the cells, singlets, the checks, pK, nExp, the singlets' counts, the CPU types and kernel, wall time and peak RSS; totals by check and by arm. The same per object as TSV: [`doubletfinder_m7.tsv`](doubletfinder_m7.tsv). |
| `tolerances.toml`, PLAN and README updated | **Yes.** `[doubletfinder.full_panel]` settles the full panel; `[doubletfinder.saved_outputs]` and `[doubletfinder.singlets]` record it. PLAN B4, B5 and B6 add the singlets' counts to the full panel. The manifest lists the 145 new exports (2,465 entries; no existing entry changed). |
| `cargo test --workspace` passes | **Yes,** all 160 tests, without warnings. On main one test failed in debug builds: an unsigned underflow in the 10x writer that release builds wrap around. It is fixed, without changing any output ([Findings](#findings)). |

In all: 1,350 of 1,350 checks identical (nine per object), 3,295,570,596 values compared, none differing.

## The panel

**The objects.** The workflow runs DoubletFinder once per sample and arm: 30 samples (the samplesheet's 31, less `cteleta`, which `config/config.yaml` excludes) by 5 arms. The step's input and saved output, in the pipeline repo:

| arm | input | saved output |
|---|---|---|
| soupx_tenx, soupx_emptydrops | `results/soupx/seurat_soupx_<tenx\|emptydrops>_<sample>.rds` | `results/doubletfinder/seurat_doubletfinder_<arm>_<sample>.rds` |
| decontx_raw_tenx, decontx_raw_emptydrops | `results/decontx_raw/seurat_decontx_raw_<tenx\|emptydrops>_<sample>.rds` | as above |
| cellbender_fromraw | `results/cellbender_fromraw/seurat_cellbender_fromraw_<sample>.rds` | as above |

[`validation/doubletfinder_panel.tsv`](../../validation/doubletfinder_panel.tsv) lists the 150 objects with these paths and the job that wrote each saved output. [`validation/slurm/doubletfinder_panel.sh`](../../validation/slurm/doubletfinder_panel.sh) writes the table from the Snakemake logs, sacct and scontrol, and reproduces it byte for byte.

The panel spans 1,120 to 79,409 cells and 10,172 to 99,452 features. 98 objects have more than 10,000 cells, so the sweep subsamples them. k runs from 15 to 20,349. The chosen pK takes 32 of the grid's 33 values, all but 5e-04.

**The jobs that wrote the saved outputs.** They come from the workflow's rerun #2, made of two Snakemake runs:

- **The run started on 8 November at 21:41** (log `.snakemake/log/2026-11-08T214133.605174.snakemake.log`, master job 59050095). It ran 148 of the 150 DoubletFinder jobs. It ended failed, because the steps that build two of the inputs failed: `soupx_emptydrops` for soupx_emptydrops_sscovelli, `cellbender2seurat` for cellbender_fromraw_ghirsutum2.
- **Its follow-up, started on 9 November at 22:02** (log `2026-11-09T220250.045638.snakemake.log`, master job 59391964, completed on 10 November). It rebuilt those two inputs and ran their DoubletFinder jobs, 59409415 and 59406621.

How each saved output was traced to its job:

- **The Snakemake logs** give each DoubletFinder job's SLURM id ("submitted with SLURM jobid"). There are 150 jobs, one per object across the two runs. None failed, so none was retried.
- **sacct.** All 150 are `COMPLETED` with exit code 0, one CPU each, in the `sapphire` partition, on 38 nodes. Step `.0`, the `srun` that runs `doubletfinder.R`, holds the job's peak RSS.
- **The files.** Every input `.rds` was last written before its job started. Every saved `.rds` was last written between 3.5 minutes after its job started and 93 s after the job's end in sacct, and every log during its job. No Snakemake run came later. So each saved output is its job's, and each input the one that job read. Each export also records the sha256 of the `.rds` it read, and the singlets' export checks the saved object's against the first export's.

**CPU type and OpenBLAS kernel.** All 38 nodes carry the SLURM features `intel,holyndr,sapphirerapids,avx,avx2,avx512`. The workflow sets no `OPENBLAS_CORETYPE`: neither `profiles/slurm/config.yaml` nor the rules or scripts do. So OpenBLAS detected the CPU, which gives `SapphireRapids` on these nodes, as M2 measured. Each Rust run therefore ran on a `sapphire` node with `OPENBLAS_CORETYPE` unset (`OPENBLAS_CORETYPE=auto source env.sh`). Before every run, the wrapper checked the node's features against the object's `workflow_cpu`, and the core OpenBLAS reported against `SapphireRapids`.

## How the comparison works

- **Exports** ([`validation/slurm/m7_export.sbatch`](../../validation/slurm/m7_export.sbatch), one array task per object). For each of the 145 new objects, the task runs M2's and M6b's export wrappers, which now also find objects in the panel table:
  - `export_doubletfinder.R`: the input object's RNA counts as stored, `seurat_clusters` and its level order; the saved output's kept barcodes with their pANN, the doublet set, and the column names with pK and nExp;
  - `export_doubletfinder_singlets.R`: the saved output's RNA counts.

  Both read their files back and check them bit for bit. The singlets' export also checks that its cells are the kept cells, in order, and records whether the saved counts equal the input's restricted to the kept cells, slot for slot. The exports read `.rds` files but call no BLAS, so they ran on any node. Job 60160712 made 143 of them; the trial job 60159280 made the other 2, with the same wrapper. Each took 1 CPU, up to 9.7 minutes and 8.6 GiB per export, and the 145 objects' exports hold 20 GB. The five P objects kept their M2 and M6b exports.
- **Rust runs** ([`validation/slurm/m7_panel.sbatch`](../../validation/slurm/m7_panel.sbatch), one array task per object). `scrust doubletfinder` on the export's counts and cluster labels, without stage dumps, on 8 threads, then `scrust validate doubletfinder --saved`. Label `m7-b8c8a8c`, array job 60163542, binary `data/doubletfinder/logs/bin/scrust-m7-b8c8a8c` (sha256 `808f630007c04a19`), built from commit b8c8a8c: main plus the writer fix ([Findings](#findings)).
- **The checks.** Nine per object, each against the export of the saved output:
  - the pANN and classification column names, which carry pN, pK and nExp as R wrote them;
  - the pK label and nExp;
  - the doublet set: the input's barcodes the workflow did not keep, against Rust's `Doublet` calls;
  - the kept cells' pANN: barcode and pANN text, in order, against R's 17 significant digits, so bitwise;
  - the singlets' counts: the barcodes in order, the feature lines, and the matrix (its header field, dimensions, every stored entry, explicit zeros included, and every value, bitwise).

  An object passes only when all nine checks are made and identical. That is stricter than the validator's exit status, which fails on a difference but not on a check it could not make. The wrapper and the report both require nine of nine.
- **A second binary.** Main's own binary ran the whole panel too, before the fix: label `m7-d95f7f7`, array job 60160715, `scrust-m7-d95f7f7` (sha256 `dfd3b9c697dcf263`). It is byte-identical to M6b's validated `scrust-m6b-13592a6`. It passed all 150 objects in the same way, and its outputs are identical to the label of record's, file for file ([below](#against-the-runs-of-label-m7-d95f7f7)).

## Results

Validation label `m7-b8c8a8c` (`data/doubletfinder/<object>/rust/m7-b8c8a8c/production/`), array job 60163542, 17 Sapphire Rapids nodes (Intel Xeon Platinum 8480+ and 8480CL), 8 threads per run, on 12 November. The tables come from [`validation/R/report_doubletfinder_m7.R`](../../validation/R/report_doubletfinder_m7.R), which reads plain files only. The per-object table is last.

*Cost* compares the Rust runs with the workflow's jobs for context only. Those jobs ran on one CPU and include `doubletfinder.R`'s Seurat processing after DoubletFinder, so they are not a benchmark of the tool; M6 benchmarks Rust against the direct R run instead.

In the per-object table, `= (n)` means identical over n elements: the doublets, the kept cells, or the stored entries of the singlets' counts. *wall* and *peak RSS* are the Rust run's, from `/usr/bin/time -v`, writing `singlets/` included.

### Totals

| objects | compared | passed | checks identical | values compared | values differing | not made or missing |
|---|---|---|---|---|---|---|
| 150 | 150 | 150 | 1,350 of 1,350 | 3,295,570,596 | 0 | 0 |

### By check

| check | objects | identical | values compared |
|---|---|---|---|
| pANN_column | 150 | 150 | 150 |
| classification_column | 150 | 150 | 150 |
| pK_label | 150 | 150 | 150 |
| nExp | 150 | 150 | 150 |
| doublet set | 150 | 150 | 353,762 |
| singlet pANN | 150 | 150 | 2,217,328 |
| singlets barcodes | 150 | 150 | 2,217,328 |
| singlets features | 150 | 150 | 5,733,960 |
| singlets matrix | 150 | 150 | 3,285,047,618 |

### By arm

| arm | objects | passed | cells | cells, total | doublets, total | singlets, total | singlets' stored entries |
|---|---|---|---|---|---|---|---|
| soupx_tenx | 30 | 30 | 2,676 to 34,421 | 411,215 | 57,731 | 353,484 | 666,317,716 |
| soupx_emptydrops | 30 | 30 | 1,122 to 79,409 | 620,233 | 84,530 | 535,703 | 714,364,146 |
| decontx_raw_tenx | 30 | 30 | 2,673 to 34,419 | 410,031 | 57,745 | 352,286 | 603,901,258 |
| decontx_raw_emptydrops | 30 | 30 | 1,120 to 78,687 | 619,341 | 83,401 | 535,940 | 658,713,598 |
| cellbender_fromraw | 30 | 30 | 4,314 to 60,329 | 510,270 | 70,355 | 439,915 | 641,750,900 |
| all | 150 | 150 | 1,120 to 79,409 | 2,571,090 | 353,762 | 2,217,328 | 3,285,047,618 |

### CPU types and OpenBLAS kernels

| side | objects | nodes | CPU type | CPU model | OpenBLAS kernel |
|---|---|---|---|---|---|
| workflow jobs (sacct node, SLURM features) | 150 | 38 | sapphirerapids (150) |  | not recorded; OPENBLAS_CORETYPE unset (production profile) |
| Rust runs (node's SLURM features; scrust fingerprint) | 150 | 17 | sapphirerapids (150) | 8480+ (109), 8480CL (41) | SapphireRapids (150) |

Rust runs on the CPU type of the object's workflow job: 150 of 150.

### Cost

| side | wall, total | wall, median | wall, largest | CPU time, total | peak RSS GiB, median | peak RSS GiB, largest |
|---|---|---|---|---|---|---|
| Rust (scrust doubletfinder, 8 threads, singlets/ included) | 1:53:49 | 43.6 s | 2:21.7 (cellbender_fromraw_ghirsutum2) | 8:18:41 | 1.15 | 6.64 (soupx_emptydrops_amellifera) |
| the workflow's jobs (DoubletFinder plus doubletfinder.R's Seurat processing; 1 CPU) | 32:42:41 | 12:58.0 | 44:08.0 (cellbender_fromraw_ghirsutum2) |  | 20.29 | 137.09 (cellbender_fromraw_ghirsutum2) |

### The saved objects' RNA counts

Equal to the input's counts restricted to the kept cells, slot for slot: 150 of 150 objects. Explicit zeros among the singlets' stored entries: 120,288,482.

### Against the runs of label m7-d95f7f7

| objects run with both binaries | files identical | result sections identical | objects identical |
|---|---|---|---|
| 150 | 900 of 900 | 150 of 150 | 150 of 150 |

### Per object

| object | cells | singlets | doublet set | singlet pANN | pK | nExp | singlets' counts | workflow job: node, CPU type | Rust: node, CPU, kernel | wall | peak RSS GiB | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_4500871lane | 10,492 | 9,004 | = (1,488) | = (9,004) | 0.08 | 1,488 | = (32,373,270) | holy8a32507, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 59.0 s | 2.60 | yes |
| soupx_emptydrops_4500871lane | 26,598 | 23,239 | = (3,359) | = (23,239) | 0.02 | 3,359 | = (35,907,669) | holy8a32606, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 56.6 s | 2.42 | yes |
| decontx_raw_tenx_4500871lane | 10,491 | 9,003 | = (1,488) | = (9,003) | 0.08 | 1,488 | = (30,829,177) | holy8a30309, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 57.4 s | 2.43 | yes |
| decontx_raw_emptydrops_4500871lane | 26,598 | 23,192 | = (3,406) | = (23,192) | 0.02 | 3,406 | = (34,158,798) | holy8a30507, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 54.5 s | 2.26 | yes |
| cellbender_fromraw_4500871lane | 11,530 | 9,887 | = (1,643) | = (9,887) | 0.005 | 1,643 | = (29,427,871) | holy8a30507, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 65.8 s | 2.32 | yes |
| soupx_tenx_563792olfmobaob | 10,895 | 9,340 | = (1,555) | = (9,340) | 0.06 | 1,555 | = (34,604,134) | holy8a30604, sapphirerapids | holy8a24107, 8480+, SapphireRapids | 61.1 s | 2.63 | yes |
| soupx_emptydrops_563792olfmobaob | 12,885 | 11,056 | = (1,829) | = (11,056) | 0.04 | 1,829 | = (34,564,360) | holy8a30302, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 60.0 s | 2.38 | yes |
| decontx_raw_tenx_563792olfmobaob | 10,889 | 9,336 | = (1,553) | = (9,336) | 0.27 | 1,553 | = (32,401,403) | holy8a32606, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 58.5 s | 2.46 | yes |
| decontx_raw_emptydrops_563792olfmobaob | 12,885 | 11,055 | = (1,830) | = (11,055) | 0.09 | 1,830 | = (33,173,572) | holy8a32608, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 57.7 s | 2.24 | yes |
| cellbender_fromraw_563792olfmobaob | 11,062 | 9,490 | = (1,572) | = (9,490) | 0.04 | 1,572 | = (32,792,072) | holy8a32603, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 59.5 s | 2.50 | yes |
| soupx_tenx_SRR18256555 | 3,108 | 2,688 | = (420) | = (2,688) | 0.005 | 420 | = (2,461,479) | holy8a30604, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 13.4 s | 0.42 | yes |
| soupx_emptydrops_SRR18256555 | 3,880 | 3,339 | = (541) | = (3,339) | 0.02 | 541 | = (2,051,247) | holy8a30606, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 14.9 s | 0.51 | yes |
| decontx_raw_tenx_SRR18256555 | 3,107 | 2,686 | = (421) | = (2,686) | 0.02 | 421 | = (2,372,841) | holy8a30412, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 13.2 s | 0.41 | yes |
| decontx_raw_emptydrops_SRR18256555 | 3,880 | 3,348 | = (532) | = (3,348) | 0.02 | 532 | = (1,941,932) | holy8a30607, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 15.2 s | 0.48 | yes |
| cellbender_fromraw_SRR18256555 | 10,602 | 9,140 | = (1,462) | = (9,140) | 0.21 | 1,462 | = (3,285,016) | holy8a32603, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 25.9 s | 0.64 | yes |
| soupx_tenx_SRR18256506 | 3,581 | 3,083 | = (498) | = (3,083) | 0.03 | 498 | = (7,670,095) | holy8a30607, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 23.3 s | 0.78 | yes |
| soupx_emptydrops_SRR18256506 | 4,053 | 3,486 | = (567) | = (3,486) | 0.03 | 567 | = (7,639,229) | holy8a30309, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 22.1 s | 0.81 | yes |
| decontx_raw_tenx_SRR18256506 | 3,576 | 3,078 | = (498) | = (3,078) | 0.03 | 498 | = (7,431,666) | holy8a32606, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 20.6 s | 0.78 | yes |
| decontx_raw_emptydrops_SRR18256506 | 4,035 | 3,467 | = (568) | = (3,467) | 0.005 | 568 | = (7,066,029) | holy8a30407, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 23.1 s | 0.83 | yes |
| cellbender_fromraw_SRR18256506 | 6,173 | 5,311 | = (862) | = (5,311) | 0.15 | 862 | = (6,673,701) | holy8a32604, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 26.8 s | 0.80 | yes |
| soupx_tenx_obOLsample1 | 22,175 | 19,017 | = (3,158) | = (19,017) | 0.17 | 3,158 | = (38,752,116) | holy8a32603, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 82.0 s | 2.42 | yes |
| soupx_emptydrops_obOLsample1 | 26,925 | 23,120 | = (3,805) | = (23,120) | 0.29 | 3,805 | = (37,870,855) | holy8a32508, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 53.9 s | 2.53 | yes |
| decontx_raw_tenx_obOLsample1 | 22,157 | 18,999 | = (3,158) | = (18,999) | 0.24 | 3,158 | = (39,490,095) | holy8a32607, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 51.6 s | 2.30 | yes |
| decontx_raw_emptydrops_obOLsample1 | 26,922 | 23,080 | = (3,842) | = (23,080) | 0.04 | 3,842 | = (35,706,163) | holy8a32607, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 68.7 s | 2.40 | yes |
| cellbender_fromraw_obOLsample1 | 22,918 | 19,650 | = (3,268) | = (19,650) | 0.17 | 3,268 | = (33,013,120) | holy8a32507, sapphirerapids | holy8a28607, 8480CL, SapphireRapids | 63.3 s | 2.19 | yes |
| soupx_tenx_obOLsample2 | 17,493 | 15,010 | = (2,483) | = (15,010) | 0.3 | 2,483 | = (26,715,435) | holy8a32603, sapphirerapids | holy8a24110, 8480+, SapphireRapids | 43.2 s | 1.66 | yes |
| soupx_emptydrops_obOLsample2 | 18,717 | 16,066 | = (2,651) | = (16,066) | 0.14 | 2,651 | = (26,438,641) | holy8a32507, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 43.3 s | 1.68 | yes |
| decontx_raw_tenx_obOLsample2 | 17,471 | 15,004 | = (2,467) | = (15,004) | 0.25 | 2,467 | = (26,252,301) | holy8a32508, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 42.7 s | 1.61 | yes |
| decontx_raw_emptydrops_obOLsample2 | 18,717 | 16,057 | = (2,660) | = (16,057) | 0.29 | 2,660 | = (27,354,740) | holy8a32602, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 43.0 s | 1.63 | yes |
| cellbender_fromraw_obOLsample2 | 19,301 | 16,543 | = (2,758) | = (16,543) | 0.26 | 2,758 | = (24,544,338) | holy8a32505, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 42.0 s | 1.57 | yes |
| soupx_tenx_chickenamygdala | 8,654 | 7,453 | = (1,201) | = (7,453) | 0.29 | 1,201 | = (7,585,456) | holy8a32503, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 30.5 s | 0.70 | yes |
| soupx_emptydrops_chickenamygdala | 47,621 | 41,457 | = (6,164) | = (41,457) | 0.04 | 6,164 | = (18,498,174) | holy8a32506, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 48.4 s | 2.35 | yes |
| decontx_raw_tenx_chickenamygdala | 8,654 | 7,453 | = (1,201) | = (7,453) | 0.005 | 1,201 | = (6,553,580) | holy8a32604, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 29.0 s | 0.67 | yes |
| decontx_raw_emptydrops_chickenamygdala | 47,621 | 41,962 | = (5,659) | = (41,962) | 0.11 | 5,659 | = (16,930,436) | holy8a32507, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 50.5 s | 2.32 | yes |
| cellbender_fromraw_chickenamygdala | 9,780 | 8,446 | = (1,334) | = (8,446) | 0.22 | 1,334 | = (5,897,924) | holy8a32603, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 29.6 s | 0.69 | yes |
| soupx_tenx_zmaysGSM5552419 | 2,676 | 2,307 | = (369) | = (2,307) | 0.14 | 369 | = (4,141,321) | holy8a32502, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 16.8 s | 0.55 | yes |
| soupx_emptydrops_zmaysGSM5552419 | 3,022 | 2,605 | = (417) | = (2,605) | 0.18 | 417 | = (4,408,726) | holy8a32607, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 18.0 s | 0.57 | yes |
| decontx_raw_tenx_zmaysGSM5552419 | 2,673 | 2,301 | = (372) | = (2,301) | 0.03 | 372 | = (3,980,486) | holy8a32508, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 16.6 s | 0.49 | yes |
| decontx_raw_emptydrops_zmaysGSM5552419 | 3,022 | 2,608 | = (414) | = (2,608) | 0.005 | 414 | = (3,720,706) | holy8a30504, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 17.9 s | 0.49 | yes |
| cellbender_fromraw_zmaysGSM5552419 | 4,685 | 4,035 | = (650) | = (4,035) | 0.06 | 650 | = (4,495,842) | holy8a32607, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 25.6 s | 0.68 | yes |
| soupx_tenx_pvitticeps | 7,880 | 6,787 | = (1,093) | = (6,787) | 0.1 | 1,093 | = (15,263,978) | holy8a32503, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 42.2 s | 1.25 | yes |
| soupx_emptydrops_pvitticeps | 17,476 | 15,047 | = (2,429) | = (15,047) | 0.04 | 2,429 | = (14,584,851) | holy8a32601, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 40.7 s | 1.10 | yes |
| decontx_raw_tenx_pvitticeps | 7,879 | 6,787 | = (1,092) | = (6,787) | 0.09 | 1,092 | = (14,529,118) | holy8a32508, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 36.2 s | 1.19 | yes |
| decontx_raw_emptydrops_pvitticeps | 17,476 | 15,053 | = (2,423) | = (15,053) | 0.04 | 2,423 | = (13,953,622) | holy8a32605, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 39.9 s | 1.09 | yes |
| cellbender_fromraw_pvitticeps | 17,939 | 15,431 | = (2,508) | = (15,431) | 0.03 | 2,508 | = (13,171,450) | holy8a32502, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 41.7 s | 1.10 | yes |
| soupx_tenx_10xmousebrain | 11,831 | 10,185 | = (1,646) | = (10,185) | 0.07 | 1,646 | = (36,923,827) | holy8a32604, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 61.1 s | 2.55 | yes |
| soupx_emptydrops_10xmousebrain | 12,289 | 10,565 | = (1,724) | = (10,565) | 0.15 | 1,724 | = (36,614,504) | holy8a32607, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 65.1 s | 2.47 | yes |
| decontx_raw_tenx_10xmousebrain | 11,799 | 10,159 | = (1,640) | = (10,159) | 0.06 | 1,640 | = (35,393,797) | holy8a32608, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 63.5 s | 2.47 | yes |
| decontx_raw_emptydrops_10xmousebrain | 12,285 | 10,561 | = (1,724) | = (10,561) | 0.03 | 1,724 | = (36,499,172) | holy8a32606, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 57.0 s | 2.40 | yes |
| cellbender_fromraw_10xmousebrain | 12,203 | 10,485 | = (1,718) | = (10,485) | 0.3 | 1,718 | = (36,431,656) | holy8a30605, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 64.6 s | 2.45 | yes |
| soupx_tenx_sc108 | 11,205 | 9,604 | = (1,601) | = (9,604) | 0.05 | 1,601 | = (16,498,489) | holy8a32607, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 48.6 s | 1.44 | yes |
| soupx_emptydrops_sc108 | 1,943 | 1,683 | = (260) | = (1,683) | 0.1 | 260 | = (9,100,919) | holy8a30606, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 19.3 s | 0.78 | yes |
| decontx_raw_tenx_sc108 | 11,204 | 9,612 | = (1,592) | = (9,612) | 0.24 | 1,592 | = (14,828,318) | holy8a32604, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 47.0 s | 1.38 | yes |
| decontx_raw_emptydrops_sc108 | 1,943 | 1,683 | = (260) | = (1,683) | 0.09 | 260 | = (9,027,313) | holy8a30508, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 19.1 s | 0.75 | yes |
| cellbender_fromraw_sc108 | 13,643 | 11,691 | = (1,952) | = (11,691) | 0.15 | 1,952 | = (13,053,876) | holy8a32608, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 45.7 s | 1.13 | yes |
| soupx_tenx_sc109 | 9,386 | 8,057 | = (1,329) | = (8,057) | 0.2 | 1,329 | = (17,876,123) | holy8a32602, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 64.3 s | 1.51 | yes |
| soupx_emptydrops_sc109 | 1,122 | 975 | = (147) | = (975) | 0.07 | 147 | = (7,350,106) | holy8a32502, sapphirerapids | holy8a24109, 8480+, SapphireRapids | 17.2 s | 0.62 | yes |
| decontx_raw_tenx_sc109 | 9,373 | 8,056 | = (1,317) | = (8,056) | 0.13 | 1,317 | = (15,993,358) | holy8a32604, sapphirerapids | holy8a24110, 8480+, SapphireRapids | 46.0 s | 1.42 | yes |
| decontx_raw_emptydrops_sc109 | 1,120 | 975 | = (145) | = (975) | 0.02 | 145 | = (7,189,144) | holy8a32502, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 14.7 s | 0.63 | yes |
| cellbender_fromraw_sc109 | 9,853 | 8,461 | = (1,392) | = (8,461) | 0.02 | 1,392 | = (12,875,076) | holy8a30506, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 46.6 s | 1.35 | yes |
| soupx_tenx_amillepora | 15,558 | 13,363 | = (2,195) | = (13,363) | 0.07 | 2,195 | = (14,892,898) | holy8a30606, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 58.8 s | 1.03 | yes |
| soupx_emptydrops_amillepora | 16,775 | 14,405 | = (2,370) | = (14,405) | 0.08 | 2,370 | = (15,101,057) | holy8a30407, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 43.5 s | 1.07 | yes |
| decontx_raw_tenx_amillepora | 15,558 | 13,367 | = (2,191) | = (13,367) | 0.08 | 2,191 | = (14,797,902) | holy8a32605, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 47.1 s | 1.00 | yes |
| decontx_raw_emptydrops_amillepora | 16,774 | 14,401 | = (2,373) | = (14,401) | 0.09 | 2,373 | = (14,822,734) | holy8a32608, sapphirerapids | holy8a28607, 8480CL, SapphireRapids | 46.4 s | 1.06 | yes |
| cellbender_fromraw_amillepora | 16,218 | 13,918 | = (2,300) | = (13,918) | 0.04 | 2,300 | = (14,569,619) | holy8a32602, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 41.4 s | 1.03 | yes |
| soupx_tenx_amellifera | 29,235 | 25,047 | = (4,188) | = (25,047) | 0.21 | 4,188 | = (67,139,676) | holy8a30502, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 65.8 s | 4.08 | yes |
| soupx_emptydrops_amellifera | 72,702 | 63,004 | = (9,698) | = (63,004) | 0.16 | 9,698 | = (100,345,292) | holy8a30602, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 91.1 s | 6.64 | yes |
| decontx_raw_tenx_amellifera | 29,233 | 25,028 | = (4,205) | = (25,028) | 0.3 | 4,205 | = (59,057,260) | holy8a30607, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 57.8 s | 3.32 | yes |
| decontx_raw_emptydrops_amellifera | 72,676 | 62,891 | = (9,785) | = (62,891) | 0.21 | 9,785 | = (90,239,217) | holy8a30606, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 85.2 s | 6.01 | yes |
| cellbender_fromraw_amellifera | 28,082 | 24,039 | = (4,043) | = (24,039) | 0.08 | 4,043 | = (44,473,879) | holy8a30507, sapphirerapids | holy8a24107, 8480+, SapphireRapids | 55.1 s | 3.01 | yes |
| soupx_tenx_cintestinalis | 14,579 | 12,491 | = (2,088) | = (12,491) | 0.22 | 2,088 | = (14,205,986) | holy8a32605, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 34.6 s | 0.94 | yes |
| soupx_emptydrops_cintestinalis | 18,983 | 16,264 | = (2,719) | = (16,264) | 0.2 | 2,719 | = (13,857,318) | holy8a30504, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 34.9 s | 1.13 | yes |
| decontx_raw_tenx_cintestinalis | 14,577 | 12,487 | = (2,090) | = (12,487) | 0.06 | 2,090 | = (13,348,691) | holy8a32604, sapphirerapids | holy8a24109, 8480+, SapphireRapids | 43.0 s | 0.92 | yes |
| decontx_raw_emptydrops_cintestinalis | 18,975 | 16,252 | = (2,723) | = (16,252) | 0.06 | 2,723 | = (13,744,837) | holy8a32608, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 34.5 s | 1.11 | yes |
| cellbender_fromraw_cintestinalis | 8,982 | 7,696 | = (1,286) | = (7,696) | 0.3 | 1,286 | = (9,825,212) | holy8a32606, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 30.6 s | 0.83 | yes |
| soupx_tenx_eberryi | 34,124 | 29,328 | = (4,796) | = (29,328) | 0.01 | 4,796 | = (27,046,912) | holy8a32606, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 46.6 s | 1.98 | yes |
| soupx_emptydrops_eberryi | 51,873 | 44,410 | = (7,463) | = (44,410) | 0.25 | 7,463 | = (32,147,334) | holy8a30606, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 54.5 s | 2.82 | yes |
| decontx_raw_tenx_eberryi | 34,124 | 29,348 | = (4,776) | = (29,348) | 0.11 | 4,776 | = (20,678,898) | holy8a32605, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 42.2 s | 1.87 | yes |
| decontx_raw_emptydrops_eberryi | 51,870 | 44,468 | = (7,402) | = (44,468) | 0.001 | 7,402 | = (27,439,212) | holy8a30602, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 51.0 s | 2.72 | yes |
| cellbender_fromraw_eberryi | 60,329 | 51,776 | = (8,553) | = (51,776) | 0.05 | 8,553 | = (33,767,365) | holy8a30604, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 58.4 s | 3.18 | yes |
| soupx_tenx_nvectensis | 3,264 | 2,811 | = (453) | = (2,811) | 0.23 | 453 | = (1,612,819) | holy8a24105, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 13.8 s | 0.37 | yes |
| soupx_emptydrops_nvectensis | 4,255 | 3,669 | = (586) | = (3,669) | 0.11 | 586 | = (1,732,341) | holy8a28105, sapphirerapids | holy8a24110, 8480+, SapphireRapids | 16.2 s | 0.47 | yes |
| decontx_raw_tenx_nvectensis | 3,264 | 2,811 | = (453) | = (2,811) | 0.005 | 453 | = (1,545,145) | holy8a28310, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 13.5 s | 0.40 | yes |
| decontx_raw_emptydrops_nvectensis | 4,255 | 3,671 | = (584) | = (3,671) | 0.04 | 584 | = (1,723,617) | holy8a30604, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 16.7 s | 0.49 | yes |
| cellbender_fromraw_nvectensis | 4,421 | 3,821 | = (600) | = (3,821) | 0.12 | 600 | = (1,704,943) | holy8a32603, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 16.3 s | 0.54 | yes |
| soupx_tenx_opatagonica | 15,636 | 13,424 | = (2,212) | = (13,424) | 0.05 | 2,212 | = (9,911,388) | holy8a32508, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 39.2 s | 0.90 | yes |
| soupx_emptydrops_opatagonica | 20,096 | 17,247 | = (2,849) | = (17,247) | 0.001 | 2,849 | = (10,253,709) | holy8a32607, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 40.2 s | 1.10 | yes |
| decontx_raw_tenx_opatagonica | 15,330 | 13,159 | = (2,171) | = (13,159) | 0.04 | 2,171 | = (8,735,103) | holy8a32603, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 37.5 s | 0.87 | yes |
| decontx_raw_emptydrops_opatagonica | 20,092 | 17,234 | = (2,858) | = (17,234) | 0.005 | 2,858 | = (9,668,978) | holy8a32607, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 39.8 s | 1.08 | yes |
| cellbender_fromraw_opatagonica | 20,475 | 17,588 | = (2,887) | = (17,588) | 0.04 | 2,887 | = (10,015,019) | holy8a32606, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 40.4 s | 1.10 | yes |
| soupx_tenx_raegyptiacus | 34,421 | 29,572 | = (4,849) | = (29,572) | 0.1 | 4,849 | = (76,036,391) | holy8a32505, sapphirerapids | holy8a24109, 8480+, SapphireRapids | 97.3 s | 4.91 | yes |
| soupx_emptydrops_raegyptiacus | 36,863 | 31,636 | = (5,227) | = (31,636) | 0.13 | 5,227 | = (77,243,595) | holy8a30502, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 82.3 s | 4.97 | yes |
| decontx_raw_tenx_raegyptiacus | 34,419 | 29,570 | = (4,849) | = (29,570) | 0.18 | 4,849 | = (76,378,251) | holy8a32503, sapphirerapids | holy8a24110, 8480+, SapphireRapids | 78.9 s | 4.77 | yes |
| decontx_raw_emptydrops_raegyptiacus | 36,856 | 31,627 | = (5,229) | = (31,627) | 0.25 | 5,229 | = (70,640,601) | holy8a30602, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 80.3 s | 4.81 | yes |
| cellbender_fromraw_raegyptiacus | 35,768 | 30,728 | = (5,040) | = (30,728) | 0.27 | 5,040 | = (71,320,787) | holy8a30604, sapphirerapids | holy8a24107, 8480+, SapphireRapids | 80.8 s | 4.78 | yes |
| soupx_tenx_slycopersicum | 23,149 | 19,927 | = (3,222) | = (19,927) | 0.04 | 3,222 | = (67,039,705) | holy8a30507, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 83.3 s | 4.01 | yes |
| soupx_emptydrops_slycopersicum | 22,237 | 19,145 | = (3,092) | = (19,145) | 0.24 | 3,092 | = (63,312,423) | holy8a30604, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 75.6 s | 3.86 | yes |
| decontx_raw_tenx_slycopersicum | 23,117 | 19,901 | = (3,216) | = (19,901) | 0.06 | 3,216 | = (58,450,732) | holy8a30601, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 70.4 s | 3.55 | yes |
| decontx_raw_emptydrops_slycopersicum | 22,167 | 19,089 | = (3,078) | = (19,089) | 0.07 | 3,078 | = (54,969,297) | holy8a30605, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 69.1 s | 3.38 | yes |
| cellbender_fromraw_slycopersicum | 20,923 | 18,035 | = (2,888) | = (18,035) | 0.13 | 2,888 | = (58,577,631) | holy8a32604, sapphirerapids | holy8a24106, 8480+, SapphireRapids | 82.5 s | 3.54 | yes |
| soupx_tenx_smansoni | 15,094 | 12,968 | = (2,126) | = (12,968) | 0.02 | 2,126 | = (17,409,162) | holy8a32606, sapphirerapids | holy8a24109, 8480+, SapphireRapids | 42.7 s | 1.13 | yes |
| soupx_emptydrops_smansoni | 12,865 | 11,029 | = (1,836) | = (11,029) | 0.03 | 1,836 | = (14,843,595) | holy8a32505, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 32.5 s | 1.06 | yes |
| decontx_raw_tenx_smansoni | 15,094 | 12,964 | = (2,130) | = (12,964) | 0.03 | 2,130 | = (15,497,786) | holy8a32603, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 32.2 s | 1.03 | yes |
| decontx_raw_emptydrops_smansoni | 12,864 | 11,027 | = (1,837) | = (11,027) | 0.005 | 1,837 | = (13,779,044) | holy8a32505, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 31.0 s | 0.96 | yes |
| cellbender_fromraw_smansoni | 14,886 | 12,783 | = (2,103) | = (12,783) | 0.01 | 2,103 | = (15,005,522) | holy8a32607, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 31.3 s | 1.00 | yes |
| soupx_tenx_smediterranea | 11,318 | 9,765 | = (1,553) | = (9,765) | 0.23 | 1,553 | = (22,457,642) | holy8a32603, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 55.0 s | 1.58 | yes |
| soupx_emptydrops_smediterranea | 10,432 | 8,988 | = (1,444) | = (8,988) | 0.16 | 1,444 | = (20,108,773) | holy8a30402, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 54.0 s | 1.57 | yes |
| decontx_raw_tenx_smediterranea | 11,315 | 9,763 | = (1,552) | = (9,763) | 0.001 | 1,552 | = (20,096,339) | holy8a32607, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 63.0 s | 1.53 | yes |
| decontx_raw_emptydrops_smediterranea | 10,432 | 8,993 | = (1,439) | = (8,993) | 0.27 | 1,439 | = (20,515,756) | holy8a32508, sapphirerapids | holy8a26107, 8480+, SapphireRapids | 53.2 s | 1.52 | yes |
| cellbender_fromraw_smediterranea | 11,656 | 10,054 | = (1,602) | = (10,054) | 0.3 | 1,602 | = (21,354,407) | holy8a30508, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 56.9 s | 1.54 | yes |
| soupx_tenx_spurpuratus | 19,307 | 16,765 | = (2,542) | = (16,765) | 0.05 | 2,542 | = (10,547,820) | holy8a32502, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 34.4 s | 1.05 | yes |
| soupx_emptydrops_spurpuratus | 9,139 | 7,933 | = (1,206) | = (7,933) | 0.02 | 1,206 | = (5,043,688) | holy8a32502, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 29.6 s | 0.64 | yes |
| decontx_raw_tenx_spurpuratus | 19,307 | 16,613 | = (2,694) | = (16,613) | 0.06 | 2,694 | = (6,591,814) | holy8a32503, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 30.4 s | 0.98 | yes |
| decontx_raw_emptydrops_spurpuratus | 9,139 | 7,942 | = (1,197) | = (7,942) | 0.24 | 1,197 | = (3,996,045) | holy8a30410, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 26.4 s | 0.61 | yes |
| cellbender_fromraw_spurpuratus | 9,902 | 8,601 | = (1,301) | = (8,601) | 0.08 | 1,301 | = (3,593,950) | holy8a30412, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 27.7 s | 0.63 | yes |
| soupx_tenx_spistillata | 18,175 | 15,577 | = (2,598) | = (15,577) | 0.12 | 2,598 | = (12,059,059) | holy8a32605, sapphirerapids | holy8a28101, 8480CL, SapphireRapids | 48.8 s | 1.06 | yes |
| soupx_emptydrops_spistillata | 21,064 | 18,055 | = (3,009) | = (18,055) | 0.001 | 3,009 | = (12,176,927) | holy8a32503, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 46.0 s | 1.19 | yes |
| decontx_raw_tenx_spistillata | 18,175 | 15,575 | = (2,600) | = (15,575) | 0.03 | 2,600 | = (11,390,846) | holy8a30504, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 43.9 s | 1.04 | yes |
| decontx_raw_emptydrops_spistillata | 21,064 | 18,050 | = (3,014) | = (18,050) | 0.001 | 3,014 | = (11,834,392) | holy8a30603, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 44.1 s | 1.17 | yes |
| cellbender_fromraw_spistillata | 20,880 | 17,890 | = (2,990) | = (17,890) | 0.02 | 2,990 | = (11,604,229) | holy8a30604, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 57.2 s | 1.17 | yes |
| soupx_tenx_sscovelli | 19,107 | 16,389 | = (2,718) | = (16,389) | 0.07 | 2,718 | = (31,343,110) | holy8a32503, sapphirerapids | holy8a28607, 8480CL, SapphireRapids | 63.2 s | 2.22 | yes |
| soupx_emptydrops_sscovelli | 79,409 | 69,500 | = (9,909) | = (69,500) | 0.19 | 9,909 | = (62,191,630) | holy8a24105, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 82.8 s | 4.43 | yes |
| decontx_raw_tenx_sscovelli | 18,446 | 15,832 | = (2,614) | = (15,832) | 0.04 | 2,614 | = (26,961,757) | holy8a30607, sapphirerapids | holy8a30102, 8480CL, SapphireRapids | 54.3 s | 1.97 | yes |
| decontx_raw_emptydrops_sscovelli | 78,687 | 69,639 | = (9,048) | = (69,639) | 0.17 | 9,048 | = (58,309,011) | holy8a24209, sapphirerapids | holy8a24108, 8480+, SapphireRapids | 87.8 s | 4.32 | yes |
| cellbender_fromraw_sscovelli | 19,342 | 16,583 | = (2,759) | = (16,583) | 0.18 | 2,759 | = (31,637,558) | holy8a32605, sapphirerapids | holy8a30103, 8480CL, SapphireRapids | 53.9 s | 1.93 | yes |
| soupx_tenx_drerio | 8,322 | 7,164 | = (1,158) | = (7,164) | 0.06 | 1,158 | = (5,585,390) | holy8a30601, sapphirerapids | holy8a30103, 8480CL, SapphireRapids | 32.0 s | 0.65 | yes |
| soupx_emptydrops_drerio | 32,676 | 28,141 | = (4,535) | = (28,141) | 0.07 | 4,535 | = (8,391,001) | holy8a32608, sapphirerapids | holy8a30103, 8480CL, SapphireRapids | 40.8 s | 1.57 | yes |
| decontx_raw_tenx_drerio | 8,322 | 7,158 | = (1,164) | = (7,158) | 0.23 | 1,164 | = (4,916,680) | holy8a32604, sapphirerapids | holy8a30106, 8480CL, SapphireRapids | 33.8 s | 0.63 | yes |
| decontx_raw_emptydrops_drerio | 32,676 | 28,030 | = (4,646) | = (28,030) | 0.3 | 4,646 | = (9,030,615) | holy8a32606, sapphirerapids | holy8a24110, 8480+, SapphireRapids | 38.1 s | 1.58 | yes |
| cellbender_fromraw_drerio | 11,356 | 9,734 | = (1,622) | = (9,734) | 0.22 | 1,622 | = (3,991,732) | holy8a30606, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 32.7 s | 0.68 | yes |
| soupx_tenx_aaegypti | 11,317 | 9,800 | = (1,517) | = (9,800) | 0.04 | 1,517 | = (3,870,107) | holy8a32508, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 23.5 s | 0.65 | yes |
| soupx_emptydrops_aaegypti | 19,483 | 16,838 | = (2,645) | = (16,838) | 0.12 | 2,645 | = (5,705,718) | holy8a32502, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 25.1 s | 0.96 | yes |
| decontx_raw_tenx_aaegypti | 11,309 | 9,760 | = (1,549) | = (9,760) | 0.005 | 1,549 | = (3,983,059) | holy8a30607, sapphirerapids | holy8a26107, 8480+, SapphireRapids | 23.8 s | 0.64 | yes |
| decontx_raw_emptydrops_aaegypti | 19,483 | 16,818 | = (2,665) | = (16,818) | 0.14 | 2,665 | = (5,209,620) | holy8a32605, sapphirerapids | holy8a26107, 8480+, SapphireRapids | 24.7 s | 0.95 | yes |
| cellbender_fromraw_aaegypti | 7,717 | 6,654 | = (1,063) | = (6,654) | 0.005 | 1,063 | = (2,929,253) | holy8a32603, sapphirerapids | holy8a26107, 8480+, SapphireRapids | 19.9 s | 0.55 | yes |
| soupx_tenx_ghirsutum1 | 6,821 | 5,889 | = (932) | = (5,889) | 0.01 | 932 | = (11,367,236) | holy8a32602, sapphirerapids | holy8a26107, 8480+, SapphireRapids | 56.5 s | 1.13 | yes |
| soupx_emptydrops_ghirsutum1 | 6,416 | 5,541 | = (875) | = (5,541) | 0.2 | 875 | = (12,898,194) | holy8a30506, sapphirerapids | holy8a30106, 8480CL, SapphireRapids | 52.5 s | 1.10 | yes |
| decontx_raw_tenx_ghirsutum1 | 6,798 | 5,854 | = (944) | = (5,854) | 0.28 | 944 | = (8,220,276) | holy8a32605, sapphirerapids | holy8a30106, 8480CL, SapphireRapids | 50.6 s | 0.96 | yes |
| decontx_raw_emptydrops_ghirsutum1 | 6,398 | 5,514 | = (884) | = (5,514) | 0.02 | 884 | = (7,302,120) | holy8a30504, sapphirerapids | holy8a30106, 8480CL, SapphireRapids | 50.9 s | 0.95 | yes |
| cellbender_fromraw_ghirsutum1 | 6,007 | 5,192 | = (815) | = (5,192) | 0.17 | 815 | = (10,025,692) | holy8a30104, sapphirerapids | holy8a30107, 8480CL, SapphireRapids | 43.7 s | 0.84 | yes |
| soupx_tenx_ghirsutum2 | 6,777 | 5,838 | = (939) | = (5,838) | 0.19 | 939 | = (17,586,640) | holy8a32505, sapphirerapids | holy8a24111, 8480+, SapphireRapids | 52.2 s | 1.50 | yes |
| soupx_emptydrops_ghirsutum2 | 6,895 | 5,933 | = (962) | = (5,933) | 0.03 | 962 | = (14,275,069) | holy8a32502, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 63.2 s | 1.50 | yes |
| decontx_raw_tenx_ghirsutum2 | 6,772 | 5,825 | = (947) | = (5,825) | 0.02 | 947 | = (9,078,025) | holy8a32605, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 56.8 s | 1.09 | yes |
| decontx_raw_emptydrops_ghirsutum2 | 6,890 | 5,925 | = (965) | = (5,925) | 0.03 | 965 | = (9,440,786) | holy8a32607, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 56.1 s | 1.11 | yes |
| cellbender_fromraw_ghirsutum2 | 59,323 | 52,554 | = (6,769) | = (52,554) | 0.03 | 6,769 | = (67,180,706) | holy8a32603, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 2:21.7 | 4.04 | yes |
| soupx_tenx_athaliana | 5,635 | 4,831 | = (804) | = (4,831) | 0.04 | 804 | = (15,340,052) | holy8a32606, sapphirerapids | holy8a30103, 8480CL, SapphireRapids | 37.3 s | 1.38 | yes |
| soupx_emptydrops_athaliana | 1,539 | 1,327 | = (212) | = (1,327) | 0.19 | 212 | = (9,707,201) | holy8a30606, sapphirerapids | holy8a24112, 8480+, SapphireRapids | 17.9 s | 0.77 | yes |
| decontx_raw_tenx_athaliana | 5,598 | 4,797 | = (801) | = (4,797) | 0.005 | 801 | = (14,116,554) | holy8a32508, sapphirerapids | holy8a24206, 8480+, SapphireRapids | 38.9 s | 1.28 | yes |
| decontx_raw_emptydrops_athaliana | 1,539 | 1,328 | = (211) | = (1,328) | 0.19 | 211 | = (9,326,089) | holy8a30608, sapphirerapids | holy8a24207, 8480+, SapphireRapids | 17.3 s | 0.77 | yes |
| cellbender_fromraw_athaliana | 4,314 | 3,699 | = (615) | = (3,699) | 0.03 | 615 | = (14,511,454) | holy8a30606, sapphirerapids | holy8a28103, 8480CL, SapphireRapids | 35.3 s | 1.24 | yes |

## Findings

- **No difference anywhere.** 150 of 150 objects, 1,350 of 1,350 checks, 3.3 billion values. What M6 and M6b established on ten objects holds on all 150: the five arms and the 30 samples, 1,120 to 79,409 cells, 10,172 to 99,452 features, k from 15 to 20,349, and 32 different pK values chosen (every value of the grid but 5e-04).
- **The saved counts are the input's restricted to the kept cells, on all 150.** As on M6b's ten objects, R's `subset()` keeps every stored entry as it is. That includes the explicit zeros: 120 million among the singlets, all in the 60 SoupX objects (decontX and CellBender store none). So the input without the doublets, every entry as stored, which is what `singlets/` writes, is exactly R's result on the whole panel.
- **Memory follows the stored entries, not the cells.** The panel's peak is soupx_emptydrops_amellifera: 6.64 GiB on 8 threads, with 72,702 cells. That is above sscovelli's 4.43 GiB with 79,409 cells, because amellifera stores 126 million entries (1,733 per cell over 12,398 features) against sscovelli's 66 million. Both peak in the final run's SCTransform, as M6 described. The workflow's job for amellifera peaked at 117.3 GiB. sscovelli and ghirsutum2 peak where M6's benchmark put them (4.43 and 4.04 GiB here, 4.41 to 4.51 and 4.02 to 4.08 GiB there).
- **Time.** The 150 Rust runs took 1:53:49 of wall time added together, on 8 threads each (median 43.6 s, longest 2:21.7 for cellbender_fromraw_ghirsutum2), and 8:18:41 of CPU time. The workflow's 150 jobs took 32:42:41 added together, on one CPU each. They include the Seurat processing that follows DoubletFinder in `doubletfinder.R`, so this is context, not a benchmark.
- **The 10x writer's test failed in debug builds.** On main (d95f7f7), `cargo test --workspace` failed one test, deterministically: `scio::tenx::tests::same_bytes_for_any_thread_count` panicked with "attempt to subtract with overflow" at `crates/scio/src/tenx.rs:230`. M6b's report records the suite passing; on main, this test fails in the default (debug) profile and passes with `--release`.
  - **The cause.** `format_piece` computed the input index of an output entry as `m.p[j] - starts[c] + lo`. When the columns are written out of input order, as that test writes them (columns 4, 0, 2, 3), `m.p[j]` can be smaller than `starts[c]`, and the unsigned subtraction underflows. Debug builds trap it.
  - **No output was affected.** Release builds wrap around and back, so the index was right. And `scrust doubletfinder` writes the singlets in input order, where the difference never goes below zero.
  - **The fix.** Commit b8c8a8c computes `m.p[j] + (lo - starts[c])`. All 160 tests pass, without warnings. Before the fix, `cargo test` stopped at the failing test binary after 141 tests.
  - **The binary.** The fix changes the binary (sha256 `808f6300` against `dfd3b9c6`) but not its outputs: both binaries ran the whole panel, and their outputs are identical, 900 of 900 files and 150 of 150 `result` sections.

## Plan corrections

- **B4 and B6 (this milestone's brief):** the full panel now includes the singlets' counts, with their one-time export. B5's singlets row covers the full panel.
- **B4's range holds:** 1,120 cells (decontx_raw_emptydrops_sc109) to 79,409 (soupx_emptydrops_sscovelli).
- **A4's CPU types hold** for the jobs that wrote the saved outputs: all 150 on Sapphire Rapids, 38 nodes. A4 and B4 now point to the panel table, which records them per object.

## Deviations

- **8 threads per run**, where the workflow's jobs had one CPU. The results do not depend on the thread count (M6, M6b).
- **The production profile only.** The saved outputs were made under the kernel OpenBLAS detects on Sapphire Rapids; the portable profile has no reference for these objects.
- **Two Rust runs per object.** The label of record is `m7-b8c8a8c`, the branch's code. `m7-d95f7f7`, main's binary, ran first and stays as a record, identical. So does the trial `m7test-d95f7f7` (job 60159281: soupx_tenx_SRR18256555, decontx_raw_tenx_SRR18256555 and decontx_raw_tenx_nvectensis), also identical to the label of record, file for file.
- **The P objects were compared again**, under M7's labels and without stage dumps. Their M6 and M6b results stand.

## Reproducing

From the repository's top directory, after `source env.sh`. The panel table, from the pipeline's Snakemake logs (paths in the pipeline repo), sacct and scontrol (no `.rds`; a login node is enough):

```bash
bash validation/slurm/doubletfinder_panel.sh .snakemake/log/2026-11-08T214133.605174.snakemake.log \
    .snakemake/log/2026-11-09T220250.045638.snakemake.log > validation/doubletfinder_panel.tsv
```

The exports, then the Rust runs, each object's run once its own export task has succeeded. Existing exports are kept, so on this checkout the export tasks only check:

```bash
cargo build --release -p scrust
c=$(git rev-parse --short HEAD)
cp target/release/scrust data/doubletfinder/logs/bin/scrust-m7-$c
mkdir -p data/doubletfinder/logs/m7
E=$(sbatch --parsable --array=1-150%50 --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_export.sbatch)
SCRUST_BIN=$PWD/data/doubletfinder/logs/bin/scrust-m7-$c sbatch --export=ALL --array=1-150%50 --dependency=aftercorr:$E \
    --output=data/doubletfinder/logs/m7/%x-%A_%a.sbatch.log validation/slurm/m7_panel.sbatch m7-$c
```

The report's tables and the TSV come from [`validation/R/report_doubletfinder_m7.R`](../../validation/R/report_doubletfinder_m7.R), with the production R env and the hermetic flags. It reads only plain files, but `--same-as` hashes about 20 GB of outputs, so it ran as a batch job (10 minutes):

```bash
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 "$R_ENV/bin/Rscript" validation/R/report_doubletfinder_m7.R \
    m7-b8c8a8c --tsv docs/validation/doubletfinder_m7.tsv --same-as m7-d95f7f7 > tables.md
```

The manifest: `validation/R/manifest_doubletfinder.R`, in a batch job (13 minutes). The SLURM logs of every M7 job, and the batch scripts of the report and the manifest, are in `data/doubletfinder/logs/m7/`.
