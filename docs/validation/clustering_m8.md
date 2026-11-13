# Clustering reference data (milestone M8)

M8 builds the R reference data that Step 2's clustering crate (milestones M9–M17) is checked against (docs/PLAN.md C4 and C6). It ports nothing to Rust. This report gives the checks behind each "done when" criterion. Recipes: [`validation/R/README.md`](../../validation/R/README.md) and [`validation/slurm/README.md`](../../validation/slurm/README.md); every file is listed in [`validation/manifest.toml`](../../validation/manifest.toml).

## Outcome

All of plan C6's M8 criteria are met, with the first read per profile (below).

| Criterion | Result |
|---|---|
| Every direct run reproduces its saved labels and levels, silhouette widths, purities and SCT layers | **Yes, in the production profile, the saved objects' own:** all 104 production-profile runs of the 34 objects (27 T, 7 P; direct, instrumented, without RunUMAP, and through the v5 assay path) reproduce their saved object in every file: labels and levels, silhouette widths and purities bit for bit, `nCount_SCT` and `nFeature_SCT`, the SCT layers and model, the PCA and both graphs. The 118 portable-profile runs (Haswell kernel) keep every label, level and graph, the corrected counts everywhere and, on tier T, every SCT layer; the PCA's last bits move, and with them the silhouette widths' last bits and, in 28 runs, a few purities (R-vs-R floor). |
| The run without UMAP gives identical graphs and labels | **Yes, 34 of 34 objects:** every result file identical, not only the graphs and labels; only the digests of `.Random.seed` after RunUMAP differ. |
| The marker reference of every T and P object settled | **Yes:** the workflow's marker steps on each saved `.rds` reproduce its saved `<prefix>_markergenes.csv` byte for byte, 34 of 34, once the clusters are combined in the saved table's order (the workflow's is the filesystem's, not the level order). Every reference is the saved table. |
| The instrumented optimizer reproduces R's labels | **Yes, 34 of 34:** on each saved SNN graph its partition equals R's optimizer result (before GroupSingletons), and the saved labels themselves for the 30 objects without singletons; its printed nodes, edges, maximum modularity and communities equal R's and, for the 7 production objects, the workflow job's log. It also gives every start's modularity as a double. |
| The floor measured | **Yes:** the same kernel on two CPU types changes nothing (27 T objects and the small P object, every file, stage array and RNG checkpoint). The Haswell and SapphireRapids kernels on one CPU type (all 34 objects) never change a label, a neighbour id or a graph entry; the PCA's first 30 dimensions move by at most 1.5e-10, which moves 121 float item values (116 by one float ulp, at most 9.5e-7), 10 cells' neighbour distances, every silhouette width in the last bits, 19 purities (by up to 0.019), and on tier P scale.data in the last bits. |
| R's memory and time for the chain measured | **Yes,** on the 7 P objects with and without RunUMAP: the chain takes 0.4 to 14.6 minutes and peaks at 2.5 to 104.1 GiB, in SCTransform (RunPCA on the two largest objects); RunUMAP adds 0.2 to 2.0 minutes and no memory. |
| Fixtures in the manifest | **Yes:** every export, run, check, K0 vector and committed fixture of M8 is listed in `validation/manifest.toml` with its sha256: 5,629 file entries, which record 7,847 byte-identical copies in other runs, and 359 run tables. |

**How the first criterion was read.** A saved object's silhouette widths come from its PCA, and the PCA's last bits depend on OpenBLAS's kernel (A5b), so only the saved objects' own kernel (SapphireRapids, detected on Sapphire Rapids) can reproduce them bit for bit. Plan K8 says the same ("identical under the reference's kernel"). The criterion was therefore checked on every production-profile run, which all pass in every file; the portable-profile runs, which are the references for portable-profile Rust runs, are compared with the saved objects as part of the floor. This reading was put to the coordinating session as a question during M8; no answer had arrived when this report was written.

## What was run

**Reference.** The workflow's clustering chain as its scripts call it (plan C1), seed 12345, in the production R env (R 4.5.3, Seurat 5.5.1, SeuratObject 5.4.0, sctransform 0.4.3, glmGamPoi 1.22.0, irlba 2.4.1, uwot 0.2.5, RcppAnnoy 0.0.23, bluster 1.20.0, BiocNeighbors 2.4.0; OpenBLAS 0.3.34, one thread): `SCTransform(seed.use = 12345)`, the workflow's PCA guard, `RunPCA(seed.use = 12345)`, `RunUMAP(dims = 1:30, seed.use = 12345)`, `FindNeighbors(dims = 1:30)`, `FindClusters(random.seed = 12345)`, the workflow's `add_silhouette_to_metadata()` and `write_cluster_metadata()`. The calls live in one place, [`validation/R/lib/clustering.R`](../../validation/R/lib/clustering.R), copied from the pipeline's scripts at commit `ea053fa`; those scripts last changed on 2026-10-30 (`8d34de1`), before both the test-reference run and the production run.

**Inputs.** The chain's input is the object's RNA `counts` layer as saved (C1). `export_clustering.R` exported it from each saved object, or linked an earlier export of the same counts and checked them equal slot for slot (M2's `input/` of the 5 + 5 DoubletFinder inputs, M6b's `saved/singlets/` of the five tier-T DoubletFinder outputs). The reference runs rebuild the identical `dgCMatrix` and wrap it with `CreateSeuratObject()`: a v5 `Assay5`, or, for the emptyDrops objects, a v3 `Assay`, as `emptydrops.R`'s `as.Seurat()` made it. Seurat passes `meta.data` to `vst()` as `cell_attr`, but vst uses only columns named `umi`, `gene`, `log_umi`, `log_gene`, `umi_per_gene` or `log_umi_per_gene`; the exports check that no saved object has one, so the counts are the whole input. Every later step reads only what the chain made, and R's stream is reseeded inside every step that draws from it (C2), so the scripts' work before the chain (emptyDrops, SoupX, DoubletFinder, `subset()`, ...) does not matter.

**Objects.** Tier T: all 27 objects of the workflow's test references (`tests/reference_outputs/`). Tier P: the 7 objects of plan C4. Paths, roles and dump levels: [`validation/clustering_datasets.tsv`](../../validation/clustering_datasets.tsv).

| object | tier | role | assay | features | cells | entries | explicit_zeros | sct_features | variable_features | clusters | singletons | counts |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| filtered_seurat_tenx_test | T | base_tenx | v5 | 2,000 | 300 | 64,934 | 0 | 821 | 821 | 3 | 0 | exported |
| filtered_seurat_emptydrops_test | T | base_emptydrops | v3 | 2,000 | 267 | 62,981 | 0 | 821 | 821 | 4 | 0 | exported |
| soupx_tenx_test | T | doubletfinder_input | v5 | 2,000 | 300 | 64,934 | 2,686 | 815 | 815 | 3 | 0 | reused (soupx_tenx_test/input) |
| soupx_emptydrops_test | T | doubletfinder_input | v5 | 2,000 | 267 | 62,981 | 809 | 818 | 818 | 4 | 0 | reused (soupx_emptydrops_test/input) |
| decontx_raw_tenx_test | T | doubletfinder_input | v5 | 2,000 | 300 | 57,708 | 0 | 782 | 782 | 3 | 0 | reused (decontx_raw_tenx_test/input) |
| decontx_raw_emptydrops_test | T | doubletfinder_input | v5 | 2,000 | 267 | 58,001 | 0 | 794 | 794 | 4 | 0 | reused (decontx_raw_emptydrops_test/input) |
| cellbender_fromraw_test | T | doubletfinder_input | v5 | 2,000 | 278 | 63,782 | 0 | 821 | 821 | 4 | 0 | reused (cellbender_fromraw_test/input) |
| doubletfinder_soupx_tenx_test | T | doubletfinder_output | v5 | 2,000 | 276 | 59,103 | 2,453 | 802 | 802 | 3 | 0 | reused (soupx_tenx_test/saved/singlets) |
| doubletfinder_soupx_emptydrops_test | T | doubletfinder_output | v5 | 2,000 | 237 | 54,893 | 684 | 799 | 799 | 3 | 0 | reused (soupx_emptydrops_test/saved/singlets) |
| doubletfinder_decontx_raw_tenx_test | T | doubletfinder_output | v5 | 2,000 | 273 | 52,332 | 0 | 776 | 776 | 3 | 0 | reused (decontx_raw_tenx_test/saved/singlets) |
| doubletfinder_decontx_raw_emptydrops_test | T | doubletfinder_output | v5 | 2,000 | 237 | 50,896 | 0 | 786 | 786 | 4 | 0 | reused (decontx_raw_emptydrops_test/saved/singlets) |
| doubletfinder_cellbender_fromraw_test | T | doubletfinder_output | v5 | 2,000 | 247 | 54,936 | 0 | 796 | 796 | 3 | 0 | reused (cellbender_fromraw_test/saved/singlets) |
| scdblfinder_soupx_tenx_test | T | scdblfinder_output | v5 | 2,000 | 292 | 63,286 | 2,615 | 812 | 812 | 3 | 0 | exported |
| scdblfinder_soupx_emptydrops_test | T | scdblfinder_output | v5 | 2,000 | 264 | 62,278 | 800 | 818 | 818 | 4 | 0 | exported |
| scdblfinder_decontx_raw_tenx_test | T | scdblfinder_output | v5 | 2,000 | 296 | 57,138 | 0 | 781 | 781 | 3 | 0 | exported |
| scdblfinder_decontx_raw_emptydrops_test | T | scdblfinder_output | v5 | 2,000 | 262 | 57,086 | 0 | 791 | 791 | 4 | 0 | exported |
| scdblfinder_cellbender_fromraw_test | T | scdblfinder_output | v5 | 2,000 | 278 | 63,782 | 0 | 821 | 821 | 4 | 0 | exported |
| posthocfilt_threshold_doubletfinder_soupx_tenx_test | T | posthoc | v5 | 2,000 | 142 | 39,797 | 1,854 | 768 | 768 | 2 | 0 | exported |
| posthocfilt_threshold_doubletfinder_soupx_emptydrops_test | T | posthoc | v5 | 2,000 | 137 | 38,303 | 509 | 766 | 766 | 2 | 0 | exported |
| posthocfilt_threshold_doubletfinder_decontx_raw_tenx_test | T | posthoc | v5 | 2,000 | 80 | 22,821 | 0 | 694 | 694 | 2 | 0 | exported |
| posthocfilt_threshold_doubletfinder_decontx_raw_emptydrops_test | T | posthoc | v5 | 2,000 | 92 | 26,072 | 0 | 718 | 718 | 3 | 0 | exported |
| posthocfilt_threshold_doubletfinder_cellbender_fromraw_test | T | posthoc | v5 | 2,000 | 135 | 37,370 | 0 | 756 | 756 | 2 | 0 | exported |
| posthocfilt_threshold_scdblfinder_soupx_tenx_test | T | posthoc | v5 | 2,000 | 154 | 43,566 | 2,000 | 776 | 776 | 2 | 0 | exported |
| posthocfilt_threshold_scdblfinder_soupx_emptydrops_test | T | posthoc | v5 | 2,000 | 161 | 45,107 | 622 | 787 | 787 | 3 | 0 | exported |
| posthocfilt_threshold_scdblfinder_decontx_raw_tenx_test | T | posthoc | v5 | 2,000 | 87 | 24,617 | 0 | 702 | 702 | 2 | 0 | exported |
| posthocfilt_threshold_scdblfinder_decontx_raw_emptydrops_test | T | posthoc | v5 | 2,000 | 112 | 31,462 | 0 | 738 | 738 | 3 | 0 | exported |
| posthocfilt_threshold_scdblfinder_cellbender_fromraw_test | T | posthoc | v5 | 2,000 | 163 | 45,591 | 0 | 793 | 793 | 2 | 0 | exported |
| soupx_tenx_SRR18256555 | P | small | v5 | 32,285 | 3,108 | 2,817,955 | 74,901 | 15,256 | 3,000 | 15 | 0 | reused (soupx_tenx_SRR18256555/input) |
| soupx_tenx_10xmousebrain | P | medium | v5 | 32,285 | 11,831 | 42,440,187 | 445,249 | 20,926 | 3,000 | 21 | 0 | reused (soupx_tenx_10xmousebrain/input) |
| filtered_seurat_tenx_10xmousebrain | P | base_tenx | v5 | 32,285 | 11,831 | 42,440,187 | 0 | 20,948 | 3,000 | 21 | 0 | exported |
| soupx_emptydrops_4500871lane | P | medium2 | v5 | 32,285 | 26,598 | 43,674,465 | 1,467,924 | 22,831 | 3,000 | 31 | 1 | reused (soupx_emptydrops_4500871lane/input) |
| soupx_emptydrops_sscovelli | P | large | v5 | 28,754 | 79,409 | 66,036,718 | 3,602,702 | 22,014 | 3,000 | 36 | 18 | reused (soupx_emptydrops_sscovelli/input) |
| filtered_seurat_emptydrops_sscovelli | P | base_emptydrops | v3 | 28,754 | 79,409 | 66,036,718 | 0 | 22,025 | 3,000 | 37 | 22 | exported |
| cellbender_fromraw_ghirsutum2 | P | many_singletons | v5 | 99,452 | 59,323 | 73,602,177 | 0 | 56,106 | 3,000 | 22 | 40 | reused (cellbender_fromraw_ghirsutum2/input) |

entries: stored entries of the RNA counts; explicit zeros: stored entries equal to 0 (the SoupX objects have them); SCT features: the genes SCTransform keeps (detected in at least 5 cells), all of them variable features on tier T; singletons: GroupSingletons' count in the direct runs (C4's tier-P table lists the same).

**Profiles (A5b).** *Portable*: `OPENBLAS_CORETYPE=Haswell`. *Production*: `OPENBLAS_CORETYPE` unset on the CPU type of the object's workflow job. Every job that made a saved object, the 840 production jobs and the 27 test-reference jobs, ran on Sapphire Rapids, so every production-profile run ran on Sapphire Rapids, where OpenBLAS detected `SapphireRapids` (checked with `scrust fingerprint` before each run). All of M8's runs on Sapphire Rapids used Intel Xeon Platinum 8480CL nodes (`test` partition), the CPU model of M2's runs on the `sapphire` partition, where the production jobs ran.

**Runs.** Run directories are `data/clustering/<object>/runs/<mode>-<profile>-<cpu arch>[-noumap][-v5]/`.

| run | P | T |
|---|---|---|
| direct-portable-cascadelake | 1 | 27 |
| direct-portable-sapphirerapids | 1 | 27 |
| direct-production-sapphirerapids | 7 | 27 |
| direct-production-sapphirerapids-noumap | 7 | 27 |
| instrumented-portable-cascadelake | 1 | 27 |
| instrumented-portable-sapphirerapids | 7 | 27 |
| instrumented-production-sapphirerapids | 7 | 27 |
| instrumented-production-sapphirerapids-v5 | 1 | 1 |

- **Tier T,** all 27 objects: direct and instrumented runs in both profiles on Sapphire Rapids, the portable ones also on Cascade Lake (`shared`, Xeon Platinum 8268), and a production direct run without RunUMAP; the emptyDrops object also through the v5 path.
- **Tier P,** all 7 objects: direct production runs with and without RunUMAP (the cross-check and the memory and time baseline), instrumented runs in both profiles; for the small object also direct portable runs and both portable runs on Cascade Lake; the emptyDrops object also through the v5 path. The portable instrumented runs stand in for portable direct runs, as in M2: instrumented and direct runs agree file for file in all 63 pairs (next sections).
- **Dump levels:** full on tier T (every array, the clipped residual matrix, purity's neighbour lists), standard on the small and medium P objects (scale.data; the small object also purity's lists), light on the three large ones (no residual-sized matrix; digests of the log1p data and scale.data).

**Queues.** During M8 the `sapphire` partition was held up by another run under the same account (more than 100 jobs pending). Every Sapphire Rapids job of M8 therefore ran on the `test` partition (Sapphire Rapids nodes, Xeon Platinum 8480CL, at most 5 jobs per user), several configurations or objects per job; jobs that need no particular CPU went to `shared` and `serial_requeue`.

## The full panel's table

[`validation/clustering_panel.tsv`](../../validation/clustering_panel.tsv) lists all 840 objects the production run clustered (30 samples × 28 objects), written by [`validation/R/clustering_panel.R`](../../validation/R/clustering_panel.R) from the two Snakemake logs of rerun #2 (2026-11-08T214133 and 2026-11-09T220250), sacct, scontrol and the jobs' own logs. For each object: its kind, rule and assay class, the job's data inputs, the saved object, its log and its combined marker table; the job that wrote it (id, Snakemake run, node, CPU type, partition, start, end, elapsed time, MaxRSS of step `.0`); and the numbers the job's log printed (nodes, edges, maximum modularity of the 10 starts, communities, singletons, final clusters).

- **Jobs.** 843 clustering jobs were submitted for the 840 objects. Three objects (`soupx_emptydrops_sscovelli`, `cellbender_fromraw_ghirsutum2`, `filtered_seurat_cellbendercells_ghirsutum2`) had a failed job in the first run and were written by the relaunch; the table lists the failed jobs as superseded. Snakemake touches each output when it sees its job end, so a saved object's modification time comes 1 second to 5.3 minutes after its job's end (median 32 s), never inside the job: the table takes each object's last completed job and checks that.
- **CPU types.** All 840 jobs ran on Sapphire Rapids nodes (`sapphirerapids`).
- **Singletons.** 36 logs report singletons, up to 40 (`cellbender_fromraw_ghirsutum2`), as plan C1 says. Every log holds exactly one optimizer block.

## Where the test references were made

The tier-T references were made by the test-workflow run that started at 14:36 on 2026-11-06 (M2). Its 27 clustering jobs (the 2 base objects, the 5 DoubletFinder inputs, the 5 DoubletFinder and 5 scDblFinder outputs, the 10 post-hoc objects) ran on `sapphire` nodes (holy8a32608, 32607, 32606, 32601, 32505, 32503, 30606, 30604; all `sapphirerapids`), 45 s to 1:53 min and at most 1.1 GiB each (sacct). Their job logs were overwritten by the second test run at 15:01, so for tier T the optimizer's numbers come from the reference runs, which reproduce the saved objects. The test run's conda env differs from the production env only in xgboost (M2), which the chain does not use.

## Reference runs against the saved objects

Each run against its saved object, file by file (`data/clustering/compare/<object>/<run>_vs_saved.json`, by [`compare_clustering.R`](../../validation/R/compare_clustering.R)): the labels and both columns' levels, the silhouette widths and purities as written (17 significant digits, so bitwise), `nCount_SCT` and `nFeature_SCT`; the SCT assay's corrected counts, their log1p (in full at level full, otherwise as a digest of the doubles), scale.data (in full at levels full and standard, otherwise as a digest), the variable features and the SCTModel's slots; the PCA's embeddings, loadings and standard deviations; both graphs, every stored entry and value.

Objects whose run reproduces the saved object, per check:

| tier | run | objects | every file identical | labels and levels | silhouette | purity | SCT counts, data, scale.data | SCT model | PCA | nn, snn |
|---|---|---|---|---|---|---|---|---|---|---|
| T | direct-portable-cascadelake | 27 | 0 | 27 | 0 | 21 | 27 | 27 | 0 | 27 |
| T | direct-portable-sapphirerapids | 27 | 0 | 27 | 0 | 21 | 27 | 27 | 0 | 27 |
| T | direct-production-sapphirerapids | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 |
| T | direct-production-sapphirerapids-noumap | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 |
| T | instrumented-portable-cascadelake | 27 | 0 | 27 | 0 | 21 | 27 | 27 | 0 | 27 |
| T | instrumented-portable-sapphirerapids | 27 | 0 | 27 | 0 | 21 | 27 | 27 | 0 | 27 |
| T | instrumented-production-sapphirerapids | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 | 27 |
| T | instrumented-production-sapphirerapids-v5 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| P | direct-portable-cascadelake | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| P | direct-portable-sapphirerapids | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| P | direct-production-sapphirerapids | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 |
| P | direct-production-sapphirerapids-noumap | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 |
| P | instrumented-portable-cascadelake | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| P | instrumented-portable-sapphirerapids | 7 | 0 | 7 | 0 | 3 | 0 | 0 | 0 | 7 |
| P | instrumented-production-sapphirerapids | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 | 7 |
| P | instrumented-production-sapphirerapids-v5 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

Every one of the 104 production-profile runs reproduces its saved object in every file.

Portable-profile runs (118): labels and levels identical in 118, SCT layers and model in 108, both graphs in 118, purities in 90; silhouette widths differ in the last bits in 118, the PCA in 118 (see the floor).

On tier P, the portable runs' SCT corrected counts are identical to the saved ones in every object; what differs is scale.data and the model's residual statistics, in the last bits, as M2 saw on the larger DoubletFinder objects (S2d). Per-comparison details: [`clustering_m8.tsv`](clustering_m8.tsv) (one row per comparison) and `data/clustering/compare/`.

## Runs against runs

| comparison | tier | pairs | every result file identical | summary result identical | stage arrays identical | RNG checkpoints identical | Annoy index identical |
|---|---|---|---|---|---|---|---|
| instr_vs_direct | T | 54 | 54 | 54 |  |  |  |
| instr_vs_direct | P | 9 | 9 | 9 |  |  |  |
| noumap_vs_umap | T | 27 | 27 | 0 |  |  |  |
| noumap_vs_umap | P | 7 | 7 | 0 |  |  |  |
| v5_vs_v3 | T | 1 | 1 | 0 | 43 / 43 | 1 | 1 |
| v5_vs_v3 | P | 1 | 1 | 0 | 41 / 41 | 1 | 1 |
| cascadelake_vs_sapphirerapids | T | 54 | 54 | 54 | 1161 / 1161 | 27 | 27 |
| cascadelake_vs_sapphirerapids | P | 2 | 2 | 2 | 44 / 44 | 1 | 1 |

- **Instrumented against direct** (63 pairs, both profiles, every object with both): every result file and the `result` section of `summary.json` identical, which includes the digests of `.Random.seed` after every step. The hooks change nothing.
- **Without RunUMAP** (34 objects, production direct runs): every result file identical. The `result` sections differ only in the digests of `.Random.seed` after RunUMAP and later steps; GroupSingletons reseeds before it draws, so even objects with singletons are unaffected.
- **v5 against v3** (the two emptyDrops objects, production instrumented runs): every result file, every stage array and RNG checkpoint, every stage scalar and the Annoy index identical; the `result` sections differ only in the recorded assay class. C3i's check holds.
- **Cascade Lake against Sapphire Rapids,** portable profile: identical in every respect, the `result` sections included.

## R-vs-R floor

**Same kernel, different CPU types.** OpenBLAS's Haswell kernel on Cascade Lake (Xeon Platinum 8268) and Sapphire Rapids (8480CL), all 27 tier-T objects and soupx_tenx_SRR18256555, direct and instrumented: every file, all 1,205 stage arrays, every RNG checkpoint and the Annoy index identical (previous section). The portable references can be regenerated on another CPU type and checked by checksum, as A5b intends; as in M2, an AMD machine remains untested.

**Different kernels, same CPU type.** Haswell (portable) against SapphireRapids (production) on Sapphire Rapids, instrumented runs of all 34 objects, block by block (`data/clustering/compare/<object>/portable_vs_production.json`):

Tier T, 27 objects: identical across kernels in every object that dumps it: K10_genes, K5_nn_index, K7_final_ids, K7_optimizer_ids, K8_silhouette_other, K9_maximum, K9_weights, S2a_cells_step1, S2a_cells, S2a_gene_stats, S2a_genes_step1_pool, S2a_genes_step1, S2a_genes, S2a_log_gmean_step1, S2b_glm, S2b_model_pars_raw, S2b_model_pars, S2c_model_pars_fit, S2c_outliers, S2c_poisson_genes, S2c_smoothing_genes, S2d_pca_features, S2d_residual_stats, S2d_scale_data.chain, S2d_top_features, S2d_variable_features, S2d_y.chain.

Tier P, 7 objects: identical across kernels in every object that dumps it: K10_genes, K5_nn_index, K7_final_ids, K7_optimizer_ids, K8_silhouette_other, K9_maximum, K9_neighbors_idx, K9_neighbors_ptr, K9_weights, S2a_cells_step1, S2a_cells, S2a_density, S2a_gene_stats, S2a_genes_step1_pool, S2a_genes_step1, S2a_genes, S2a_log_gmean_step1, S2a_sampling_prob, S2b_glm, S2b_model_pars_raw, S2b_model_pars, S2c_model_pars_fit, S2c_outliers, S2c_poisson_genes, S2c_smoothing_genes, S2d_pca_features, S2d_top_features, S2d_variable_features.

Stage arrays that differ between the kernels. Embeddings and loadings are compared after aligning each component's sign; items and centroids are not, so their counts include whole components that changed sign (the floor statistics below align them):

| tier | array | objects | objects differing | elements differing | max abs diff |
|---|---|---|---|---|---|
| T | K3_items | 27 | 27 | 74,844 | 32.6 |
| T | K5_nn_dist | 27 | 2 | 2 | 1.91e-06 |
| T | K8_centroids | 27 | 27 | 2,447 | 8.24 |
| T | K8_cluster_var | 27 | 27 | 82 | 1.88e-12 |
| T | K8_silhouette_width | 27 | 27 | 5,871 | 2.41e-14 |
| T | K9_dist50 | 27 | 27 | 5,858 | 3.2e-12 |
| T | K9_neighbor_count | 27 | 6 | 12 |    0 |
| T | K9_neighbors_idx | 27 | 6 | 0 (lengths differ in 6) |    0 |
| T | K9_neighbors_ptr | 27 | 6 | 540 |    0 |
| T | K9_purity | 27 | 6 | 11 | 0.0142 |
| T | K9_sums | 27 | 6 | 12 | 0.0179 |
| T | K9_totals | 27 | 6 | 12 | 0.0179 |
| T | S3_embeddings | 27 | 27 | 294,404 | 4.27e-06 |
| T | S3_irlba_d | 27 | 27 | 1,307 | 2.15e-12 |
| T | S3_loadings | 27 | 27 | 1,054,669 | 1.19e-07 |
| T | S3_sdev | 27 | 27 | 1,302 | 1.85e-13 |
| P | K3_items | 7 | 7 | 2,844,663 |  400 |
| P | K5_nn_dist | 7 | 3 | 9 | 9.54e-07 |
| P | K8_centroids | 7 | 7 | 5,460 |  164 |
| P | K8_cluster_var | 7 | 7 | 179 | 1.75e-09 |
| P | K8_silhouette_width | 7 | 7 | 269,054 | 3.26e-13 |
| P | K9_dist50 | 7 | 7 | 268,809 | 6.43e-11 |
| P | K9_neighbor_count | 7 | 5 | 10 |    0 |
| P | K9_purity | 7 | 4 | 8 | 0.0187 |
| P | K9_sums | 7 | 5 | 10 | 0.00314 |
| P | K9_totals | 7 | 5 | 10 | 0.00314 |
| P | S2d_residual_stats | 7 | 7 | 1,430 | 1.42e-14 |
| P | S2d_scale_data.chain | 4 | 4 | 746,447 | 2.84e-14 |
| P | S3_embeddings | 7 | 7 | 13,540,801 | 1.53e-10 |
| P | S3_irlba_d | 7 | 7 | 336 | 6.93e-12 |
| P | S3_loadings | 7 | 7 | 1,045,990 | 1.56e-12 |
| P | S3_sdev | 7 | 7 | 333 | 4.44e-14 |

Per object (instrumented runs, portable against production, on Sapphire Rapids; for tier T, sums over the objects, and for the TRUE/FALSE columns the number of objects with TRUE):

| object | tier | cells | labels_identical | ari | silhouette_differing | purity_differing | sct_counts_identical | sct_scale_data_identical | pca_embeddings_max_abs_diff | item_components_flipped | items_differing_aligned | neighbour_lists_differing | neighbour_distances_differing | purity_neighbor_counts_differing | embeddings_1_30_max_abs_diff_aligned | snn_graph_identical | optimizer_ids_differing | rng_checkpoints_identical |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all 27 tier-T objects | T | 5904 | 27 | 1 | 5871 | 11 | 27 | 27 | 4.27045538284581e-06 | 341 | 7 | 0 | 2 | 12 | 3.65649732714246e-11 | 27 | 0 | 27 |
| soupx_tenx_SRR18256555 | P | 3108 | 1 | 1 | 3028 | 0 | 1 | 0 | 1.53340729003304e-10 | 10 | 5 | 0 | 0 | 0 | 1.53340729003304e-10 | 1 | 0 | 1 |
| soupx_tenx_10xmousebrain | P | 11831 | 1 | 1 | 11578 | 0 | 1 | 0 | 2.34194885706529e-11 | 13 | 0 | 0 | 0 | 2 | 1.87014848052058e-11 | 1 | 0 | 1 |
| filtered_seurat_tenx_10xmousebrain | P | 11831 | 1 | 1 | 11425 | 2 | 1 | 0 | 1.95612415154756e-11 | 9 | 1 | 0 | 0 | 2 | 1.70796710108334e-11 | 1 | 0 | 1 |
| soupx_emptydrops_4500871lane | P | 26598 | 1 | 1 | 26332 | 0 | 1 | 0 | 1.35852218363652e-10 | 17 | 9 | 0 | 2 | 0 | 3.95061761082616e-11 | 1 | 0 | 1 |
| soupx_emptydrops_sscovelli | P | 79409 | 1 | 1 | 79128 | 2 | 1 | 0 | 1.10429887456576e-10 | 7 | 39 | 0 | 4 | 2 | 1.10429887456576e-10 | 1 | 0 | 1 |
| filtered_seurat_emptydrops_sscovelli | P | 79409 | 1 | 1 | 78628 | 2 | 1 | 0 | 4.91198193230957e-11 | 9 | 33 | 0 | 2 | 2 | 4.78337369713699e-11 | 1 | 0 | 1 |
| cellbender_fromraw_ghirsutum2 | P | 59323 | 1 | 1 | 58935 | 2 | 1 | 0 | 1.15623066676562e-10 | 14 | 27 | 0 | 0 | 2 | 2.64925859028153e-11 | 1 | 0 | 1 |

- **Identical across kernels everywhere:** SCTransform's setup, fits and regularisation (S2a–S2c), the variable features, every RNG checkpoint (so irlba's draw count), the neighbour ids, both graphs, the optimizer's result, the labels and levels, and the corrected counts. On tier T, SCTransform's residual statistics and scale.data (S2d) are identical too; on tier P they differ in the last bits (at most 2.8e-14), as in M2.
- **The PCA** differs in every object. The 30 dimensions the neighbour search, silhouette and purity read differ by at most 1.5e-10 after aligning signs (2.3e-12 of the largest value on tier T); the later dimensions of the smallest tier-T objects, where the singular values crowd together, by up to 4.3e-6. Components often change sign between kernels (341 of tier T's 810 used components), which changes no distance.
- **Annoy's float items:** after aligning signs, 7 of tier T's 5,904 cells and 114 of tier P's 271,509 have an item value that differs: 116 of these 121 values by one float ulp, 3 by two to four, and two values near zero (2e-8, 1.8e-6) by more; at most 9.5e-7. No cell's neighbour ids change, and 10 cells' neighbour distances move by about a float ulp (at most 1.9e-6).
- **Silhouette widths** differ in the last bits for nearly every cell (at most 3.3e-13).
- **Purity:** the 50th distances and the radius move in the last bits, and a cell whose distance to another lies within an ulp of the radius gains or loses it: 22 neighbour counts changed in 11 objects, 19 purities in 10 objects (by up to 0.0187).
- **So:** a portable-profile run matches a production-profile reference in every label, graph and neighbour id, but not bitwise in the PCA, silhouette, a few purities or (on larger objects) scale.data. Exact checks of those need the reference's kernel, as plan K8 and K9 assume; the full panel (M17) uses each job's CPU type.

## Memory and time of R's chain

chain: from the loaded object to the end of add_silhouette_to_metadata (VmHWM at its end); script: the whole Rscript under /usr/bin/time -v, which adds R's startup, reading the export and writing the outputs.

| object | umap | cells | chain min | chain peak RSS GiB | script min | script peak GiB | load min | SCTransform min | RunPCA min | RunUMAP min | FindNeighbors min | FindClusters min | silhouette min | purity min | workflow job | workflow MaxRSS GiB | workflow min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | yes | 3,108 | 0.7 | 2.5 | 1.2 | 2.5 | 0.1 | 0.5 | 0.1 | 0.2 | 0.0 | 0.0 | 0.0 | 0.0 | 59051782 | 2.0 | 1.6 |
| soupx_tenx_SRR18256555 | no | 3,108 | 0.4 | 2.5 | 0.8 | 2.5 | 0.1 | 0.3 | 0.1 |  | 0.0 | 0.0 | 0.0 | 0.0 | 59051782 | 2.0 | 1.6 |
| soupx_tenx_10xmousebrain | yes | 11,831 | 1.8 | 12.0 | 4.9 | 12.0 | 0.5 | 1.3 | 0.2 | 0.3 | 0.1 | 0.0 | 0.0 | 0.0 | 59057762 | 12.0 | 5.8 |
| soupx_tenx_10xmousebrain | no | 11,831 | 1.2 | 12.0 | 3.6 | 12.0 | 0.4 | 1.0 | 0.1 |  | 0.0 | 0.0 | 0.0 | 0.0 | 59057762 | 12.0 | 5.8 |
| filtered_seurat_tenx_10xmousebrain | yes | 11,831 | 1.3 | 11.9 | 3.8 | 11.9 | 0.4 | 0.9 | 0.1 | 0.2 | 0.0 | 0.0 | 0.0 | 0.0 | 59051641 | 10.4 | 4.2 |
| filtered_seurat_tenx_10xmousebrain | no | 11,831 | 1.1 | 11.9 | 3.5 | 12.1 | 0.3 | 0.9 | 0.1 |  | 0.0 | 0.0 | 0.0 | 0.0 | 59051641 | 10.4 | 4.2 |
| soupx_emptydrops_4500871lane | yes | 26,598 | 2.9 | 23.4 | 4.0 | 23.4 | 0.4 | 1.7 | 0.3 | 0.5 | 0.1 | 0.1 | 0.0 | 0.2 | 59059692 | 22.9 | 7.4 |
| soupx_emptydrops_4500871lane | no | 26,598 | 2.3 | 23.4 | 3.4 | 23.4 | 0.4 | 1.7 | 0.3 |  | 0.1 | 0.1 | 0.0 | 0.2 | 59059692 | 22.9 | 7.4 |
| soupx_emptydrops_sscovelli | yes | 79,409 | 11.6 | 62.1 | 14.7 | 62.1 | 0.7 | 4.8 | 1.2 | 1.8 | 0.5 | 0.5 | 0.0 | 2.8 | 59399440 | 60.5 | 19.9 |
| soupx_emptydrops_sscovelli | no | 79,409 | 10.4 | 62.1 | 13.5 | 62.1 | 0.7 | 4.9 | 1.2 |  | 0.5 | 0.5 | 0.0 | 3.2 | 59399440 | 60.5 | 19.9 |
| filtered_seurat_emptydrops_sscovelli | yes | 79,409 | 13.6 | 57.5 | 16.8 | 57.5 | 0.7 | 5.0 | 1.2 | 2.0 | 0.6 | 0.6 | 0.0 | 4.3 | 59051578 | 56.9 | 20.2 |
| filtered_seurat_emptydrops_sscovelli | no | 79,409 | 10.0 | 57.5 | 13.0 | 57.5 | 0.6 | 4.7 | 1.1 |  | 0.5 | 0.6 | 0.0 | 3.0 | 59051578 | 56.9 | 20.2 |
| cellbender_fromraw_ghirsutum2 | yes | 59,323 | 14.6 | 104.1 | 19.1 | 104.1 | 0.6 | 8.9 | 1.0 | 1.4 | 0.4 | 0.4 | 0.0 | 2.6 | 59399468 | 87.5 | 15.6 |
| cellbender_fromraw_ghirsutum2 | no | 59,323 | 10.7 | 104.1 | 15.1 | 104.1 | 0.6 | 8.0 | 0.8 |  | 0.3 | 0.3 | 0.0 | 1.2 | 59399468 | 87.5 | 15.6 |

- **SCTransform sets the peak** in five of the seven objects (VmHWM after it already equals the run's peak): vst holds the residual matrix of every kept gene (genes × cells, dense) and `correct()`'s output of the same size, 56,106 genes × 59,323 cells for cellbender_fromraw_ghirsutum2, 104.1 GiB. On the two 79,409-cell objects RunPCA adds 1.9 to 2.9 GiB on top of SCTransform's peak (59.2 to 62.1 GiB and 55.6 to 57.5 GiB). The direct runs' exact peaks exceed the workflow jobs' sampled MaxRSS there (104.1 against 87.5 GiB), as in M2.
- **Time:** the chain takes 10.0 to 14.6 minutes on the large objects, SCTransform about half of it; purity 1.2 to 4.3 minutes; RunPCA about 1 minute; FindNeighbors and FindClusters under a minute each; approxSilhouette seconds. RunUMAP adds 0.2 to 2.0 minutes and no memory. These jobs shared their nodes with others (test partition), so differences of a minute between runs are noise; the per-step times are in each run's `summary.json`.
- These are the baseline of C5's performance row: Rust is to beat the chain without RunUMAP (and stay under 10 GiB on the 79,409-cell objects), with purity's share reported.

## Marker references

[`markers_clustering.R`](../../validation/R/markers_clustering.R) runs the steps of the workflow's `find_markers.R` (once per cluster, as the rule `find_markers` does) and `combine_markers.R`, copied with their source cited, on each saved `.rds`, and compares the combined table with the saved `<prefix>_markergenes.csv` byte for byte. One R process reads the object once and runs the clusters in turn (FindMarkers does not change the object and, with its defaults, draws no random numbers); `set.seed(12345)` is called before each cluster, as each job did.

**The order of the clusters.** The rule `combine_markers` receives the per-cluster files in the order Snakemake 9.16.3's `glob_wildcards()` lists the checkpoint's manifest directory: `os.walk`, the filesystem's directory order, with no sort. The saved tables list their clusters in orders such as 3, 1, 2, 0, not in level order. The reproduction combines the clusters in the order of the saved table (clusters without rows can go anywhere). This corrects plan C4, which said "in level order".

Tier T: 27 of 27 identical. Tier P:

| object | tier | clusters | saved cluster order | rows | identical | reference | min |
|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | P | 15 | 7 13 6 14 4 3 9 1 11 2 0 5 12 8 10 | 22,943 | TRUE | saved | 0.4 |
| soupx_tenx_10xmousebrain | P | 21 | 15 7 18 13 6 14 4 3 17 9 1 11 2 0 5 12 20 16 8 19 10 | 81,448 | TRUE | saved | 2.7 |
| filtered_seurat_tenx_10xmousebrain | P | 21 | 15 7 18 13 6 14 4 3 17 9 1 11 2 0 5 12 20 16 8 19 10 | 81,317 | TRUE | saved | 2.7 |
| soupx_emptydrops_4500871lane | P | 31 | 27 15 29 7 21 18 13 6 14 4 30 3 24 17 9 1 11 22 2 0 25 5 23 28 12 20 16 8 19 10 26 | 35,901 | TRUE | saved | 1.0 |
| soupx_emptydrops_sscovelli | P | 36 | 27 15 29 7 21 35 34 18 13 6 14 4 31 30 3 24 17 9 1 11 22 33 2 0 25 5 23 32 28 12 20 16 8 19 10 26 | 57,493 | TRUE | saved | 3.7 |
| filtered_seurat_emptydrops_sscovelli | P | 37 | 27 15 29 7 21 35 34 18 13 6 14 36 4 31 30 3 24 17 9 1 11 22 33 2 0 25 5 23 32 28 12 20 16 8 19 10 26 | 58,776 | TRUE | saved | 3.6 |
| cellbender_fromraw_ghirsutum2 | P | 22 | 15 7 21 18 13 6 14 4 3 17 9 1 11 2 0 5 12 20 16 8 19 10 | 44,222 | TRUE | saved | 3.6 |

rows: the combined table's rows; min: the whole check, reading the object included (FindMarkers through presto). The saved cluster orders are those of the files' directory listing.

## Instrumented optimizer and K0 vectors

**The optimizer.** [`validation/cpp/build.sh`](../../validation/cpp/build.sh) compiles Seurat 5.5.1's `src/ModularityOptimizer.cpp` in its standalone mode (`-DSTANDALONE`, the file's own `main()`), with the R env's compiler and the flags of its `Makeconf` (`R CMD config`: `x86_64-conda-linux-gnu-c++ -std=gnu++17 -O2 -march=nocona -mtune=haswell ...`), into `data/clustering/bin/`. [`validation/cpp/ModularityOptimizer.patch`](../../validation/cpp/ModularityOptimizer.patch) adds nine lines that print, after each random start, its iterations, its modularity with `%.17g` and its partition; nothing else changes. The R env's compiler is GCC 16.2.0, while `Seurat.so` was built with GCC 14.3.0 (conda-forge); with the same flags (`-march=nocona`: no FMA, no AVX) the arithmetic is the same, which the results confirm.

[`optimizer_clustering.R`](../../validation/R/optimizer_clustering.R) writes each saved SNN graph's strictly lower triangle as the edge list `main()` reads (`col<TAB>row<TAB>weight`, 0-based, in stored order: the order `RunModularityClusteringCpp` walks the matrix, so `readInputFile()` builds R's network; the weights are copied as written with `%.17g`, which `std::stod` reads exactly), runs the optimizer with FindClusters' settings (modularity 1, resolution 0.8, Louvain, 10 starts, 10 iterations, seed 12345) and compares its final partition with R's optimizer result (the production instrumented run's `K7_optimizer_ids`, before GroupSingletons), with the saved labels (for objects without singletons) and its printed numbers with R's and, for tier P, with the workflow job's log.

Tier T: 27 of 27 objects reproduce R's optimizer ids, R's printed numbers and the saved labels (no tier-T object has singletons); in 20, the random starts end at more than one modularity. Tier P (= saved labels is FALSE where GroupSingletons merged singletons after the optimizer):

| object | tier | nodes | edges | Q | communities | distinct start modularities | best start | iterations | = R's optimizer ids | = R's printed numbers | = saved labels | = workflow log |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soupx_tenx_SRR18256555 | P | 3,108 | 109,910 | 0.8777 | 15 | 10 | 10 | 3-3 | TRUE | TRUE | TRUE | TRUE |
| soupx_tenx_10xmousebrain | P | 11,831 | 361,793 | 0.8920 | 21 | 10 | 2 | 3-3 | TRUE | TRUE | TRUE | TRUE |
| filtered_seurat_tenx_10xmousebrain | P | 11,831 | 362,146 | 0.8898 | 21 | 10 | 2 | 3-4 | TRUE | TRUE | TRUE | TRUE |
| soupx_emptydrops_4500871lane | P | 26,598 | 728,347 | 0.9015 | 32 | 10 | 1 | 3-3 | TRUE | TRUE | FALSE | TRUE |
| soupx_emptydrops_sscovelli | P | 79,409 | 1,564,680 | 0.8818 | 54 | 10 | 8 | 3-4 | TRUE | TRUE | FALSE | TRUE |
| filtered_seurat_emptydrops_sscovelli | P | 79,409 | 1,480,322 | 0.8722 | 59 | 10 | 9 | 3-4 | TRUE | TRUE | FALSE | TRUE |
| cellbender_fromraw_ghirsutum2 | P | 59,323 | 1,046,981 | 0.8567 | 62 | 10 | 4 | 3-4 | TRUE | TRUE | FALSE | TRUE |

distinct start modularities: how many different values the 10 random starts reached (as `%.17g`); best start: the first start at the maximum, which the optimizer keeps; iterations: the range of Louvain iterations per start. On tier T, 13 objects have starts that end within 1e-12 of the best start's modularity without equalling it: all 72 such starts found the best start's partition, numbered differently, so the modularity's sums ran in another order. No tier-P object has such a near-tie; their 10 starts end at 10 different partitions.

**K0 vectors** (`tests/fixtures/clustering/k0/`, 245 KB): `javarandom_vectors` links the unmodified `ModularityOptimizer.cpp` and prints `JavaRandom.nextInt(n)` from fresh generators (seeds as the `int` FindClusters passes, negative ones included; n from 1 to 2^31 − 1, powers of two included) and successive `generateRandomPermutation` draws from one generator; `new Random(0).nextInt(10)` gives 0, 8, 9, 7, 5, 3, 1, 1, 9, 4, as Java does. `kiss64_vectors` includes RcppAnnoy 0.0.23's installed `kissrandom.h` (identical to the source tarball's) and prints `kiss()`, `index(n)`, `flip()` and one stream through mixed calls, from the default seed every Seurat build uses (1234567890987654321) and others. [`k0_vectors.R`](../../validation/R/k0_vectors.R) gives R's `set.seed(1); sample(x, 1)` for 1 to 100,000 tied names, the draw GroupSingletons makes, with `.Random.seed[2]` afterwards: one to three uniforms, depending on n (rejection sampling: three for n = 65–67 and 257 after `set.seed(1)`).

## Instrumented runs

[`instrumented_clustering.R`](../../validation/R/instrumented_clustering.R) runs the same code as the direct run with passive `trace()` hooks ([`lib/clustering_hooks.R`](../../validation/R/lib/clustering_hooks.R)). The SCTransform and RunPCA hooks are M2's, unchanged (`lib/doubletfinder_hooks.R`); the new ones are on `sctransform::correct`, `RunUMAP.Seurat`, Seurat's `AnnoyBuildIndex`, `AnnoySearch`, `RunModularityClustering` and `GroupSingletons`, and bluster's `approxSilhouette` and `neighborPurity`. They read values and write files; they draw nothing from R's stream and change nothing. Two hooks compute something: the Annoy hook builds a second index from the same rows (a new object, `addItem` in cell order, `build(50)`) and saves it, because `save()` reloads the file and re-derives the roots, which the index FindNeighbors searches must not do; and the purity hook calls `findDistance()` on neighborPurity's own index with neighborPurity's arguments, to get the per-cell 50th distances whose median neighborPurity keeps. Every instrumented run also checks that each hook fired as often as expected (`stages.json`: `hook_counts`, `problems`).

**Proof that the hooks change nothing.** For every pair of instrumented and direct runs in the same profile, the final files are byte-identical and so is the `result` section of `summary.json`, which includes a digest of `.Random.seed` after every step (table above).

**What is dumped** (`stages/index.tsv` maps every array to its stage; one context, `chain`):

| Stage | Arrays |
|---|---|
| S2a–S2d (K1) | M2's SCTransform dumps at seed 12345: cell attributes, filtered genes, gene statistics, fitting cells and genes, every glmGamPoi fit, model_pars before and after the floor, model_pars_fit, outliers, Poisson and smoothing genes, residual statistics, the top and variable features, the PCA input (levels full and standard), the clipped residual matrix (level full); `vst/`: SCTransform.default's result in `scrust sctransform`'s layout |
| S3 (K2) | embeddings, loadings, sdev, irlba's d; iterations and products |
| K3 | `K3_items`: the items Annoy holds, `getItemsVector(i)` for every cell |
| K4 | `annoy.ann`: the second build's index, as `save()` writes it |
| K5 | `K5_nn_index` (1-based), `K5_nn_dist` (Annoy's floats as doubles) |
| K7 | `K7_optimizer_ids` (before GroupSingletons, numbered by size), `K7_final_ids` |
| K8 | silhouette width and nearest other cluster, centroids, within-cluster variances |
| K9 | each cell's 50th distance, the neighbour counts and (level full, and the small P object) lists, weights, per-cluster sums, totals, purity, maximum |
| K10 | the genes `correct()` corrects; its median log_umi and the residuals beyond the clip range (`stages.json`) |
| RNG | `.Random.seed` at vst's entry, after the cell sampling, at each fitting bin, at `correct()`'s entry, at vst's and SCTransform's exits, around irlba and RunPCA, around RunUMAP, after the Annoy build and search, around the optimizer and GroupSingletons, at the chain's end |

Facts every production instrumented run confirmed (`stages.json`): Annoy's items equal R's own float rounding of the embeddings' first 30 columns (`writeBin(size = 4)`); the reloaded index returns the same neighbours and distances for every cell; every cell has itself among its neighbours, and is its own first neighbour except for the 6 duplicated cells of cellbender_fromraw_ghirsutum2; neighborPurity builds a `KmknnIndex` (BiocNeighbors' default for `BNPARAM = NULL`), and its radius is the median of the 50th distances, bitwise; `correct()` regresses to R's `median()` of log_umi; and in every object some of the residuals `correct()` reads lie beyond vst's clip range (50 to 5,455).

## Findings for the port

- **The corrected counts come from `sctransform::correct()`, not `correct_counts()`.** With Seurat's default path, `vst()` itself calls `correct()` (`return_corrected_umi = TRUE`, from `do.correct.umi`) on its own residual matrix, before it clips it; `correct_counts()`, which recomputes the residuals unfloored, serves only `conserve.memory`. So the residuals `correct()` uses are floored at `min_variance` and not clipped, and every object has some beyond the clip range (87 of `soupx_tenx_test`'s; `K10` in `stages.json`). The corrected value is μ + r·√(μ + μ²/θ) at the median log_umi, with no floor on that variance; then `round(x, 0)`, negatives to 0, zeros dropped. `data` is `log1p` of it, bitwise, in every object, and `nCount_SCT` and `nFeature_SCT` are its column sums and nonzero counts.
- **`scale.data` keeps its rows in the SCT assay's gene order,** not in `VariableFeatures()` order, in every object and through both assay paths; RunPCA reorders them, and the loadings' rows follow `VariableFeatures()`. A loader rebuilding the assay has to keep that order to compare equal (C3g, K11).
- **The graphs** (all 34 objects): `SCT_snn` is symmetric, with rows sorted within each column and the whole diagonal stored (value 1); `SCT_nn` has exactly k = 20 entries per row, one of them the cell itself. The SNN's values are the Jaccard values s/(40 − s) for s = 3 to 20 shared neighbours (s ≤ 2 is pruned): all 18 occur in every tier-P object, 16 to 18 in the tier-T ones.
- **Annoy.** The items are exactly R's float rounding of the embeddings' first 30 columns, and an index saved with `save()` and reloaded searches identically, in every object. Every cell is its own first neighbour except in cellbender_fromraw_ghirsutum2, where 11 cells in 5 groups have identical items (cells with 1 to 4 UMIs, as M2 found in its DoubletFinder PCA): the 6 with a lower-numbered twin get the twin first (equal distance 0, lower id first, C3i) and still have themselves in their list. M9 has a real case of that tie order.
- **The optimizer's near-ties are real.** In 13 of the 27 tier-T objects, some random starts end within 1e-12 of the best start's modularity without equalling it. Every such start (72, tier T and P) found the best start's partition, numbered differently, so `calcQualityFunction`'s sums ran in another order (`soupx_tenx_test`: 0.33019971131680448 and 0.33019971131680459; the first start at the maximum wins). M10 has to reproduce every sum's order, as C3c says; the per-start values are in `optimizer/starts.tsv`.
- **The v3 path equals the v5 path.** On both emptyDrops objects, `SCTransform.Assay` called directly (v3, as the workflow did) and through `SCTransform.StdAssay` (v5) give identical stage arrays, RNG checkpoints and outputs; only the recorded assay class differs.
- **RunUMAP changes nothing downstream:** without it, every file of every object is identical; only the digests of `.Random.seed` after RunUMAP and later steps differ (GroupSingletons reseeds before its draws).
- **Purity across kernels.** The radius and the 50th distances move in the last bits with the PCA, and a cell whose distance to another lies within an ulp of the radius can gain or lose that neighbour: 22 neighbour counts in 11 of the 34 objects, and 19 purities (by up to 0.0187) in 10. Exact purity checks therefore need the reference's kernel, like silhouette's, and M12's certificate has real near-boundary pairs to find.
- **Marker tables** list their clusters in the filesystem's order (Snakemake's `glob_wildcards()`), not in level order.
- **Snakemake touches outputs** after their jobs end, so files' modification times say little about which job wrote them.
- **Test-reference logs** were overwritten by the second test run; tier T's optimizer numbers come from the reference runs.

## Corrections to the plan

Corrected in [`docs/PLAN.md`](../PLAN.md) on this branch:

- **`correct()`, not `correct_counts()`** (A3, A7b, C1, C2, C3, C3f, C3i, C4, C6's M8 and M13): Seurat's corrected counts come from `sctransform::correct()`, called inside `vst()` on vst's residuals before clipping; C3f now describes that port (the residuals floored at `min_variance` but not clipped; the median-depth variance not floored), and C3i and K10's hazards follow.
- **The marker table's cluster order** (C3g, C4): the workflow combines the per-cluster files in the order Snakemake's `glob_wildcards()` lists them, the filesystem's directory order, not in level order; the loader's test takes the order from the saved table.
- **`meta.data` reaches `vst()` as `cell_attr`** (C1, new): harmless for the workflow's objects, which lack vst's attribute names, but it is part of the code path.
- **`scale.data`'s row order** (C3g, new): Seurat keeps the rows in the SCT assay's gene order, not in `VariableFeatures()` order; the Seurat loader has to keep both orders.

Plan facts the data confirmed as written: all 840 production objects and all 27 test objects were clustered on Sapphire Rapids; 36 production logs report singletons, up to 40; the tier-P objects' sizes and singletons (C4's table); the emptyDrops objects' v3 assay; the SNN graph's layout and diagonal; Annoy's node size (16 + 4f bytes); RunUMAP's lack of downstream effect; GroupSingletons' `set.seed(1)`; neighborPurity's KMKNN index and median radius; Seurat 5.5.1 and RcppAnnoy built with GCC 14.3 and `-march=nocona`; Seurat built against RcppEigen 0.3.4.0.2 and RcppAnnoy 0.0.23.

## Fixtures

**Committed** (`tests/fixtures/clustering/`, 1.90 MB; [its README](../../tests/fixtures/clustering/README.md)), copied from `data/clustering/` by [`validation/slurm/fixtures_clustering.sh`](../../validation/slurm/fixtures_clustering.sh), byte for byte or gzipped (`gzip -n -9`, read through `scio::gz`):

- `k0/`: the K0 vectors (245 KB);
- for the seven objects cargo test covers (plan C4: the two base objects and the five DoubletFinder inputs, whose counts are M2's fixtures): the expected labels and scores of the portable profile (`reference/`) and of the saved object (`saved/`), and the instrumented optimizer's per-start results;
- for the two base objects also their counts and the stage references M9 and M10 test from R's own embeddings (the embeddings, `annoy.ann`, every cell's neighbours, both graphs, the optimizer's, silhouette's and purity's arrays).

`tests/fixtures/` now holds 23.07 MB (21.17 MB before M8, after Step 1). That is above A1's "about 20 MB": Step 1 had already passed it, and M8 adds 1.90 MB, less than one of M2's tier-T DoubletFinder objects (3.1–3.4 MB each). The stage references of the five DoubletFinder inputs, and everything for the other 20 tier-T objects, stay in `data/` (each about 0.5 MB compressed); committing them would add about 2.5 MB. Moving the base objects' stage references to `data/` too would save 1.4 MB.

**Data** (`data/clustering/`, gitignored, 19.6 GB: 18.6 GB for the seven tier-P objects, most of it scale.data and the residual-sized dumps at level standard, about 1 GB for the rest): every export, run, marker check, optimizer run and comparison, listed in `validation/manifest.toml` (5,629 file entries with 7,847 byte-identical copies recorded in `identical_on`, and 359 `[[run]]` tables; the manifest grows from 3.2 to 6.0 MB). The manifest's clustering section uses a compact form, described in its header: one `[[run]]` table per run with what its files share, and one short `[[file]]` table per distinct file, with `identical_on` naming the other runs whose copy has the same sha256.

**Tolerances.** `validation/tolerances.toml` gains the clustering entries M8's own checks need: the saved objects (exact in the production profile), the marker references (byte-identical), the instrumented optimizer (exact), and the floor as measurements. The stage criteria K0–K13 are settled at their milestones.

## Reproducing

From the repository's top directory, after `source env.sh`, `cargo build --release -p scrust` (for `scrust fingerprint`), `third_party/fetch_sources.sh` and `bash validation/cpp/build.sh`:

```bash
sbatch validation/slurm/export_clustering.sbatch soupx_tenx_test
sbatch validation/slurm/reference_clustering.sbatch production direct soupx_tenx_test
sbatch validation/slurm/reference_clustering.sbatch production instrumented soupx_tenx_test
sbatch validation/slurm/reference_clustering.sbatch portable instrumented soupx_tenx_test
sbatch validation/slurm/markers_clustering.sbatch soupx_tenx_test
bash validation/slurm/optimizer_clustering.sbatch soupx_tenx_test
bash validation/slurm/compare_clustering.sbatch soupx_tenx_test
```

Memory per object is in `validation/slurm/README.md`. The tables in this report come from `Rscript validation/R/report_clustering_m8.R docs/validation/clustering_m8.tsv`, the panel table from `validation/R/clustering_panel.R`, the manifest's clustering section from `validation/R/manifest_clustering.R` and the committed fixtures from `validation/slurm/fixtures_clustering.sh`; all read plain files only. SLURM logs of every M8 job are in `data/clustering/logs/`.
