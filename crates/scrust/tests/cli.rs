//! The `scrust` binary: `--version`, `fingerprint` (environment fingerprint, M2),
//! `doubletfinder`, `validate` and `knn` (M3), `sctransform` (M4), `pca` (M5), and the whole
//! `doubletfinder` from the counts (M6), `neighbors` (M9), `cluster` (M10), and the clustering
//! chain from the counts and from each stage, its parameters, and kernel-gated validation (M11).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn scrust(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_scrust"))
        .args(args)
        .output()
        .expect("failed to run the scrust binary")
}

#[test]
fn version() {
    let out = scrust(&["--version"]);
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        format!("scrust {}\n", env!("CARGO_PKG_VERSION"))
    );
}

/// The whole tool from the counts and the cluster labels (milestone M6), on a tier-T fixture:
/// the same files as R's direct reference run, byte for byte, the same `result` (the digest of
/// R's final `.Random.seed` included), every stage array, RNG checkpoint and scalar of R's
/// instrumented run, and the saved workflow output; `singlets/` (M6b) holds the saved output's
/// kept cells. The products are OpenBLAS calls and the fixtures are portable runs: the test first
/// checks that the kernel in use is the fixture's (Haswell, env.sh's default).
#[test]
fn doubletfinder_runs_every_stage_from_the_counts() {
    let f = fixture("decontx_raw_emptydrops_test");
    rcompat::blas::require_fixture_kernel(&f, "scrust doubletfinder from the counts");
    let dir = tmp("decontx_raw_emptydrops_test_m6");
    let out = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--out", &dir.display().to_string(), "--dump-stages", "light", "--threads", "2",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
        assert!(std::fs::read(f.join("reference").join(file)).unwrap() == std::fs::read(dir.join(file)).unwrap(), "{file}");
    }
    let summary = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    let want = std::fs::read_to_string(f.join("reference/summary.json")).unwrap();
    let digest = |s: &str| s.lines().find(|l| l.contains("\"rng_end_sha256\"")).map(str::trim).map(String::from);
    assert_eq!(digest(&summary), digest(&want), "the final .Random.seed differs from R's");
    let report = dir.join("validate.json");
    let v = scrust(&[
        "validate", "doubletfinder", "--reference", &f.join("reference").display().to_string(),
        "--reference-stages", &f.join("stages").display().to_string(), "--saved", &f.join("saved").display().to_string(),
        "--candidate", &dir.display().to_string(), "--json", &report.display().to_string(),
    ]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    let r = scio::Json::parse(&std::fs::read_to_string(&report).unwrap()).unwrap();
    assert_eq!(r.get("pass"), Some(&scio::Json::Bool(true)));
    let Some(scio::Json::List(by_stage)) = r.get("by_stage") else { panic!("no by_stage") };
    for stage in ["S1", "S2a", "S2b", "S2c", "S2d", "S3", "S4", "S5", "S6", "S7", "saved"] {
        let e = by_stage.iter().find(|e| e.get("stage").and_then(|v| v.as_str()) == Some(stage)).unwrap_or_else(|| panic!("no {stage} checks"));
        assert!(e.get("checks_identical").and_then(|v| v.as_f64()).unwrap() > 0.0, "{stage}");
        assert_eq!(e.get("checks_differing").and_then(|v| v.as_f64()), Some(0.0), "{stage}");
    }
    assert_eq!(r.get("grading").and_then(|g| g.get("mode")).and_then(|m| m.as_str()), Some("exact"));
    // the same run, recorded as if under another kernel: tolerance mode, said so, and still a
    // pass (nothing differs); the kernel check no longer fails the report
    let other = tmp("decontx_raw_emptydrops_test_m6_other_kernel");
    std::fs::create_dir_all(&other).unwrap();
    for f in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
        std::fs::copy(dir.join(f), other.join(f)).unwrap();
    }
    let copy_dir = |from: &Path, to: &Path| {
        std::fs::create_dir_all(to).unwrap();
        for e in std::fs::read_dir(from).unwrap() {
            let e = e.unwrap();
            std::fs::copy(e.path(), to.join(e.file_name())).unwrap();
        }
    };
    copy_dir(&dir.join("stages"), &other.join("stages"));
    std::fs::write(other.join("summary.json"), summary.replace("\"openblas_corename\": \"Haswell\"", "\"openblas_corename\": \"SapphireRapids\"")).unwrap();
    let report2 = other.join("validate.json");
    let v = scrust(&[
        "validate", "doubletfinder", "--reference", &f.join("reference").display().to_string(),
        "--reference-stages", &f.join("stages").display().to_string(), "--candidate", &other.display().to_string(), "--json", &report2.display().to_string(),
    ]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    let text = std::fs::read_to_string(&report2).unwrap();
    assert!(text.contains("\"mode\": \"tolerance\""), "tolerance mode");
    assert!(text.contains("the reference recorded OpenBLAS kernel Haswell, the candidate SapphireRapids"), "the reason");
    assert!(text.contains("\"class\": \"BLAS-dependent\"") && text.contains("\"class\": \"BLAS-free\""), "classes");
    assert!(String::from_utf8_lossy(&v.stderr).contains("grading: tolerance mode"));
    // singlets/: the saved output's kept cells, a 10x directory any reader takes
    let kept: Vec<String> = std::fs::read_to_string(f.join("saved/cells.tsv")).unwrap().lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect();
    let singlets = dir.join("singlets");
    assert_eq!(scio::tenx::read_barcodes(&singlets).unwrap(), kept);
    let m = scio::sparse::read_mtx(&singlets).unwrap();
    assert_eq!((m.nrow, m.ncol), (2000, kept.len()));
    assert!(summary.contains(&format!("\"cells\": {}", kept.len())), "summary.json's singlets record");
}

fn fixture(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder").join(object)
}

fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli").join(name);
    let _ = std::fs::remove_dir_all(&d);
    d
}

/// Teacher-forced run on a tier-T fixture: the same files as R's reference run, and
/// `scrust validate` finds nothing that differs.
#[test]
fn doubletfinder_reproduces_a_reference_run() {
    let f = fixture("soupx_emptydrops_test");
    let dir = tmp("soupx_emptydrops_test");
    let out = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--from-stages", &f.join("stages").display().to_string(),
        "--out", &dir.display().to_string(), "--dump-stages", "light", "--threads", "3",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
        assert!(std::fs::read(f.join("reference").join(file)).unwrap() == std::fs::read(dir.join(file)).unwrap(), "{file}");
    }
    let summary = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    assert!(summary.contains("\"pANN_column\": \"pANN_0.25_0.28_30\""));
    let report = dir.join("validate.json");
    let v = scrust(&[
        "validate", "doubletfinder", "--reference", &f.join("reference").display().to_string(),
        "--reference-stages", &f.join("stages").display().to_string(),
        "--candidate", &dir.display().to_string(), "--json", &report.display().to_string(),
    ]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    assert!(std::fs::read_to_string(&report).unwrap().contains("\"pass\": true"));
    // a second run into the same directory is refused
    let again = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--from-stages", &f.join("stages").display().to_string(), "--out", &dir.display().to_string(),
    ]);
    assert_eq!(again.status.code(), Some(2));
}

/// `scrust knn --output count` on the final run's embeddings gives pANN's numerators.
#[test]
fn knn_counts_give_the_final_pann() {
    let f = fixture("decontx_raw_tenx_test");
    let dir = tmp("knn");
    std::fs::create_dir_all(&dir).unwrap();
    let counts = dir.join("counts.npy");
    // 300 real cells, k = round(400 * 0.09) = 36
    let out = scrust(&[
        "knn", "--from-stages", &f.join("stages").display().to_string(), "--context", "final",
        "--queries", "300", "--neighbours", "37", "--output", "count", "--out", &counts.display().to_string(),
        "--threads", "2",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let c = scio::Npy::read(&counts).unwrap();
    let pann = scio::Npy::read(&f.join("stages/S5_pann_final.npy")).unwrap();
    for i in 0..300 {
        assert_eq!((c.get(i, 0) / 36.0).to_bits(), pann.get(i, 0).to_bits(), "cell {i}");
    }
}

#[test]
fn fingerprint_reports_the_pinned_openblas_kernel() {
    // OpenBLAS reads OPENBLAS_CORETYPE when it loads, so set it for the child process: the
    // portable profile's kernel, which every x86-64 CPU with AVX2 can run (docs/PLAN.md A5b).
    let out = Command::new(env!("CARGO_BIN_EXE_scrust"))
        .arg("fingerprint")
        .env("OPENBLAS_CORETYPE", "Haswell")
        .env("OPENBLAS_NUM_THREADS", "1")
        .output()
        .expect("failed to run the scrust binary");
    assert!(out.status.success());
    let json = String::from_utf8_lossy(&out.stdout);
    for line in [
        "\"openblas_corename\": \"Haswell\"",
        "\"openblas_num_threads\": 1",
        "\"OPENBLAS_CORETYPE\": \"Haswell\"",
    ] {
        assert!(json.contains(line), "missing {line} in:\n{json}");
    }
    assert!(json.contains("\"glibc\": \"2."), "no glibc version in:\n{json}");
    assert!(json.contains("libopenblas"), "OpenBLAS library not found in:\n{json}");
}

#[test]
fn unknown_command_is_a_usage_error() {
    let out = scrust(&["no-such-tool"]);
    assert_eq!(out.status.code(), Some(2));
}

/// SHA-256, for comparing a large output with the digest R recorded.
fn sha256_hex(data: &[u8]) -> String {
    scio::sha256::hex(data)
}

/// `scrust sctransform` on a tier-T object's counts gives R's SCTransform() outputs byte for
/// byte (validation/R/sctransform_cli_reference.R; scale.data by its SHA-256). Its products are
/// OpenBLAS calls, so it first checks that the kernel in use is the reference's (Haswell).
#[test]
fn sctransform_matches_r_byte_for_byte() {
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/sctransform/cli_soupx_tenx_test");
    rcompat::blas::require_fixture_kernel(&reference, "scrust sctransform against R");
    let out = tmp("sctransform");
    let run = scrust(&[
        "sctransform", "--counts", &fixture("soupx_tenx_test").join("input").display().to_string(),
        "--out", &out.display().to_string(), "--threads", "3",
    ]);
    assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
    for f in ["cell_attr.tsv", "gene_attr.tsv", "model_pars.tsv", "model_pars_fit.tsv", "variable_features.tsv"] {
        let got = std::fs::read(out.join(f)).unwrap();
        let want = std::fs::read(reference.join(f)).unwrap();
        assert!(got == want, "{f} differs from R's");
    }
    let digest = std::fs::read_to_string(reference.join("scale_data.npy.sha256")).unwrap();
    assert_eq!(sha256_hex(&std::fs::read(out.join("scale_data.npy")).unwrap()), digest.trim(), "scale_data.npy differs from R's");
    assert!(out.join("summary.json").is_file());
}

/// `scrust pca` on `scrust sctransform`'s output for a tier-T object gives R's RunPCA() outputs
/// (Seurat's default npcs = 50) byte for byte (validation/R/runpca_reference.R on R's
/// SCTransform output; large files by their SHA-256). Its products are OpenBLAS calls, so it first
/// checks that the kernel in use is the reference's (Haswell).
#[test]
fn pca_matches_r_byte_for_byte() {
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/pca/cli_soupx_tenx_test");
    rcompat::blas::require_fixture_kernel(&reference, "scrust pca against R");
    let sct = tmp("pca_sct");
    let run = scrust(&[
        "sctransform", "--counts", &fixture("soupx_tenx_test").join("input").display().to_string(),
        "--out", &sct.display().to_string(), "--threads", "2",
    ]);
    assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
    for threads in ["1", "3"] {
        let out = tmp(&format!("pca_t{threads}"));
        let run = scrust(&["pca", "--sct", &sct.display().to_string(), "--out", &out.display().to_string(), "--threads", threads]);
        assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
        let got = std::fs::read(out.join("stdev.tsv")).unwrap();
        assert!(got == std::fs::read(reference.join("stdev.tsv")).unwrap(), "stdev.tsv differs from R's");
        for f in ["embeddings.tsv", "loadings.tsv", "random_seed.txt"] {
            let digest = std::fs::read_to_string(reference.join(format!("{f}.sha256"))).unwrap();
            assert_eq!(sha256_hex(&std::fs::read(out.join(f)).unwrap()), digest.trim(), "{f} differs from R's ({threads} threads)");
        }
        assert!(out.join("summary.json").is_file());
    }
}

/// `scrust neighbors` on a base object's PCA embeddings (milestone M9's committed fixture, the
/// portable instrumented R run of M8) at the workflow's settings, then `scrust validate
/// neighbors` against that R run: the parameters, Annoy's items and index file (K3, K4), every
/// cell's neighbours (K5) and both graphs (K6) identical. No BLAS: any kernel, any machine.
#[test]
fn neighbors_match_r_and_validate() {
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/filtered_seurat_emptydrops_test/reference");
    for threads in ["1", "4"] {
        let out = tmp(&format!("neighbors_t{threads}"));
        let run = scrust(&[
            "neighbors", "--pca", &reference.join("pca").display().to_string(), "--out", &out.display().to_string(),
            "--dims", "1:30", "--annoy", "yes", "--threads", threads,
        ]);
        assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
        let report = out.join("validate.json");
        let v = scrust(&[
            "validate", "neighbors", "--reference", &reference.display().to_string(), "--candidate", &out.display().to_string(),
            "--json", &report.display().to_string(),
        ]);
        assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
        let r = scio::Json::parse(&std::fs::read_to_string(&report).unwrap()).unwrap();
        assert_eq!(r.get("pass"), Some(&scio::Json::Bool(true)));
        let text = std::fs::read_to_string(&report).unwrap();
        for stage in ["\"K3\"", "\"K4\"", "\"K5\"", "\"K6\""] {
            assert!(text.contains(stage), "{stage} not checked");
        }
        let ids = std::fs::read_to_string(out.join("knn_index.tsv")).unwrap();
        assert!(ids.starts_with("barcode\tnn_1\t"));
        assert_eq!(ids.lines().count(), 268);
    }
    // the same files at one and at four threads
    for f in ["knn_index.tsv", "knn_distance.tsv", "nn.mtx.gz", "snn.mtx.gz", "annoy.ann"] {
        let a = std::fs::read(PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli/neighbors_t1").join(f)).unwrap();
        let b = std::fs::read(PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli/neighbors_t4").join(f)).unwrap();
        assert!(a == b, "{f} depends on the thread count");
    }
    // refusals: another code path, and dims beyond the PCA
    let out = tmp("neighbors_refused");
    let pca = reference.join("pca").display().to_string();
    let r = scrust(&["neighbors", "--pca", &pca, "--out", &out.display().to_string(), "--metric", "hamming"]);
    assert!(!r.status.success());
    let r = scrust(&["neighbors", "--pca", &pca, "--out", &out.display().to_string(), "--dims", "1:51"]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("More dimensions specified in dims"));
}

/// `scrust cluster --snn` on a base object's SNN graph (M8's committed fixture, the portable
/// instrumented R run) at the workflow's settings, then `scrust validate cluster` against that R
/// run and the instrumented optimizer (K7): the optimizer's ids, the final ids, the labels and
/// levels, the workflow's files, every start's iterations, modularity and partition. The same
/// files at one and at four threads, with five resolutions. Then `scrust neighbors`' own output
/// as the input (the chain's layout), and the refusals: Leiden, another starting point,
/// resolutions that would overwrite each other. No BLAS: any kernel, any machine.
#[test]
fn cluster_matches_r_and_validates() {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/filtered_seurat_tenx_test");
    let reference = base.join("reference");
    let out = tmp("cluster_workflow");
    let run = scrust(&[
        "cluster", "--snn", &reference.join("graphs").display().to_string(), "--out", &out.display().to_string(),
        "--seed", "12345", "--keep-stages",
    ]);
    assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
    let err = String::from_utf8_lossy(&run.stderr);
    assert!(err.contains("Maximum modularity in 10 random starts: "), "prints as Seurat does");
    // --stats asks for silhouette and purity by default: from the SNN graph alone neither can be
    // computed (both need the PCA), and the command says so
    assert!(err.contains("silhouette asked for, but not computed"), "{err}");
    assert!(err.contains("purity asked for, but not computed"), "{err}");
    let report = out.join("validate.json");
    let v = scrust(&[
        "validate", "cluster", "--reference", &reference.display().to_string(), "--candidate", &out.display().to_string(),
        "--optimizer", &base.join("optimizer").display().to_string(), "--json", &report.display().to_string(),
    ]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    let text = std::fs::read_to_string(&report).unwrap();
    assert!(text.contains("\"pass\": true"));
    for check in ["optimizer ids", "final ids", "labels", "levels", "every start (iterations, modularity)", "every start's partition", "cluster_ids.txt", "nclusters.txt"] {
        assert!(text.contains(check), "{check} not checked");
    }
    for f in ["seurat_clusters.tsv", "seurat_clusters_levels.txt", "cluster_ids.txt", "nclusters.txt", "SCT_snn_res.0.8_levels.txt"] {
        assert!(out.join(f).is_file(), "{f}");
    }
    let sc = std::fs::read_to_string(out.join("seurat_clusters.tsv")).unwrap();
    assert!(sc.starts_with("barcode\tseurat_clusters\n"));
    assert_eq!(sc.lines().count(), 301);
    // five resolutions, at one and at four threads: the same files
    for t in ["1", "4"] {
        let o = tmp(&format!("cluster_res_t{t}"));
        let r = scrust(&[
            "cluster", "--snn", &reference.join("graphs").display().to_string(), "--out", &o.display().to_string(),
            "--seed", "12345", "--resolution", "0.2,0.4,0.8,1.2,2.0", "--threads", t, "--stats", "none", "--keep-stages",
        ]);
        assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
        assert!(!String::from_utf8_lossy(&r.stderr).contains("not computed"), "--stats none");
    }
    let dir = |t: &str| PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("cli/cluster_res_t{t}"));
    for f in ["clusters.tsv", "metadata.tsv", "optimizer_ids.tsv", "seurat_clusters.tsv", "seurat_clusters_levels.txt", "SCT_snn_res.2_levels.txt", "stages/SCT_snn_res.0.2_partitions.npy"] {
        assert!(std::fs::read(dir("1").join(f)).unwrap() == std::fs::read(dir("4").join(f)).unwrap(), "{f} depends on the thread count");
    }
    let result = |t: &str| {
        let s = scio::Json::parse(&std::fs::read_to_string(dir(t).join("summary.json")).unwrap()).unwrap();
        scio::json::to_file_text(s.get("result").unwrap())
    };
    assert_eq!(result("1"), result("4"));
    // scrust neighbors' output as the input: the same labels
    let nn = tmp("cluster_from_neighbors");
    let r = scrust(&["neighbors", "--pca", &reference.join("pca").display().to_string(), "--out", &nn.display().to_string(), "--dims", "1:30"]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let o = tmp("cluster_from_neighbors_out");
    let r = scrust(&["cluster", "--snn", &nn.display().to_string(), "--out", &o.display().to_string(), "--seed", "12345", "--stats", "none"]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    assert_eq!(std::fs::read(o.join("metadata.tsv")).unwrap(), std::fs::read(out.join("metadata.tsv")).unwrap());
    // refusals
    let snn = reference.join("graphs").display().to_string();
    let refused = |extra: &[&str], needle: &str| {
        let o = tmp("cluster_refused");
        let mut a = vec!["cluster", "--snn", snn.as_str(), "--out", o.to_str().unwrap()];
        a.extend_from_slice(extra);
        let r = scrust(&a);
        assert!(!r.status.success(), "{extra:?} accepted");
        let e = String::from_utf8_lossy(&r.stderr);
        assert!(e.contains(needle), "{extra:?}: {e}");
    };
    refused(&["--algorithm", "4"], "Leiden");
    refused(&["--counts", "x"], "give one starting point");
    refused(&["--export", "yes"], "milestone M13");
    refused(&["--resolution", "0.8,0.8"], "overwrite");
    refused(&["--modularity", "2", "--resolution", "1.5"], "resolution<1 for alternative modularity");
    refused(&["--starts", "0"], "Have to have at least one start");
    refused(&["--seed", "1.5"], "not an integer");
}

/// Every file of a directory tree, by relative path (summary.json, which records the run, left
/// out), for comparing two runs.
fn tree(dir: &Path, skip: &[&str]) -> std::collections::BTreeMap<String, Vec<u8>> {
    let mut out = std::collections::BTreeMap::new();
    let mut todo = vec![dir.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            let rel = p.strip_prefix(dir).unwrap().display().to_string();
            if skip.iter().any(|s| rel == *s || rel.starts_with(&format!("{s}/"))) {
                continue;
            }
            if p.is_dir() {
                todo.push(p);
            } else {
                out.insert(rel, std::fs::read(&p).unwrap());
            }
        }
    }
    out
}

/// The clustering chain (milestones M11 and M12) from the counts in one process, and from each
/// stage's output (`scrust sctransform`, `scrust pca`, `scrust neighbors`): the same files, byte
/// for byte, by every route (the purities among them), and the same clustering and scores in
/// `summary.json`. Rust against Rust under one kernel, so any kernel.
#[test]
fn cluster_routes_agree_byte_for_byte() {
    let counts = fixture("soupx_tenx_test").join("input");
    let c = counts.display().to_string();
    let run = |args: &[&str]| {
        let out = scrust(args);
        assert!(out.status.success(), "{args:?}: {}", String::from_utf8_lossy(&out.stderr));
        out
    };
    let one = tmp("chain_one");
    let o = run(&["cluster", "--counts", &c, "--preset", "workflow", "--out", one.to_str().unwrap(), "--keep-stages", "--annoy", "yes", "--threads", "2"]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("neighborPurity: radius") && !err.contains("asked for, but not computed"), "{err}");
    let (s, p, n) = (tmp("chain_sct"), tmp("chain_pca"), tmp("chain_nn"));
    run(&["sctransform", "--counts", &c, "--out", s.to_str().unwrap(), "--seed", "12345", "--threads", "2"]);
    run(&["pca", "--sct", s.to_str().unwrap(), "--out", p.to_str().unwrap(), "--seed", "12345"]);
    run(&["neighbors", "--pca", p.to_str().unwrap(), "--out", n.to_str().unwrap(), "--dims", "1:30", "--annoy", "yes"]);
    assert!(p.join("cell_metadata.tsv").is_file() && n.join("cell_metadata.tsv").is_file(), "the counts' columns are carried forward");
    let (from_sct, from_pca, from_snn) = (tmp("chain_from_sct"), tmp("chain_from_pca"), tmp("chain_from_snn"));
    run(&["cluster", "--sct", s.to_str().unwrap(), "--preset", "workflow", "--out", from_sct.to_str().unwrap(), "--keep-stages", "--annoy", "yes"]);
    run(&["cluster", "--pca", p.to_str().unwrap(), "--preset", "workflow", "--out", from_pca.to_str().unwrap(), "--annoy", "yes"]);
    run(&["cluster", "--snn", n.to_str().unwrap(), "--pca", p.to_str().unwrap(), "--preset", "workflow", "--out", from_snn.to_str().unwrap()]);
    let want = tree(&one, &["summary.json", "stages"]);
    assert!(want.contains_key("metadata.tsv") && want.contains_key("pca/embeddings.tsv") && want.contains_key("neighbors/snn.mtx.gz") && want.contains_key("sct/scale_data.npy"));
    let meta = String::from_utf8(want["metadata.tsv"].clone()).unwrap();
    assert!(meta.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tSCT_snn_res.0.8\tseurat_clusters\tsilhouette_width\tsilhouette_other\tneighborhood_purity\tpurity_maximum\n"), "{}", &meta[..160]);
    let (got_sct, got_pca, got_snn) = (tree(&from_sct, &["summary.json", "stages"]), tree(&from_pca, &["summary.json", "stages"]), tree(&from_snn, &["summary.json", "stages"]));
    assert_eq!(got_sct.keys().collect::<Vec<_>>(), want.keys().collect::<Vec<_>>());
    for (k, v) in &want {
        assert!(got_sct[k] == *v, "--sct: {k} differs");
        if !k.starts_with("sct/") {
            assert!(got_pca.get(k) == Some(v), "--pca: {k} differs or is missing");
            assert!(got_snn.get(k) == Some(v), "--snn: {k} differs or is missing");
        }
    }
    // the clustering and the scores in summary.json's result section, by every route
    let result = |d: &Path| scio::Json::parse(&std::fs::read_to_string(d.join("summary.json")).unwrap()).unwrap().get("result").cloned().unwrap();
    let r1 = result(&one);
    for d in [&from_sct, &from_pca, &from_snn] {
        let r = result(d);
        for key in ["n_cells", "nodes", "edges", "total_edge_weight", "resolutions", "seurat_clusters_levels", "silhouette", "purity"] {
            assert_eq!(r1.get(key), r.get(key), "{}: result.{key}", d.display());
        }
    }
    // the stages each route made itself are recorded too
    assert!(r1.get("rng_after_sha256").and_then(|x| x.get("RunPCA")).and_then(|x| x.as_str()).is_some());
    // the run record's phases, in order, on one clock
    for (d, want) in [(&one, "load sctransform pca neighbors cluster silhouette purity write"), (&from_sct, "load pca neighbors cluster silhouette purity write"), (&from_pca, "load neighbors cluster silhouette purity write"), (&from_snn, "load cluster silhouette purity write")] {
        let s = scio::Json::parse(&std::fs::read_to_string(d.join("summary.json")).unwrap()).unwrap();
        let Some(scio::Json::Object(phases)) = s.get("run").and_then(|r| r.get("phases")) else { panic!("{}: no run.phases", d.display()) };
        let names: Vec<&str> = phases.iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(names.join(" "), want, "{}", d.display());
        let t: Vec<f64> = phases.iter().map(|(_, v)| v.get("elapsed_s").and_then(|x| x.as_f64()).unwrap()).collect();
        assert!(t.windows(2).all(|w| w[0] <= w[1]), "{}: the phases go back in time: {t:?}", d.display());
    }
    assert_eq!(result(&from_pca).get("pca"), r1.get("pca"), "--pca copies scrust pca's result");
}

/// The chain's parameters (milestone M11): a TOML file under the preset, then flags; the values
/// recorded in summary.json; and the code paths that are refused.
#[test]
fn cluster_parameters_and_refusals() {
    let counts = fixture("decontx_raw_tenx_test").join("input").display().to_string();
    let dir = tmp("chain_params");
    std::fs::create_dir_all(&dir).unwrap();
    let toml = dir.join("params.toml");
    std::fs::write(&toml, "preset = \"workflow\"\n[neighbors]\nk = 15\nprune = \"1/15\"\n[cluster]\nresolution = [0.4, 0.8]\n[stats]\ncompute = \"silhouette\"\n").unwrap();
    let out = dir.join("out");
    let r = scrust(&["cluster", "--counts", &counts, "--params", toml.to_str().unwrap(), "--trees", "30", "--out", out.to_str().unwrap()]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    assert!(!String::from_utf8_lossy(&r.stderr).contains("not computed"), "purity not asked for");
    let s = scio::Json::parse(&std::fs::read_to_string(out.join("summary.json")).unwrap()).unwrap();
    let p = s.get("run").unwrap().get("parameters").unwrap();
    assert_eq!(p.get("preset").and_then(|v| v.as_str()), Some("workflow"));
    let nb = p.get("neighbors").unwrap();
    assert_eq!(nb.get("k_param").and_then(|v| v.as_f64()), Some(15.0));
    assert_eq!(nb.get("n_trees").and_then(|v| v.as_f64()), Some(30.0), "the flag after the file");
    assert_eq!(p.get("sctransform").unwrap().get("seed_use").and_then(|v| v.as_f64()), Some(12345.0));
    let clusters = std::fs::read_to_string(out.join("clusters.tsv")).unwrap();
    assert!(clusters.starts_with("barcode\tSCT_snn_res.0.4\tSCT_snn_res.0.8\n"));
    // refusals: other code paths, a stage dump without the counts, two starting points
    let refused = |extra: &[&str], needle: &str| {
        let o = tmp("chain_refused");
        let mut a = vec!["cluster", "--counts", counts.as_str(), "--out", o.to_str().unwrap()];
        a.extend_from_slice(extra);
        let r = scrust(&a);
        assert!(!r.status.success(), "{extra:?} accepted");
        let e = String::from_utf8_lossy(&r.stderr);
        assert!(e.contains(needle), "{extra:?}: {e}");
    };
    refused(&["--nn-method", "rann"], "only Annoy");
    refused(&["--l2-norm", "true"], "l2.norm");
    refused(&["--vars-to-regress", "percent.mt"], "vars.to.regress");
    refused(&["--algorithm", "4"], "Leiden");
    refused(&["--preset", "fast"], "unknown preset");
    refused(&["--pca", "x"], "give one starting point");
    std::fs::write(&toml, "[neighbors]\nl2_norm = true\n").unwrap();
    refused(&["--params", toml.to_str().unwrap()], "l2.norm");
    // the workflow's guard: 50 cells or fewer stop the run with its message
    let few = tmp("chain_guard_input");
    std::fs::create_dir_all(&few).unwrap();
    let m = scio::sparse::read_mtx(Path::new(&counts)).unwrap();
    let barcodes = scio::tenx::read_barcodes(Path::new(&counts)).unwrap();
    let features: Vec<String> = scio::tenx::read_lines(&scio::tenx::find(Path::new(&counts), "features.tsv").unwrap()).unwrap();
    scio::tenx::write_columns(&few, &m, &(0..50).collect::<Vec<_>>(), &barcodes[..50], &features, "scrust test", 1).unwrap();
    let o = tmp("chain_guard");
    let r = scrust(&["cluster", "--counts", few.to_str().unwrap(), "--preset", "workflow", "--out", o.to_str().unwrap()]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("[LOW_QUALITY_SAMPLE] only 50 cell(s) remain"), "{}", String::from_utf8_lossy(&r.stderr));
}

/// `scrust validate cluster` on the chain's run against a committed reference that records no
/// OpenBLAS kernel (the fixtures keep no fingerprint): tolerance mode, never exact, said so, and
/// every check identical, the purity certificate's flagged pairs reported. No fixture value is compared with a BLAS-dependent result unless the
/// kernel is the fixture's, so the run is made under the reference's kernel first.
#[test]
fn cluster_validate_without_a_recorded_kernel_is_not_exact() {
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/filtered_seurat_emptydrops_test/reference");
    rcompat::blas::require_fixture_kernel(&reference, "the clustering chain against a portable reference");
    let input = reference.parent().unwrap().join("input");
    let out = tmp("chain_validate");
    let r = scrust(&["cluster", "--counts", input.to_str().unwrap(), "--preset", "workflow", "--out", out.to_str().unwrap(), "--threads", "2"]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let report = out.join("validate.json");
    let v = scrust(&["validate", "cluster", "--reference", reference.to_str().unwrap(), "--candidate", out.to_str().unwrap(), "--json", report.to_str().unwrap()]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    let text = std::fs::read_to_string(&report).unwrap();
    assert!(text.contains("\"mode\": \"tolerance\"") && text.contains("the reference records no OpenBLAS kernel"));
    assert!(!text.contains("\"status\": \"within floor\"") && !text.contains("\"status\": \"differ\""), "everything identical");
    // an odd number of cells (267): the median cell and its 50th neighbour lie on the radius, and
    // the certificate's flags are reported, not failed (M12's decision), labelled as that pair
    assert!(text.contains("\"status\": \"flagged, reported\"") && text.contains("2 pair(s) flagged (2 the median pair, 0 other"), "the certificate's report");
    for check in ["seurat_clusters: labels", "SCT_snn_res.0.8: levels", "metadata.tsv: silhouette_width", "approxSilhouette width (R's embeddings and labels)", "pca/ embeddings, dims 1:30", "metadata.tsv: neighborhood_purity", "neighborPurity purity (R's embeddings and labels)", "certificate of the candidate's purity search"] {
        assert!(text.contains(check), "{check} not checked");
    }
}
