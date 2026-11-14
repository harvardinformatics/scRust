//! The `scrust` binary: `--version`, `fingerprint` (environment fingerprint, M2),
//! `doubletfinder`, `validate` and `knn` (M3), `sctransform` (M4), `pca` (M5), and the whole
//! `doubletfinder` from the counts (M6).

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
/// kept cells. The products are OpenBLAS calls and the fixtures are portable runs, so this needs
/// `OPENBLAS_CORETYPE=Haswell` (env.sh's default).
#[test]
fn doubletfinder_runs_every_stage_from_the_counts() {
    if std::env::var("OPENBLAS_CORETYPE").ok().as_deref() != Some("Haswell") {
        eprintln!("skipped: needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
        return;
    }
    let f = fixture("decontx_raw_emptydrops_test");
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
/// OpenBLAS calls, so it needs the reference's kernel (env.sh's OPENBLAS_CORETYPE=Haswell).
#[test]
fn sctransform_matches_r_byte_for_byte() {
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    if std::env::var("OPENBLAS_CORETYPE").ok().as_deref() != Some("Haswell") {
        eprintln!("skipped: needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
        return;
    }
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/sctransform/cli_soupx_tenx_test");
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
/// SCTransform output; large files by their SHA-256). Its products are OpenBLAS calls, so it needs
/// the reference's kernel (env.sh's OPENBLAS_CORETYPE=Haswell).
#[test]
fn pca_matches_r_byte_for_byte() {
    if std::env::var("OPENBLAS_CORETYPE").ok().as_deref() != Some("Haswell") {
        eprintln!("skipped: needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
        return;
    }
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/pca/cli_soupx_tenx_test");
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
