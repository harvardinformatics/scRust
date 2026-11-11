//! The `scrust` binary: `--version`, `fingerprint` (environment fingerprint, M2), and
//! `doubletfinder`, `validate` and `knn` (M3).

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

#[test]
fn doubletfinder_needs_stages_until_sctransform_and_pca_are_ported() {
    let f = fixture("soupx_tenx_test");
    let out = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--out", &tmp("no_stages").display().to_string(),
    ]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("--from-stages is required"));
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
