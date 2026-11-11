//! The `scrust` binary: `--version`, `fingerprint` (environment fingerprint, M2) and the
//! `doubletfinder` placeholder.

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
fn doubletfinder_is_not_implemented() {
    let out = scrust(&["doubletfinder"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("not implemented"));
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
