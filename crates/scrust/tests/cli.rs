//! The `scrust` binary at milestone M0: `--version` works and `doubletfinder` is a placeholder.

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
fn unknown_command_is_a_usage_error() {
    let out = scrust(&["no-such-tool"]);
    assert_eq!(out.status.code(), Some(2));
}
