//! `scrust`, the single command-line entry point to scRust's tools (docs/PLAN.md A1).
//!
//! Each tool is a library crate under `crates/` plus a subcommand here, so the repository ships
//! as one binary. At milestone M0 only `--version` and `--help` work; `doubletfinder` is a
//! placeholder.

use std::process::ExitCode;

const USAGE: &str = "\
Usage: scrust <COMMAND> [ARGS]...

Commands:
  doubletfinder  Doublet detection reproducing DoubletFinder (not implemented yet)

Options:
  -h, --help     Print help
  -V, --version  Print version
";

fn main() -> ExitCode {
    let Some(command) = std::env::args_os().nth(1) else {
        eprint!("{USAGE}");
        return ExitCode::from(2);
    };
    match command.to_str() {
        Some("-V" | "--version") => {
            println!("scrust {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        Some("-h" | "--help" | "help") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("doubletfinder") => {
            eprintln!(
                "scrust doubletfinder: not implemented yet; it will reproduce {}. \
                 See docs/PLAN.md, Part B.",
                doubletfinder::REFERENCE
            );
            ExitCode::FAILURE
        }
        _ => {
            eprint!(
                "scrust: unrecognized command '{}'\n\n{USAGE}",
                command.to_string_lossy()
            );
            ExitCode::from(2)
        }
    }
}
