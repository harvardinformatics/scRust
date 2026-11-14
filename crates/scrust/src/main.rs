//! `scrust`, the single command-line entry point to scRust's tools (docs/PLAN.md A1).
//!
//! Each tool is a library crate under `crates/` plus a subcommand here, so the repository ships
//! as one binary: `doubletfinder` (milestone M6: the whole tool; M3: on R's own embeddings), `knn` (the shared
//! neighbour stage on its own), `sctransform` (SCTransform v2, milestone M4), `pca` (RunPCA,
//! milestone M5), `neighbors` (FindNeighbors: Annoy and the SNN graph, milestone M9),
//! `validate` (Rust against R, plan A4) and `fingerprint` (the environment fingerprint of plan
//! A5b, which the R reference runs record too).

mod cli;
mod cmd_doubletfinder;
mod cmd_knn;
mod cmd_neighbors;
mod cmd_pca;
mod cmd_sctransform;
mod cmd_validate;
mod fingerprint;

use std::process::ExitCode;

const USAGE: &str = "\
Usage: scrust <COMMAND> [ARGS]...

Commands:
  doubletfinder  Doublet detection reproducing DoubletFinder 2.0.6, from the counts and the
                 cluster labels (milestone M6; --from-stages: on an instrumented R run's
                 embeddings, M3)
  knn            Exact nearest neighbours as DoubletFinder ranks them
  sctransform    SCTransform v2 as Seurat 5.5.1 runs it (sctransform 0.4.3, glmGamPoi 1.22.0):
                 Pearson residuals, variable features and scale.data (milestone M4)
  pca            RunPCA as Seurat 5.5.1 runs it (irlba 2.4.1) on an SCTransform output
                 (milestone M5)
  neighbors      FindNeighbors as Seurat 5.5.1 runs it (RcppAnnoy 0.0.23's Annoy, the SNN
                 graph) on cell embeddings (milestone M9)
  validate       Compare a Rust run with an R reference run
                 (validate doubletfinder|sctransform|pca|neighbors)
  fingerprint    Print the environment fingerprint as JSON: node, CPU, OpenBLAS kernel
                 and threads, glibc (docs/PLAN.md A5b)

Run 'scrust <COMMAND> --help' for a command's options.

Options:
  -h, --help     Print help
  -V, --version  Print version
";

fn main() -> ExitCode {
    // R's reference runs use one OpenBLAS thread; the ports' products must too (rcompat::blas)
    rcompat::blas::pin_single_thread();
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
        Some("fingerprint") => {
            print!("{}", fingerprint::json());
            ExitCode::SUCCESS
        }
        Some(cmd @ ("doubletfinder" | "knn" | "neighbors" | "pca" | "sctransform" | "validate")) => {
            let rest: Vec<String> = std::env::args().skip(2).collect();
            let usage = match cmd {
                "doubletfinder" => cmd_doubletfinder::USAGE,
                "knn" => cmd_knn::USAGE,
                "neighbors" => cmd_neighbors::USAGE,
                "pca" => cmd_pca::USAGE,
                "sctransform" => cmd_sctransform::USAGE,
                _ => cmd_validate::USAGE,
            };
            if rest.iter().any(|a| a == "-h" || a == "--help") {
                print!("{usage}");
                return ExitCode::SUCCESS;
            }
            let result = match cmd {
                "doubletfinder" => cmd_doubletfinder::main(&rest).map(|_| true),
                "knn" => cmd_knn::main(&rest).map(|_| true),
                "neighbors" => cmd_neighbors::main(&rest).map(|_| true),
                "pca" => cmd_pca::main(&rest).map(|_| true),
                "sctransform" => cmd_sctransform::main(&rest).map(|_| true),
                _ => cmd_validate::main(&rest),
            };
            match result {
                Ok(true) => ExitCode::SUCCESS,
                Ok(false) => ExitCode::FAILURE,
                Err(e) => {
                    eprintln!("scrust {cmd}: {e}\n\n{usage}");
                    ExitCode::from(2)
                }
            }
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
