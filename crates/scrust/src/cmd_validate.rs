//! `scrust validate <tool>`: compares a Rust run with an R reference run on plain files and
//! writes a JSON and a markdown report (docs/PLAN.md A4).

use crate::cli::{self, Args};
use std::path::PathBuf;

pub const USAGE: &str = "\
Usage: scrust validate doubletfinder --reference <R run dir> --candidate <Rust output dir>
                                     [--reference-stages <dir>] [--json <file>]
                                     [--markdown <file>] [--threads N]

Compares cells.tsv, sweep_stats.tsv and bcmvn.tsv byte for byte, the result section of
summary.json field by field, every stage array block (S4-S7) the two runs both dumped bit for
bit, and S6 and S7 recomputed from R's own S5 (teacher forcing). Exit status 0 when nothing
differs, 1 otherwise.
";

pub fn main(argv: &[String]) -> Result<bool, String> {
    let Some(tool) = argv.first() else {
        return Err("which tool? (doubletfinder)".into());
    };
    if tool != "doubletfinder" {
        return Err(format!("no validation for '{tool}' yet"));
    }
    let args = Args::parse(&argv[1..])?;
    args.only(&["reference", "reference-stages", "candidate", "json", "markdown", "threads"])?;
    let reference = args.path("reference")?;
    let reference_stages = args.get("reference-stages").map(PathBuf::from).unwrap_or_else(|| reference.join("stages"));
    let candidate = args.path("candidate")?;
    let threads = cli::threads(&args)?;
    let report = doubletfinder::validate::compare(&reference, &reference_stages, &candidate, threads).map_err(|e| e.to_string())?;
    if let Some(p) = args.get("json").map(PathBuf::from) {
        std::fs::write(&p, scio::json::to_file_text(&report.to_json())).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    if let Some(p) = args.get("markdown").map(PathBuf::from) {
        std::fs::write(&p, report.to_markdown()).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    for (stage, ok, bad, n, nd) in report.by_stage() {
        eprintln!("{stage:>7}: {ok} checks identical, {bad} differing; {nd} of {n} elements differ");
    }
    for c in report.checks.iter().filter(|c| c.status != doubletfinder::validate::Status::Identical) {
        eprintln!("  {} {}: {} {}", c.stage, c.name, c.status.as_str(), c.detail);
    }
    eprintln!("{}", if report.pass() { "pass: everything compared is identical" } else { "FAIL" });
    Ok(report.pass())
}
