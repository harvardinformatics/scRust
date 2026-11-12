//! `scrust validate <tool>`: compares a Rust run with an R reference run on plain files and
//! writes a JSON and a markdown report (docs/PLAN.md A4).

use crate::cli::{self, Args};
use std::path::PathBuf;

pub const USAGE: &str = "\
Usage: scrust validate doubletfinder --candidate <Rust output dir> [--reference <R run dir>]
                                     [--reference-stages <dir>] [--saved <dir>] [--json <file>]
                                     [--markdown <file>] [--threads N]
       scrust validate sctransform --input <10x dir> --stages <R stage dir>
                                   [--mode both|teacher|chained] [--contexts <ctx,ctx,...>]
                                   [--json <file>] [--markdown <file>] [--threads N]
       scrust validate pca --input <10x dir> --stages <R stage dir>
                           [--mode both|teacher|chained] [--contexts <ctx,ctx,...>]
                           [--json <file>] [--markdown <file>] [--threads N]

doubletfinder: compares a Rust run with an R reference run (--reference): cells.tsv,
sweep_stats.tsv and bcmvn.tsv byte for byte and the result section of summary.json field by
field; against an instrumented run (with stages/, or --reference-stages), also every block of
every stage array R dumped (S1-S7) bit for bit, every .Random.seed checkpoint by context and
label, every per-context scalar, and S6 and S7 recomputed from R's own S5. With --saved (an
export's saved/ directory), also the saved workflow output: the doublet set, the kept cells'
pANN as R wrote them, the column names (pK, nExp).

sctransform: runs Rust's SCTransform on the SCT input of every context of an instrumented
DoubletFinder run (the counts in --input, R's S1 draws from --stages) and compares stage S2
with R's dumps bit for bit, plus the RNG checkpoints and scalars: each sub-stage teacher-forced
on R's own inputs (S2a-S2d; --mode teacher), and the whole chain on Rust's own results
(S2-chain; --mode chained); both by default. S2d's products are OpenBLAS calls: run under the
reference run's OpenBLAS kernel (OPENBLAS_CORETYPE).

pca: runs Rust's RunPCA (npcs = 10) in every context of an instrumented DoubletFinder run and
compares stage S3 with R's dumps bit for bit (irlba's d, sdev, embeddings, loadings, the PCA
features, irlba's iterations and products, total.variance) and the RNG checkpoints around it:
teacher-forced on R's PCA input (or, where R did not dump it, Rust's S2d on R's S2c
parameters; --mode teacher), and after Rust's own SCTransform (S3-chain; --mode chained); both
by default. Run under the reference run's OpenBLAS kernel.

Exit status 0 when nothing differs, 1 otherwise.
";

pub fn main(argv: &[String]) -> Result<bool, String> {
    let Some(tool) = argv.first() else {
        return Err("which tool? (doubletfinder, sctransform, pca)".into());
    };
    let report = match tool.as_str() {
        "doubletfinder" => doubletfinder_report(&argv[1..])?,
        "sctransform" => sctransform_report(&argv[1..])?,
        "pca" => pca_report(&argv[1..])?,
        _ => return Err(format!("no validation for '{tool}' yet")),
    };
    let (report, args) = report;
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

fn sctransform_report(argv: &[String]) -> Result<(doubletfinder::validate::Report, Args), String> {
    let args = Args::parse(argv)?;
    args.only(&["input", "stages", "mode", "contexts", "json", "markdown", "threads"])?;
    let input = args.path("input")?;
    let stages = args.path("stages")?;
    let only: Option<Vec<String>> = args.get("contexts").map(|s| s.split(',').map(String::from).collect());
    use doubletfinder::sct::Which;
    let which = match args.get("mode").unwrap_or("both") {
        "both" => Which::ALL,
        "teacher" => Which::TEACHER,
        "chained" => Which::CHAIN,
        m => return Err(format!("--mode: expected both, teacher or chained, not {m}")),
    };
    let params = sctransform::SctParams::default();
    let threads = cli::threads(&args)?;
    let report = doubletfinder::sct::validate(&input, &stages, &params, which, only.as_deref(), threads).map_err(|e| e.to_string())?;
    Ok((report, args))
}

fn pca_report(argv: &[String]) -> Result<(doubletfinder::validate::Report, Args), String> {
    let args = Args::parse(argv)?;
    args.only(&["input", "stages", "mode", "contexts", "json", "markdown", "threads"])?;
    let input = args.path("input")?;
    let stages = args.path("stages")?;
    let only: Option<Vec<String>> = args.get("contexts").map(|s| s.split(',').map(String::from).collect());
    use doubletfinder::pca::Which;
    let which = match args.get("mode").unwrap_or("both") {
        "both" => Which::ALL,
        "teacher" => Which::TEACHER,
        "chained" => Which::CHAIN,
        m => return Err(format!("--mode: expected both, teacher or chained, not {m}")),
    };
    let threads = cli::threads(&args)?;
    let (report, sources) = doubletfinder::pca::validate(&input, &stages, which, only.as_deref(), threads).map_err(|e| e.to_string())?;
    for (ctx, source) in &sources {
        eprintln!("{ctx}: teacher-forced PCA input from {source}");
    }
    Ok((report, args))
}

fn doubletfinder_report(argv: &[String]) -> Result<(doubletfinder::validate::Report, Args), String> {
    let args = Args::parse(argv)?;
    args.only(&["reference", "reference-stages", "candidate", "saved", "json", "markdown", "threads"])?;
    let candidate = args.path("candidate")?;
    let threads = cli::threads(&args)?;
    let mut report = match args.get("reference").map(PathBuf::from) {
        Some(reference) => {
            let reference_stages = args.get("reference-stages").map(PathBuf::from).unwrap_or_else(|| reference.join("stages"));
            Some(doubletfinder::validate::compare(&reference, &reference_stages, &candidate, threads).map_err(|e| e.to_string())?)
        }
        None => None,
    };
    if let Some(saved) = args.get("saved").map(PathBuf::from) {
        let s = doubletfinder::validate::compare_saved(&saved, &candidate).map_err(|e| e.to_string())?;
        report = Some(match report {
            Some(mut r) => {
                r.checks.extend(s.checks);
                r
            }
            None => s,
        });
    }
    let report = report.ok_or("--reference or --saved is required")?;
    Ok((report, args))
}
