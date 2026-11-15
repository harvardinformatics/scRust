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
       scrust validate neighbors --reference <R run dir> --candidate <Rust output dir>
                                 [--json <file>] [--markdown <file>]
       scrust validate cluster --reference <R run dir> --candidate <Rust output dir>
                               [--optimizer <dir>] [--panel <tsv> --object <name>]
                               [--json <file>] [--markdown <file>]

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

neighbors: compares a `scrust neighbors --annoy yes` run with an R run of FindNeighbors on the
same embeddings (an instrumented clustering run, or a neighbour reference run in its layout:
stages/ and graphs/): the parameters R recorded; Annoy's items (K3_items.npy against the
candidate's annoy.ann) and the index file byte for byte (K3, K4); every cell's neighbour ids and
float distances (K5); SCT_nn and SCT_snn, every stored entry in order and every value bitwise
(K6). No BLAS is involved, so the comparison is exact on any machine.

cluster: compares a `scrust cluster --snn` run with an R run of FindClusters on the same SNN
graph (an instrumented clustering run of M8, a direct run, or a FindClusters reference run in
the same layout), resolution by resolution (K7): the settings; the optimizer's ids before
GroupSingletons (stages/K7_optimizer_ids.npy) and the final ids (K7_final_ids.npy); each
column's labels (metadata.tsv) and levels (<column>_levels.txt, byte for byte); every line
FindClusters printed but the elapsed time (nodes, edges, the maximum modularity as %.4f,
communities) and GroupSingletons' messages; the singletons in order; R's .Random.seed after
GroupSingletons (stages/rng.npy); seurat_clusters with its levels, idents_levels.txt,
cluster_ids.txt and nclusters.txt. With --optimizer (the instrumented optimizer's output: a
directory holding starts.tsv, or one per column or res.<r>), every start's iterations and
modularity as a double, the best start, ids.txt, and each start's partition when the candidate
ran with --keep-stages. With --panel (validation/clustering_panel.tsv) and --object, the
numbers the workflow job logged. No BLAS is involved, so the comparison is exact on any
machine.

Exit status 0 when nothing differs, 1 otherwise.
";

pub fn main(argv: &[String]) -> Result<bool, String> {
    let Some(tool) = argv.first() else {
        return Err("which tool? (doubletfinder, sctransform, pca, neighbors, cluster)".into());
    };
    if tool == "neighbors" {
        return neighbors_report(&argv[1..]);
    }
    if tool == "cluster" {
        return cluster_report(&argv[1..]);
    }
    if tool == "floor" {
        return floor_report(&argv[1..]);
    }
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
    if let Some((m, r)) = report.mode_text() {
        eprintln!("grading: {m} mode: {r}");
    }
    for (stage, ok, within, bad, n, nd) in report.by_stage() {
        if within > 0 {
            eprintln!("{stage:>7}: {ok} checks identical, {within} within the floor, {bad} differing; {nd} of {n} elements differ");
        } else {
            eprintln!("{stage:>7}: {ok} checks identical, {bad} differing; {nd} of {n} elements differ");
        }
    }
    use doubletfinder::validate::Status;
    for c in report.checks.iter().filter(|c| matches!(c.status, Status::Differ | Status::Missing)) {
        eprintln!("  {} {}: {} {}", c.stage, c.name, c.status.as_str(), c.detail);
    }
    let na = report.checks.iter().filter(|c| c.status == Status::NotApplicable).count();
    if na > 0 {
        eprintln!("  {na} checks not applicable (R's structural checks of its own containers; stage dumps of a direct run; the RNG digest under teacher forcing)");
    }
    let within = report.checks.iter().filter(|c| c.status == Status::WithinFloor).count();
    eprintln!(
        "{}",
        if !report.pass() {
            "FAIL".to_string()
        } else if within > 0 {
            format!("pass: {within} BLAS-dependent checks differ within the R-vs-R floor (tolerance mode), everything else is identical")
        } else {
            "pass: everything compared is identical".to_string()
        }
    );
    Ok(report.pass())
}

fn neighbors_report(argv: &[String]) -> Result<bool, String> {
    let args = Args::parse(argv)?;
    args.only(&["reference", "candidate", "json", "markdown"])?;
    let reference = args.path("reference")?;
    let candidate = args.path("candidate")?;
    let report = neighbors::validate::compare(&reference, &candidate)?;
    if let Some(p) = args.get("json").map(PathBuf::from) {
        std::fs::write(&p, scio::json::to_file_text(&report.to_json())).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    if let Some(p) = args.get("markdown").map(PathBuf::from) {
        std::fs::write(&p, report.to_markdown()).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    for c in &report.checks {
        eprintln!("{:>6} {:<36} {:<9} {} of {} differ; {}", c.stage, c.name, c.status.as_str(), c.n_differ, c.n, c.detail);
    }
    eprintln!("{}", if report.pass() { "pass: everything compared is identical" } else { "FAIL" });
    Ok(report.pass())
}

fn cluster_report(argv: &[String]) -> Result<bool, String> {
    let args = Args::parse(argv)?;
    args.only(&["reference", "candidate", "optimizer", "panel", "object", "saved", "saved-kernel", "json", "markdown"])?;
    let reference = args.path("reference")?;
    let candidate = args.path("candidate")?;
    let panel = match (args.get("panel"), args.get("object")) {
        (Some(p), Some(o)) => Some((PathBuf::from(p), o.to_string())),
        (None, None) => None,
        _ => return Err("--panel and --object go together".into()),
    };
    let opts = cluster::validate::CompareOptions {
        optimizer: args.get("optimizer").map(PathBuf::from),
        panel,
        saved: args.get("saved").map(PathBuf::from),
        saved_kernel: args.get("saved-kernel").map(String::from),
    };
    let report = cluster::validate::compare(&reference, &candidate, &opts)?;
    if let Some(p) = args.get("json").map(PathBuf::from) {
        std::fs::write(&p, scio::json::to_file_text(&report.to_json())).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    if let Some(p) = args.get("markdown").map(PathBuf::from) {
        std::fs::write(&p, report.to_markdown()).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    eprintln!("grading: {} mode: {}", report.mode.name(), report.mode.reason());
    if let Some(m) = &report.saved_mode {
        eprintln!("against the saved object: {} mode: {}", m.name(), m.reason());
    }
    for c in &report.checks {
        eprintln!("{:>11} {:<62} {:<14} {:<12} {} of {} differ; {}", c.stage, c.name, c.class.name(), c.status.as_str(), c.n_differ, c.n, c.detail);
    }
    use cluster::validate::Status as S;
    let within = report.checks.iter().filter(|c| c.status == S::WithinFloor && c.stage != "saved").count();
    let within_saved = report.checks.iter().filter(|c| c.status == S::WithinFloor && c.stage == "saved").count();
    let mut verdict = if !report.pass() {
        "FAIL".to_string()
    } else if within > 0 {
        format!("pass: {within} BLAS-dependent checks differ within the R-vs-R floor (tolerance mode), everything else compared with the R run is identical")
    } else {
        "pass: everything compared with the R run is identical".to_string()
    };
    if report.pass() && within_saved > 0 {
        verdict.push_str(&format!("; against the saved object (made under another kernel), {within_saved} checks within the floor"));
    }
    eprintln!("{verdict}");
    Ok(report.pass())
}

/// `scrust validate floor`: the R-vs-R floor between pairs of runs that differ only in their
/// OpenBLAS kernel, measured with the statistics tolerance mode grades with (scio::grade).
fn floor_report(argv: &[String]) -> Result<bool, String> {
    let args = Args::parse(argv)?;
    args.only(&["pairs", "out", "summary"])?;
    let pairs = scio::Table::read(&args.path("pairs")?).map_err(|e| format!("--pairs: {e}"))?;
    let (lc, ac, bc) = (pairs.col("label").ok_or("--pairs: no label column")?, pairs.col("a").ok_or("--pairs: no a column")?, pairs.col("b").ok_or("--pairs: no b column")?);
    let mut rows: Vec<Vec<String>> = Vec::new();
    // key -> (blocks, blocks differing, values differing, rows differing, max abs, max scaled, max rel, max fraction)
    let mut agg: Vec<(String, [f64; 8])> = Vec::new();
    for p in &pairs.rows {
        let (label, a, b) = (&p[lc], PathBuf::from(&p[ac]), PathBuf::from(&p[bc]));
        let m = scio::grade::measure_stages(&a, &b).map_err(|e| format!("{label}: {e}"))?;
        for r in m {
            let s = &r.stats;
            rows.push(vec![
                label.clone(),
                r.key.clone(),
                r.ctx.clone(),
                s.n.to_string(),
                s.n_differ.to_string(),
                s.rows_differ.to_string(),
                scio::fmt17(s.max_abs),
                scio::fmt17(s.max_scaled),
                scio::fmt17(s.max_rel),
                s.flipped.to_string(),
                s.incomparable.to_string(),
            ]);
            let i = match agg.iter().position(|x| x.0 == r.key) {
                Some(i) => i,
                None => {
                    agg.push((r.key.clone(), [0.0; 8]));
                    agg.len() - 1
                }
            };
            let v = &mut agg[i].1;
            v[0] += 1.0;
            v[1] += (s.n_differ > 0) as usize as f64;
            v[2] += s.n_differ as f64;
            v[3] += s.rows_differ as f64;
            v[4] = v[4].max(s.max_abs);
            v[5] = v[5].max(s.max_scaled);
            v[6] = v[6].max(s.max_rel);
            v[7] = v[7].max(s.fraction());
        }
    }
    let header = ["pair", "key", "ctx", "n", "n_differ", "rows_differ", "max_abs", "max_scaled", "max_rel", "flipped", "incomparable"];
    scio::Table { header: header.iter().map(|s| s.to_string()).collect(), rows }.write(&args.path("out")?).map_err(|e| e.to_string())?;
    let mut md = String::from("| key | blocks | blocks differing | values differing | rows differing | max abs | max scaled | max rel | max fraction differing |\n|---|---|---|---|---|---|---|---|---|\n");
    agg.sort_by(|x, y| x.0.cmp(&y.0));
    for (k, v) in &agg {
        md.push_str(&format!("| {k} | {} | {} | {} | {} | {:.3e} | {:.3e} | {:.3e} | {:.4} |\n", v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]));
    }
    match args.get("summary") {
        Some(p) => std::fs::write(p, &md).map_err(|e| format!("{p}: {e}"))?,
        None => print!("{md}"),
    }
    Ok(true)
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
                // the saved output's checks read the calls, which come after every BLAS product:
                // exact (the floor found S4-S7 identical across kernels)
                for _ in &s.checks {
                    r.grading.push((scio::grade::Class::BlasDependent, "exact (the floor: identical across kernels)".into()));
                }
                r.checks.extend(s.checks);
                r
            }
            None => s,
        });
    }
    let report = report.ok_or("--reference or --saved is required")?;
    Ok((report, args))
}
