//! Comparison of a `scrust cluster` run with an R run of FindClusters on the same SNN graph
//! (docs/PLAN.md C5, stage K7), on plain files. M11 extends it to the whole chain.
//!
//! The R side is an instrumented clustering run (`validation/R/instrumented_clustering.R`, M8),
//! a direct run (labels and printed numbers only), or a FindClusters reference run in the same
//! layout (`validation/R/cluster_references.R`, M10), optionally with the instrumented
//! optimizer's per-start results (`--optimizer`). For each resolution, in call order:
//!
//! - **the optimizer's ids** (`stages/K7_optimizer_ids.npy`, before GroupSingletons, numbered
//!   by size) against the candidate's `optimizer_ids.tsv`;
//! - **the final labels** (`stages/K7_final_ids.npy`, and the column in `metadata.tsv`) and the
//!   column's levels (`<column>_levels.txt`, byte for byte);
//! - **the printed numbers**: every line the optimizer printed but the elapsed time (nodes,
//!   edges, the maximum modularity as `%.4f`, communities), and GroupSingletons' message
//!   (singletons and final clusters); the singletons in order (`stages.json`);
//! - **every start** against the instrumented optimizer (`starts.tsv`: iterations, and the
//!   modularity as a double, bitwise; `ids.txt`; and each start's partition,
//!   `partitions.npy`, when the candidate kept them with `--keep-stages`);
//! - **R's `.Random.seed` after GroupSingletons** (`stages/rng.npy`, label
//!   `groupsingletons_exit`) where GroupSingletons drew; where it did not, that it did not;
//!
//! then `seurat_clusters` (the column, its levels, `idents_levels.txt`, `cluster_ids.txt`,
//! `nclusters.txt`) and, with `--panel`, the numbers the workflow job logged. Every criterion
//! is exact (`validation/tolerances.toml`, K7): no BLAS is involved.

use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

/// The outcome of one check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Identical,
    Differ,
    /// The reference or the candidate lacks the file.
    Missing,
    /// The reference does not hold what the check needs (a direct run has no stage dumps).
    NotApplicable,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Identical => "identical",
            Status::Differ => "differ",
            Status::Missing => "missing",
            Status::NotApplicable => "n/a",
        }
    }
}

/// One check.
#[derive(Debug, Clone)]
pub struct Check {
    pub stage: &'static str,
    pub name: String,
    pub status: Status,
    /// Elements compared and elements that differ.
    pub n: usize,
    pub n_differ: usize,
    pub detail: String,
}

/// All checks of one comparison.
#[derive(Debug, Clone)]
pub struct Report {
    pub reference: PathBuf,
    pub candidate: PathBuf,
    pub checks: Vec<Check>,
}

impl Report {
    /// True when every applicable check is identical (a missing file fails).
    pub fn pass(&self) -> bool {
        self.checks.iter().all(|c| matches!(c.status, Status::Identical | Status::NotApplicable))
    }

    pub fn to_json(&self) -> Json {
        Json::obj(vec![
            ("reference", Json::Str(self.reference.display().to_string())),
            ("candidate", Json::Str(self.candidate.display().to_string())),
            ("pass", Json::Bool(self.pass())),
            (
                "checks",
                Json::List(
                    self.checks
                        .iter()
                        .map(|c| {
                            Json::obj(vec![
                                ("stage", Json::Str(c.stage.into())),
                                ("name", Json::Str(c.name.clone())),
                                ("status", Json::Str(c.status.as_str().into())),
                                ("n", Json::Int(c.n as i64)),
                                ("n_differ", Json::Int(c.n_differ as i64)),
                                ("detail", Json::Str(c.detail.clone())),
                            ])
                        })
                        .collect(),
                ),
            ),
        ])
    }

    pub fn to_markdown(&self) -> String {
        let mut s = format!(
            "# scrust validate cluster\n\nreference: `{}`\ncandidate: `{}`\n\n**{}**\n\n| stage | check | status | compared | differing | detail |\n|---|---|---|---|---|---|\n",
            self.reference.display(),
            self.candidate.display(),
            if self.pass() { "pass: everything compared is identical" } else { "FAIL" }
        );
        for c in &self.checks {
            s.push_str(&format!("| {} | {} | {} | {} | {} | {} |\n", c.stage, c.name, c.status.as_str(), c.n, c.n_differ, c.detail.replace('|', "\\|")));
        }
        s
    }
}

fn check(stage: &'static str, name: &str, n: usize, n_differ: usize, detail: String) -> Check {
    let status = if n_differ == 0 { Status::Identical } else { Status::Differ };
    Check { stage, name: name.into(), status, n, n_differ, detail }
}

fn missing(stage: &'static str, name: &str, detail: String) -> Check {
    Check { stage, name: name.into(), status: Status::Missing, n: 0, n_differ: 0, detail }
}

fn not_applicable(stage: &'static str, name: &str, detail: String) -> Check {
    Check { stage, name: name.into(), status: Status::NotApplicable, n: 0, n_differ: 0, detail }
}

/// A file's bytes, from `path` or `path.gz` (the committed fixtures are gzipped), decompressed.
fn read_any(path: &Path) -> Result<Vec<u8>, String> {
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let p = if path.exists() || !gz.exists() { path } else { gz.as_path() };
    scio::tenx::read_maybe_gz(p).map_err(|e| e.to_string())
}

fn exists_any(path: &Path) -> bool {
    path.exists() || PathBuf::from(format!("{}.gz", path.display())).exists()
}

fn read_text(path: &Path) -> Result<String, String> {
    String::from_utf8(read_any(path)?).map_err(|_| format!("{}: not text", path.display()))
}

fn read_table(path: &Path) -> Result<Table, String> {
    Table::parse(&read_text(path)?).ok_or_else(|| format!("{}: empty", path.display()))
}

fn read_json(path: &Path) -> Result<Json, String> {
    Json::parse(&read_text(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_npy(path: &Path) -> Result<Npy, String> {
    Npy::parse(&read_any(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn str_list(j: Option<&Json>) -> Vec<String> {
    match j {
        Some(Json::Vector(v) | Json::List(v)) => v.iter().map(|x| x.as_str().map(String::from).unwrap_or_else(|| x.to_text(0))).collect(),
        Some(Json::Str(s)) => vec![s.clone()],
        _ => Vec::new(),
    }
}

fn num_list(j: Option<&Json>) -> Vec<f64> {
    match j {
        Some(Json::Vector(v) | Json::List(v)) => v.iter().filter_map(Json::as_f64).collect(),
        Some(x) => x.as_f64().into_iter().collect(),
        None => Vec::new(),
    }
}

fn jnum(j: &Json, key: &str) -> Option<f64> {
    j.get(key).and_then(Json::as_f64)
}

/// FindClusters' settings as an R run records them.
#[derive(Debug, Clone, PartialEq)]
struct Settings {
    modularity_fxn: i64,
    resolutions: Vec<f64>,
    algorithm: i64,
    n_start: i64,
    n_iter: i64,
    random_seed: i64,
    group_singletons: bool,
    source: String,
}

impl Settings {
    fn workflow(source: &str) -> Settings {
        Settings { modularity_fxn: 1, resolutions: vec![0.8], algorithm: 1, n_start: 10, n_iter: 10, random_seed: 12345, group_singletons: true, source: source.into() }
    }

    fn from_params(p: &Json, source: &str) -> Option<Settings> {
        Some(Settings {
            modularity_fxn: jnum(p, "modularity_fxn")? as i64,
            resolutions: num_list(p.get("resolution")),
            algorithm: jnum(p, "algorithm")? as i64,
            n_start: jnum(p, "n_start")? as i64,
            n_iter: jnum(p, "n_iter")? as i64,
            random_seed: jnum(p, "random_seed")? as i64,
            group_singletons: matches!(p.get("group_singletons"), Some(Json::Bool(true))),
            source: source.into(),
        })
    }

    fn describe(&self) -> Vec<(&'static str, String)> {
        vec![
            ("modularity.fxn", self.modularity_fxn.to_string()),
            ("resolution", self.resolutions.iter().map(|r| scio::fmt17(*r)).collect::<Vec<_>>().join(",")),
            ("algorithm", self.algorithm.to_string()),
            ("n.start", self.n_start.to_string()),
            ("n.iter", self.n_iter.to_string()),
            ("random.seed", self.random_seed.to_string()),
            ("group.singletons", self.group_singletons.to_string()),
        ]
    }
}

/// The R run's settings: `summary.json`'s `parameters` (M10's reference runs), else the
/// optimizer settings an instrumented run's hooks recorded, else the workflow's (a direct run of
/// the workflow's chain).
fn reference_settings(reference: &Path, summary: Option<&Json>, scalars: Option<&Json>) -> Settings {
    if let Some(p) = summary.and_then(|s| s.get("parameters")) {
        if let Some(s) = Settings::from_params(p, "summary.json parameters") {
            return s;
        }
    }
    if let Some(sc) = scalars {
        if let (Some(m), Some(a), Some(st), Some(it), Some(seed)) = (
            jnum(sc, "optimizer_modularity_function"),
            jnum(sc, "optimizer_algorithm"),
            jnum(sc, "optimizer_n_start"),
            jnum(sc, "optimizer_n_iter"),
            jnum(sc, "optimizer_random_seed"),
        ) {
            return Settings {
                modularity_fxn: m as i64,
                resolutions: num_list(sc.get("optimizer_resolution")),
                algorithm: a as i64,
                n_start: st as i64,
                n_iter: it as i64,
                random_seed: seed as i64,
                group_singletons: true,
                source: "stages.json (the hooks of an instrumented run)".into(),
            };
        }
    }
    let _ = reference;
    Settings::workflow("the workflow's call (a direct run records no settings)")
}

/// R's printed output split into one block per optimizer run (each starts with the banner).
fn printed_blocks(lines: &[String]) -> Vec<Vec<String>> {
    let mut out: Vec<Vec<String>> = Vec::new();
    for l in lines {
        if l.starts_with("Modularity Optimizer version") {
            out.push(Vec::new());
        }
        if let Some(b) = out.last_mut() {
            b.push(l.clone());
        }
    }
    out
}

/// The `n`-th row of R's RNG checkpoints labelled `label` (stages/rng.npy and rng.tsv).
fn rng_rows(reference: &Path, label: &str) -> Result<Vec<Vec<i32>>, String> {
    let tsv = read_table(&reference.join("stages/rng.tsv"))?;
    let a = read_npy(&reference.join("stages/rng.npy"))?;
    let lc = tsv.col("label").ok_or("rng.tsv: no label column")?;
    let rc = tsv.col("row").ok_or("rng.tsv: no row column")?;
    let mut out = Vec::new();
    for r in &tsv.rows {
        if r.get(lc).map(String::as_str) == Some(label) {
            let row: usize = r[rc].parse().map_err(|_| "rng.tsv: bad row".to_string())?;
            out.push((0..a.ncol()).map(|j| a.get_i32(row, j)).collect());
        }
    }
    Ok(out)
}

/// Options of the comparison.
#[derive(Debug, Clone, Default)]
pub struct CompareOptions {
    /// The instrumented optimizer's output: `starts.tsv`, `ids.txt` and `partitions.npy` in
    /// `<dir>/<column>/`, `<dir>/<res.r>/`, or (one resolution) `<dir>/`.
    pub optimizer: Option<PathBuf>,
    /// The panel table (`validation/clustering_panel.tsv`) and the object's name: the numbers
    /// the workflow job logged.
    pub panel: Option<(PathBuf, String)>,
}

/// Compares the candidate directory (`scrust cluster --out`) with the R run directory
/// `reference`.
pub fn compare(reference: &Path, candidate: &Path, opts: &CompareOptions) -> Result<Report, String> {
    let mut checks = Vec::new();
    let r_summary = if exists_any(&reference.join("summary.json")) { Some(read_json(&reference.join("summary.json"))?) } else { None };
    let has_stages = exists_any(&reference.join("stages/stages.json"));
    let r_scalars = if has_stages {
        read_json(&reference.join("stages/stages.json"))?.get("scalars").and_then(|s| s.get("chain")).cloned()
    } else {
        None
    };
    let c_summary = read_json(&candidate.join("summary.json"))?;
    let c_result = c_summary.get("result").cloned().ok_or("candidate summary.json: no result section")?;
    let c_params = c_summary.get("run").and_then(|r| r.get("parameters")).cloned().ok_or("candidate summary.json: no run.parameters")?;

    // settings
    let rs = reference_settings(reference, r_summary.as_ref(), r_scalars.as_ref());
    let cs = Settings::from_params(&c_params, "candidate").ok_or("candidate summary.json: incomplete parameters")?;
    let (ra, ca) = (rs.describe(), cs.describe());
    let bad: Vec<String> = ra.iter().zip(&ca).filter(|(a, b)| a.1 != b.1).map(|(a, b)| format!("{}: R {}, Rust {}", a.0, a.1, b.1)).collect();
    checks.push(check("params", "FindClusters' settings", ra.len(), bad.len(), format!("R's from {}{}", rs.source, if bad.is_empty() { String::new() } else { format!("; {}", bad.join("; ")) })));

    // the cells
    let r_meta = read_table(&reference.join("metadata.tsv"))?;
    let c_meta = read_table(&candidate.join("metadata.tsv"))?;
    let rb = r_meta.col("barcode").ok_or("R's metadata.tsv: no barcode column")?;
    let cb = c_meta.col("barcode").ok_or("candidate metadata.tsv: no barcode column")?;
    let n = r_meta.rows.len();
    let bc_differ = if c_meta.rows.len() != n { n.max(c_meta.rows.len()) } else { (0..n).filter(|&i| r_meta.rows[i].get(rb) != c_meta.rows[i].get(cb)).count() };
    checks.push(check("cells", "barcodes", n, bc_differ, format!("{n} cells in R's metadata.tsv, {} in Rust's", c_meta.rows.len())));

    let columns = str_list(c_result.get("columns"));
    let c_res = match c_result.get("resolutions") {
        Some(Json::List(v) | Json::Vector(v)) => v.clone(),
        _ => Vec::new(),
    };
    let r_blocks = printed_blocks(&str_list(r_summary.as_ref().and_then(|s| s.get("run")).and_then(|r| r.get("printed"))));
    let r_messages: Vec<String> = str_list(r_summary.as_ref().and_then(|s| s.get("run")).and_then(|r| r.get("messages"))).iter().map(|m| m.trim().to_string()).filter(|m| m.contains("singletons identified")).collect();
    let r_opt_ids = if has_stages && exists_any(&reference.join("stages/K7_optimizer_ids.npy")) { Some(read_npy(&reference.join("stages/K7_optimizer_ids.npy"))?) } else { None };
    let r_final_ids = if has_stages && exists_any(&reference.join("stages/K7_final_ids.npy")) { Some(read_npy(&reference.join("stages/K7_final_ids.npy"))?) } else { None };
    let c_opt = read_table(&candidate.join("optimizer_ids.tsv"))?;
    let gs_entry = if has_stages { rng_rows(reference, "groupsingletons_entry").ok() } else { None };
    let gs_exit = if has_stages { rng_rows(reference, "groupsingletons_exit").ok() } else { None };

    for (k, col) in columns.iter().enumerate() {
        let cr = c_res.get(k).cloned().unwrap_or(Json::Null);
        let res_name = cr.get("name").and_then(Json::as_str).unwrap_or("").to_string();

        // K7: the optimizer's ids
        let c_ids: Vec<Option<i64>> = match c_opt.col(col) {
            Some(j) => c_opt.rows.iter().map(|r| r.get(j).and_then(|v| v.parse().ok())).collect(),
            None => Vec::new(),
        };
        match &r_opt_ids {
            Some(a) if a.ncol() > k => {
                let differ = if c_ids.len() != a.nrow() { a.nrow().max(c_ids.len()) } else { (0..a.nrow()).filter(|&i| c_ids[i] != Some(a.get_i32(i, k) as i64)).count() };
                checks.push(check("K7", &format!("{col}: optimizer ids"), a.nrow(), differ, format!("K7_optimizer_ids.npy column {} against optimizer_ids.tsv", k + 1)));
            }
            Some(_) => checks.push(missing("K7", &format!("{col}: optimizer ids"), format!("K7_optimizer_ids.npy has no column {}", k + 1))),
            None => checks.push(not_applicable("K7", &format!("{col}: optimizer ids"), "no stage dumps in the reference".into())),
        }

        // K7: the final labels
        let rc = r_meta.col(col);
        let cc = c_meta.col(col);
        match (rc, cc) {
            (Some(rc), Some(cc)) if c_meta.rows.len() == n => {
                let mut differ = 0;
                let mut first = None;
                for i in 0..n {
                    if r_meta.rows[i].get(rc) != c_meta.rows[i].get(cc) {
                        differ += 1;
                        if first.is_none() {
                            first = Some(format!("first at cell {}: R {:?}, Rust {:?}", i + 1, r_meta.rows[i].get(rc), c_meta.rows[i].get(cc)));
                        }
                    }
                }
                checks.push(check("K7", &format!("{col}: labels"), n, differ, first.unwrap_or_else(|| "metadata.tsv".into())));
            }
            _ => checks.push(missing("K7", &format!("{col}: labels"), format!("column {col} missing from R's or Rust's metadata.tsv"))),
        }
        if let Some(a) = &r_final_ids {
            if a.ncol() > k && cc.is_some() && c_meta.rows.len() == n {
                let cc = cc.unwrap();
                let differ = (0..n)
                    .filter(|&i| {
                        let want = a.get_i32(i, k);
                        let got = c_meta.rows[i].get(cc).map(String::as_str);
                        let want_s = if want < 0 { "singleton".to_string() } else { want.to_string() };
                        got != Some(want_s.as_str())
                    })
                    .count();
                checks.push(check("K7", &format!("{col}: final ids"), n, differ, format!("K7_final_ids.npy column {}", k + 1)));
            }
        }
        // the column's levels
        let lv = format!("{col}_levels.txt");
        match (read_any(&reference.join(&lv)), read_any(&candidate.join(&lv))) {
            (Ok(a), Ok(b)) => {
                let (la, lb) = (String::from_utf8_lossy(&a).lines().count(), String::from_utf8_lossy(&b).lines().count());
                checks.push(check("K7", &format!("{col}: levels"), la.max(lb), (a != b) as usize, format!("{la} levels in R's, {lb} in Rust's{}", if a == b { String::new() } else { format!(": R {:?}, Rust {:?}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b)) })));
            }
            (Err(e), _) | (_, Err(e)) => checks.push(missing("K7", &format!("{col}: levels"), e)),
        }

        // the printed numbers
        let c_printed = str_list(cr.get("printed"));
        match r_blocks.get(k) {
            Some(rb) => {
                let ra: Vec<&String> = rb.iter().filter(|l| !l.starts_with("Elapsed time")).collect();
                let ca: Vec<&String> = c_printed.iter().filter(|l| !l.starts_with("Elapsed time")).collect();
                let m = ra.len().max(ca.len());
                let differ = (0..m).filter(|&i| ra.get(i) != ca.get(i)).count();
                let mut detail = rb.iter().filter(|l| l.starts_with("Number of") || l.contains("odularity:") || l.starts_with("Maximum")).cloned().collect::<Vec<_>>().join("; ");
                if differ > 0 {
                    detail = format!("R {:?}; Rust {:?}", ra, ca);
                }
                checks.push(check("K7", &format!("{col}: printed numbers"), m, differ, detail));
            }
            None => checks.push(not_applicable("K7", &format!("{col}: printed numbers"), "R's printed output not recorded (no summary.json)".into())),
        }

        // the singletons
        let c_sing = str_list(cr.get("singletons"));
        if let Some(sc) = &r_scalars {
            // M10's references record each call's singletons; M8's one call's
            let r_sing = match sc.get("singletons_by_call") {
                Some(Json::List(v) | Json::Vector(v)) => v.get(k).map(|x| str_list(Some(x))),
                _ if columns.len() == 1 && sc.get("singletons").is_some() => Some(str_list(sc.get("singletons"))),
                _ => None,
            };
            if let Some(r_sing) = r_sing {
                let differ = if r_sing == c_sing { 0 } else { r_sing.len().max(c_sing.len()).max(1) };
                let shown = |v: &Vec<String>| if v.len() > 12 { format!("{:?} ... ({} in all)", &v[..12], v.len()) } else { format!("{v:?}") };
                checks.push(check("K7", &format!("{col}: singletons in order"), r_sing.len().max(c_sing.len()), differ, format!("R {}, Rust {}", shown(&r_sing), shown(&c_sing))));
            }
        }
        // per-start results, against the instrumented optimizer
        if let Some(dir) = &opts.optimizer {
            let d = [dir.join(col), dir.join(&res_name), dir.clone()].into_iter().find(|d| d.join("starts.tsv").exists() && (d != dir || columns.len() == 1));
            match d {
                Some(d) => checks.extend(compare_starts(&d, col, &cr, &c_ids, candidate)?),
                None => checks.push(missing("K7", &format!("{col}: every start"), format!("{}: no starts.tsv for {col}", dir.display()))),
            }
        }
        // R's .Random.seed after GroupSingletons
        match (&gs_entry, &gs_exit) {
            (Some(en), Some(ex)) if en.len() > k && ex.len() > k => {
                let drew = en[k] != ex[k];
                let c_rng: Option<Vec<i64>> = match cr.get("rng_after_groupsingletons") {
                    Some(Json::Vector(v) | Json::List(v)) => Some(v.iter().filter_map(Json::as_f64).map(|x| x as i64).collect()),
                    _ => None,
                };
                let name = format!("{col}: .Random.seed after GroupSingletons");
                if drew {
                    match c_rng {
                        Some(v) => {
                            let differ = if v.len() != ex[k].len() { ex[k].len() } else { v.iter().zip(&ex[k]).filter(|(a, b)| **a != **b as i64).count() };
                            checks.push(check("K7", &name, ex[k].len(), differ, format!("R drew (position {} after set.seed(1))", ex[k][1])));
                        }
                        None => checks.push(check("K7", &name, ex[k].len(), ex[k].len(), "R drew, Rust did not".into())),
                    }
                } else {
                    checks.push(check("K7", &name, 1, c_rng.is_some() as usize, "no draw in R (no singletons, or group.singletons = FALSE)".into()));
                }
            }
            _ => checks.push(not_applicable("K7", &format!("{col}: .Random.seed after GroupSingletons"), "no RNG checkpoints in the reference".into())),
        }
    }
    // GroupSingletons' messages, all resolutions together
    if r_summary.is_some() {
        let c_msgs: Vec<String> = c_res.iter().filter_map(|r| r.get("message").and_then(Json::as_str).map(String::from)).collect();
        let differ = if c_msgs == r_messages { 0 } else { r_messages.len().max(c_msgs.len()).max(1) };
        checks.push(check("K7", "GroupSingletons' messages", r_messages.len().max(c_msgs.len()), differ, format!("R {:?}, Rust {:?}", r_messages, c_msgs)));
    }

    // seurat_clusters and the workflow's files
    match (r_meta.col("seurat_clusters"), c_meta.col("seurat_clusters")) {
        (Some(a), Some(b)) if c_meta.rows.len() == n => {
            let differ = (0..n).filter(|&i| r_meta.rows[i].get(a) != c_meta.rows[i].get(b)).count();
            checks.push(check("K7", "seurat_clusters: labels", n, differ, "metadata.tsv".into()));
        }
        _ => checks.push(missing("K7", "seurat_clusters: labels", "seurat_clusters missing from a metadata.tsv".into())),
    }
    for (rname, cname) in [
        ("seurat_clusters_levels.txt", "seurat_clusters_levels.txt"),
        ("idents_levels.txt", "seurat_clusters_levels.txt"),
        ("cluster_ids.txt", "cluster_ids.txt"),
        ("nclusters.txt", "nclusters.txt"),
    ] {
        let label = if rname == cname { rname.to_string() } else { format!("{rname} (R) = {cname} (Rust)") };
        match (read_any(&reference.join(rname)), read_any(&candidate.join(cname))) {
            (Ok(a), Ok(b)) => checks.push(check("K7", &label, 1, (a != b) as usize, if a == b { format!("{} lines", String::from_utf8_lossy(&a).lines().count()) } else { format!("R {:?}, Rust {:?}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b)) })),
            (Err(_), _) if rname == "idents_levels.txt" => {}
            (Err(e), _) | (_, Err(e)) => checks.push(missing("K7", &label, e)),
        }
    }
    // the seurat_clusters.tsv that scrust doubletfinder reads
    match (read_table(&candidate.join("seurat_clusters.tsv")), c_meta.col("seurat_clusters")) {
        (Ok(t), Some(b)) => {
            let differ = if t.rows.len() != c_meta.rows.len() { n } else { (0..t.rows.len()).filter(|&i| t.rows[i].get(1) != c_meta.rows[i].get(b) || t.rows[i].first() != c_meta.rows[i].get(cb)).count() };
            checks.push(check("K7", "seurat_clusters.tsv = metadata.tsv", t.rows.len(), differ, format!("header {:?}", t.header)));
        }
        (Err(e), _) => checks.push(missing("K7", "seurat_clusters.tsv", e)),
        _ => {}
    }

    // the workflow job's log
    if let Some((panel, object)) = &opts.panel {
        checks.push(compare_panel(panel, object, &c_res)?);
    }
    Ok(Report { reference: reference.to_path_buf(), candidate: candidate.to_path_buf(), checks })
}

/// The candidate's starts against the instrumented optimizer's output in `dir`.
fn compare_starts(dir: &Path, col: &str, cr: &Json, c_ids: &[Option<i64>], candidate: &Path) -> Result<Vec<Check>, String> {
    let mut out = Vec::new();
    let t = read_table(&dir.join("starts.tsv"))?;
    let (ci, cm) = (t.col("iterations").ok_or("starts.tsv: no iterations")?, t.col("modularity").ok_or("starts.tsv: no modularity")?);
    let want_it: Vec<i64> = t.rows.iter().map(|r| r[ci].parse().unwrap_or(-1)).collect();
    let want_q: Vec<f64> = t.rows.iter().map(|r| scio::parse17(&r[cm]).unwrap_or(f64::NAN)).collect();
    let got_it: Vec<i64> = num_list(cr.get("start_iterations")).iter().map(|&v| v as i64).collect();
    let got_q = num_list(cr.get("start_modularity"));
    let m = want_q.len().max(got_q.len());
    let mut differ = 0;
    let mut first = None;
    for s in 0..m {
        let (a, b) = (want_q.get(s), got_q.get(s));
        let (ia, ib) = (want_it.get(s), got_it.get(s));
        if a.map(|v| v.to_bits()) != b.map(|v| v.to_bits()) || ia != ib {
            differ += 1;
            if first.is_none() {
                first = Some(format!("first at start {}: optimizer {:?} after {:?} iterations, Rust {:?} after {:?}", s + 1, a, ia, b, ib));
            }
        }
    }
    let distinct = {
        let mut v: Vec<u64> = want_q.iter().map(|q| q.to_bits()).collect();
        v.sort_unstable();
        v.dedup();
        v.len()
    };
    let best = want_q.iter().enumerate().fold((0usize, f64::NEG_INFINITY), |acc, (i, &q)| if q > acc.1 { (i, q) } else { acc });
    let c_best = jnum(cr, "best_start").map(|v| v as usize);
    let best_ok = c_best == Some(best.0 + 1);
    let detail = first.unwrap_or_else(|| format!("{m} starts, {distinct} distinct modularities; best start {} ({})", best.0 + 1, scio::fmt17(best.1)));
    out.push(check("K7", &format!("{col}: every start (iterations, modularity)"), m, differ, detail));
    out.push(check("K7", &format!("{col}: best start"), 1, (!best_ok) as usize, format!("the first start at the maximum: optimizer {}, Rust {:?}", best.0 + 1, c_best)));
    // the optimizer's final ids
    if let Ok(text) = read_text(&dir.join("ids.txt")) {
        let want: Vec<Option<i64>> = text.lines().map(|l| l.trim().parse().ok()).collect();
        let differ = if want.len() != c_ids.len() { want.len().max(c_ids.len()) } else { want.iter().zip(c_ids).filter(|(a, b)| a != b).count() };
        out.push(check("K7", &format!("{col}: optimizer ids = instrumented optimizer's ids.txt"), want.len(), differ, String::new()));
    }
    // each start's partition, when the candidate kept them
    let cp = candidate.join("stages").join(format!("{col}_partitions.npy"));
    if cp.exists() && exists_any(&dir.join("partitions.npy")) {
        let a = read_npy(&dir.join("partitions.npy"))?;
        let b = read_npy(&cp)?;
        let differ = if (a.nrow(), a.ncol()) != (b.nrow(), b.ncol()) {
            a.nrow() * a.ncol()
        } else {
            let mut d = 0;
            for s in 0..a.nrow() {
                for c in 0..a.ncol() {
                    d += (a.get_i32(s, c) != b.get_i32(s, c)) as usize;
                }
            }
            d
        };
        out.push(check("K7", &format!("{col}: every start's partition"), a.nrow() * a.ncol(), differ, format!("{} starts x {} cells", a.nrow(), a.ncol())));
    }
    Ok(out)
}

/// The candidate's numbers against the workflow job's log (`clustering_panel.tsv`).
fn compare_panel(panel: &Path, object: &str, c_res: &[Json]) -> Result<Check, String> {
    let t = read_table(panel)?;
    let oc = t.col("object").ok_or("panel: no object column")?;
    let row = t.rows.iter().find(|r| r.get(oc).map(String::as_str) == Some(object)).ok_or_else(|| format!("{}: no row for {object}", panel.display()))?;
    let get = |name: &str| t.col(name).and_then(|c| row.get(c)).cloned().unwrap_or_default();
    let Some(cr) = c_res.first() else {
        return Ok(missing("log", "workflow log", "no resolution in the candidate".into()));
    };
    let int_of = |key: &str| jnum(cr, key).map(|v| (v as i64).to_string()).unwrap_or_default();
    // nodes and edges as the candidate printed them
    let printed = str_list(cr.get("printed"));
    let pnum = |prefix: &str| printed.iter().find_map(|l| l.strip_prefix(prefix)).unwrap_or("").to_string();
    let all: Vec<(&str, String, String)> = vec![
        ("nodes", get("log_nodes"), pnum("Number of nodes: ")),
        ("edges", get("log_edges"), pnum("Number of edges: ")),
        ("max modularity", get("log_max_modularity"), cr.get("max_modularity_printed").and_then(Json::as_str).unwrap_or("").to_string()),
        ("communities", get("log_communities"), int_of("communities")),
        ("singletons", get("log_singletons"), int_of("n_singletons")),
        ("final clusters", get("log_final_clusters"), int_of("final_clusters")),
    ];
    let bad: Vec<String> = all.iter().filter(|(_, a, b)| a != b).map(|(n, a, b)| format!("{n}: log {a}, Rust {b}")).collect();
    let detail = if bad.is_empty() { all.iter().map(|(n, a, _)| format!("{n} {a}")).collect::<Vec<_>>().join(", ") } else { bad.join("; ") };
    Ok(check("log", "the workflow job's logged numbers", all.len(), bad.len(), detail))
}
