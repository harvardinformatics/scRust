//! Comparison of a `scrust cluster` run with an R run of the workflow's chain (docs/PLAN.md C5:
//! K1-K9; milestones M10, M11 and M12), on plain files, with kernel-gated grading (A5b
//! mechanism 3, [`scio::grade`]).
//!
//! The R side is an instrumented clustering run (`validation/R/instrumented_clustering.R`, M8,
//! with its stage dumps), a direct run (its final files), or a FindClusters reference run in the
//! same layout (`validation/R/cluster_references.R`, M10). What is compared depends on what both
//! sides hold:
//!
//! - **K1, K2, K3, K5, K7, K8, K9: the stage dumps**, when both runs have `stages/`: every block
//!   of every array R dumped (S2a-S2d: SCTransform; S3: RunPCA; K3: Annoy's items; K5: the
//!   neighbours; K7: the optimizer's and the final ids; K8: approxSilhouette's widths, closest
//!   clusters, centroids and within-cluster variances; K9: neighborPurity's k-th distances,
//!   neighbour counts, neighbour lists (cell by cell, where R kept them), weights, sums, totals,
//!   purities and maxima), every `.Random.seed` checkpoint both recorded, and every scalar both
//!   recorded (the radius among them). For runs that kept `scale.data` only as a digest (level
//!   light), its digest (`result.scale_data_sha256`). Arrays of later milestones (K10: the
//!   corrected counts, M13) are listed as not compared;
//! - **K4: Annoy's index** (`stages/annoy.ann`), byte for byte;
//! - **K6: the graphs**, R's `graphs/` against the candidate's `neighbors/`: dimensions, every
//!   stored entry in order, every value;
//! - **K7: FindClusters**, resolution by resolution (M10): the settings, the optimizer's ids
//!   (`stages/K7_optimizer_ids.npy` against `optimizer_ids.tsv`), the final ids, the labels and
//!   levels, what FindClusters printed, GroupSingletons' singletons, messages and R's
//!   `.Random.seed` after it, every start against the instrumented optimizer (`--optimizer`),
//!   `seurat_clusters` and the workflow's files, the workflow job's logged numbers (`--panel`);
//! - **K8 and K9: the files**, against R's final files: `metadata.tsv`'s columns (`nCount_RNA`,
//!   `nFeature_RNA`, the labels, `silhouette_width`, `neighborhood_purity`, and the candidate's
//!   `silhouette_other` and `purity_maximum` against R's `K8_silhouette_other` and `K9_maximum`);
//!   `pca/` (embeddings, loadings and their features, standard deviations, total variance); the
//!   `result` section's numbers (cells, clusters and their sizes, `scale.data`'s shape and digest,
//!   the variable features, the optimizer's numbers, the digests of `.Random.seed` after
//!   SCTransform and RunPCA). Columns of later milestones (`nCount_SCT`, `nFeature_SCT`: M13) are
//!   left out;
//! - **K9: the certificate** of the candidate's purity search (`result.purity.certificate`): the
//!   pairs close enough to a threshold for KMKNN's floating-point bounds to treat them
//!   differently, each labelled by kind (the median pair, on the radius by construction, or
//!   other). Flagged pairs are reported, not failed ([`Status::Reported`]): the user's decision at
//!   M12 keeps the exhaustive result, and M17's saved purities check every flagged pair;
//! - **the saved object** (`--saved`): the labels and levels, and the silhouette widths and
//!   purities, graded under the saved object's kernel (`--saved-kernel`, SapphireRapids: M8 found
//!   it on every production job);
//! - **silhouette and purity, teacher-forced**: approxSilhouette and neighborPurity on R's own
//!   embeddings and labels, against R's `K8_*` and `K9_*` arrays, the radius, and the
//!   `silhouette_width` and `neighborhood_purity` columns, with the certificate of that search.
//!   BLAS-free, so exact on any machine.
//!
//! **Grading.** Each check is BLAS-free or BLAS-dependent ([`Class`]): BLAS-dependent when any of
//! its inputs passed through a BLAS product in the candidate's run. From the counts, that is
//! every check from SCTransform's fits on (S2b's Cox-Reid term and nlminb call BLAS), except
//! SCTransform's random-number states and the counts' own columns; from R's own embeddings or
//! SNN graph, none. When both runs recorded the same OpenBLAS kernel, every check is exact.
//! Otherwise ([`Mode::Tolerance`]) BLAS-dependent checks are graded against
//! `[clustering.floor.kernel]` (`validation/tolerances.toml`, embedded at build time), and the
//! report says so.

use crate::find_clusters::Factor;
use crate::purity::neighbor_purity;
use crate::silhouette::approx_silhouette;
use scio::grade::{self, Class, Floor, Mode, Outcome, Stats};
use scio::stages::{RngCheckpoints, StageDir};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

/// `validation/tolerances.toml`, as built: the floors tolerance mode grades against.
pub const TOLERANCES: &str = include_str!("../../../validation/tolerances.toml");

/// The floor of the clustering chain's checks.
pub fn floor() -> Floor {
    Floor::from_toml(TOLERANCES, "clustering.floor.kernel").expect("validation/tolerances.toml has [clustering.floor.kernel]")
}

/// The outcome of one check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Identical,
    /// It differs, within the R-vs-R floor (tolerance mode, BLAS-dependent checks only).
    WithinFloor,
    Differ,
    /// The reference or the candidate lacks what the check needs.
    Missing,
    /// The check does not apply (a direct run has no stage dumps; a later milestone's output).
    NotApplicable,
    /// Flagged and reported, not failed: the purity certificate's flagged pairs (M12's decision).
    Reported,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Identical => "identical",
            Status::WithinFloor => "within floor",
            Status::Differ => "differ",
            Status::Missing => "missing",
            Status::NotApplicable => "n/a",
            Status::Reported => "flagged, reported",
        }
    }
}

/// One check.
#[derive(Debug, Clone)]
pub struct Check {
    pub stage: &'static str,
    pub name: String,
    pub class: Class,
    /// `exact`, or the floor entry it was graded against.
    pub criterion: String,
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
    pub mode: Mode,
    /// The mode of the checks against the saved object, under its own kernel.
    pub saved_mode: Option<Mode>,
    pub checks: Vec<Check>,
}

impl Report {
    /// True when every check is identical, within the floor, not applicable, or reported (the
    /// purity certificate's flagged pairs); a missing file fails.
    pub fn pass(&self) -> bool {
        self.checks.iter().all(|c| matches!(c.status, Status::Identical | Status::WithinFloor | Status::NotApplicable | Status::Reported))
    }

    /// Checks per stage: (stage, identical, within floor, differing, missing, n/a, reported).
    pub fn by_stage(&self) -> Vec<(&'static str, [usize; 6])> {
        let mut out: Vec<(&'static str, [usize; 6])> = Vec::new();
        for c in &self.checks {
            let i = match out.iter().position(|o| o.0 == c.stage) {
                Some(i) => i,
                None => {
                    out.push((c.stage, [0; 6]));
                    out.len() - 1
                }
            };
            let k = match c.status {
                Status::Identical => 0,
                Status::WithinFloor => 1,
                Status::Differ => 2,
                Status::Missing => 3,
                Status::NotApplicable => 4,
                Status::Reported => 5,
            };
            out[i].1[k] += 1;
        }
        out
    }

    pub fn to_json(&self) -> Json {
        let mode = |m: &Mode| Json::obj(vec![("mode", Json::Str(m.name().into())), ("reason", Json::Str(m.reason()))]);
        Json::obj(vec![
            ("reference", Json::Str(self.reference.display().to_string())),
            ("candidate", Json::Str(self.candidate.display().to_string())),
            ("pass", Json::Bool(self.pass())),
            ("grading", mode(&self.mode)),
            ("saved_grading", self.saved_mode.as_ref().map(mode).unwrap_or(Json::Null)),
            (
                "by_stage",
                Json::List(
                    self.by_stage()
                        .into_iter()
                        .map(|(s, k)| {
                            Json::obj(vec![
                                ("stage", Json::Str(s.into())),
                                ("identical", Json::Int(k[0] as i64)),
                                ("within_floor", Json::Int(k[1] as i64)),
                                ("differing", Json::Int(k[2] as i64)),
                                ("missing", Json::Int(k[3] as i64)),
                                ("not_applicable", Json::Int(k[4] as i64)),
                                ("reported", Json::Int(k[5] as i64)),
                            ])
                        })
                        .collect(),
                ),
            ),
            (
                "checks",
                Json::List(
                    self.checks
                        .iter()
                        .map(|c| {
                            Json::obj(vec![
                                ("stage", Json::Str(c.stage.into())),
                                ("name", Json::Str(c.name.clone())),
                                ("class", Json::Str(c.class.name().into())),
                                ("criterion", Json::Str(c.criterion.clone())),
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
            "# scrust validate cluster\n\nreference: `{}`\ncandidate: `{}`\n\n**{}**\n\nGrading: **{} mode**: {}.\n",
            self.reference.display(),
            self.candidate.display(),
            if self.pass() { "pass" } else { "FAIL" },
            self.mode.name(),
            self.mode.reason()
        );
        if let Some(m) = &self.saved_mode {
            s.push_str(&format!("Against the saved object: **{} mode**: {}.\n", m.name(), m.reason()));
        }
        s.push_str("\n| stage | identical | within floor | differing | missing | n/a | flagged, reported |\n|---|---|---|---|---|---|---|\n");
        for (st, k) in self.by_stage() {
            s.push_str(&format!("| {st} | {} | {} | {} | {} | {} | {} |\n", k[0], k[1], k[2], k[3], k[4], k[5]));
        }
        s.push_str("\n| stage | check | class | criterion | status | compared | differing | detail |\n|---|---|---|---|---|---|---|---|\n");
        for c in &self.checks {
            s.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} | {} | {} |\n",
                c.stage,
                c.name,
                c.class.name(),
                c.criterion.replace('|', "\\|"),
                c.status.as_str(),
                c.n,
                c.n_differ,
                c.detail.replace('|', "\\|")
            ));
        }
        s
    }
}

/// A file's bytes, from `path` or `path.gz` (the committed fixtures are gzipped), decompressed.
fn read_any(path: &Path) -> Result<Vec<u8>, String> {
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let p = if path.exists() || !gz.exists() { path } else { gz.as_path() };
    scio::tenx::read_maybe_gz(p).map_err(|e| format!("{}: {e}", p.display()))
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

/// Builds checks, classified and graded for the comparison's mode.
struct Grader {
    mode: Mode,
    floor: Floor,
    /// Where the candidate's chain started: counts, sct, pca or snn.
    start: String,
    checks: Vec<Check>,
}

/// The stage a check is attached to, for its class: SCTransform's setup and random-number states
/// (no BLAS), SCTransform's fits on, RunPCA, or what follows the PCA.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Inputs {
    /// The counts alone, or R's own data (settings, barcodes, the counts' columns, S2a, the
    /// RNG states SCTransform draws).
    Free,
    /// SCTransform's fits and residuals (S2b-S2d).
    Sct,
    /// RunPCA (S3) and everything after it.
    Pca,
}

impl Grader {
    fn class(&self, inputs: Inputs) -> Class {
        let blas = match (inputs, self.start.as_str()) {
            (Inputs::Free, _) => false,
            (Inputs::Sct, "counts") => true,
            (Inputs::Pca, "counts" | "sct") => true,
            _ => false,
        };
        if blas {
            Class::BlasDependent
        } else {
            Class::BlasFree
        }
    }

    fn push(&mut self, stage: &'static str, name: &str, class: Class, criterion: String, status: Status, n: usize, n_differ: usize, detail: String) {
        self.checks.push(Check { stage, name: name.into(), class, criterion, status, n, n_differ, detail });
    }

    /// An exact check (labels, ids, text): identical or not, in any mode.
    fn exact(&mut self, stage: &'static str, name: &str, inputs: Inputs, n: usize, n_differ: usize, detail: String) {
        let class = self.class(inputs);
        let status = if n_differ == 0 { Status::Identical } else { Status::Differ };
        let criterion = if class == Class::BlasDependent && !self.mode.is_exact() { "exact (the floor: identical across kernels)".into() } else { "exact".into() };
        self.push(stage, name, class, criterion, status, n, n_differ, detail);
    }

    /// A numeric check, graded on its statistics against the floor entry `key`.
    fn numeric(&mut self, stage: &'static str, name: &str, inputs: Inputs, key: &str, s: &Stats, detail: &str) {
        let class = self.class(inputs);
        let outcome = grade::grade(s, class, &self.mode, &self.floor, key);
        let criterion = if class == Class::BlasDependent && !self.mode.is_exact() { self.floor.describe(key) } else { "exact".into() };
        let (status, why) = match outcome {
            Outcome::Identical => (Status::Identical, String::new()),
            Outcome::WithinFloor => (Status::WithinFloor, String::new()),
            Outcome::Differ(m) => (Status::Differ, format!("; {m}")),
        };
        let d = if detail.is_empty() { format!("{}{why}", s.describe()) } else { format!("{detail}: {}{why}", s.describe()) };
        self.push(stage, name, class, criterion, status, s.n, s.n_differ, d);
    }

    fn missing(&mut self, stage: &'static str, name: &str, inputs: Inputs, detail: String) {
        let class = self.class(inputs);
        self.push(stage, name, class, "exact".into(), Status::Missing, 0, 0, detail);
    }

    fn not_applicable(&mut self, stage: &'static str, name: &str, inputs: Inputs, detail: String) {
        let class = self.class(inputs);
        self.push(stage, name, class, "-".into(), Status::NotApplicable, 0, 0, detail);
    }
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
fn reference_settings(summary: Option<&Json>, scalars: Option<&Json>) -> Settings {
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

/// The rows of R's RNG checkpoints labelled `label` (stages/rng.npy and rng.tsv).
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
    /// The saved object's export (`data/clustering/<object>/saved`): its labels, levels and
    /// silhouette widths.
    pub saved: Option<PathBuf>,
    /// The OpenBLAS kernel the saved objects were made with (default SapphireRapids, M8).
    pub saved_kernel: Option<String>,
}

/// The kernel the saved objects were made with: every production job ran on Sapphire Rapids,
/// where OpenBLAS detected SapphireRapids (M8).
pub const SAVED_KERNEL: &str = "SapphireRapids";

/// Compares the candidate directory (`scrust cluster --out`) with the R run directory
/// `reference`.
pub fn compare(reference: &Path, candidate: &Path, opts: &CompareOptions) -> Result<Report, String> {
    let r_summary = if exists_any(&reference.join("summary.json")) { Some(read_json(&reference.join("summary.json"))?) } else { None };
    // R's stage dumps (an instrumented run), and its scalars (also in the committed fixtures, which
    // keep stages.json without the index)
    let has_stages = exists_any(&reference.join("stages/index.tsv"));
    let r_scalars = if exists_any(&reference.join("stages/stages.json")) { read_json(&reference.join("stages/stages.json"))?.get("scalars").and_then(|s| s.get("chain")).cloned() } else { None };
    let c_summary = read_json(&candidate.join("summary.json"))?;
    let c_result = c_summary.get("result").cloned().ok_or("candidate summary.json: no result section")?;
    let c_run = c_summary.get("run").cloned().ok_or("candidate summary.json: no run section")?;
    let c_params_all = c_run.get("parameters").cloned().ok_or("candidate summary.json: no run.parameters")?;
    let c_params = c_params_all.get("cluster").cloned().unwrap_or(c_params_all.clone());
    let start = c_run.get("start").and_then(Json::as_str).unwrap_or("snn").to_string();
    let mode = grade::mode(grade::recorded_kernel(r_summary.as_ref()).as_deref(), grade::recorded_kernel(Some(&c_summary)).as_deref());
    let mut g = Grader { mode, floor: floor(), start, checks: Vec::new() };

    // the environment: which mode, and why
    {
        let rk = grade::recorded_kernel(r_summary.as_ref()).unwrap_or_else(|| "none recorded".into());
        let ck = grade::recorded_kernel(Some(&c_summary)).unwrap_or_else(|| "none recorded".into());
        let detail = format!("R {rk}, Rust {ck}; {} mode: {}", g.mode.name(), g.mode.reason());
        g.push("environment", "OpenBLAS kernel", Class::BlasFree, "-".into(), Status::NotApplicable, 1, 0, detail);
    }

    // settings
    let rs = reference_settings(r_summary.as_ref(), r_scalars.as_ref());
    let cs = Settings::from_params(&c_params, "candidate").ok_or("candidate summary.json: incomplete FindClusters parameters")?;
    let (ra, ca) = (rs.describe(), cs.describe());
    let bad: Vec<String> = ra.iter().zip(&ca).filter(|(a, b)| a.1 != b.1).map(|(a, b)| format!("{}: R {}, Rust {}", a.0, a.1, b.1)).collect();
    g.exact("params", "FindClusters' settings", Inputs::Free, ra.len(), bad.len(), format!("R's from {}{}", rs.source, if bad.is_empty() { String::new() } else { format!("; {}", bad.join("; ")) }));

    // the cells
    let r_meta = read_table(&reference.join("metadata.tsv"))?;
    let c_meta = read_table(&candidate.join("metadata.tsv"))?;
    let rb = r_meta.col("barcode").ok_or("R's metadata.tsv: no barcode column")?;
    let cb = c_meta.col("barcode").ok_or("candidate metadata.tsv: no barcode column")?;
    let n = r_meta.rows.len();
    let same_n = c_meta.rows.len() == n;
    let bc_differ = if !same_n { n.max(c_meta.rows.len()) } else { (0..n).filter(|&i| r_meta.rows[i].get(rb) != c_meta.rows[i].get(cb)).count() };
    g.exact("cells", "barcodes", Inputs::Free, n, bc_differ, format!("{n} cells in R's metadata.tsv, {} in Rust's", c_meta.rows.len()));

    compare_k7(&mut g, reference, candidate, opts, &r_summary, &r_scalars, has_stages, &c_result, &r_meta, &c_meta)?;
    if has_stages {
        compare_stages(&mut g, reference, candidate, &r_summary, &c_result)?;
    } else {
        g.not_applicable("K1-K8", "stage dumps", Inputs::Free, "the reference has none (a direct run): its final files are compared".into());
    }
    compare_files(&mut g, reference, candidate, &r_summary, &c_result, &c_run, &r_meta, &c_meta)?;
    let saved_mode = match &opts.saved {
        Some(s) => Some(compare_saved(&mut g, s, opts.saved_kernel.as_deref().unwrap_or(SAVED_KERNEL), &c_summary, &c_meta, candidate)?),
        None => None,
    };
    teacher_forced_silhouette(&mut g, reference, &r_meta, r_scalars.as_ref())?;
    teacher_forced_purity(&mut g, reference, &r_meta, r_scalars.as_ref())?;
    Ok(Report { reference: reference.to_path_buf(), candidate: candidate.to_path_buf(), mode: g.mode, saved_mode, checks: g.checks })
}

/// K7 (M10): FindClusters, resolution by resolution, and the workflow's files.
#[allow(clippy::too_many_arguments)]
fn compare_k7(g: &mut Grader, reference: &Path, candidate: &Path, opts: &CompareOptions, r_summary: &Option<Json>, r_scalars: &Option<Json>, has_stages: bool, c_result: &Json, r_meta: &Table, c_meta: &Table) -> Result<(), String> {
    let n = r_meta.rows.len();
    let cb = c_meta.col("barcode").unwrap_or(0);
    let columns = str_list(c_result.get("columns"));
    let c_res = match c_result.get("resolutions") {
        Some(Json::List(v) | Json::Vector(v)) => v.clone(),
        _ => Vec::new(),
    };
    let r_blocks = printed_blocks(&str_list(r_summary.as_ref().and_then(|s| s.get("run")).and_then(|r| r.get("printed"))));
    let r_messages: Vec<String> = str_list(r_summary.as_ref().and_then(|s| s.get("run")).and_then(|r| r.get("messages"))).iter().map(|m| m.trim().to_string()).filter(|m| m.contains("singletons identified")).collect();
    let _ = has_stages;
    let r_opt_ids = if exists_any(&reference.join("stages/K7_optimizer_ids.npy")) { Some(read_npy(&reference.join("stages/K7_optimizer_ids.npy"))?) } else { None };
    let r_final_ids = if exists_any(&reference.join("stages/K7_final_ids.npy")) { Some(read_npy(&reference.join("stages/K7_final_ids.npy"))?) } else { None };
    let c_opt = read_table(&candidate.join("optimizer_ids.tsv"))?;
    let gs_entry = if exists_any(&reference.join("stages/rng.tsv")) { rng_rows(reference, "groupsingletons_entry").ok() } else { None };
    let gs_exit = if exists_any(&reference.join("stages/rng.tsv")) { rng_rows(reference, "groupsingletons_exit").ok() } else { None };
    let k7 = Inputs::Pca;

    for (k, col) in columns.iter().enumerate() {
        let cr = c_res.get(k).cloned().unwrap_or(Json::Null);
        let res_name = cr.get("name").and_then(Json::as_str).unwrap_or("").to_string();

        // the optimizer's ids
        let c_ids: Vec<Option<i64>> = match c_opt.col(col) {
            Some(j) => c_opt.rows.iter().map(|r| r.get(j).and_then(|v| v.parse().ok())).collect(),
            None => Vec::new(),
        };
        match &r_opt_ids {
            Some(a) if a.ncol() > k => {
                let differ = if c_ids.len() != a.nrow() { a.nrow().max(c_ids.len()) } else { (0..a.nrow()).filter(|&i| c_ids[i] != Some(a.get_i32(i, k) as i64)).count() };
                g.exact("K7", &format!("{col}: optimizer ids"), k7, a.nrow(), differ, format!("K7_optimizer_ids.npy column {} against optimizer_ids.tsv", k + 1));
            }
            Some(_) => g.missing("K7", &format!("{col}: optimizer ids"), k7, format!("K7_optimizer_ids.npy has no column {}", k + 1)),
            None => g.not_applicable("K7", &format!("{col}: optimizer ids"), k7, "no stage dumps in the reference".into()),
        }

        // the final labels
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
                g.exact("K7", &format!("{col}: labels"), k7, n, differ, first.unwrap_or_else(|| "metadata.tsv".into()));
            }
            _ => g.missing("K7", &format!("{col}: labels"), k7, format!("column {col} missing from R's or Rust's metadata.tsv")),
        }
        if let (Some(a), Some(cc)) = (&r_final_ids, cc) {
            if a.ncol() > k && c_meta.rows.len() == n {
                let differ = (0..n)
                    .filter(|&i| {
                        let want = a.get_i32(i, k);
                        let got = c_meta.rows[i].get(cc).map(String::as_str);
                        let want_s = if want < 0 { "singleton".to_string() } else { want.to_string() };
                        got != Some(want_s.as_str())
                    })
                    .count();
                g.exact("K7", &format!("{col}: final ids"), k7, n, differ, format!("K7_final_ids.npy column {}", k + 1));
            }
        }
        // the column's levels
        let lv = format!("{col}_levels.txt");
        match (read_any(&reference.join(&lv)), read_any(&candidate.join(&lv))) {
            (Ok(a), Ok(b)) => {
                let (la, lb) = (String::from_utf8_lossy(&a).lines().count(), String::from_utf8_lossy(&b).lines().count());
                let detail = format!("{la} levels in R's, {lb} in Rust's{}", if a == b { String::new() } else { format!(": R {:?}, Rust {:?}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b)) });
                g.exact("K7", &format!("{col}: levels"), k7, la.max(lb), (a != b) as usize, detail);
            }
            (Err(e), _) | (_, Err(e)) => g.missing("K7", &format!("{col}: levels"), k7, e),
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
                g.exact("K7", &format!("{col}: printed numbers"), k7, m, differ, detail);
            }
            None => g.not_applicable("K7", &format!("{col}: printed numbers"), k7, "R's printed output not recorded (no summary.json)".into()),
        }

        // the singletons
        let c_sing = str_list(cr.get("singletons"));
        if let Some(sc) = r_scalars {
            let r_sing = match sc.get("singletons_by_call") {
                Some(Json::List(v) | Json::Vector(v)) => v.get(k).map(|x| str_list(Some(x))),
                _ if columns.len() == 1 && sc.get("singletons").is_some() => Some(str_list(sc.get("singletons"))),
                _ => None,
            };
            if let Some(r_sing) = r_sing {
                let differ = if r_sing == c_sing { 0 } else { r_sing.len().max(c_sing.len()).max(1) };
                let shown = |v: &Vec<String>| if v.len() > 12 { format!("{:?} ... ({} in all)", &v[..12], v.len()) } else { format!("{v:?}") };
                g.exact("K7", &format!("{col}: singletons in order"), k7, r_sing.len().max(c_sing.len()), differ, format!("R {}, Rust {}", shown(&r_sing), shown(&c_sing)));
            }
        }
        // per-start results, against the instrumented optimizer
        if let Some(dir) = &opts.optimizer {
            let d = [dir.join(col), dir.join(&res_name), dir.clone()].into_iter().find(|d| d.join("starts.tsv").exists() && (d != dir || columns.len() == 1));
            match d {
                Some(d) => compare_starts(g, &d, col, &cr, &c_ids, candidate)?,
                None => g.missing("K7", &format!("{col}: every start"), k7, format!("{}: no starts.tsv for {col}", dir.display())),
            }
        }
        // R's .Random.seed after GroupSingletons. Whether R drew cannot be read from its entry
        // and exit states: a call that draws after an earlier call drew can end where it started
        // (set.seed(1) and as many uniforms either way). So: where Rust drew, R's state after the
        // call must be Rust's; where Rust did not, R's state must not have changed.
        match (&gs_entry, &gs_exit) {
            (Some(en), Some(ex)) if en.len() > k && ex.len() > k => {
                let c_rng: Option<Vec<i64>> = match cr.get("rng_after_groupsingletons") {
                    Some(Json::Vector(v) | Json::List(v)) => Some(v.iter().filter_map(Json::as_f64).map(|x| x as i64).collect()),
                    _ => None,
                };
                let name = format!("{col}: .Random.seed after GroupSingletons");
                match c_rng {
                    Some(v) => {
                        let differ = if v.len() != ex[k].len() { ex[k].len() } else { v.iter().zip(&ex[k]).filter(|(a, b)| **a != **b as i64).count() };
                        g.exact("K7", &name, k7, ex[k].len(), differ, format!("Rust drew; R's state after the call (position {} after set.seed(1))", ex[k][1]));
                    }
                    None => g.exact("K7", &name, k7, 1, (en[k] != ex[k]) as usize, "Rust drew nothing (no singletons, or group.singletons = FALSE); R's state unchanged by the call".into()),
                }
            }
            _ => g.not_applicable("K7", &format!("{col}: .Random.seed after GroupSingletons"), k7, "no RNG checkpoints in the reference".into()),
        }
    }
    // GroupSingletons' messages, all resolutions together
    if r_summary.is_some() {
        let c_msgs: Vec<String> = c_res.iter().filter_map(|r| r.get("message").and_then(Json::as_str).map(String::from)).collect();
        let differ = if c_msgs == r_messages { 0 } else { r_messages.len().max(c_msgs.len()).max(1) };
        g.exact("K7", "GroupSingletons' messages", k7, r_messages.len().max(c_msgs.len()), differ, format!("R {:?}, Rust {:?}", r_messages, c_msgs));
    }

    // seurat_clusters and the workflow's files
    match (r_meta.col("seurat_clusters"), c_meta.col("seurat_clusters")) {
        (Some(a), Some(b)) if c_meta.rows.len() == n => {
            let differ = (0..n).filter(|&i| r_meta.rows[i].get(a) != c_meta.rows[i].get(b)).count();
            g.exact("K7", "seurat_clusters: labels", k7, n, differ, "metadata.tsv".into());
        }
        _ => g.missing("K7", "seurat_clusters: labels", k7, "seurat_clusters missing from a metadata.tsv".into()),
    }
    for (rname, cname) in [
        ("seurat_clusters_levels.txt", "seurat_clusters_levels.txt"),
        ("idents_levels.txt", "seurat_clusters_levels.txt"),
        ("cluster_ids.txt", "cluster_ids.txt"),
        ("nclusters.txt", "nclusters.txt"),
    ] {
        let label = if rname == cname { rname.to_string() } else { format!("{rname} (R) = {cname} (Rust)") };
        match (read_any(&reference.join(rname)), read_any(&candidate.join(cname))) {
            (Ok(a), Ok(b)) => {
                let detail = if a == b { format!("{} lines", String::from_utf8_lossy(&a).lines().count()) } else { format!("R {:?}, Rust {:?}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b)) };
                g.exact("K7", &label, k7, 1, (a != b) as usize, detail);
            }
            (Err(_), _) if rname == "idents_levels.txt" => {}
            (Err(e), _) | (_, Err(e)) => g.missing("K7", &label, k7, e),
        }
    }
    // the seurat_clusters.tsv that scrust doubletfinder reads
    match (read_table(&candidate.join("seurat_clusters.tsv")), c_meta.col("seurat_clusters")) {
        (Ok(t), Some(b)) => {
            let differ = if t.rows.len() != c_meta.rows.len() { n } else { (0..t.rows.len()).filter(|&i| t.rows[i].get(1) != c_meta.rows[i].get(b) || t.rows[i].first() != c_meta.rows[i].get(cb)).count() };
            g.exact("K7", "seurat_clusters.tsv = metadata.tsv", k7, t.rows.len(), differ, format!("header {:?}", t.header));
        }
        (Err(e), _) => g.missing("K7", "seurat_clusters.tsv", k7, e),
        _ => {}
    }
    if let Some((panel, object)) = &opts.panel {
        compare_panel(g, panel, object, &c_res)?;
    }
    Ok(())
}

/// The candidate's starts against the instrumented optimizer's output in `dir`.
fn compare_starts(g: &mut Grader, dir: &Path, col: &str, cr: &Json, c_ids: &[Option<i64>], candidate: &Path) -> Result<(), String> {
    let k7 = Inputs::Pca;
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
    g.exact("K7", &format!("{col}: every start (iterations, modularity)"), k7, m, differ, detail);
    g.exact("K7", &format!("{col}: best start"), k7, 1, (!best_ok) as usize, format!("the first start at the maximum: optimizer {}, Rust {:?}", best.0 + 1, c_best));
    if let Ok(text) = read_text(&dir.join("ids.txt")) {
        let want: Vec<Option<i64>> = text.lines().map(|l| l.trim().parse().ok()).collect();
        let differ = if want.len() != c_ids.len() { want.len().max(c_ids.len()) } else { want.iter().zip(c_ids).filter(|(a, b)| a != b).count() };
        g.exact("K7", &format!("{col}: optimizer ids = instrumented optimizer's ids.txt"), k7, want.len(), differ, String::new());
    }
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
        g.exact("K7", &format!("{col}: every start's partition"), k7, a.nrow() * a.ncol(), differ, format!("{} starts x {} cells", a.nrow(), a.ncol()));
    }
    Ok(())
}

/// The candidate's numbers against the workflow job's log (`clustering_panel.tsv`).
fn compare_panel(g: &mut Grader, panel: &Path, object: &str, c_res: &[Json]) -> Result<(), String> {
    let t = read_table(panel)?;
    let oc = t.col("object").ok_or("panel: no object column")?;
    let row = t.rows.iter().find(|r| r.get(oc).map(String::as_str) == Some(object)).ok_or_else(|| format!("{}: no row for {object}", panel.display()))?;
    let get = |name: &str| t.col(name).and_then(|c| row.get(c)).cloned().unwrap_or_default();
    let Some(cr) = c_res.first() else {
        g.missing("log", "workflow log", Inputs::Pca, "no resolution in the candidate".into());
        return Ok(());
    };
    let int_of = |key: &str| jnum(cr, key).map(|v| (v as i64).to_string()).unwrap_or_default();
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
    g.exact("log", "the workflow job's logged numbers", Inputs::Pca, all.len(), bad.len(), detail);
    Ok(())
}

/// The stage, inputs, floor key and sign alignment of a stage array R dumps.
fn array_kind(array: &str) -> Option<(&'static str, Inputs, String, bool)> {
    let stage = match array.split('_').next().unwrap_or("") {
        "S2a" => "K1",
        "S2b" | "S2c" => "K1",
        "S2d" if array == "S2d_pca_features" => "K2",
        "S2d" => "K1",
        "S3" => "K2",
        "K3" => "K3",
        "K5" => "K5",
        "K7" => "K7",
        "K8" => "K8",
        "K9" => "K9",
        _ => return None,
    };
    let inputs = if array.starts_with("S2a") { Inputs::Free } else if stage == "K1" { Inputs::Sct } else { Inputs::Pca };
    let key = array.to_string();
    let align = matches!(array, "S3_embeddings" | "S3_loadings" | "K3_items" | "K8_centroids");
    Some((stage, inputs, key, align))
}

/// The later milestones' arrays, and what supplies them.
fn later_milestone(array: &str) -> Option<&'static str> {
    if array.starts_with("K10") {
        Some("the corrected counts, milestone M13")
    } else {
        None
    }
}

/// K1-K9 from the stage dumps: every block of every array R dumped against the candidate's,
/// every RNG checkpoint and scalar both recorded, Annoy's index, and (light references) the
/// digest of scale.data.
fn compare_stages(g: &mut Grader, reference: &Path, candidate: &Path, r_summary: &Option<Json>, c_result: &Json) -> Result<(), String> {
    let rs = StageDir::open(&reference.join("stages")).map_err(|e| e.to_string())?;
    let csd = candidate.join("stages");
    if !csd.join("index.tsv").exists() {
        g.missing("K1-K8", "stage dumps", Inputs::Free, "the candidate has none (run scrust cluster --counts with --dump-stages)".into());
        return Ok(());
    }
    let cs = StageDir::open(&csd).map_err(|e| e.to_string())?;
    let mut files: Vec<String> = rs.blocks.iter().map(|b| b.file.clone()).collect();
    files.sort();
    files.dedup();
    let mut later: Vec<String> = Vec::new();
    for file in &files {
        for rb in rs.blocks_of(file) {
            let array = file.trim_end_matches(".npy").trim_end_matches(&format!(".{}", rb.ctx)).to_string();
            if let Some(m) = later_milestone(&array) {
                later.push(format!("{array} ({m})"));
                continue;
            }
            if array == "K9_neighbors_ptr" || array == "K9_neighbors_idx" {
                // compared cell by cell, below
                continue;
            }
            let Some((stage, inputs, key, align)) = array_kind(&array) else {
                g.not_applicable("stages", &array, Inputs::Free, "not an array of the chain's stages".into());
                continue;
            };
            let name = format!("{}: {array}", rb.ctx);
            if !rs.has(file, &rb.ctx) {
                g.not_applicable(stage, &name, inputs, "listed in R's index, not held (committed fixtures leave large arrays out)".into());
                continue;
            }
            let Some(cb) = cs.block(file, &rb.ctx).cloned() else {
                g.missing(stage, &name, inputs, "not in the candidate's dumps".into());
                continue;
            };
            if cb.dtype != rb.dtype || (cb.nrow, cb.ncol) != (rb.nrow, rb.ncol) {
                g.exact(stage, &name, inputs, rb.nrow * rb.ncol, rb.nrow * rb.ncol, format!("dtype {} and shape {}x{} against R's {} and {}x{}", cb.dtype, cb.nrow, cb.ncol, rb.dtype, rb.nrow, rb.ncol));
                continue;
            }
            let ra = rs.read(file).map_err(|e| e.to_string())?;
            let ca = cs.read(file).map_err(|e| e.to_string())?;
            if rb.dtype.contains('i') {
                if g.floor.get(&key).is_some() {
                    // integers the floor bounds (K9's neighbour counts): graded as numbers
                    let a: Vec<f64> = (0..rb.ncol).flat_map(|c| (0..rb.nrow).map(move |r| (c, r))).map(|(c, r)| ra.get_i32(rb.row_start + r, c) as f64).collect();
                    let b: Vec<f64> = (0..cb.ncol).flat_map(|c| (0..cb.nrow).map(move |r| (c, r))).map(|(c, r)| ca.get_i32(cb.row_start + r, c) as f64).collect();
                    let s = grade::stats(&a, &b, rb.nrow, rb.ncol, false);
                    g.numeric(stage, &name, inputs, &key, &s, &format!("{} x {}", rb.nrow, rb.ncol));
                    continue;
                }
                let mut nd = 0usize;
                for c in 0..rb.ncol {
                    for r in 0..rb.nrow {
                        nd += (ra.get_i32(rb.row_start + r, c) != ca.get_i32(cb.row_start + r, c)) as usize;
                    }
                }
                g.exact(stage, &name, inputs, rb.nrow * rb.ncol, nd, format!("{} x {}", rb.nrow, rb.ncol));
                continue;
            }
            let a: Vec<f64> = (0..rb.ncol).flat_map(|c| ra.col_f64(c, rb.row_start, rb.nrow)).collect();
            let b: Vec<f64> = (0..cb.ncol).flat_map(|c| ca.col_f64(c, cb.row_start, cb.nrow)).collect();
            if array == "S3_embeddings" && rb.ncol > 30 {
                // the floor separates the 30 dimensions the later stages read from the rest
                let split = 30 * rb.nrow;
                let s1 = grade::stats(&a[..split], &b[..split], rb.nrow, 30, true);
                let s2 = grade::stats(&a[split..], &b[split..], rb.nrow, rb.ncol - 30, true);
                g.numeric(stage, &format!("{name}, dims 1:30"), inputs, "S3_embeddings_dims_1_30", &s1, "");
                g.numeric(stage, &format!("{name}, dims 31:{}", rb.ncol), inputs, "S3_embeddings_dims_31_50", &s2, "");
                continue;
            }
            let fkey = if array == "S3_embeddings" { "S3_embeddings_dims_1_30".to_string() } else { key };
            let s = grade::stats(&a, &b, rb.nrow, rb.ncol, align);
            g.numeric(stage, &name, inputs, &fkey, &s, &format!("{} x {}", rb.nrow, rb.ncol));
        }
    }
    if !later.is_empty() {
        later.dedup();
        g.not_applicable("later", "arrays of later milestones", Inputs::Free, format!("not compared: {}", later.join(", ")));
    }
    // K9: purity's neighbour lists, cell by cell (where R kept them)
    if rs.has("K9_neighbors_ptr.npy", "chain") && rs.has("K9_neighbors_idx.npy", "chain") {
        let r_lists = read_lists(&rs)?;
        let c_lists = if cs.has("K9_neighbors_ptr.npy", "chain") && cs.has("K9_neighbors_idx.npy", "chain") { Some(read_lists(&cs)?) } else { None };
        match c_lists {
            Some(c) => compare_lists(g, "K9", "chain: K9_neighbors (cell by cell)", Inputs::Pca, &r_lists, &c),
            None => g.missing("K9", "chain: K9_neighbors (cell by cell)", Inputs::Pca, "not in the candidate's dumps (scrust cluster --purity-lists yes)".into()),
        }
    }
    // scale.data as a digest where R kept no copy (level light)
    if !rs.has("S2d_scale_data.chain.npy", "chain") {
        let rd = r_summary.as_ref().and_then(|s| s.get("result")).and_then(|r| r.get("scale_data_sha256")).and_then(Json::as_str);
        let cd = c_result.get("scale_data_sha256").and_then(Json::as_str);
        match (rd, cd) {
            (Some(a), Some(b)) => {
                let class = g.class(Inputs::Sct);
                if a == b {
                    g.exact("K1", "scale.data (digest; R kept no copy)", Inputs::Sct, 1, 0, format!("sha256 {a}"));
                } else if class == Class::BlasDependent && !g.mode.is_exact() {
                    g.push("K1", "scale.data (digest; R kept no copy)", class, g.floor.describe("S2d_scale_data"), Status::NotApplicable, 1, 1, "the digests differ; a digest cannot be graded against the floor's bounds (compare a run that kept scale.data)".into());
                } else {
                    g.exact("K1", "scale.data (digest; R kept no copy)", Inputs::Sct, 1, 1, format!("R {a}, Rust {b}"));
                }
            }
            _ => g.missing("K1", "scale.data (digest; R kept no copy)", Inputs::Sct, "no scale_data_sha256 in one of the summaries".into()),
        }
    }
    // the RNG checkpoints
    let rr = RngCheckpoints::open(&rs).map_err(|e| e.to_string())?;
    let cr = RngCheckpoints::open(&cs).map_err(|e| e.to_string())?;
    let mut na: Vec<String> = Vec::new();
    for (ctx, stage, label, state) in rr.all() {
        let (st, inputs, why_na): (&'static str, Inputs, Option<&str>) = match label.as_str() {
            "vst_entry" | "after_cell_sampling" | "fits_entry" | "fits_exit" | "vst_exit" | "sctransform_exit" => ("K1", Inputs::Free, None),
            l if l.starts_with("bin") => ("K1", Inputs::Free, None),
            "pca_entry" | "irlba_entry" => ("K2", Inputs::Free, None),
            "irlba_exit" | "pca_exit" => ("K2", Inputs::Pca, None),
            "groupsingletons_exit" => ("K7", Inputs::Pca, Some("compared under K7")),
            "correct_entry" => ("later", Inputs::Free, Some("correct(): milestone M13")),
            "chain_exit" => ("K9", Inputs::Free, Some("after neighborPurity, which draws nothing from R's stream (Rcpp's rng = false; KMKNN has its own generator): R's state after GroupSingletons where it drew, compared under K7, else after RunUMAP")),
            _ => ("K2", Inputs::Free, Some("after RunUMAP, which draws from R's stream and is not ported; nothing later reads these states (GroupSingletons reseeds)")),
        };
        if let Some(w) = why_na {
            na.push(format!("{label} ({w})"));
            continue;
        }
        let _ = stage;
        match cr.state(&ctx, &label) {
            Some(c) => g.exact(st, &format!("{ctx}: RNG {label}"), inputs, 1, (c != state) as usize, String::new()),
            None => g.missing(st, &format!("{ctx}: RNG {label}"), inputs, "not among the candidate's checkpoints".into()),
        }
    }
    if !na.is_empty() {
        g.not_applicable("RNG", "checkpoints not compared here", Inputs::Free, na.join("; "));
    }
    // the scalars both recorded
    let read_sc = |d: &StageDir| -> Option<Json> { Json::parse(&String::from_utf8(read_any(&d.dir.join("stages.json")).ok()?).ok()?).ok()?.get("scalars")?.get("chain").cloned() };
    if let (Some(Json::Object(r)), Some(c)) = (read_sc(&rs), read_sc(&cs)) {
        let mut r_only: Vec<String> = Vec::new();
        for (key, rv) in &r {
            let Some(cv) = c.get(key) else {
                r_only.push(key.clone());
                continue;
            };
            let inputs = match key.as_str() {
                "vst_n_cells" | "n_genes" | "n_genes_step1" | "n_cells_step1" | "n_genes_step1_pool" | "gene_sampling" | "n_cells" | "density_bw" | "density_n" | "min_variance" | "res_clip_range" | "clip_range" | "n_bins_fit" => Inputs::Free,
                "n_theta_below_1e7" | "n_theta_inf" | "reg_bw" | "mean_cell_sum" | "n_outliers" | "n_all_poisson_genes" | "n_smoothing_genes" => Inputs::Sct,
                _ => Inputs::Pca,
            };
            let st: &'static str = match inputs {
                Inputs::Free | Inputs::Sct => "K1",
                Inputs::Pca if key.starts_with("pca") || key.starts_with("irlba") || key == "total_variance" => "K2",
                Inputs::Pca if key.starts_with("annoy") => "K4",
                Inputs::Pca if key.starts_with("knn") => "K5",
                Inputs::Pca if key.starts_with("optimizer") || key.contains("singleton") => "K7",
                Inputs::Pca if key.starts_with("purity") => "K9",
                _ => "K8",
            };
            // numbers, graded (the floor lists irlba's iterations and products); the rest as text
            if let (Some(x), Some(y)) = (atoms_num(rv), atoms_num(cv)) {
                if x.len() == y.len() {
                    let s = grade::stats(&x, &y, x.len(), 1, false);
                    g.numeric(st, &format!("chain: {key}"), inputs, key, &s, &atoms_text(rv));
                    continue;
                }
            }
            // as atoms: R writes a vector of length one as a scalar
            let (a, b) = (atoms_text(rv), atoms_text(cv));
            g.exact(st, &format!("chain: {key}"), inputs, 1, (a != b) as usize, if a == b { a } else { format!("R {a}, Rust {b}") });
        }
        if !r_only.is_empty() {
            g.not_applicable("scalars", "R's own scalars", Inputs::Free, format!("recorded by R only (its structural checks of its own containers, RunUMAP's, later milestones'): {}", r_only.join(", ")));
        }
    }
    // K4: Annoy's index, byte for byte
    match (read_any(&reference.join("stages/annoy.ann")), read_any(&candidate.join("stages/annoy.ann"))) {
        (Ok(a), Ok(b)) => {
            let nd = if a.len() != b.len() { a.len().max(b.len()) } else { a.iter().zip(&b).filter(|(x, y)| x != y).count() };
            let s = Stats { n: a.len(), n_differ: nd, rows_differ: nd.min(1), ..Stats::default() };
            g.numeric("K4", "annoy.ann", Inputs::Pca, "K4_annoy_index", &s, &format!("{} bytes in R's, {} in Rust's", a.len(), b.len()));
        }
        (Err(e), _) | (_, Err(e)) => g.missing("K4", "annoy.ann", Inputs::Pca, e),
    }
    Ok(())
}

/// A JSON value's scalars as text, in order (a vector of length one reads as its scalar; numbers
/// with 17 significant digits, so equal text is equal bits).
fn atoms_text(v: &Json) -> String {
    fn walk(v: &Json, out: &mut Vec<String>) {
        match v {
            Json::Vector(items) | Json::List(items) => items.iter().for_each(|x| walk(x, out)),
            Json::Int(_) | Json::Num(_) => out.push(scio::fmt17(v.as_f64().unwrap())),
            other => out.push(other.to_text(0)),
        }
    }
    let mut out = Vec::new();
    walk(v, &mut out);
    out.join(", ")
}

/// A JSON value's numbers, in order (`None` if any atom is not a number).
fn atoms_num(v: &Json) -> Option<Vec<f64>> {
    match v {
        Json::Int(_) | Json::Num(_) => Some(vec![v.as_f64()?]),
        Json::Vector(items) | Json::List(items) => {
            let mut out = Vec::new();
            for x in items {
                out.extend(atoms_num(x)?);
            }
            Some(out)
        }
        _ => None,
    }
}

/// The rows of the candidate's TSV of numbers (`pca/embeddings.tsv`: barcode then values), as a
/// column-major matrix.
fn tsv_matrix(t: &Table, skip: usize) -> Result<(Vec<f64>, usize, usize), String> {
    let nrow = t.rows.len();
    let ncol = t.header.len().saturating_sub(skip);
    let mut m = vec![0.0; nrow * ncol];
    for (i, r) in t.rows.iter().enumerate() {
        for c in 0..ncol {
            m[c * nrow + i] = r.get(c + skip).and_then(|v| scio::parse17(v)).ok_or_else(|| format!("row {}: bad value", i + 1))?;
        }
    }
    Ok((m, nrow, ncol))
}

/// One column of numbers from two metadata tables, by name.
fn meta_numbers(t: &Table, col: usize) -> Vec<f64> {
    t.rows.iter().map(|r| r.get(col).and_then(|v| scio::parse17(v)).unwrap_or(f64::NAN)).collect()
}

/// K8: the candidate's final files against R's.
#[allow(clippy::too_many_arguments)]
fn compare_files(g: &mut Grader, reference: &Path, candidate: &Path, r_summary: &Option<Json>, c_result: &Json, c_run: &Json, r_meta: &Table, c_meta: &Table) -> Result<(), String> {
    let n = r_meta.rows.len();
    let same_n = c_meta.rows.len() == n;
    let stats_computed = str_list(c_run.get("stats").and_then(|s| s.get("computed")));
    // metadata.tsv: the counts' columns and the silhouette widths (the labels are K7's)
    for (col, inputs) in [("nCount_RNA", Inputs::Free), ("nFeature_RNA", Inputs::Free)] {
        match (r_meta.col(col), c_meta.col(col)) {
            (Some(a), Some(b)) if same_n => {
                let nd = (0..n).filter(|&i| r_meta.rows[i].get(a) != c_meta.rows[i].get(b)).count();
                g.exact("files", &format!("metadata.tsv: {col}"), inputs, n, nd, "as text".into());
            }
            (Some(_), None) => g.not_applicable("files", &format!("metadata.tsv: {col}"), inputs, "not in the candidate's metadata.tsv (a staged start without the counts' columns)".into()),
            _ => {}
        }
    }
    match (r_meta.col("silhouette_width"), c_meta.col("silhouette_width")) {
        (Some(a), Some(b)) if same_n => {
            let s = grade::stats(&meta_numbers(r_meta, a), &meta_numbers(c_meta, b), n, 1, false);
            g.numeric("K8", "metadata.tsv: silhouette_width", Inputs::Pca, "K8_silhouette_width", &s, "");
        }
        (Some(_), _) if !stats_computed.iter().any(|s| s == "silhouette") => {
            g.not_applicable("K8", "metadata.tsv: silhouette_width", Inputs::Pca, "the candidate did not compute silhouette widths".into())
        }
        (Some(_), _) => g.missing("K8", "metadata.tsv: silhouette_width", Inputs::Pca, "not in the candidate's metadata.tsv".into()),
        _ => {}
    }
    // silhouette_other against R's K8_silhouette_other and the levels approxSilhouette saw
    if let (Ok(o), Some(cc)) = (read_npy(&reference.join("stages/K8_silhouette_other.npy")), c_meta.col("silhouette_other")) {
        let levels: Vec<String> = read_json(&reference.join("stages/stages.json"))
            .ok()
            .and_then(|j| j.get("scalars").and_then(|s| s.get("chain")).map(|c| str_list(c.get("silhouette_levels"))))
            .unwrap_or_default();
        if same_n && o.nrow() == n && !levels.is_empty() {
            let nd = (0..n).filter(|&i| c_meta.rows[i].get(cc).map(String::as_str) != levels.get(o.get_i32(i, 0) as usize).map(String::as_str)).count();
            g.exact("K8", "metadata.tsv: silhouette_other = K8_silhouette_other", Inputs::Pca, n, nd, format!("{} levels", levels.len()));
        }
    }
    match (r_meta.col("neighborhood_purity"), c_meta.col("neighborhood_purity")) {
        (Some(a), Some(b)) if same_n => {
            let s = grade::stats(&meta_numbers(r_meta, a), &meta_numbers(c_meta, b), n, 1, false);
            g.numeric("K9", "metadata.tsv: neighborhood_purity", Inputs::Pca, "K9_purity", &s, "");
        }
        (Some(_), _) if !stats_computed.iter().any(|s| s == "purity") => {
            g.not_applicable("K9", "metadata.tsv: neighborhood_purity", Inputs::Pca, "the candidate did not compute purities".into())
        }
        (Some(_), _) => g.missing("K9", "metadata.tsv: neighborhood_purity", Inputs::Pca, "not in the candidate's metadata.tsv".into()),
        _ => {}
    }
    // purity_maximum against R's K9_maximum and the levels neighborPurity saw
    if let (Ok(o), Some(cc)) = (read_npy(&reference.join("stages/K9_maximum.npy")), c_meta.col("purity_maximum")) {
        let levels: Vec<String> = read_json(&reference.join("stages/stages.json"))
            .ok()
            .and_then(|j| j.get("scalars").and_then(|s| s.get("chain")).map(|c| str_list(c.get("purity_levels"))))
            .unwrap_or_default();
        if same_n && o.nrow() == n && !levels.is_empty() {
            let nd = (0..n).filter(|&i| c_meta.rows[i].get(cc).map(String::as_str) != levels.get(o.get_i32(i, 0) as usize).map(String::as_str)).count();
            // the floor found the maxima identical across kernels
            g.exact("K9", "metadata.tsv: purity_maximum = K9_maximum", Inputs::Pca, n, nd, format!("{} levels", levels.len()));
        }
    }
    // the candidate's certificate (C5: no flagged pair)
    certificate_check(g, "K9", "certificate of the candidate's purity search", Inputs::Free, c_result.get("purity"));
    let later: Vec<&str> = ["nCount_SCT", "nFeature_SCT"].into_iter().filter(|c| r_meta.col(c).is_some()).collect();
    if !later.is_empty() {
        g.not_applicable("later", "metadata.tsv: columns of later milestones", Inputs::Free, format!("not compared: {} (M13)", later.join(", ")));
    }

    // pca/: R's embeddings.npy, loadings.npy, loadings_features.tsv, stdev.npy, pca.json against
    // the candidate's scrust pca layout
    let rp = reference.join("pca");
    let cp = candidate.join("pca");
    if exists_any(&rp.join("embeddings.npy")) && cp.join("embeddings.tsv").exists() {
        let a = read_npy(&rp.join("embeddings.npy"))?;
        let (b, bn, bc) = tsv_matrix(&read_table(&cp.join("embeddings.tsv"))?, 1).map_err(|e| format!("{}: {e}", cp.join("embeddings.tsv").display()))?;
        let (an, ac) = (a.nrow(), a.ncol());
        if (an, ac) != (bn, bc) {
            g.exact("K2", "pca/ embeddings", Inputs::Pca, an * ac, an * ac, format!("{an} x {ac} in R's, {bn} x {bc} in Rust's"));
        } else {
            let av = a.block_f64(0, an);
            let k = ac.min(30);
            let s1 = grade::stats(&av[..k * an], &b[..k * an], an, k, true);
            g.numeric("K2", "pca/ embeddings, dims 1:30", Inputs::Pca, "S3_embeddings_dims_1_30", &s1, "");
            if ac > 30 {
                let s2 = grade::stats(&av[k * an..], &b[k * an..], an, ac - k, true);
                g.numeric("K2", &format!("pca/ embeddings, dims 31:{ac}"), Inputs::Pca, "S3_embeddings_dims_31_50", &s2, "");
            }
        }
        if let (Ok(la), Ok(lt)) = (read_npy(&rp.join("loadings.npy")), read_table(&cp.join("loadings.tsv"))) {
            let rf: Vec<String> = read_text(&rp.join("loadings_features.tsv")).unwrap_or_default().lines().map(String::from).collect();
            let nc = lt.col("name").unwrap_or(1);
            let cf: Vec<String> = lt.rows.iter().map(|r| r.get(nc).cloned().unwrap_or_default()).collect();
            let nd = if rf.len() != cf.len() { rf.len().max(cf.len()) } else { rf.iter().zip(&cf).filter(|(x, y)| x != y).count() };
            g.exact("K2", "pca/ loadings' features, in order", Inputs::Pca, rf.len(), nd, String::new());
            if let Ok((lb, ln, lc)) = tsv_matrix(&lt, 2) {
                if (la.nrow(), la.ncol()) == (ln, lc) {
                    let s = grade::stats(&la.block_f64(0, ln), &lb, ln, lc, true);
                    g.numeric("K2", "pca/ loadings", Inputs::Pca, "S3_loadings", &s, "");
                } else {
                    g.exact("K2", "pca/ loadings", Inputs::Pca, la.nrow() * la.ncol(), la.nrow() * la.ncol(), format!("{} x {} in R's, {ln} x {lc} in Rust's", la.nrow(), la.ncol()));
                }
            }
        }
        if let (Ok(sa), Ok(st)) = (read_npy(&rp.join("stdev.npy")), read_table(&cp.join("stdev.tsv"))) {
            let a = sa.f64s();
            let b: Vec<f64> = st.rows.iter().map(|r| r.get(1).and_then(|v| scio::parse17(v)).unwrap_or(f64::NAN)).collect();
            let s = grade::stats(&a, &b, a.len(), 1, false);
            g.numeric("K2", "pca/ stdev", Inputs::Pca, "S3_sdev", &s, "");
        }
        if let (Ok(pj), Some(tv)) = (read_json(&rp.join("pca.json")), c_result.get("pca").and_then(|p| p.get("total_variance")).and_then(Json::as_f64)) {
            if let Some(rv) = pj.get("total_variance").and_then(Json::as_f64) {
                let s = grade::stats(&[rv], &[tv], 1, 1, false);
                g.numeric("K2", "total variance", Inputs::Pca, "total_variance", &s, "");
            }
        }
    } else if exists_any(&rp.join("embeddings.npy")) {
        g.not_applicable("K2", "pca/", Inputs::Pca, "the candidate has no pca/embeddings.tsv (a run started from the SNN graph without --pca)".into());
    }

    // K6: the graphs
    for (name, file) in [("SCT_nn", "nn.mtx.gz"), ("SCT_snn", "snn.mtx.gz")] {
        let rpth = reference.join("graphs").join(file);
        let cpth = candidate.join("neighbors").join(file);
        if !exists_any(&rpth) {
            continue;
        }
        let (a, b) = match (scio::mtx::read_entries(&rpth), scio::mtx::read_entries(&cpth)) {
            (Ok(a), Ok(b)) => (a, b),
            (Err(e), _) | (_, Err(e)) => {
                g.missing("K6", name, Inputs::Pca, e.to_string());
                continue;
            }
        };
        let m = a.vals.len().min(b.vals.len());
        let mut differ = 0usize;
        let mut first = None;
        for t in 0..m {
            if a.rows[t] != b.rows[t] || a.cols[t] != b.cols[t] || a.vals[t].to_bits() != b.vals[t].to_bits() {
                differ += 1;
                if first.is_none() {
                    first = Some(format!("first at entry {}: R ({}, {}) {}, Rust ({}, {}) {}", t + 1, a.rows[t] + 1, a.cols[t] + 1, scio::fmt17(a.vals[t]), b.rows[t] + 1, b.cols[t] + 1, scio::fmt17(b.vals[t])));
                }
            }
        }
        differ += a.vals.len().max(b.vals.len()) - m;
        if (a.nrow, a.ncol) != (b.nrow, b.ncol) {
            differ = differ.max(1);
        }
        let mut detail = format!("{} x {}, {} entries in R's, {} x {}, {} in Rust's", a.nrow, a.ncol, a.vals.len(), b.nrow, b.ncol, b.vals.len());
        if let Some(x) = first {
            detail.push_str("; ");
            detail.push_str(&x);
        }
        g.exact("K6", name, Inputs::Pca, a.vals.len().max(b.vals.len()), differ, detail);
    }

    // the result section's numbers
    if let Some(rr) = r_summary.as_ref().and_then(|s| s.get("result")) {
        let sizes = {
            let col = c_meta.col("seurat_clusters");
            let levels: Vec<String> = read_text(&candidate.join("seurat_clusters_levels.txt")).unwrap_or_default().lines().map(String::from).collect();
            col.map(|c| levels.iter().map(|l| c_meta.rows.iter().filter(|r| r.get(c) == Some(l)).count() as f64).collect::<Vec<f64>>())
        };
        let c_sct = c_result.get("sctransform");
        let int_text = |v: Option<f64>| v.map(|x| (x as i64).to_string()).unwrap_or_else(|| "-".into());
        let mut pairs: Vec<(&str, Inputs, String, String)> = vec![
            ("n_cells", Inputs::Free, int_text(jnum(rr, "n_cells")), int_text(jnum(c_result, "n_cells"))),
            ("n_clusters", Inputs::Pca, int_text(jnum(rr, "n_clusters")), int_text(jnum(c_result, "n_clusters"))),
            (
                "cluster_sizes",
                Inputs::Pca,
                num_list(rr.get("cluster_sizes")).iter().map(|v| (*v as i64).to_string()).collect::<Vec<_>>().join(","),
                sizes.unwrap_or_default().iter().map(|v| (*v as i64).to_string()).collect::<Vec<_>>().join(","),
            ),
        ];
        if let Some(cs) = c_sct.filter(|s| !matches!(s, Json::Null)) {
            pairs.push(("n_variable_features", Inputs::Sct, int_text(jnum(rr, "n_variable_features")), int_text(jnum(cs, "n_variable_features"))));
            pairs.push((
                "scale_data_dim",
                Inputs::Sct,
                num_list(rr.get("scale_data_dim")).iter().map(|v| (*v as i64).to_string()).collect::<Vec<_>>().join(","),
                format!("{},{}", int_text(jnum(cs, "n_variable_features")), int_text(jnum(cs, "n_cells"))),
            ));
        }
        if let Some(o) = rr.get("optimizer") {
            let cr0 = match c_result.get("resolutions") {
                Some(Json::List(v) | Json::Vector(v)) => v.first().cloned(),
                _ => None,
            };
            if let Some(cr0) = cr0 {
                let printed = str_list(cr0.get("printed"));
                let pnum = |p: &str| printed.iter().find_map(|l| l.strip_prefix(p)).unwrap_or("-").to_string();
                pairs.push(("optimizer: nodes", Inputs::Pca, int_text(jnum(o, "nodes")), pnum("Number of nodes: ")));
                pairs.push(("optimizer: edges", Inputs::Pca, int_text(jnum(o, "edges")), pnum("Number of edges: ")));
                pairs.push(("optimizer: max_modularity", Inputs::Pca, o.get("max_modularity").and_then(Json::as_str).unwrap_or("-").to_string(), cr0.get("max_modularity_printed").and_then(Json::as_str).unwrap_or("-").to_string()));
                pairs.push(("optimizer: communities", Inputs::Pca, int_text(jnum(o, "communities")), int_text(jnum(&cr0, "communities"))));
                pairs.push(("optimizer: singletons", Inputs::Pca, int_text(jnum(o, "singletons")), int_text(jnum(&cr0, "n_singletons"))));
                pairs.push(("optimizer: final_clusters", Inputs::Pca, int_text(jnum(o, "final_clusters")), int_text(jnum(&cr0, "final_clusters"))));
            }
        }
        for (key, inputs, a, b) in pairs {
            g.exact("files", &format!("result: {key}"), inputs, 1, (a != b) as usize, if a == b { a } else { format!("R {a}, Rust {b}") });
        }
        // R's digests: scale.data, and .Random.seed after SCTransform and RunPCA
        if let Some(rd) = rr.get("scale_data_sha256").and_then(Json::as_str) {
            match c_result.get("scale_data_sha256").and_then(Json::as_str) {
                Some(cd) if cd == rd => g.exact("K1", "result: scale_data_sha256", Inputs::Sct, 1, 0, rd.to_string()),
                Some(cd) if g.class(Inputs::Sct) == Class::BlasDependent && !g.mode.is_exact() => {
                    let crit = g.floor.describe("S2d_scale_data");
                    g.push("K1", "result: scale_data_sha256", Class::BlasDependent, crit, Status::NotApplicable, 1, 1, format!("the digests differ (R {rd}, Rust {cd}); a digest cannot be graded against the floor's bounds"))
                }
                Some(cd) => g.exact("K1", "result: scale_data_sha256", Inputs::Sct, 1, 1, format!("R {rd}, Rust {cd}")),
                None => g.not_applicable("K1", "result: scale_data_sha256", Inputs::Sct, "the candidate did not run SCTransform".into()),
            }
        }
        let rng_after = rr.get("rng_after");
        for (key, label, inputs) in [("SCTransform", "K1", Inputs::Free), ("RunPCA", "K2", Inputs::Pca)] {
            let a = rng_after.and_then(|x| x.get(key)).and_then(Json::as_str);
            let b = c_result.get("rng_after_sha256").and_then(|x| x.get(key)).and_then(Json::as_str);
            match (a, b) {
                (Some(a), Some(b)) => g.exact(if label == "K1" { "K1" } else { "K2" }, &format!("result: .Random.seed after {key} (digest)"), inputs, 1, (a != b) as usize, if a == b { a.to_string() } else { format!("R {a}, Rust {b}") }),
                (Some(_), None) => g.not_applicable(if label == "K1" { "K1" } else { "K2" }, &format!("result: .Random.seed after {key} (digest)"), inputs, format!("the candidate did not run {key}")),
                _ => {}
            }
        }
    }
    Ok(())
}

/// The labels, levels and silhouette widths against the saved object's export; graded under the
/// saved object's kernel.
fn compare_saved(g: &mut Grader, saved: &Path, saved_kernel: &str, c_summary: &Json, c_meta: &Table, candidate: &Path) -> Result<Mode, String> {
    let s_meta = read_table(&saved.join("metadata.tsv"))?;
    let n = s_meta.rows.len();
    let mode = grade::mode(Some(saved_kernel), grade::recorded_kernel(Some(c_summary)).as_deref());
    let same_n = c_meta.rows.len() == n;
    let sb = s_meta.col("barcode").unwrap_or(0);
    let cb = c_meta.col("barcode").unwrap_or(0);
    let nd = if same_n { (0..n).filter(|&i| s_meta.rows[i].get(sb) != c_meta.rows[i].get(cb)).count() } else { n };
    g.exact("saved", "barcodes", Inputs::Free, n, nd, String::new());
    for col in ["SCT_snn_res.0.8", "seurat_clusters"] {
        match (s_meta.col(col), c_meta.col(col)) {
            (Some(a), Some(b)) if same_n => {
                let nd = (0..n).filter(|&i| s_meta.rows[i].get(a) != c_meta.rows[i].get(b)).count();
                g.exact("saved", &format!("{col}: labels"), Inputs::Pca, n, nd, "metadata.tsv".into());
            }
            _ => g.missing("saved", &format!("{col}: labels"), Inputs::Pca, format!("{col} missing")),
        }
        let f = format!("{col}_levels.txt");
        match (read_any(&saved.join(&f)), read_any(&candidate.join(&f))) {
            (Ok(a), Ok(b)) => g.exact("saved", &format!("{col}: levels"), Inputs::Pca, 1, (a != b) as usize, format!("{} levels", String::from_utf8_lossy(&a).lines().count())),
            (Err(e), _) | (_, Err(e)) => g.missing("saved", &format!("{col}: levels"), Inputs::Pca, e),
        }
    }
    match (s_meta.col("silhouette_width"), c_meta.col("silhouette_width")) {
        (Some(a), Some(b)) if same_n => {
            let s = grade::stats(&meta_numbers(&s_meta, a), &meta_numbers(c_meta, b), n, 1, false);
            // graded under the saved object's kernel
            let keep = std::mem::replace(&mut g.mode, mode.clone());
            g.numeric("saved", "silhouette_width", Inputs::Pca, "K8_silhouette_width", &s, &format!("saved under {saved_kernel}"));
            g.mode = keep;
        }
        (Some(_), _) => g.not_applicable("saved", "silhouette_width", Inputs::Pca, "the candidate did not compute silhouette widths".into()),
        _ => {}
    }
    match (s_meta.col("neighborhood_purity"), c_meta.col("neighborhood_purity")) {
        (Some(a), Some(b)) if same_n => {
            let s = grade::stats(&meta_numbers(&s_meta, a), &meta_numbers(c_meta, b), n, 1, false);
            let keep = std::mem::replace(&mut g.mode, mode.clone());
            g.numeric("saved", "neighborhood_purity", Inputs::Pca, "K9_purity", &s, &format!("saved under {saved_kernel}"));
            g.mode = keep;
        }
        (Some(_), _) => g.not_applicable("saved", "neighborhood_purity", Inputs::Pca, "the candidate did not compute purities".into()),
        _ => {}
    }
    let later: Vec<&str> = ["nCount_SCT", "nFeature_SCT"].into_iter().filter(|c| s_meta.col(c).is_some()).collect();
    if !later.is_empty() {
        g.not_applicable("later", "saved: columns of later milestones", Inputs::Free, format!("not compared: {} (M13)", later.join(", ")));
    }
    Ok(mode)
}

/// approxSilhouette on R's own embeddings and labels (teacher forcing): exact on any machine.
fn teacher_forced_silhouette(g: &mut Grader, reference: &Path, r_meta: &Table, r_scalars: Option<&Json>) -> Result<(), String> {
    let emb_path = reference.join("pca/embeddings.npy");
    let w_path = reference.join("stages/K8_silhouette_width.npy");
    if !exists_any(&emb_path) || !exists_any(&w_path) {
        return Ok(());
    }
    let emb = read_npy(&emb_path)?;
    let (n, ncol) = (emb.nrow(), emb.ncol());
    let e = emb.block_f64(0, n);
    let col = r_meta.col("seurat_clusters").ok_or("R's metadata.tsv: no seurat_clusters")?;
    let levels: Vec<String> = read_text(&reference.join("seurat_clusters_levels.txt"))?.lines().map(String::from).collect();
    let mut codes = Vec::with_capacity(n);
    for r in &r_meta.rows {
        let l = r.get(col).ok_or("seurat_clusters: short row")?;
        codes.push(levels.iter().position(|v| v == l).ok_or_else(|| format!("label {l} is not a level"))? as u32);
    }
    let f = Factor { name: "seurat_clusters".into(), levels, codes };
    let nd = r_scalars.and_then(|s| jnum(s, "silhouette_dims")).map(|v| v as usize).unwrap_or(30);
    let dims: Vec<usize> = (1..=nd).collect();
    let s = approx_silhouette(&e, n, ncol, &dims, &f)?;
    let free = Inputs::Free;
    let want = read_npy(&w_path)?.f64s();
    let st = grade::stats(&want, &s.width, n, 1, false);
    g.numeric("teacher", "approxSilhouette width (R's embeddings and labels)", free, "K8_silhouette_width", &st, "");
    if let Ok(o) = read_npy(&reference.join("stages/K8_silhouette_other.npy")) {
        let want: Vec<i64> = (0..o.nrow()).map(|i| o.get_i32(i, 0) as i64).collect();
        let got: Vec<i64> = s.other.iter().map(|x| x.map_or(-1, |v| v as i64)).collect();
        let nd = if want.len() != got.len() { want.len().max(got.len()) } else { want.iter().zip(&got).filter(|(a, b)| a != b).count() };
        g.exact("teacher", "approxSilhouette other (R's embeddings and labels)", free, want.len(), nd, String::new());
    }
    if let Ok(c) = read_npy(&reference.join("stages/K8_centroids.npy")) {
        let (k, d) = (c.nrow(), c.ncol());
        let mut got = vec![0.0; k * d];
        if s.centroids.len() == k * d {
            for a in 0..k {
                for b in 0..d {
                    got[a + b * k] = s.centroids[a * d + b];
                }
            }
        }
        let st = grade::stats(&c.block_f64(0, k), &got, k, d, false);
        g.numeric("teacher", "approxSilhouette centroids (R's embeddings and labels)", free, "K8_centroids", &st, "");
    }
    if let Ok(v) = read_npy(&reference.join("stages/K8_cluster_var.npy")) {
        let st = grade::stats(&v.f64s(), &s.cluster_var, v.nrow(), 1, false);
        g.numeric("teacher", "approxSilhouette cluster_var (R's embeddings and labels)", free, "K8_cluster_var", &st, "");
    }
    if let Some(c) = r_meta.col("silhouette_width") {
        let st = grade::stats(&meta_numbers(r_meta, c), &s.width, n, 1, false);
        g.numeric("teacher", "approxSilhouette width = R's metadata.tsv silhouette_width", free, "K8_silhouette_width", &st, "");
    }
    Ok(())
}

/// Neighbour lists from a stage directory's `K9_neighbors_ptr` and `K9_neighbors_idx` (1-based
/// indices, as R writes them): one list per cell.
fn read_lists(d: &StageDir) -> Result<Vec<Vec<i32>>, String> {
    let ptr = d.read("K9_neighbors_ptr.npy").map_err(|e| e.to_string())?;
    let idx = d.read("K9_neighbors_idx.npy").map_err(|e| e.to_string())?;
    split_lists(&ptr, &idx).map_err(|e| format!("{}: {e}", d.dir.display()))
}

/// One list per cell from a pointer array and the indices.
fn split_lists(ptr: &Npy, idx: &Npy) -> Result<Vec<Vec<i32>>, String> {
    let n = ptr.nrow().saturating_sub(1);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = (ptr.get_i32(i, 0) as usize, ptr.get_i32(i + 1, 0) as usize);
        if a > b || b > idx.nrow() {
            return Err(format!("K9_neighbors_ptr is not a pointer array at cell {}", i + 1));
        }
        out.push((a..b).map(|r| idx.get_i32(r, 0)).collect());
    }
    Ok(out)
}

/// Compares two sets of neighbour lists cell by cell, in order: graded on the cells whose list
/// differs (floor key `K9_neighbors`). A list that holds the same cells in another order is
/// counted apart in the detail: BiocNeighbors orders equal raw distances by KMKNN's internal
/// index, the port by the cell's index.
fn compare_lists(g: &mut Grader, stage: &'static str, name: &str, inputs: Inputs, a: &[Vec<i32>], b: &[Vec<i32>]) {
    if a.len() != b.len() {
        g.exact(stage, name, inputs, a.len().max(b.len()), a.len().max(b.len()), format!("{} cells in R's lists, {} in the candidate's", a.len(), b.len()));
        return;
    }
    let (mut differ, mut reordered, mut first) = (0usize, 0usize, None);
    let mut entries = 0usize;
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        entries += x.len();
        if x == y {
            continue;
        }
        differ += 1;
        let (mut xs, mut ys) = (x.clone(), y.clone());
        xs.sort_unstable();
        ys.sort_unstable();
        if xs == ys {
            reordered += 1;
        }
        if first.is_none() {
            first = Some(format!("first at cell {}: {} neighbours in R's, {} in the candidate's", i + 1, x.len(), y.len()));
        }
    }
    let s = Stats { n: a.len(), n_differ: differ, rows_differ: differ, ..Stats::default() };
    let mut detail = format!("{} cells, {entries} neighbours in R's lists", a.len());
    if differ > 0 {
        detail.push_str(&format!("; {differ} cells' lists differ ({reordered} only in the order of their cells)"));
    }
    if let Some(f) = first {
        detail.push_str("; ");
        detail.push_str(&f);
    }
    g.numeric(stage, name, inputs, "K9_neighbors", &s, &detail);
}

/// The certificate of a purity search (`result.purity` of a summary): identical when nothing is
/// flagged; otherwise the flagged pairs are reported, with their kinds, and the check passes
/// (C5's K9 criterion as M12's decision settled it).
fn certificate_check(g: &mut Grader, stage: &'static str, name: &str, inputs: Inputs, purity: Option<&Json>) {
    let Some(p) = purity.filter(|p| !matches!(p, Json::Null)) else {
        return;
    };
    if !matches!(p.get("computed"), Some(Json::Bool(true))) {
        g.not_applicable(stage, name, inputs, "purity not computed (the workflow's guard: fewer than 2 cells or clusters)".into());
        return;
    }
    let Some(c) = p.get("certificate") else {
        g.missing(stage, name, inputs, "no certificate in the candidate's summary.json".into());
        return;
    };
    let flagged = jnum(c, "flagged_pairs").unwrap_or(f64::NAN) as usize;
    let mut detail = format!(
        "{} pair(s) flagged ({} the median pair, {} other; {} near a k-th distance, {} near the radius); step 3's band {} ulps of the radius squared",
        flagged,
        jnum(c, "flagged_median_pair").unwrap_or(0.0) as usize,
        jnum(c, "flagged_other").unwrap_or(0.0) as usize,
        jnum(c, "flagged_kth_distance").unwrap_or(0.0) as usize,
        jnum(c, "flagged_radius").unwrap_or(0.0) as usize,
        jnum(c, "radius_band_ulps").map(|v| format!("{v:.0}")).unwrap_or_else(|| "-".into())
    );
    if let Some(Json::List(pairs) | Json::Vector(pairs)) = c.get("pairs") {
        let shown: Vec<String> = pairs
            .iter()
            .take(6)
            .map(|q| {
                format!(
                    "{} {} cells {}-{} at {} ulps",
                    q.get("kind").and_then(Json::as_str).unwrap_or("?"),
                    q.get("step").and_then(Json::as_str).unwrap_or("?"),
                    jnum(q, "cell").unwrap_or(0.0) as i64,
                    jnum(q, "other").unwrap_or(0.0) as i64,
                    jnum(q, "ulps_from_threshold").map(|v| format!("{v:.0}")).unwrap_or_else(|| "?".into())
                )
            })
            .collect();
        if !shown.is_empty() {
            detail.push_str(&format!(": {}{}", shown.join(", "), if pairs.len() > 6 { ", ..." } else { "" }));
        }
    }
    // a property of the search, not a comparison: the same in any mode
    let status = if flagged == 0 { Status::Identical } else { Status::Reported };
    g.push(stage, name, Class::BlasFree, "no flagged pair; flagged pairs are reported, not failed (M12's decision)".into(), status, 1, (flagged > 0) as usize, detail);
}

/// neighborPurity on R's own embeddings and labels (teacher forcing): exact on any machine.
fn teacher_forced_purity(g: &mut Grader, reference: &Path, r_meta: &Table, r_scalars: Option<&Json>) -> Result<(), String> {
    let emb_path = reference.join("pca/embeddings.npy");
    let d_path = reference.join("stages/K9_dist50.npy");
    if !exists_any(&emb_path) || !exists_any(&d_path) {
        return Ok(());
    }
    let emb = read_npy(&emb_path)?;
    let (n, ncol) = (emb.nrow(), emb.ncol());
    let e = emb.block_f64(0, n);
    let col = r_meta.col("seurat_clusters").ok_or("R's metadata.tsv: no seurat_clusters")?;
    let levels: Vec<String> = read_text(&reference.join("seurat_clusters_levels.txt"))?.lines().map(String::from).collect();
    let mut codes = Vec::with_capacity(n);
    for r in &r_meta.rows {
        let l = r.get(col).ok_or("seurat_clusters: short row")?;
        codes.push(levels.iter().position(|v| v == l).ok_or_else(|| format!("label {l} is not a level"))? as u32);
    }
    let f = Factor { name: "seurat_clusters".into(), levels, codes };
    let nd = r_scalars.and_then(|s| jnum(s, "silhouette_dims")).map(|v| v as usize).unwrap_or(30);
    let dims: Vec<usize> = (1..=nd).collect();
    let k = r_scalars.and_then(|s| jnum(s, "purity_k")).map(|v| v as usize).unwrap_or(50);
    let has_lists = exists_any(&reference.join("stages/K9_neighbors_ptr.npy"));
    let threads = std::thread::available_parallelism().map(|v| v.get()).unwrap_or(1);
    let p = neighbor_purity(&e, n, ncol, &dims, &f, k, true, has_lists, threads)?;
    let free = Inputs::Free;
    let st = |name: &str| -> Result<Npy, String> { read_npy(&reference.join(format!("stages/{name}.npy"))) };
    let s = grade::stats(&st("K9_dist50")?.f64s(), &p.dist_k, n, 1, false);
    g.numeric("teacher", "neighborPurity k-th distances (R's embeddings)", free, "K9_dist50", &s, "");
    if let Some(r) = r_scalars.and_then(|sc| jnum(sc, "purity_radius")) {
        let s = grade::stats(&[r], &[p.radius], 1, 1, false);
        g.numeric("teacher", "neighborPurity radius (R's embeddings)", free, "purity_radius", &s, &scio::fmt17(r));
    }
    if let Ok(c) = st("K9_neighbor_count") {
        let a: Vec<f64> = (0..c.nrow()).map(|i| c.get_i32(i, 0) as f64).collect();
        let b: Vec<f64> = p.neighbor_count.iter().map(|&v| v as f64).collect();
        let s = grade::stats(&a, &b, n, 1, false);
        g.numeric("teacher", "neighborPurity neighbour counts (R's embeddings)", free, "K9_neighbor_count", &s, "");
    }
    if has_lists {
        let r_lists = split_lists(&st("K9_neighbors_ptr")?, &st("K9_neighbors_idx")?)?;
        let (ptr, idx) = p.neighbor_lists.clone().unwrap_or_default();
        let c_lists: Vec<Vec<i32>> = (0..n).map(|i| idx[ptr[i] as usize..ptr[i + 1] as usize].iter().map(|&j| j as i32 + 1).collect()).collect();
        compare_lists(g, "teacher", "neighborPurity neighbour lists (R's embeddings)", free, &r_lists, &c_lists);
    }
    for (name, got) in [("K9_weights", &p.weights), ("K9_totals", &p.totals), ("K9_purity", &p.purity)] {
        if let Ok(a) = st(name) {
            let s = grade::stats(&a.f64s(), got, n, 1, false);
            g.numeric("teacher", &format!("neighborPurity {} (R's embeddings and labels)", name.trim_start_matches("K9_")), free, name, &s, "");
        }
    }
    if let Ok(a) = st("K9_sums") {
        let nk = p.uclust.len();
        let mut got = vec![0.0; n * nk];
        if p.sums.len() == n * nk {
            for i in 0..n {
                for c in 0..nk {
                    got[i + c * n] = p.sums[i * nk + c];
                }
            }
        }
        let s = grade::stats(&a.block_f64(0, a.nrow()), &got, a.nrow(), a.ncol(), false);
        g.numeric("teacher", "neighborPurity sums (R's embeddings and labels)", free, "K9_sums", &s, &format!("{} x {}", a.nrow(), a.ncol()));
    }
    if let Ok(o) = st("K9_maximum") {
        let want: Vec<i64> = (0..o.nrow()).map(|i| o.get_i32(i, 0) as i64).collect();
        let got: Vec<i64> = p.maximum.iter().map(|x| x.map_or(-1, |v| v as i64)).collect();
        let nd = if want.len() != got.len() { want.len().max(got.len()) } else { want.iter().zip(&got).filter(|(a, b)| a != b).count() };
        g.exact("teacher", "neighborPurity maximum (R's embeddings and labels)", free, want.len(), nd, String::new());
    }
    if let Some(c) = r_meta.col("neighborhood_purity") {
        let s = grade::stats(&meta_numbers(r_meta, c), &p.purity, n, 1, false);
        g.numeric("teacher", "neighborPurity purity = R's metadata.tsv neighborhood_purity", free, "K9_purity", &s, "");
    }
    let levels: Vec<String> = p.uclust.iter().map(|&c| f.levels[c as usize].clone()).collect();
    let pj = crate::chain_output::purity_json(&p, &levels, &[]);
    certificate_check(g, "teacher", "certificate of the search on R's embeddings", free, Some(&pj));
    Ok(())
}
