//! Comparison of a `scrust neighbors` run with an R run of FindNeighbors on the same
//! embeddings (docs/PLAN.md C5, stages K3 to K6), on plain files:
//!
//! - **parameters**: the metric, the number of dimensions, trees and neighbours and search.k
//!   that R recorded (`stages/stages.json`) against the candidate's `summary.json`;
//! - **K3, the items**: R's `getItemsVector` for every cell (`stages/K3_items.npy`) against the
//!   first n nodes of the candidate's `annoy.ann`, every float bitwise;
//! - **K4, the index**: `annoy.ann` against R's saved index (`stages/annoy.ann`), byte for byte;
//!   on a difference, the first differing node and field, and how many nodes differ;
//! - **K5, the neighbours**: every cell's k ids (`stages/K5_nn_index.npy`, 1-based) and float
//!   distances (`stages/K5_nn_dist.npy`) against `knn_index.tsv` and `knn_distance.tsv`, bitwise;
//!   a NaN distance (only overflowing inputs make one) matches a NaN, because the text spells
//!   every NaN alike;
//! - **K6, the graphs**: `graphs/nn.mtx.gz` and `graphs/snn.mtx.gz` (`SCT_nn` and `SCT_snn` as
//!   R stored them) against the candidate's: dimensions, every stored entry in order, values
//!   bitwise.
//!
//! The R side is an instrumented clustering run (`validation/R/instrumented_clustering.R`, M8)
//! or a neighbour reference run in the same layout (`validation/R/neighbors_references.R`,
//! M9). Every criterion is exact (`validation/tolerances.toml`, K3 to K6).

use crate::annoy::{AnnoyIndex, Metric};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

/// The outcome of one check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Identical,
    Differ,
    /// The reference or the candidate lacks the file.
    Missing,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Identical => "identical",
            Status::Differ => "differ",
            Status::Missing => "missing",
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
    /// True when every check is identical (a missing file fails).
    pub fn pass(&self) -> bool {
        self.checks.iter().all(|c| c.status == Status::Identical)
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
            "# scrust validate neighbors\n\nreference: `{}`\ncandidate: `{}`\n\n**{}**\n\n| stage | check | status | compared | differing | detail |\n|---|---|---|---|---|---|\n",
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

/// A file's bytes, from `path` or `path.gz` (the committed fixtures are gzipped), decompressed.
fn read_any(path: &Path) -> Result<Vec<u8>, String> {
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let p = if path.exists() || !gz.exists() { path } else { gz.as_path() };
    scio::tenx::read_maybe_gz(p).map_err(|e| e.to_string())
}

fn read_npy(path: &Path) -> Result<Npy, String> {
    Npy::parse(&read_any(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

/// The scalars of R's `stages.json` (context `chain`).
fn r_scalars(reference: &Path) -> Result<Json, String> {
    let p = reference.join("stages/stages.json");
    let text = String::from_utf8(read_any(&p)?).map_err(|_| format!("{}: not text", p.display()))?;
    let j = Json::parse(&text).map_err(|e| format!("{}: {e}", p.display()))?;
    j.get("scalars").and_then(|s| s.get("chain")).cloned().ok_or_else(|| format!("{}: no scalars for context chain", p.display()))
}

/// A knn table written by `scrust neighbors`: barcode, then k values per row.
fn read_knn_table(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>), String> {
    let t = Table::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let barcodes = t.rows.iter().map(|r| r[0].clone()).collect();
    let vals = t.rows.iter().map(|r| r[1..].to_vec()).collect();
    Ok((barcodes, vals))
}

fn json_num(j: &Json, key: &str) -> Option<f64> {
    j.get(key).and_then(Json::as_f64)
}

/// Compares the candidate directory (`scrust neighbors --out`, written with `--annoy yes`) with
/// the R run directory `reference`.
pub fn compare(reference: &Path, candidate: &Path) -> Result<Report, String> {
    let mut checks = Vec::new();
    let sc = r_scalars(reference)?;
    let metric = Metric::parse(sc.get("annoy_metric").and_then(Json::as_str).ok_or("stages.json: no annoy_metric")?)?;
    let f = json_num(&sc, "annoy_f").ok_or("stages.json: no annoy_f")? as usize;
    let k = json_num(&sc, "knn_k").ok_or("stages.json: no knn_k")? as usize;

    // parameters
    let sp = candidate.join("summary.json");
    let summary = Json::parse(&std::fs::read_to_string(&sp).map_err(|e| format!("{}: {e}", sp.display()))?).map_err(|e| format!("{}: {e}", sp.display()))?;
    let par = summary.get("parameters").cloned().unwrap_or(Json::Null);
    let res = summary.get("result").cloned().unwrap_or(Json::Null);
    let pairs: Vec<(&str, String, String)> = vec![
        ("metric", metric.seurat_name().to_string(), par.get("metric").and_then(Json::as_str).unwrap_or("?").to_string()),
        ("dims", f.to_string(), par.get("dims").map(|d| if let Json::List(v) | Json::Vector(v) = d { v.len().to_string() } else { "?".into() }).unwrap_or_default()),
        ("n_trees", json_num(&sc, "annoy_n_trees_arg").map(|v| v.to_string()).unwrap_or_default(), json_num(&par, "n_trees").map(|v| v.to_string()).unwrap_or_default()),
        ("k", k.to_string(), json_num(&res, "k").map(|v| v.to_string()).unwrap_or_default()),
        ("search_k", json_num(&sc, "knn_search_k").map(|v| v.to_string()).unwrap_or_default(), json_num(&par, "search_k").map(|v| v.to_string()).unwrap_or_default()),
    ];
    let bad: Vec<String> = pairs.iter().filter(|(_, a, b)| a != b).map(|(n, a, b)| format!("{n}: R {a}, Rust {b}")).collect();
    checks.push(check("params", "metric, dims, n.trees, k, search.k", pairs.len(), bad.len(), bad.join("; ")));

    // K4 and K3: the index and its items
    let r_ann = reference.join("stages/annoy.ann");
    let c_ann = candidate.join("annoy.ann");
    let mut n_cells = None;
    match (read_any(&r_ann), std::fs::read(&c_ann).map_err(|e| e.to_string())) {
        (Ok(a), Ok(b)) => {
            let ri = AnnoyIndex::from_bytes(metric, f, &a).map_err(|e| format!("{}: {e}", r_ann.display()))?;
            let s = ri.node_bytes();
            let (na, nb) = (a.len() / s, b.len() / s);
            let mut nodes_differ = 0usize;
            let mut first = None;
            for t in 0..na.min(nb) {
                let (x, y) = (&a[t * s..(t + 1) * s], &b[t * s..(t + 1) * s]);
                if x != y {
                    nodes_differ += 1;
                    if first.is_none() {
                        let byte = x.iter().zip(y).position(|(u, v)| u != v).unwrap();
                        first = Some(format!("first at node {t} ({}), field {}", ri.describe_node(t), ri.field_at(byte)));
                    }
                }
            }
            nodes_differ += na.max(nb) - na.min(nb);
            let mut detail = format!("{na} nodes of {s} bytes in R's, {nb} in Rust's");
            if let Some(fd) = first {
                detail.push_str("; ");
                detail.push_str(&fd);
            }
            checks.push(check("K4", "annoy.ann", na.max(nb), nodes_differ, detail));
            n_cells = Some(ri.n_items());
            // K3: R's getItemsVector (or, where it was not dumped, the items of R's index file)
            // against the candidate's items
            let (items, source): (Vec<f64>, &str) = match read_npy(&reference.join("stages/K3_items.npy")) {
                Ok(a) if a.ncol() == f => {
                    let n = a.nrow();
                    ((0..n).flat_map(|i| (0..f).map(move |z| (i, z))).map(|(i, z)| a.get(i, z)).collect(), "K3_items.npy")
                }
                _ => ((0..ri.n_items()).flat_map(|i| ri.get_item(i).iter().map(|&v| v as f64).collect::<Vec<_>>()).collect(), "R's annoy.ann"),
            };
            match AnnoyIndex::from_bytes(metric, f, &b) {
                Ok(ci) => {
                    let n = items.len() / f;
                    let mut differ = 0usize;
                    if ci.n_items() < n {
                        differ = n * f;
                    } else {
                        for i in 0..n {
                            let v = ci.get_item(i);
                            for z in 0..f {
                                if (v[z] as f64).to_bits() != items[i * f + z].to_bits() {
                                    differ += 1;
                                }
                            }
                        }
                    }
                    checks.push(check("K3", "items", n * f, differ, format!("{n} cells x {f} dims, against {source}")));
                }
                Err(e) => checks.push(missing("K3", "items", format!("candidate annoy.ann: {e}"))),
            }
        }
        (Err(e), _) => {
            checks.push(missing("K4", "annoy.ann", format!("{}: {e}", r_ann.display())));
            checks.push(missing("K3", "items", "no reference index".into()));
        }
        (_, Err(e)) => {
            checks.push(missing("K4", "annoy.ann", format!("{}: {e} (run scrust neighbors with --annoy yes)", c_ann.display())));
            checks.push(missing("K3", "items", "no candidate index".into()));
        }
    }

    // K5: ids and distances
    let ri = read_npy(&reference.join("stages/K5_nn_index.npy"))?;
    let rd = read_npy(&reference.join("stages/K5_nn_dist.npy"))?;
    let n = ri.nrow();
    if let Some(nc) = n_cells {
        if nc != n {
            return Err(format!("R's index holds {nc} items but its neighbour table {n} cells"));
        }
    }
    let (_, ci) = read_knn_table(&candidate.join("knn_index.tsv"))?;
    let (_, cd) = read_knn_table(&candidate.join("knn_distance.tsv"))?;
    let kk = ri.ncol();
    let (mut id_cells, mut id_vals, mut d_cells, mut d_vals, mut d_nan) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut first_id = None;
    let mut first_d = None;
    for i in 0..n {
        let (mut ci_bad, mut cd_bad) = (false, false);
        for r in 0..kk {
            let want = ri.get_i32(i, r);
            let got = ci.get(i).and_then(|row| row.get(r)).and_then(|v| v.parse::<i32>().ok());
            if got != Some(want) {
                id_vals += 1;
                ci_bad = true;
                if first_id.is_none() {
                    first_id = Some(format!("first at cell {}, position {}: R {want}, Rust {got:?}", i + 1, r + 1));
                }
            }
            let want = rd.get(i, r);
            let got = cd.get(i).and_then(|row| row.get(r)).and_then(|v| scio::parse17(v));
            // the text spells every NaN "NaN", without its sign or payload (cargo test compares
            // NaN distances bitwise, through the arrays themselves)
            let both_nan = want.is_nan() && got.is_some_and(f64::is_nan);
            d_nan += both_nan as usize;
            if !both_nan && got.map(f64::to_bits) != Some(want.to_bits()) {
                d_vals += 1;
                cd_bad = true;
                if first_d.is_none() {
                    first_d = Some(format!("first at cell {}, position {}: R {want:e}, Rust {got:?}", i + 1, r + 1));
                }
            }
        }
        id_cells += ci_bad as usize;
        d_cells += cd_bad as usize;
    }
    let shape = |rows: &Vec<Vec<String>>| rows.len() != n || rows.iter().any(|r| r.len() != kk);
    let mut d1 = format!("{n} cells x {kk}; {id_cells} cells differ");
    if shape(&ci) {
        d1.push_str(&format!("; candidate table is {} rows", ci.len()));
        id_vals = id_vals.max(1);
    }
    if let Some(x) = first_id {
        d1.push_str("; ");
        d1.push_str(&x);
    }
    checks.push(check("K5", "nn_index", n * kk, id_vals, d1));
    let mut d2 = format!("{n} cells x {kk}; {d_cells} cells differ");
    if d_nan > 0 {
        d2.push_str(&format!("; {d_nan} distances NaN in both (compared as NaN: the text keeps no NaN sign)"));
    }
    if shape(&cd) {
        d2.push_str(&format!("; candidate table is {} rows", cd.len()));
        d_vals = d_vals.max(1);
    }
    if let Some(x) = first_d {
        d2.push_str("; ");
        d2.push_str(&x);
    }
    checks.push(check("K5", "nn_dist", n * kk, d_vals, d2));

    // K6: the graphs
    for (name, file) in [("SCT_nn", "nn.mtx.gz"), ("SCT_snn", "snn.mtx.gz")] {
        let rp = reference.join("graphs").join(file);
        let cp = candidate.join(file);
        let (a, b) = match (scio::mtx::read_entries(&rp), scio::mtx::read_entries(&cp)) {
            (Ok(a), Ok(b)) => (a, b),
            (Err(e), _) | (_, Err(e)) => {
                checks.push(missing("K6", name, e.to_string()));
                continue;
            }
        };
        let mut differ = 0usize;
        let mut first = None;
        let dims_ok = (a.nrow, a.ncol) == (b.nrow, b.ncol);
        let m = a.vals.len().min(b.vals.len());
        for t in 0..m {
            if a.rows[t] != b.rows[t] || a.cols[t] != b.cols[t] || a.vals[t].to_bits() != b.vals[t].to_bits() {
                differ += 1;
                if first.is_none() {
                    first = Some(format!(
                        "first at entry {}: R ({}, {}) {}, Rust ({}, {}) {}",
                        t + 1,
                        a.rows[t] + 1,
                        a.cols[t] + 1,
                        scio::fmt17(a.vals[t]),
                        b.rows[t] + 1,
                        b.cols[t] + 1,
                        scio::fmt17(b.vals[t])
                    ));
                }
            }
        }
        differ += a.vals.len().max(b.vals.len()) - m;
        if !dims_ok {
            differ = differ.max(1);
        }
        let mut detail = format!("{} x {}, {} entries in R's, {} x {}, {} in Rust's", a.nrow, a.ncol, a.vals.len(), b.nrow, b.ncol, b.vals.len());
        if let Some(x) = first {
            detail.push_str("; ");
            detail.push_str(&x);
        }
        checks.push(check("K6", name, a.vals.len().max(b.vals.len()), differ, detail));
    }
    Ok(Report { reference: reference.to_path_buf(), candidate: candidate.to_path_buf(), checks })
}
