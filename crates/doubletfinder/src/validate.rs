//! Comparison of a Rust run with an R reference run (docs/PLAN.md A4 and B5), on plain files:
//!
//! - the final files `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv`, byte for byte (and, when
//!   they differ, how many rows, pANN values and calls differ);
//! - the `result` section of `summary.json`, field by field as text (`rng_end_sha256`, the
//!   digest of R's final `.Random.seed`, is not applicable when Rust drew no random numbers,
//!   under teacher forcing);
//!
//! and, against an instrumented R run (with `stages/`):
//!
//! - every block of every stage array R dumped, S1 to S7 (one block per array and context),
//!   against the candidate's: dtype, shape and every element, bit for bit;
//! - every `.Random.seed` checkpoint, by context and label, and their sequence;
//! - every scalar R recorded per context (`stages.json`), bit for bit; R's structural checks of
//!   its own containers do not apply to Rust;
//! - S6 and S7 teacher-forced on R's own S5: BCreal and the values passed to
//!   `bimodality_coefficient` recomputed from R's sweep pANN, and the calls from R's final pANN
//!   and nExp, so that a difference in S6 or S7 cannot hide behind one in S5;
//! - S1 from R's own random-number states: the sweep subsample and every doublet draw, redrawn
//!   with rcompat from R's `.Random.seed` before them.
//!
//! [`compare_saved`] compares a run with the saved workflow output instead: the doublet set, the
//! singlets' pANN as R wrote them, the column names (so pK and nExp), and, where they were
//! exported (M6b), the singlets' counts: barcodes, feature lines and every stored entry.
//!
//! Every criterion is exact (`validation/tolerances.toml`).

use crate::pann::classify;
use crate::stages::StageDir;
use crate::summary::bc_of_pann;
use crate::DfError;
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

/// The outcome of one check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Identical,
    Differ,
    /// The reference or the candidate lacks it (e.g. a light dump has no S4 positions).
    Missing,
    /// Not produced by this kind of run (e.g. the RNG digest under teacher forcing).
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

/// One check: a file, a summary field or an array block.
#[derive(Debug, Clone)]
pub struct Check {
    pub stage: String,
    pub name: String,
    pub status: Status,
    /// Elements compared and elements that differ (0 when not counted).
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
    /// True when nothing differs (missing and not-applicable checks are listed, not failed).
    pub fn pass(&self) -> bool {
        self.checks.iter().all(|c| c.status != Status::Differ)
    }

    /// Per stage: (checks identical, checks differing, elements compared, elements differing).
    pub fn by_stage(&self) -> Vec<(String, usize, usize, usize, usize)> {
        let mut out: Vec<(String, usize, usize, usize, usize)> = Vec::new();
        for c in &self.checks {
            if !matches!(c.status, Status::Identical | Status::Differ) {
                continue;
            }
            let i = match out.iter().position(|o| o.0 == c.stage) {
                Some(i) => i,
                None => {
                    out.push((c.stage.clone(), 0, 0, 0, 0));
                    out.len() - 1
                }
            };
            if c.status == Status::Identical {
                out[i].1 += 1;
            } else {
                out[i].2 += 1;
            }
            out[i].3 += c.n;
            out[i].4 += c.n_differ;
        }
        out
    }

    pub fn to_json(&self) -> Json {
        Json::obj(vec![
            ("reference", Json::Str(self.reference.display().to_string())),
            ("candidate", Json::Str(self.candidate.display().to_string())),
            ("pass", Json::Bool(self.pass())),
            (
                "by_stage",
                Json::List(
                    self.by_stage()
                        .into_iter()
                        .map(|(s, ok, bad, n, nd)| {
                            Json::obj(vec![
                                ("stage", Json::Str(s)),
                                ("checks_identical", Json::Int(ok as i64)),
                                ("checks_differing", Json::Int(bad as i64)),
                                ("elements", Json::Int(n as i64)),
                                ("elements_differing", Json::Int(nd as i64)),
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
                                ("stage", Json::Str(c.stage.clone())),
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
            "# DoubletFinder: Rust against R\n\n- reference: `{}`\n- candidate: `{}`\n- result: **{}**\n\n",
            self.reference.display(),
            self.candidate.display(),
            if self.pass() { "pass (everything compared is identical)" } else { "FAIL" }
        );
        s.push_str("| stage | checks identical | checks differing | elements | elements differing |\n|---|---|---|---|---|\n");
        for (st, ok, bad, n, nd) in self.by_stage() {
            s.push_str(&format!("| {st} | {ok} | {bad} | {n} | {nd} |\n"));
        }
        s.push_str("\n| stage | check | status | elements | differing | detail |\n|---|---|---|---|---|---|\n");
        for c in &self.checks {
            s.push_str(&format!("| {} | {} | {} | {} | {} | {} |\n", c.stage, c.name, c.status.as_str(), c.n, c.n_differ, c.detail));
        }
        s
    }
}

fn check(stage: &str, name: &str, status: Status, n: usize, n_differ: usize, detail: impl Into<String>) -> Check {
    Check { stage: stage.into(), name: name.into(), status, n, n_differ, detail: detail.into() }
}

fn same_bits(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// Final files, byte for byte, with row counts when they differ.
fn compare_files(reference: &Path, candidate: &Path, checks: &mut Vec<Check>) {
    for (file, stage) in [("cells.tsv", "S5/S7"), ("sweep_stats.tsv", "S6"), ("bcmvn.tsv", "S6")] {
        let (r, c) = (std::fs::read(reference.join(file)), std::fs::read(candidate.join(file)));
        let (Ok(r), Ok(c)) = (r, c) else {
            checks.push(check(stage, file, Status::Missing, 0, 0, "file absent in one of the runs"));
            continue;
        };
        if r == c {
            let rows = r.iter().filter(|&&b| b == b'\n').count().saturating_sub(1);
            checks.push(check(stage, file, Status::Identical, rows, 0, "byte-identical"));
            continue;
        }
        let (rt, ct) = (Table::parse(&String::from_utf8_lossy(&r)), Table::parse(&String::from_utf8_lossy(&c)));
        let detail = match (rt, ct) {
            (Some(rt), Some(ct)) if rt.header == ct.header && rt.rows.len() == ct.rows.len() => {
                let mut per_col = vec![0usize; rt.header.len()];
                for (a, b) in rt.rows.iter().zip(&ct.rows) {
                    for (j, (x, y)) in a.iter().zip(b).enumerate() {
                        per_col[j] += (x != y) as usize;
                    }
                }
                let n_rows = rt.rows.iter().zip(&ct.rows).filter(|(a, b)| a != b).count();
                let cols: Vec<String> = rt.header.iter().zip(&per_col).filter(|(_, &n)| n > 0).map(|(h, n)| format!("{h}: {n}")).collect();
                checks.push(check(stage, file, Status::Differ, rt.rows.len(), n_rows, format!("rows differ; by column {}", cols.join(", "))));
                continue;
            }
            _ => "different headers or row counts".to_string(),
        };
        checks.push(check(stage, file, Status::Differ, 0, 0, detail));
    }
}

/// The `result` sections of two `summary.json` files, field by field.
fn compare_results(reference: &Path, candidate: &Path, checks: &mut Vec<Check>) -> Option<Json> {
    let read = |p: &Path| -> Option<Json> { Json::parse(&std::fs::read_to_string(p.join("summary.json")).ok()?).ok()?.get("result").cloned() };
    let (Some(r), Some(c)) = (read(reference), read(candidate)) else {
        checks.push(check("result", "summary.json", Status::Missing, 0, 0, "no result section in one of the runs"));
        return None;
    };
    let Json::Object(fields) = &r else { return None };
    for (key, rv) in fields {
        let cv = c.get(key);
        let (status, detail) = match cv {
            None => (Status::Missing, "absent from the candidate".to_string()),
            Some(Json::Null) if key == "rng_end_sha256" => {
                (Status::NotApplicable, "no random numbers drawn (teacher forcing)".to_string())
            }
            Some(cv) => {
                let (rt, ct) = (rv.to_text(0), cv.to_text(0));
                if rt == ct {
                    (Status::Identical, rt)
                } else {
                    (Status::Differ, format!("R {rt}, Rust {ct}"))
                }
            }
        };
        checks.push(check("result", key, status, 1, (status == Status::Differ) as usize, detail));
    }
    Some(r)
}

/// The array name of a stage file: `S2b_glm.npy` -> `S2b_glm`, `S2d_scale_data.final.npy` ->
/// `S2d_scale_data`.
fn array_name(file: &str, ctx: &str) -> String {
    let f = file.trim_end_matches(".npy");
    f.strip_suffix(&format!(".{ctx}")).unwrap_or(f).to_string()
}

/// How one block of a stage array compares.
enum Outcome {
    /// Elements compared, elements differing.
    Compared(usize, usize),
    /// One of the runs does not hold it.
    Missing(&'static str),
    /// Different dtype or shape; the block's size.
    Mismatch(String, usize),
}

fn compare_block(rb: &crate::stages::Block, ra: Option<&Npy>, cb: Option<&crate::stages::Block>, ca: Option<&Npy>) -> Outcome {
    let Some(ra) = ra else { return Outcome::Missing("not in the reference's directory") };
    let (Some(cb), Some(ca)) = (cb, ca) else { return Outcome::Missing("not in the candidate's stages") };
    let size = rb.nrow * rb.ncol;
    if cb.dtype != rb.dtype {
        return Outcome::Mismatch(format!("dtype {} against R's {}", cb.dtype, rb.dtype), size);
    }
    if (cb.nrow, cb.ncol) != (rb.nrow, rb.ncol) {
        return Outcome::Mismatch(format!("shape {}x{} against R's {}x{}", cb.nrow, cb.ncol, rb.nrow, rb.ncol), size);
    }
    let mut nd = 0usize;
    for col in 0..rb.ncol {
        for row in 0..rb.nrow {
            nd += !same_bits(ra.get(rb.row_start + row, col), ca.get(cb.row_start + row, col)) as usize;
        }
    }
    Outcome::Compared(size, nd)
}

/// The summary contexts of one array (one per pN and pK), gathered into one check.
#[derive(Default)]
struct Gathered {
    stage: String,
    name: String,
    blocks: usize,
    missing: usize,
    missing_why: &'static str,
    n: usize,
    n_differ: usize,
    bad: Vec<String>,
}

/// Every block of every stage array the reference run dumped (S1 to S7), against the
/// candidate's block of the same array and context: dtype, shape and every element, bit for
/// bit. One check per array and context; S6's summary contexts (one per pN and pK) are gathered
/// in one check per array. Files the reference lists but does not hold (the tier-T fixtures
/// leave out the per-context matrices) are reported missing, not compared.
fn compare_stage_blocks(rs: &StageDir, cs: &StageDir, checks: &mut Vec<Check>) -> Result<(), DfError> {
    let mut files: Vec<String> = rs.blocks.iter().map(|b| b.file.clone()).collect();
    files.sort();
    files.dedup();
    for file in &files {
        let ra = if rs.dir.join(file).is_file() { Some(rs.read(file)?) } else { None };
        let ca = if cs.dir.join(file).is_file() { Some(cs.read(file)?) } else { None };
        let mut agg: Option<Gathered> = None;
        for rb in rs.blocks_of(file) {
            let array = array_name(file, &rb.ctx);
            let outcome = compare_block(rb, ra.as_ref(), cs.block(file, &rb.ctx), ca.as_ref());
            if rb.ctx.starts_with("summary:") {
                let a = agg.get_or_insert_with(|| Gathered { stage: rb.stage.clone(), name: array.clone(), ..Gathered::default() });
                a.blocks += 1;
                match outcome {
                    Outcome::Compared(n, nd) => {
                        a.n += n;
                        a.n_differ += nd;
                        if nd > 0 {
                            a.bad.push(format!("{} ({nd})", rb.ctx));
                        }
                    }
                    Outcome::Missing(why) => {
                        a.missing += 1;
                        a.missing_why = why;
                    }
                    Outcome::Mismatch(why, n) => {
                        a.n += n;
                        a.n_differ += n;
                        a.bad.push(format!("{} ({why})", rb.ctx));
                    }
                }
                continue;
            }
            let name = format!("{}: {array}", rb.ctx);
            checks.push(match outcome {
                Outcome::Compared(n, 0) => check(&rb.stage, &name, Status::Identical, n, 0, ""),
                Outcome::Compared(n, nd) => check(&rb.stage, &name, Status::Differ, n, nd, format!("{nd} of {n} elements differ")),
                Outcome::Missing(why) => check(&rb.stage, &name, Status::Missing, 0, 0, why),
                Outcome::Mismatch(why, n) => check(&rb.stage, &name, Status::Differ, n, n, why),
            });
        }
        if let Some(a) = agg {
            let status = if a.n_differ > 0 {
                Status::Differ
            } else if a.missing == a.blocks {
                Status::Missing
            } else {
                Status::Identical
            };
            let mut detail = format!("{} (pN, pK) blocks compared", a.blocks - a.missing);
            if a.missing > 0 {
                detail.push_str(&format!("; {} {}", a.missing, a.missing_why));
            }
            if !a.bad.is_empty() {
                detail.push_str(&format!("; differing: {}", a.bad.join(", ")));
            }
            checks.push(check(&a.stage, &a.name, status, a.n, a.n_differ, detail));
        }
    }
    Ok(())
}

/// Every `.Random.seed` checkpoint of the reference (`rng.npy`, `rng.tsv`) against the
/// candidate's checkpoint of the same context and label (one check each, under the stage R
/// recorded it in), and whether the two runs recorded the same checkpoints in the same order.
fn compare_rng(rs: &StageDir, cs: &StageDir, checks: &mut Vec<Check>) -> Result<(), DfError> {
    let read = |d: &StageDir| -> Option<(Npy, Table)> {
        let t = Table::read(&d.dir.join("rng.tsv")).ok()?;
        Some((d.read("rng.npy").ok()?, t))
    };
    let Some((rn, rt)) = read(rs) else { return Ok(()) };
    let Some((cn, ct)) = read(cs) else {
        checks.push(check("RNG", "RNG checkpoints", Status::Missing, 0, 0, "the candidate recorded none"));
        return Ok(());
    };
    let state = |a: &Npy, row: usize| -> Vec<i32> { (0..a.ncol()).map(|c| a.get_i32(row, c)).collect() };
    let row_of = |t: &Table, ctx: &str, label: &str| -> Option<usize> { t.rows.iter().find(|r| r[1] == ctx && r[3] == label).and_then(|r| r[0].parse().ok()) };
    for r in &rt.rows {
        let (ctx, stage, label) = (&r[1], &r[2], &r[3]);
        let name = format!("{ctx}: RNG {label}");
        let rrow: usize = r[0].parse().map_err(|_| DfError::Input(format!("rng.tsv: bad row {}", r[0])))?;
        match row_of(&ct, ctx, label) {
            None => checks.push(check(stage, &name, Status::Missing, 0, 0, "not among the candidate's checkpoints")),
            Some(crow) => {
                let same = state(&rn, rrow) == state(&cn, crow);
                checks.push(check(stage, &name, if same { Status::Identical } else { Status::Differ }, 1, (!same) as usize, ""));
            }
        }
    }
    let seq = |t: &Table| -> Vec<(String, String)> { t.rows.iter().map(|r| (r[1].clone(), r[3].clone())).collect() };
    let (a, b) = (seq(&rt), seq(&ct));
    let same = a == b;
    let detail = if same { format!("{} checkpoints", a.len()) } else { format!("R {} checkpoints, Rust {}; first difference at {}", a.len(), b.len(), a.iter().zip(&b).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()))) };
    checks.push(check("RNG", "checkpoint sequence", if same { Status::Identical } else { Status::Differ }, 1, (!same) as usize, detail));
    Ok(())
}

/// R's structural checks of its own containers (always TRUE in the reference runs), which hold
/// by construction in Rust and are not recorded there.
pub const R_ONLY_SCALARS: [&str; 5] = [
    "vst_rows_are_input_features",
    "scale_data_rows_are_top_features",
    "variable_features_equal_top_features",
    "scale_data_equals_sctransform_default_y",
    "pca_cells_are_vst_cells",
];

/// The stage a scalar of `stages.json` belongs to.
pub fn scalar_stage(key: &str) -> &'static str {
    match key {
        "pN" | "n_real" | "subsampled" | "n_doublets" => "S1",
        "vst_n_cells" | "density_bw" | "density_n" | "n_genes" | "n_genes_step1" | "n_cells_step1" | "n_genes_step1_pool" | "gene_sampling" | "n_cells" => "S2a",
        "n_bins_fit" | "n_theta_below_1e7" | "n_theta_inf" => "S2b",
        "reg_bw" | "mean_cell_sum" | "n_outliers" | "n_all_poisson_genes" | "n_smoothing_genes" => "S2c",
        "min_variance" | "res_clip_range" | "clip_range" => "S2d",
        "pca_n_features" | "total_variance" | "irlba_iter" | "irlba_mprod" => "S3",
        "nCells" | "ind" | "k" | "first_position_is_self" => "S4",
        "pK" | "pk_temp" => "S5",
        "nExp" => "S7",
        _ => "scalars",
    }
}

/// A JSON value's scalars in order (a scalar is a vector of length one, as in R's `to_json`):
/// numbers as doubles, the rest as text.
fn atoms(v: &Json) -> Vec<Result<f64, String>> {
    match v {
        Json::Vector(items) | Json::List(items) => items.iter().flat_map(atoms).collect(),
        Json::Int(_) | Json::Num(_) => vec![Ok(v.as_f64().unwrap())],
        other => vec![Err(other.to_text(0))],
    }
}

/// Every scalar of the reference's `stages.json`, per context, against the candidate's, bit for
/// bit (R wrote them with `%.17g`, which reads back exactly).
fn compare_scalars(rs: &StageDir, cs: &StageDir, checks: &mut Vec<Check>) -> Result<(), DfError> {
    let read = |d: &StageDir| -> Option<Json> { Json::parse(&std::fs::read_to_string(d.dir.join("stages.json")).ok()?).ok()?.get("scalars").cloned() };
    let Some(Json::Object(r)) = read(rs) else { return Ok(()) };
    let Some(c) = read(cs) else {
        checks.push(check("scalars", "stages.json scalars", Status::Missing, 0, 0, "the candidate recorded none"));
        return Ok(());
    };
    for (ctx, members) in &r {
        let Json::Object(members) = members else { continue };
        for (key, rv) in members {
            let name = format!("{ctx}: {key}");
            let stage = scalar_stage(key);
            if R_ONLY_SCALARS.contains(&key.as_str()) {
                checks.push(check(stage, &name, Status::NotApplicable, 0, 0, format!("R's own structural check ({})", rv.to_text(0))));
                continue;
            }
            match c.get(ctx).and_then(|m| m.get(key)) {
                None => checks.push(check(stage, &name, Status::Missing, 0, 0, "not recorded by the candidate")),
                Some(cv) => {
                    let (a, b) = (atoms(rv), atoms(cv));
                    let nd = if a.len() != b.len() {
                        a.len().max(1)
                    } else {
                        a.iter()
                            .zip(&b)
                            .filter(|(x, y)| match (x, y) {
                                (Ok(x), Ok(y)) => !same_bits(*x, *y),
                                (Err(x), Err(y)) => x != y,
                                _ => true,
                            })
                            .count()
                    };
                    let detail = if nd == 0 { String::new() } else { format!("R {}, Rust {}", rv.to_text(0), cv.to_text(0)) };
                    checks.push(check(stage, &name, if nd == 0 { Status::Identical } else { Status::Differ }, a.len(), nd, detail));
                }
            }
        }
    }
    Ok(())
}

/// S6 and S7 recomputed from R's own S5 arrays (teacher forcing), against R's S6 and S7.
fn teacher_forced_s6_s7(reference: &Path, rs: &StageDir, result: Option<&Json>, checks: &mut Vec<Check>, threads: usize) -> Result<(), DfError> {
    // S6: every (pN, pK) of the sweep, in sweep order
    let s5 = rs.read("S5_pann_sweep.npy").ok();
    let gk = rs.read("S6_gkde.npy").ok();
    let sweep_stats = Table::read(&reference.join("sweep_stats.tsv")).ok();
    if let (Some(s5), Some(gk), Some(stats)) = (s5, gk, sweep_stats) {
        let s5_blocks = rs.blocks_of("S5_pann_sweep.npy");
        let g_blocks = rs.blocks_of("S6_gkde.npy");
        let npk = s5.ncol();
        if g_blocks.len() != s5_blocks.len() * npk || stats.rows.len() != g_blocks.len() {
            checks.push(check("S6", "S6 (teacher-forced on R's S5)", Status::Differ, 0, 0, "sweep layout differs"));
        } else {
            let items: Vec<usize> = (0..g_blocks.len()).collect();
            let compute = |&item: &usize| -> (usize, usize, bool) {
                let b = s5_blocks[item / npk];
                let pann = s5.col_f64(item % npk, b.row_start, b.nrow);
                let gb = g_blocks[item];
                match bc_of_pann(&pann) {
                    Ok((bc, _, g)) => {
                        let want = gk.col_f64(0, gb.row_start, gb.nrow);
                        let d = g.iter().zip(&want).filter(|(a, b)| !same_bits(**a, **b)).count() + g.len().abs_diff(want.len());
                        let bc_r = scio::parse17(&stats.rows[item][2]).unwrap_or(f64::NAN);
                        (want.len(), d, same_bits(bc, bc_r))
                    }
                    Err(_) => (gb.nrow, gb.nrow, false),
                }
            };
            let results: Vec<(usize, usize, bool)> = crate::par_map(&items, threads, compute);
            let (n, nd) = results.iter().fold((0, 0), |(a, b), r| (a + r.0, b + r.1));
            let bc_bad = results.iter().filter(|r| !r.2).count();
            checks.push(check("S6", "S6_gkde (teacher-forced on R's S5)", if nd == 0 { Status::Identical } else { Status::Differ }, n, nd, format!("{} (pN, pK) pairs", results.len())));
            checks.push(check("S6", "BCreal (teacher-forced on R's S5)", if bc_bad == 0 { Status::Identical } else { Status::Differ }, results.len(), bc_bad, "against sweep_stats.tsv"));
        }
    } else {
        checks.push(check("S6", "S6 (teacher-forced on R's S5)", Status::Missing, 0, 0, "R's S5 sweep, S6_gkde or sweep_stats.tsv not found"));
    }
    // S7: the calls from R's final pANN and nExp
    let pann = rs.read("S5_pann_final.npy").ok();
    let calls = rs.read("S7_doublet.npy").ok();
    let nexp = result.and_then(|r| r.get("nExp")).and_then(Json::as_f64);
    match (pann, calls, nexp) {
        (Some(p), Some(c), Some(nexp)) => {
            let pv = p.f64s();
            let got = classify(&pv, nexp)?;
            let want = c.f64s();
            let nd = got.iter().zip(&want).filter(|(g, w)| (**g as i32 as f64) != **w).count();
            checks.push(check("S7", "S7_doublet (teacher-forced on R's S5)", if nd == 0 { Status::Identical } else { Status::Differ }, want.len(), nd, format!("nExp = {nexp}")));
        }
        _ => checks.push(check("S7", "S7_doublet (teacher-forced on R's S5)", Status::Missing, 0, 0, "R's S5_pann_final, S7_doublet or nExp not found")),
    }
    Ok(())
}

/// Stage S1 from R's own `.Random.seed` checkpoints: the sweep subsample (`sample(1:n, 10000)`
/// from the state at `paramSweep` entry, when n > 10,000) and, for every sweep iteration and the
/// final run, the two `sample(real.cells, n_doublets, replace = TRUE)` draws from the state
/// before them; the draws must equal R's (S1_subsample, S1_doublets) and leave R's state after
/// them. This checks rcompat's sampling on DoubletFinder's own draws; Rust's teacher-forced runs
/// take the subsample from R.
fn s1_from_rng(rs: &StageDir, n_cells: Option<usize>, checks: &mut Vec<Check>) -> Result<(), DfError> {
    let (Ok(rng), Ok(labels), Ok(doublets)) = (rs.read("rng.npy"), Table::read(&rs.dir.join("rng.tsv")), rs.read("S1_doublets.npy")) else {
        checks.push(check("S1", "S1 (from R's RNG checkpoints)", Status::Missing, 0, 0, "no rng.npy, rng.tsv or S1_doublets.npy"));
        return Ok(());
    };
    let n_cells = match n_cells.or_else(|| rs.block("S5_pann_final.npy", "final").map(|b| b.nrow)) {
        Some(n) => n,
        None => return Ok(()),
    };
    let state = |row: usize| -> Vec<i32> { (0..rng.ncol()).map(|c| rng.get_i32(row, c)).collect() };
    let find = |ctx: &str, label: &str| -> Option<usize> {
        labels.rows.iter().find(|r| r[1] == ctx && r[3] == label).and_then(|r| r[0].parse().ok())
    };
    let rng_at = |row: usize| rcompat::RRng::from_random_seed(&state(row)).map_err(|e| DfError::Input(format!("rng.npy row {row}: {e:?}")));
    // the sweep subsample
    let first_ctx = rs.blocks_of("S1_doublets.npy").first().map(|b| b.ctx.clone()).unwrap_or_default();
    let (entry, first) = (find("sweep", "paramSweep_entry"), find(&first_ctx, "before_doublets"));
    let mut real = n_cells;
    if let (Some(entry), Some(first)) = (entry, first) {
        let mut g = rng_at(entry)?;
        let (mut n, mut nd) = (0usize, 0usize);
        if let Some(b) = rs.block("S1_subsample.npy", "sweep") {
            let want = rs.read("S1_subsample.npy")?.col_i32(0, b.row_start, b.nrow);
            let got = rcompat::sample_int(&mut g, n_cells, want.len(), false).map_err(|e| DfError::Input(format!("{e:?}")))?;
            n = want.len();
            nd = got.iter().zip(&want).filter(|(a, b)| **a as i32 != **b).count();
            real = want.len();
        }
        let same_state = g.random_seed() == state(first);
        let status = if nd == 0 && same_state { Status::Identical } else { Status::Differ };
        let detail = format!("{}; RNG state afterwards {}", if n > 0 { format!("{n} cells drawn") } else { "no subsample (n <= 10,000)".into() }, if same_state { "identical" } else { "DIFFERS" });
        checks.push(check("S1", "S1_subsample (from R's RNG state)", status, n.max(1), nd + (!same_state) as usize, detail));
    }
    // the doublet pairs of every context
    let (mut n, mut nd, mut states_bad, mut blocks) = (0usize, 0usize, 0usize, 0usize);
    for b in rs.blocks_of("S1_doublets.npy") {
        let (Some(before), Some(after)) = (find(&b.ctx, "before_doublets"), find(&b.ctx, "after_doublets")) else { continue };
        let n_real = if b.ctx == "final" { n_cells } else { real };
        let mut g = rng_at(before)?;
        let d1 = rcompat::sample_int(&mut g, n_real, b.nrow, true).map_err(|e| DfError::Input(format!("{e:?}")))?;
        let d2 = rcompat::sample_int(&mut g, n_real, b.nrow, true).map_err(|e| DfError::Input(format!("{e:?}")))?;
        let (w1, w2) = (doublets.col_i32(0, b.row_start, b.nrow), doublets.col_i32(1, b.row_start, b.nrow));
        n += 2 * b.nrow;
        nd += d1.iter().zip(&w1).chain(d2.iter().zip(&w2)).filter(|(a, b)| **a as i32 != **b).count();
        states_bad += (g.random_seed() != state(after)) as usize;
        blocks += 1;
    }
    if blocks > 0 {
        let status = if nd == 0 && states_bad == 0 { Status::Identical } else { Status::Differ };
        checks.push(check("S1", "S1_doublets (from R's RNG states)", status, n, nd + states_bad,
            format!("{blocks} contexts; RNG state after the draws differs in {states_bad}")));
    }
    Ok(())
}

/// The OpenBLAS kernel each run recorded (`run.fingerprint.openblas_corename` in `summary.json`):
/// the BLAS-relevant field of the environment fingerprint (plan A5b). The exact criteria for
/// SCTransform's residuals and the PCA hold under the same kernel only.
fn compare_kernel(reference: &Path, candidate: &Path, checks: &mut Vec<Check>) {
    let core = |p: &Path| -> Option<String> {
        let j = Json::parse(&std::fs::read_to_string(p.join("summary.json")).ok()?).ok()?;
        j.get("run")?.get("fingerprint")?.get("openblas_corename")?.as_str().map(String::from)
    };
    match (core(reference), core(candidate)) {
        (Some(r), Some(c)) => {
            let same = r == c;
            checks.push(check("environment", "OpenBLAS kernel", if same { Status::Identical } else { Status::Differ }, 1, (!same) as usize, format!("R {r}, Rust {c}")));
        }
        _ => checks.push(check("environment", "OpenBLAS kernel", Status::Missing, 0, 0, "not recorded by one of the runs")),
    }
}

/// Compares the candidate (a Rust output directory: the final files, `summary.json` and
/// `stages/`) with an R reference run: its directory (the final files and `summary.json`) and,
/// for an instrumented run, its stage directory (usually `<reference>/stages`). A direct run has
/// no stage directory: then only the files and the `result` section are compared.
pub fn compare(reference: &Path, reference_stages: &Path, candidate: &Path, threads: usize) -> Result<Report, DfError> {
    let mut checks = Vec::new();
    compare_files(reference, candidate, &mut checks);
    let result = compare_results(reference, candidate, &mut checks);
    compare_kernel(reference, candidate, &mut checks);
    let (rsd, csd) = (reference_stages.to_path_buf(), candidate.join("stages"));
    if !rsd.join("index.tsv").exists() {
        checks.push(check("stages", "stage dumps", Status::NotApplicable, 0, 0, "the reference has none (a direct run)"));
    } else if !csd.join("index.tsv").exists() {
        checks.push(check("stages", "stage dumps", Status::Missing, 0, 0, "the candidate has none"));
    } else {
        let (rs, cs) = (StageDir::open(&rsd)?, StageDir::open(&csd)?);
        compare_stage_blocks(&rs, &cs, &mut checks)?;
        compare_rng(&rs, &cs, &mut checks)?;
        compare_scalars(&rs, &cs, &mut checks)?;
        teacher_forced_s6_s7(reference, &rs, result.as_ref(), &mut checks, threads)?;
        let n_cells = result.as_ref().and_then(|r| r.get("n_cells")).and_then(Json::as_f64).map(|v| v as usize);
        s1_from_rng(&rs, n_cells, &mut checks)?;
    }
    Ok(Report { reference: reference.to_path_buf(), candidate: candidate.to_path_buf(), checks })
}

/// Compares a run's outputs (`cells.tsv` and the `result` section of `summary.json`) with the
/// saved workflow output as exported (`saved/`: `cells.tsv` with the kept cells and their
/// pANN as R wrote them, `doublets.txt`, `saved.json` with the column names, pK and nExp):
/// plan B5's row for saved outputs. The kept cells' pANN are compared as text (17 significant
/// digits, so bitwise). With `saved/singlets/` (the saved object's RNA counts, M6b), the run's
/// `singlets/` too ([`compare_singlets`]).
pub fn compare_saved(saved: &Path, candidate: &Path) -> Result<Report, DfError> {
    let stage = "saved";
    let mut checks = Vec::new();
    let read_json = |p: &Path| -> Result<Json, DfError> {
        let text = std::fs::read_to_string(p).map_err(|e| DfError::Input(format!("{}: {e}", p.display())))?;
        Json::parse(&text).map_err(|e| DfError::Input(format!("{}: {e}", p.display())))
    };
    let info = read_json(&saved.join("saved.json"))?;
    let result = read_json(&candidate.join("summary.json"))?.get("result").cloned().ok_or_else(|| DfError::Input("summary.json: no result section".into()))?;
    for (field, saved_field) in [("pANN_column", "pANN_column"), ("classification_column", "classification_column"), ("pK_label", "pK_label")] {
        let (a, b) = (result.get(field).and_then(Json::as_str), info.get(saved_field).and_then(Json::as_str));
        let same = a.is_some() && a == b;
        checks.push(check(stage, field, if same { Status::Identical } else { Status::Differ }, 1, (!same) as usize, format!("saved {}, run {}", b.unwrap_or("-"), a.unwrap_or("-"))));
    }
    let (a, b) = (result.get("nExp").and_then(Json::as_f64), info.get("nExp").and_then(Json::as_f64));
    let same = matches!((a, b), (Some(x), Some(y)) if same_bits(x, y));
    checks.push(check(stage, "nExp", if same { Status::Identical } else { Status::Differ }, 1, (!same) as usize, format!("saved {b:?}, run {a:?}")));

    let cells = Table::read(&candidate.join("cells.tsv")).map_err(DfError::Io)?;
    let (bc, pa, cl) = (cells.col("barcode"), cells.col("pANN"), cells.col("classification"));
    let (Some(bc), Some(pa), Some(cl)) = (bc, pa, cl) else {
        return Err(DfError::Input("cells.tsv: expected barcode, pANN and classification columns".into()));
    };
    let run_doublets: std::collections::BTreeSet<&str> = cells.rows.iter().filter(|r| r[cl] == "Doublet").map(|r| r[bc].as_str()).collect();
    let saved_doublets_text = std::fs::read_to_string(saved.join("doublets.txt")).map_err(DfError::Io)?;
    let saved_doublets: std::collections::BTreeSet<&str> = saved_doublets_text.lines().filter(|l| !l.is_empty()).collect();
    let nd = run_doublets.symmetric_difference(&saved_doublets).count();
    checks.push(check(stage, "doublet set", if nd == 0 { Status::Identical } else { Status::Differ }, saved_doublets.len(), nd, format!("{} saved, {} called", saved_doublets.len(), run_doublets.len())));

    let kept = Table::read(&saved.join("cells.tsv")).map_err(DfError::Io)?;
    let (Some(kb), Some(kp)) = (kept.col("barcode"), kept.col("pANN")) else {
        return Err(DfError::Input("saved cells.tsv: expected barcode and pANN columns".into()));
    };
    let run_singlets: Vec<(&str, &str)> = cells.rows.iter().filter(|r| r[cl] == "Singlet").map(|r| (r[bc].as_str(), r[pa].as_str())).collect();
    let saved_singlets: Vec<(&str, &str)> = kept.rows.iter().map(|r| (r[kb].as_str(), r[kp].as_str())).collect();
    let nd = if run_singlets.len() == saved_singlets.len() {
        run_singlets.iter().zip(&saved_singlets).filter(|(a, b)| a != b).count()
    } else {
        saved_singlets.len().max(run_singlets.len())
    };
    checks.push(check(stage, "singlet pANN", if nd == 0 { Status::Identical } else { Status::Differ }, saved_singlets.len(), nd, format!("{} kept cells, in order, barcode and pANN text", saved_singlets.len())));
    compare_singlets(&saved.join(crate::singlets::DIR), &candidate.join(crate::singlets::DIR), &mut checks)?;
    Ok(Report { reference: saved.to_path_buf(), candidate: candidate.to_path_buf(), checks })
}

/// The run's singlets (`singlets/`) against the saved object's RNA counts as exported
/// (`saved/singlets/`, `validation/R/export_doubletfinder_singlets.R`): the barcodes in order,
/// the feature lines, and the matrix: its header field, dimensions, every stored entry (explicit
/// zeros included) and every value, bit for bit. Not applicable when they were not exported.
pub fn compare_singlets(saved: &Path, candidate: &Path, checks: &mut Vec<Check>) -> Result<(), DfError> {
    use scio::tenx;
    let stage = "saved";
    if tenx::find(saved, "matrix.mtx").is_err() {
        checks.push(check(stage, "singlets", Status::NotApplicable, 0, 0, "the saved output's singlets were not exported"));
        return Ok(());
    }
    if tenx::find(candidate, "matrix.mtx").is_err() {
        checks.push(check(stage, "singlets", Status::Missing, 0, 0, "the run wrote no singlets/"));
        return Ok(());
    }
    let differ = |n: usize| if n == 0 { Status::Identical } else { Status::Differ };
    let in_order = |a: &[String], b: &[String]| if a.len() == b.len() { a.iter().zip(b).filter(|(x, y)| x != y).count() } else { a.len().max(b.len()) };
    let (rb, cb) = (tenx::read_barcodes(saved)?, tenx::read_barcodes(candidate)?);
    let nd = in_order(&rb, &cb);
    checks.push(check(stage, "singlets barcodes", differ(nd), rb.len(), nd, format!("{} saved, {} written, in order", rb.len(), cb.len())));
    let lines = |d: &Path| -> Result<Vec<String>, DfError> { Ok(tenx::read_lines(&tenx::find(d, "features.tsv")?)?) };
    let (rf, cf) = (lines(saved)?, lines(candidate)?);
    let nd = in_order(&rf, &cf);
    checks.push(check(stage, "singlets features", differ(nd), rf.len(), nd, format!("{} saved, {} written, lines in order", rf.len(), cf.len())));
    let (rh, ch) = (tenx::read_mtx_header(saved)?, tenx::read_mtx_header(candidate)?);
    let (rm, cm) = (scio::sparse::read_mtx(saved)?, scio::sparse::read_mtx(candidate)?);
    let same_entries = rm.nrow == cm.nrow && rm.ncol == cm.ncol && rm.p == cm.p && rm.i == cm.i;
    let nd = if same_entries { rm.x.iter().zip(&cm.x).filter(|(a, b)| !same_bits(**a, **b)).count() } else { rm.nnz().max(cm.nnz()) };
    let nd = nd + (rh.field != ch.field) as usize;
    let zeros = |x: &[f64]| x.iter().filter(|&&v| v == 0.0).count();
    let detail = format!(
        "saved {} x {}, {} entries ({} explicit zeros), {}; written {} x {}, {} entries ({} explicit zeros), {}{}",
        rm.nrow,
        rm.ncol,
        rm.nnz(),
        zeros(&rm.x),
        rh.field,
        cm.nrow,
        cm.ncol,
        cm.nnz(),
        zeros(&cm.x),
        ch.field,
        if same_entries { "" } else { "; the stored entries differ" }
    );
    checks.push(check(stage, "singlets matrix", differ(nd), rm.nnz(), nd, detail));
    Ok(())
}

/// Reads a `.npy` file (for callers outside this crate).
pub fn read_npy(path: &Path) -> Result<Npy, DfError> {
    Npy::read(path).map_err(|e| DfError::Input(e.to_string()))
}
