//! Comparison of a Rust run with an R reference run (docs/PLAN.md A4 and B5), on plain files:
//!
//! - the final files `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv`, byte for byte (and, when
//!   they differ, how many rows, pANN values and calls differ);
//! - the `result` section of `summary.json`, field by field as text (`rng_end_sha256`, the
//!   digest of R's final `.Random.seed`, is not applicable when Rust drew no random numbers);
//! - every block of the stage arrays S4-S7 that both runs dumped, bit for bit;
//! - S6 and S7 teacher-forced on R's own S5: BCreal and the values passed to
//!   `bimodality_coefficient` recomputed from R's sweep pANN, and the calls from R's final pANN
//!   and nExp, so that a difference in S6 or S7 cannot hide behind one in S5;
//! - S1 from R's own random-number states: the sweep subsample and every doublet draw, redrawn
//!   with rcompat from R's `.Random.seed` before them.
//!
//! Every criterion is exact (`validation/tolerances.toml`, S4-S7).

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

/// Every S4-S7 block that both stage directories have, bit for bit.
fn compare_stage_blocks(rs: &StageDir, cs: &StageDir, checks: &mut Vec<Check>) -> Result<(), DfError> {
    let mut files: Vec<&str> = vec!["S4_digest.npy", "S5_pann_sweep.npy", "S5_pann_final.npy", "S6_bkde.npy", "S6_gkde.npy", "S7_doublet.npy"];
    let positions: Vec<String> = rs.blocks.iter().filter(|b| b.file.starts_with("S4_positions.")).map(|b| b.file.clone()).collect();
    files.extend(positions.iter().map(String::as_str));
    for file in files {
        let rblocks = rs.blocks_of(file);
        if rblocks.is_empty() {
            continue;
        }
        let stage = &rblocks[0].stage;
        let (ra, ca) = (rs.read(file).ok(), cs.read(file).ok());
        let (Some(ra), Some(ca)) = (ra, ca) else {
            let which = if rs.dir.join(file).exists() { "candidate" } else { "reference" };
            checks.push(check(stage, file, Status::Missing, 0, 0, format!("not in the {which} run's stages")));
            continue;
        };
        // one check per array, with differing blocks listed
        let (mut n, mut nd) = (0usize, 0usize);
        let mut bad_blocks: Vec<String> = Vec::new();
        let mut missing: Vec<String> = Vec::new();
        for rb in &rblocks {
            let Some(cb) = cs.block(file, &rb.ctx) else {
                missing.push(rb.ctx.clone());
                continue;
            };
            if (cb.nrow, cb.ncol) != (rb.nrow, rb.ncol) {
                bad_blocks.push(format!("{} (shape {}x{} vs {}x{})", rb.ctx, cb.nrow, cb.ncol, rb.nrow, rb.ncol));
                n += rb.nrow * rb.ncol;
                nd += rb.nrow * rb.ncol;
                continue;
            }
            let mut block_diff = 0;
            for col in 0..rb.ncol {
                for row in 0..rb.nrow {
                    let (x, y) = (ra.get(rb.row_start + row, col), ca.get(cb.row_start + row, col));
                    block_diff += !same_bits(x, y) as usize;
                }
            }
            n += rb.nrow * rb.ncol;
            nd += block_diff;
            if block_diff > 0 {
                bad_blocks.push(format!("{} ({block_diff})", rb.ctx));
            }
        }
        let status = if nd > 0 { Status::Differ } else if n == 0 { Status::Missing } else { Status::Identical };
        let mut detail = format!("{} blocks", rblocks.len() - missing.len());
        if !bad_blocks.is_empty() {
            detail.push_str(&format!("; differing: {}", bad_blocks.join(", ")));
        }
        if !missing.is_empty() {
            detail.push_str(&format!("; not in the candidate: {}", missing.join(", ")));
        }
        checks.push(check(stage, file, status, n, nd, detail));
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

/// Compares the candidate (a Rust output directory: the final files, `summary.json` and
/// `stages/`) with an R reference run: its directory (the final files and `summary.json`) and
/// the stage directory of an instrumented run (usually `<reference>/stages`).
pub fn compare(reference: &Path, reference_stages: &Path, candidate: &Path, threads: usize) -> Result<Report, DfError> {
    let mut checks = Vec::new();
    compare_files(reference, candidate, &mut checks);
    let result = compare_results(reference, candidate, &mut checks);
    let (rsd, csd) = (reference_stages.to_path_buf(), candidate.join("stages"));
    if rsd.join("index.tsv").exists() && csd.join("index.tsv").exists() {
        let (rs, cs) = (StageDir::open(&rsd)?, StageDir::open(&csd)?);
        compare_stage_blocks(&rs, &cs, &mut checks)?;
        teacher_forced_s6_s7(reference, &rs, result.as_ref(), &mut checks, threads)?;
        let n_cells = result.as_ref().and_then(|r| r.get("n_cells")).and_then(Json::as_f64).map(|v| v as usize);
        s1_from_rng(&rs, n_cells, &mut checks)?;
    } else {
        checks.push(check("S4-S7", "stages", Status::Missing, 0, 0, "no stage dumps in one of the runs"));
    }
    Ok(Report { reference: reference.to_path_buf(), candidate: candidate.to_path_buf(), checks })
}

/// Reads a `.npy` file (for callers outside this crate).
pub fn read_npy(path: &Path) -> Result<Npy, DfError> {
    Npy::read(path).map_err(|e| DfError::Input(e.to_string()))
}
