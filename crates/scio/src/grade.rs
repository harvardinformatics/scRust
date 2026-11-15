//! Kernel-gated grading (docs/PLAN.md A5b, mechanism 3; milestone M11), shared by
//! `scrust validate cluster` and `scrust validate doubletfinder`.
//!
//! **The mode.** Exact checks run only when the reference's and the candidate's recorded
//! OpenBLAS kernels match (`run.fingerprint.openblas_corename` in each run's `summary.json`):
//! OpenBLAS picks its kernel per CPU, and BLAS products differ in their last bits between kernels.
//! Otherwise, and whenever a run records no kernel, the comparison runs in **tolerance mode**:
//!
//! - each **BLAS-free** check stays exact: none of its inputs passed through a BLAS product in
//!   the candidate's run, so the kernel cannot move it;
//! - each **BLAS-dependent** check is graded against the R-vs-R floor: R itself, run under the
//!   two kernels on the same CPU type (`[<tool>.floor.kernel]` in `validation/tolerances.toml`).
//!   A check whose array the floor lists passes when it stays within every bound of the entry;
//!   a check the floor does not list stays exact, because the floor found it identical across
//!   kernels.
//!
//! **Bounds.** In a floor entry, the keys starting with `max_` are bounds (`max_abs`, the largest
//! absolute difference; `max_scaled`, that difference over the reference's largest absolute value
//! in the array; `max_rel`, the largest difference relative to the reference value;
//! `max_fraction_differing`), and so are `values_differing` and `cells_differing` (totals over
//! the objects the floor measured, so they bound one object too). Other keys (`note`,
//! `fraction_differing`, `float_ulps`, `objects`, ...) are measurements, not bounds. An entry
//! with `may_differ = true` holds something the floor found differing without a numeric
//! measure (a file of bytes, such as Annoy's index); such a check passes in tolerance mode
//! whatever its bytes, and its downstream checks carry the claim.
//!
//! **Signs.** PCA components may change sign between kernels (irlba's vectors are defined up to
//! sign), which moves no distance. Arrays whose columns are components (embeddings, loadings,
//! Annoy's items, silhouette's centroids) are compared column by column after aligning each
//! column's sign, as the floor was measured.

use crate::json::Json;
use crate::toml::{self, Table, Value};

/// How a comparison grades its BLAS-dependent checks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mode {
    /// The kernels match: every check is exact.
    Exact { kernel: String },
    /// They differ, or one is not recorded: BLAS-dependent checks are graded against the floor.
    Tolerance { reason: String },
}

impl Mode {
    pub fn is_exact(&self) -> bool {
        matches!(self, Mode::Exact { .. })
    }

    pub fn name(&self) -> &'static str {
        match self {
            Mode::Exact { .. } => "exact",
            Mode::Tolerance { .. } => "tolerance",
        }
    }

    /// Why this mode, for the report.
    pub fn reason(&self) -> String {
        match self {
            Mode::Exact { kernel } => format!("the reference and the candidate both recorded OpenBLAS kernel {kernel}: every check is exact"),
            Mode::Tolerance { reason } => format!("{reason}: BLAS-dependent checks are graded against the R-vs-R floor, BLAS-free checks stay exact"),
        }
    }
}

/// The OpenBLAS kernel a run's `summary.json` records (`run.fingerprint.openblas_corename`).
pub fn recorded_kernel(summary: Option<&Json>) -> Option<String> {
    summary?.get("run")?.get("fingerprint")?.get("openblas_corename")?.as_str().map(String::from)
}

/// The mode from the two runs' recorded kernels. A run without a recorded kernel is never graded
/// as exact.
pub fn mode(reference: Option<&str>, candidate: Option<&str>) -> Mode {
    match (reference, candidate) {
        (Some(r), Some(c)) if r == c => Mode::Exact { kernel: r.to_string() },
        (Some(r), Some(c)) => Mode::Tolerance { reason: format!("the reference recorded OpenBLAS kernel {r}, the candidate {c}") },
        (None, Some(c)) => Mode::Tolerance { reason: format!("the reference records no OpenBLAS kernel (the candidate's is {c})") },
        (Some(r), None) => Mode::Tolerance { reason: format!("the candidate records no OpenBLAS kernel (the reference's is {r})") },
        (None, None) => Mode::Tolerance { reason: "neither run records an OpenBLAS kernel".into() },
    }
}

/// Whether a check's inputs passed through a BLAS product in the candidate's run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    BlasFree,
    BlasDependent,
}

impl Class {
    pub fn name(self) -> &'static str {
        match self {
            Class::BlasFree => "BLAS-free",
            Class::BlasDependent => "BLAS-dependent",
        }
    }
}

/// A floor table (`[<tool>.floor.kernel]`): array name -> its entry.
#[derive(Debug, Clone, Default)]
pub struct Floor {
    pub entries: Vec<(String, Table)>,
    /// Where it came from (`clustering.floor.kernel`).
    pub source: String,
}

impl Floor {
    /// The table at `path` (`clustering.floor.kernel`) of a TOML document.
    pub fn from_toml(text: &str, path: &str) -> Result<Floor, String> {
        let t = toml::parse(text).map_err(|e| format!("tolerances: {e}"))?;
        let tab = t.path(path).and_then(Value::as_table).ok_or_else(|| format!("tolerances: no table [{path}]"))?;
        let mut entries = Vec::new();
        for (k, v) in &tab.entries {
            if let Value::Table(e) = v {
                entries.push((k.clone(), e.clone()));
            }
        }
        Ok(Floor { entries, source: path.to_string() })
    }

    pub fn get(&self, key: &str) -> Option<&Table> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, t)| t)
    }

    /// The entry's text, for the report (`{ max_abs = 1e-12, ... }`).
    pub fn describe(&self, key: &str) -> String {
        match self.get(key) {
            None => "not in the floor: identical across kernels".into(),
            Some(t) => {
                let parts: Vec<String> = t
                    .entries
                    .iter()
                    .filter(|(k, _)| is_bound(k) || k == "may_differ")
                    .map(|(k, v)| format!("{k} = {}", v.number_text().map(String::from).unwrap_or_else(|| match v {
                        Value::Boolean(b) => b.to_string(),
                        Value::String(s) => format!("{s:?}"),
                        _ => "?".into(),
                    })))
                    .collect();
                format!("{}.{key} {{ {} }}", self.source, parts.join(", "))
            }
        }
    }
}

fn is_bound(key: &str) -> bool {
    key.starts_with("max_") || key == "values_differing" || key == "cells_differing"
}

/// The statistics of one array against the reference's.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Stats {
    /// Elements compared, and those whose bits differ (after sign alignment, where asked).
    pub n: usize,
    pub n_differ: usize,
    /// Rows with a differing element (cells, when rows are cells).
    pub rows_differ: usize,
    pub max_abs: f64,
    /// `max_abs` over the reference's largest absolute value.
    pub max_scaled: f64,
    /// Largest `|a - b| / |a|` over nonzero reference values (infinite where `a` is 0 and `b` is not).
    pub max_rel: f64,
    /// Columns whose sign was flipped to align them.
    pub flipped: usize,
    /// Elements where exactly one side is NaN, or a shape mismatch.
    pub incomparable: bool,
}

impl Stats {
    pub fn fraction(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.n_differ as f64 / self.n as f64
        }
    }

    pub fn describe(&self) -> String {
        if self.n_differ == 0 {
            return format!("{} values identical{}", self.n, if self.flipped > 0 { format!(" after flipping {} columns' signs", self.flipped) } else { String::new() });
        }
        format!(
            "{} of {} values differ ({} rows); max abs {:.3e}, max scaled {:.3e}, max rel {:.3e}{}",
            self.n_differ,
            self.n,
            self.rows_differ,
            self.max_abs,
            self.max_scaled,
            self.max_rel,
            if self.flipped > 0 { format!("; {} columns' signs flipped", self.flipped) } else { String::new() }
        )
    }
}

/// Compares two `nrow x ncol` arrays (column-major), bitwise and numerically. With `align`, each
/// column of `b` whose negation is closer to `a`'s column is negated first.
pub fn stats(a: &[f64], b: &[f64], nrow: usize, ncol: usize, align: bool) -> Stats {
    let mut s = Stats { n: a.len(), ..Stats::default() };
    if a.len() != b.len() || a.len() != nrow * ncol {
        s.incomparable = true;
        s.n_differ = a.len().max(b.len());
        return s;
    }
    let max_a = a.iter().fold(0.0f64, |m, v| if v.abs() > m { v.abs() } else { m });
    let mut row_bad = vec![false; nrow];
    for c in 0..ncol {
        let (ac, bc) = (&a[c * nrow..(c + 1) * nrow], &b[c * nrow..(c + 1) * nrow]);
        let flip = align && {
            let (mut same, mut opp) = (0.0f64, 0.0f64);
            for (x, y) in ac.iter().zip(bc) {
                same += (x - y).abs();
                opp += (x + y).abs();
            }
            opp < same
        };
        if flip {
            s.flipped += 1;
        }
        for (r, (&x, &y0)) in ac.iter().zip(bc).enumerate() {
            let y = if flip { -y0 } else { y0 };
            if x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()) {
                continue;
            }
            s.n_differ += 1;
            row_bad[r] = true;
            if x.is_nan() != y.is_nan() {
                s.incomparable = true;
                continue;
            }
            let d = (x - y).abs();
            if d > s.max_abs || d.is_nan() {
                s.max_abs = d;
            }
            let rel = if x == 0.0 { if d == 0.0 { 0.0 } else { f64::INFINITY } } else { d / x.abs() };
            if rel > s.max_rel {
                s.max_rel = rel;
            }
        }
    }
    s.rows_differ = row_bad.iter().filter(|&&b| b).count();
    s.max_scaled = if max_a > 0.0 { s.max_abs / max_a } else if s.max_abs > 0.0 { f64::INFINITY } else { 0.0 };
    s
}

/// Whether `s` stays within the floor entry's bounds; `Err` names the first bound exceeded.
pub fn within(s: &Stats, entry: &Table) -> Result<(), String> {
    if s.incomparable {
        return Err("not comparable (a shape mismatch, or a NaN on one side only)".into());
    }
    if entry.get("may_differ").and_then(Value::as_bool) == Some(true) {
        return Ok(());
    }
    for (k, v) in &entry.entries {
        if !is_bound(k) {
            continue;
        }
        let Some(bound) = v.as_f64() else { continue };
        let got = match k.as_str() {
            "max_abs" => s.max_abs,
            "max_scaled" => s.max_scaled,
            "max_rel" => s.max_rel,
            "max_fraction_differing" => s.fraction(),
            "values_differing" => s.n_differ as f64,
            "cells_differing" => s.rows_differ as f64,
            _ => continue,
        };
        if got > bound || got.is_nan() {
            return Err(format!("{k} {got:.3e} exceeds the floor's {bound:.3e}"));
        }
    }
    Ok(())
}

/// The outcome of grading one compared array.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Identical,
    /// It differs, but within the floor (tolerance mode only).
    WithinFloor,
    Differ(String),
}

/// Grades one array's statistics: identical passes; a difference passes only in tolerance mode,
/// for a BLAS-dependent check whose floor entry `floor_key` exists and holds.
pub fn grade(s: &Stats, class: Class, mode: &Mode, floor: &Floor, floor_key: &str) -> Outcome {
    if s.n_differ == 0 && !s.incomparable {
        return Outcome::Identical;
    }
    if mode.is_exact() {
        return Outcome::Differ("exact mode (the kernels match): every value must be identical".into());
    }
    if class == Class::BlasFree {
        return Outcome::Differ("BLAS-free: exact in any mode".into());
    }
    match floor.get(floor_key) {
        None => Outcome::Differ(format!("{floor_key} is not in the floor ({}): identical across kernels there, so exact", floor.source)),
        Some(e) => match within(s, e) {
            Ok(()) => Outcome::WithinFloor,
            Err(m) => Outcome::Differ(m),
        },
    }
}

/// One row of an R-vs-R floor measurement: an array (or scalar, or file) of one context,
/// compared between two runs that differ only in their OpenBLAS kernel.
#[derive(Debug, Clone)]
pub struct FloorRow {
    /// The floor key: the array's name (`S3_embeddings` split into `S3_embeddings_dims_1_30`
    /// and `S3_embeddings_dims_31_<n>`), a scalar's name, or `K4_annoy_index`.
    pub key: String,
    pub ctx: String,
    pub stats: Stats,
}

/// Whether the floor compares an array's columns after aligning their signs.
pub fn aligned(array: &str) -> bool {
    matches!(array, "S3_embeddings" | "S3_loadings" | "K3_items" | "K8_centroids")
}

/// Measures the floor between two runs' stage directories: every block of doubles both hold
/// (with each column's sign aligned where [`aligned`]), every numeric scalar both recorded, and
/// Annoy's index file where both have it. Integer arrays are compared too (any difference there
/// is a difference the floor must record).
pub fn measure_stages(a: &std::path::Path, b: &std::path::Path) -> std::io::Result<Vec<FloorRow>> {
    use crate::stages::StageDir;
    let (sa, sb) = (StageDir::open(a)?, StageDir::open(b)?);
    let mut rows = Vec::new();
    let mut files: Vec<String> = sa.blocks.iter().map(|x| x.file.clone()).collect();
    files.sort();
    files.dedup();
    for file in &files {
        for ba in sa.blocks_of(file) {
            if !sa.has(file, &ba.ctx) || !sb.has(file, &ba.ctx) {
                continue;
            }
            let Some(bb) = sb.block(file, &ba.ctx).cloned() else { continue };
            let array = file.trim_end_matches(".npy").trim_end_matches(&format!(".{}", ba.ctx)).to_string();
            let (na, nb) = (sa.read(file)?, sb.read(file)?);
            let col = |n: &crate::Npy, blk: &crate::stages::Block, c: usize| -> Vec<f64> {
                if blk.dtype.contains('i') {
                    n.col_i32(c, blk.row_start, blk.nrow).into_iter().map(|v| v as f64).collect()
                } else {
                    n.col_f64(c, blk.row_start, blk.nrow)
                }
            };
            let va: Vec<f64> = (0..ba.ncol).flat_map(|c| col(&na, ba, c)).collect();
            let vb: Vec<f64> = (0..bb.ncol).flat_map(|c| col(&nb, &bb, c)).collect();
            if (ba.nrow, ba.ncol) != (bb.nrow, bb.ncol) {
                rows.push(FloorRow { key: array.clone(), ctx: ba.ctx.clone(), stats: Stats { incomparable: true, n: va.len(), n_differ: va.len().max(vb.len()), ..Stats::default() } });
                continue;
            }
            if array == "S3_embeddings" && ba.ncol > 30 {
                let split = 30 * ba.nrow;
                rows.push(FloorRow { key: "S3_embeddings_dims_1_30".into(), ctx: ba.ctx.clone(), stats: stats(&va[..split], &vb[..split], ba.nrow, 30, true) });
                rows.push(FloorRow { key: format!("S3_embeddings_dims_31_{}", ba.ncol), ctx: ba.ctx.clone(), stats: stats(&va[split..], &vb[split..], ba.nrow, ba.ncol - 30, true) });
                continue;
            }
            let key = if array == "S3_embeddings" { "S3_embeddings_dims_1_30".to_string() } else { array.clone() };
            rows.push(FloorRow { key, ctx: ba.ctx.clone(), stats: stats(&va, &vb, ba.nrow, ba.ncol, aligned(&array)) });
        }
    }
    // numeric scalars both recorded
    let read = |d: &std::path::Path| -> Option<Json> { Json::parse(&String::from_utf8(crate::tenx::read_maybe_gz(&d.join("stages.json")).ok()?).ok()?).ok()?.get("scalars").cloned() };
    if let (Some(Json::Object(ja)), Some(jb)) = (read(a), read(b)) {
        for (ctx, members) in &ja {
            let Json::Object(m) = members else { continue };
            for (key, va) in m {
                let Some(vb) = jb.get(ctx).and_then(|c| c.get(key)) else { continue };
                let nums = |v: &Json| -> Option<Vec<f64>> {
                    match v {
                        Json::Int(_) | Json::Num(_) => Some(vec![v.as_f64()?]),
                        Json::Vector(x) | Json::List(x) => x.iter().map(Json::as_f64).collect(),
                        _ => None,
                    }
                };
                if let (Some(x), Some(y)) = (nums(va), nums(vb)) {
                    if x.len() == y.len() && !x.is_empty() {
                        rows.push(FloorRow { key: key.clone(), ctx: ctx.clone(), stats: stats(&x, &y, x.len(), 1, false) });
                    }
                }
            }
        }
    }
    // Annoy's index, byte for byte
    if let (Ok(x), Ok(y)) = (std::fs::read(a.join("annoy.ann")), std::fs::read(b.join("annoy.ann"))) {
        let nd = if x.len() != y.len() { x.len().max(y.len()) } else { x.iter().zip(&y).filter(|(p, q)| p != q).count() };
        rows.push(FloorRow { key: "K4_annoy_index".into(), ctx: "chain".into(), stats: Stats { n: x.len(), n_differ: nd, rows_differ: (nd > 0) as usize, ..Stats::default() } });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modes() {
        assert!(mode(Some("Haswell"), Some("Haswell")).is_exact());
        assert!(!mode(Some("Haswell"), Some("SapphireRapids")).is_exact());
        assert!(!mode(None, Some("Haswell")).is_exact());
        assert!(!mode(Some("Haswell"), None).is_exact());
        assert!(mode(Some("SkylakeX"), Some("Haswell")).reason().contains("SkylakeX"));
    }

    #[test]
    fn statistics_and_alignment() {
        // two columns; the second flipped in b, with one last-bit difference in the first
        let a = [1.0, 2.0, -3.0, 4.0];
        let x = 2.0f64;
        let b = [1.0, f64::from_bits(x.to_bits() + 1), 3.0, -4.0];
        let s = stats(&a, &b, 2, 2, true);
        assert_eq!((s.n, s.n_differ, s.rows_differ, s.flipped), (4, 1, 1, 1));
        assert_eq!(s.max_abs, f64::from_bits(x.to_bits() + 1) - 2.0);
        let s = stats(&a, &b, 2, 2, false);
        assert_eq!(s.n_differ, 3);
    }

    #[test]
    fn grading() {
        let floor = Floor::from_toml("[t.floor.kernel]\nA = { max_abs = 1e-12, max_scaled = 1e-13, note = \"x\" }\nB = { may_differ = true }\nC = { values_differing = 1 }\n", "t.floor.kernel").unwrap();
        let tol = mode(Some("Haswell"), Some("SapphireRapids"));
        let ex = mode(Some("Haswell"), Some("Haswell"));
        let small = stats(&[100.0], &[100.0 + 1e-13], 1, 1, false);
        assert_eq!(grade(&small, Class::BlasDependent, &tol, &floor, "A"), Outcome::WithinFloor);
        assert!(matches!(grade(&small, Class::BlasDependent, &ex, &floor, "A"), Outcome::Differ(_)));
        assert!(matches!(grade(&small, Class::BlasFree, &tol, &floor, "A"), Outcome::Differ(_)));
        assert!(matches!(grade(&small, Class::BlasDependent, &tol, &floor, "Z"), Outcome::Differ(_)));
        let big = stats(&[1.0], &[1.0 + 1e-11], 1, 1, false);
        assert!(matches!(grade(&big, Class::BlasDependent, &tol, &floor, "A"), Outcome::Differ(m) if m.contains("max_abs")));
        assert_eq!(grade(&big, Class::BlasDependent, &tol, &floor, "B"), Outcome::WithinFloor);
        let two = stats(&[1.0, 2.0], &[1.5, 2.5], 2, 1, false);
        assert!(matches!(grade(&two, Class::BlasDependent, &tol, &floor, "C"), Outcome::Differ(_)));
        assert_eq!(grade(&stats(&[1.0], &[1.0], 1, 1, false), Class::BlasFree, &ex, &floor, "A"), Outcome::Identical);
        assert!(floor.describe("A").contains("max_abs = 1e-12") && !floor.describe("A").contains("note"));
    }
}
