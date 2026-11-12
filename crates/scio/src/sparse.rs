//! A compressed sparse column matrix in R's `dgCMatrix` layout, and the 10x `matrix.mtx` reader.
//!
//! Entries are kept exactly as stored, **explicit zeros included**: the SoupX inputs store zeros,
//! and sctransform's `min_variance = (median(umi@x)/5)^2` counts every stored entry, so dropping
//! them would change results. Within each column, rows are increasing, as in a valid `dgCMatrix`.

use crate::tenx::{find, read_maybe_gz};
use std::io;
use std::path::Path;

/// `nrow x ncol`, column `j`'s entries at `p[j]..p[j + 1]` of `i` (0-based rows) and `x`.
#[derive(Debug, Clone, PartialEq)]
pub struct CscMatrix {
    pub nrow: usize,
    pub ncol: usize,
    pub p: Vec<usize>,
    pub i: Vec<u32>,
    pub x: Vec<f64>,
}

impl CscMatrix {
    /// An empty (all-zero) matrix.
    pub fn zeros(nrow: usize, ncol: usize) -> CscMatrix {
        CscMatrix { nrow, ncol, p: vec![0; ncol + 1], i: Vec::new(), x: Vec::new() }
    }

    pub fn nnz(&self) -> usize {
        self.x.len()
    }

    /// Rows and values of column `j`.
    pub fn col(&self, j: usize) -> (&[u32], &[f64]) {
        let (a, b) = (self.p[j], self.p[j + 1]);
        (&self.i[a..b], &self.x[a..b])
    }

    /// Checks the invariants of a `dgCMatrix`: `p` non-decreasing from 0 to nnz, rows in range and
    /// strictly increasing within each column.
    pub fn validate(&self) -> Result<(), String> {
        if self.p.len() != self.ncol + 1 || self.p[0] != 0 || self.p[self.ncol] != self.x.len() || self.i.len() != self.x.len() {
            return Err("inconsistent column pointers".into());
        }
        for j in 0..self.ncol {
            if self.p[j] > self.p[j + 1] {
                return Err(format!("column pointers decrease at column {j}"));
            }
            let (rows, _) = self.col(j);
            for w in rows.windows(2) {
                if w[0] >= w[1] {
                    return Err(format!("rows not increasing in column {j}"));
                }
            }
            if let Some(&r) = rows.last() {
                if r as usize >= self.nrow {
                    return Err(format!("row index out of range in column {j}"));
                }
            }
        }
        Ok(())
    }

    /// The columns `cols` in that order (R's `m[, cols]`): each column's entries are copied as
    /// stored, explicit zeros included.
    pub fn select_cols(&self, cols: &[usize]) -> CscMatrix {
        let nnz: usize = cols.iter().map(|&j| self.p[j + 1] - self.p[j]).sum();
        let mut out = CscMatrix { nrow: self.nrow, ncol: cols.len(), p: Vec::with_capacity(cols.len() + 1), i: Vec::with_capacity(nnz), x: Vec::with_capacity(nnz) };
        out.p.push(0);
        for &j in cols {
            let (r, v) = self.col(j);
            out.i.extend_from_slice(r);
            out.x.extend_from_slice(v);
            out.p.push(out.i.len());
        }
        out
    }

    /// The rows `rows` (increasing) in that order, as R's `m[rows, ]`: within each column the kept
    /// entries stay in order, explicit zeros included; row `rows[k]` becomes row `k`.
    pub fn select_rows(&self, rows: &[usize]) -> CscMatrix {
        let mut map = vec![u32::MAX; self.nrow];
        for (k, &r) in rows.iter().enumerate() {
            assert!(map[r] == u32::MAX, "select_rows: duplicate row {r}");
            if k > 0 {
                assert!(rows[k - 1] < r, "select_rows: rows must be increasing");
            }
            map[r] = k as u32;
        }
        let mut out = CscMatrix { nrow: rows.len(), ncol: self.ncol, p: Vec::with_capacity(self.ncol + 1), i: Vec::new(), x: Vec::new() };
        out.p.push(0);
        for j in 0..self.ncol {
            let (r, v) = self.col(j);
            for (&ri, &vi) in r.iter().zip(v) {
                let m = map[ri as usize];
                if m != u32::MAX {
                    out.i.push(m);
                    out.x.push(vi);
                }
            }
            out.p.push(out.i.len());
        }
        out
    }

    /// Columns of `self` followed by the columns of `other` (Matrix's `cbind`, which keeps every
    /// stored entry).
    pub fn cbind(&self, other: &CscMatrix) -> CscMatrix {
        assert_eq!(self.nrow, other.nrow, "cbind: row counts differ");
        let mut out = self.clone();
        let base = out.i.len();
        out.i.extend_from_slice(&other.i);
        out.x.extend_from_slice(&other.x);
        out.p.extend(other.p[1..].iter().map(|&q| q + base));
        out.ncol += other.ncol;
        out
    }

    /// The transpose's layout: for each row, its (column, value) entries in increasing column
    /// order, which is the order in which a column-major pass over the stored entries visits them.
    pub fn to_csr(&self) -> CsrMatrix {
        let mut counts = vec![0usize; self.nrow + 1];
        for &r in &self.i {
            counts[r as usize + 1] += 1;
        }
        for k in 0..self.nrow {
            counts[k + 1] += counts[k];
        }
        let mut next = counts.clone();
        let mut j_idx = vec![0u32; self.nnz()];
        let mut x = vec![0.0; self.nnz()];
        for j in 0..self.ncol {
            let (rows, vals) = self.col(j);
            for (&r, &v) in rows.iter().zip(vals) {
                let pos = next[r as usize];
                j_idx[pos] = j as u32;
                x[pos] = v;
                next[r as usize] += 1;
            }
        }
        CsrMatrix { nrow: self.nrow, ncol: self.ncol, p: counts, j: j_idx, x }
    }
}

/// Row-major companion of [`CscMatrix`]: row `r`'s entries at `p[r]..p[r + 1]`, in increasing
/// column order.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrMatrix {
    pub nrow: usize,
    pub ncol: usize,
    pub p: Vec<usize>,
    pub j: Vec<u32>,
    pub x: Vec<f64>,
}

impl CsrMatrix {
    pub fn row(&self, r: usize) -> (&[u32], &[f64]) {
        let (a, b) = (self.p[r], self.p[r + 1]);
        (&self.j[a..b], &self.x[a..b])
    }
}

fn invalid(m: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m)
}

/// Parses an unsigned decimal integer token.
fn parse_usize(tok: &[u8]) -> Option<usize> {
    if tok.is_empty() || tok.len() > 19 {
        return None;
    }
    let mut v: usize = 0;
    for &c in tok {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v * 10 + (c - b'0') as usize;
    }
    Some(v)
}

/// A value token: integers directly (exact), anything else through Rust's correctly rounded
/// parser. (The workflow's counts are integers, which R's `scan()` also reads exactly.)
fn parse_value(tok: &[u8]) -> Option<f64> {
    if let Some(v) = parse_usize(tok) {
        if v < (1usize << 53) {
            return Some(v as f64);
        }
    }
    std::str::from_utf8(tok).ok()?.parse::<f64>().ok()
}

/// Reads `matrix.mtx(.gz)` from a 10x directory: a `coordinate` `integer` or `real` `general`
/// Matrix Market file. Entries are kept as stored, explicit zeros included; they are sorted into
/// column-major order with increasing rows if the file has them in another order. Duplicate
/// entries are an error.
pub fn read_mtx(dir: &Path) -> io::Result<CscMatrix> {
    let path = find(dir, "matrix.mtx")?;
    let bytes = read_maybe_gz(&path)?;
    let name = path.display().to_string();
    let mut lines = bytes.split(|&b| b == b'\n');
    let banner = lines.next().unwrap_or(b"");
    let banner = String::from_utf8_lossy(banner);
    let words: Vec<String> = banner.split_whitespace().map(|w| w.to_ascii_lowercase()).collect();
    if words.len() < 5 || words[0] != "%%matrixmarket" || words[1] != "matrix" || words[2] != "coordinate"
        || !(words[3] == "integer" || words[3] == "real") || words[4] != "general"
    {
        return Err(invalid(format!("{name}: unsupported header {banner:?}")));
    }
    let size = loop {
        match lines.next() {
            None => return Err(invalid(format!("{name}: no size line"))),
            Some(l) if l.first() == Some(&b'%') => continue,
            Some(l) => break l,
        }
    };
    let dims: Vec<usize> = size.split(|c| c.is_ascii_whitespace()).filter(|t| !t.is_empty()).filter_map(parse_usize).collect();
    if dims.len() != 3 {
        return Err(invalid(format!("{name}: bad size line")));
    }
    let (nrow, ncol, nnz) = (dims[0], dims[1], dims[2]);
    if nrow > u32::MAX as usize {
        return Err(invalid(format!("{name}: too many rows")));
    }
    let mut rows: Vec<u32> = Vec::with_capacity(nnz);
    let mut cols: Vec<u32> = Vec::with_capacity(nnz);
    let mut vals: Vec<f64> = Vec::with_capacity(nnz);
    let mut sorted = true;
    let (mut prev_r, mut prev_c) = (0usize, 0usize);
    for line in lines {
        let mut toks = line.split(|c| c.is_ascii_whitespace()).filter(|t| !t.is_empty());
        let Some(t0) = toks.next() else { continue };
        let (Some(t1), Some(t2)) = (toks.next(), toks.next()) else {
            return Err(invalid(format!("{name}: short entry line")));
        };
        let r = parse_usize(t0).ok_or_else(|| invalid(format!("{name}: bad row index")))?;
        let c = parse_usize(t1).ok_or_else(|| invalid(format!("{name}: bad column index")))?;
        let v = parse_value(t2).ok_or_else(|| invalid(format!("{name}: bad value")))?;
        if r == 0 || r > nrow || c == 0 || c > ncol {
            return Err(invalid(format!("{name}: entry ({r}, {c}) out of range")));
        }
        if !rows.is_empty() && (c < prev_c || (c == prev_c && r <= prev_r)) {
            sorted = false;
        }
        prev_r = r;
        prev_c = c;
        rows.push((r - 1) as u32);
        cols.push((c - 1) as u32);
        vals.push(v);
    }
    if rows.len() != nnz {
        return Err(invalid(format!("{name}: expected {nnz} entries, read {}", rows.len())));
    }
    if !sorted {
        let mut idx: Vec<usize> = (0..nnz).collect();
        idx.sort_by_key(|&k| (cols[k], rows[k]));
        for w in idx.windows(2) {
            if (cols[w[0]], rows[w[0]]) == (cols[w[1]], rows[w[1]]) {
                return Err(invalid(format!("{name}: duplicate entry ({}, {})", rows[w[0]] + 1, cols[w[0]] + 1)));
            }
        }
        rows = idx.iter().map(|&k| rows[k]).collect();
        cols = idx.iter().map(|&k| cols[k]).collect();
        vals = idx.iter().map(|&k| vals[k]).collect();
    }
    let mut p = vec![0usize; ncol + 1];
    for &c in &cols {
        p[c as usize + 1] += 1;
    }
    for j in 0..ncol {
        p[j + 1] += p[j];
    }
    Ok(CscMatrix { nrow, ncol, p, i: rows, x: vals })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m() -> CscMatrix {
        // [[1, 0, 2], [0, 0(stored), 3]]
        CscMatrix { nrow: 2, ncol: 3, p: vec![0, 1, 2, 4], i: vec![0, 1, 0, 1], x: vec![1.0, 0.0, 2.0, 3.0] }
    }

    #[test]
    fn selection_keeps_stored_zeros() {
        let a = m();
        a.validate().unwrap();
        let s = a.select_cols(&[1, 1, 0]);
        assert_eq!((s.p.clone(), s.i.clone(), s.x.clone()), (vec![0, 1, 2, 3], vec![1, 1, 0], vec![0.0, 0.0, 1.0]));
        let r = a.select_rows(&[1]);
        assert_eq!((r.p, r.i, r.x), (vec![0, 0, 1, 2], vec![0, 0], vec![0.0, 3.0]));
        let c = a.cbind(&s);
        assert_eq!(c.ncol, 6);
        c.validate().unwrap();
        let t = a.to_csr();
        assert_eq!(t.row(1), (&[1u32, 2][..], &[0.0, 3.0][..]));
    }

    #[test]
    fn reads_a_tier_t_matrix() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder/soupx_tenx_test/input");
        let a = read_mtx(&dir).unwrap();
        a.validate().unwrap();
        assert_eq!((a.nrow, a.ncol, a.nnz()), (2000, 300, 64934));
        // soupx_tenx_test stores 2,686 explicit zeros (docs/validation/doubletfinder_m2.md)
        assert_eq!(a.x.iter().filter(|&&v| v == 0.0).count(), 2686);
    }
}
