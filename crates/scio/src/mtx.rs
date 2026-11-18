//! Plain Matrix Market files (not 10x directories): a sparse matrix as one gzipped
//! `coordinate` file with comment lines, entries in stored order (docs/PLAN.md C3, "Outputs":
//! `nn.mtx.gz`, `snn.mtx.gz`), and a reader that keeps the file's entry order, for comparisons
//! of stored order (stage K6).

use crate::deflate::{Chunk, GzWriter};
use crate::sparse::CscMatrix;
use crate::tenx::{find, push_uint, push_value, read_maybe_gz};
use std::fs;
use std::io::{self, Write};
use std::path::Path;

fn invalid(m: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m)
}

/// Writes `m` to `path` as gzipped coordinate Matrix Market: the banner, `comments` (each
/// written after `%`), the size line, then every stored entry in column-major order (1-based),
/// explicit zeros included. Values as C's `%.17g` writes them (integer-valued ones as
/// integers), as the R reference scripts write theirs (`mtx_write_plain`); the field is
/// "integer" when every value is integer-valued, "real" otherwise. Returns the field. The bytes
/// do not depend on `threads`.
pub fn write_gz(path: &Path, m: &CscMatrix, comments: &[String], threads: usize) -> io::Result<&'static str> {
    if m.p.len() != m.ncol + 1 || m.i.len() != m.x.len() || m.p[m.ncol] != m.x.len() {
        return Err(invalid(format!("{}: inconsistent matrix", path.display())));
    }
    write_gz_with(path, m.nrow, m.ncol, &m.p, &m.i, &|k| m.x[k], comments, threads)
}

/// [`write_gz`] for a matrix given by its compressed columns (`p`, `i`: rows within each column,
/// as stored) and `value(k)`, the value of stored entry `k`: the layout of matrices whose values
/// are computed as they are written (the corrected counts' integers and their `log1p`, M13).
#[allow(clippy::too_many_arguments)]
pub fn write_gz_with(path: &Path, nrow: usize, ncol: usize, p: &[usize], i: &[u32], value: &(dyn Fn(usize) -> f64 + Sync), comments: &[String], threads: usize) -> io::Result<&'static str> {
    if p.len() != ncol + 1 || p[ncol] != i.len() {
        return Err(invalid(format!("{}: inconsistent matrix", path.display())));
    }
    let nnz = i.len();
    let mut integer = true;
    for k in 0..nnz {
        let v = value(k);
        if !v.is_finite() {
            return Err(invalid(format!("{}: a non-finite value cannot be written", path.display())));
        }
        integer &= v == v.trunc();
    }
    let field = if integer { "integer" } else { "real" };
    write_gz_formatted(path, nrow, ncol, p, i, field, &|k, out| push_value(out, value(k)), comments, threads)?;
    Ok(field)
}

/// [`write_gz_with`] with each value's text given by `format(k, out)`, which appends stored entry
/// `k`'s value to `out` exactly as `%.17g` writes it (a caller that knows its values' text, such
/// as the `log1p` of small counts, formats each distinct value once), under the given `field`
/// (`integer` or `real`).
#[allow(clippy::too_many_arguments)]
pub fn write_gz_formatted(path: &Path, nrow: usize, ncol: usize, p: &[usize], i: &[u32], field: &str, format_value: &(dyn Fn(usize, &mut Vec<u8>) + Sync), comments: &[String], threads: usize) -> io::Result<()> {
    if p.len() != ncol + 1 || p[ncol] != i.len() {
        return Err(invalid(format!("{}: inconsistent matrix", path.display())));
    }
    let nnz = i.len();
    let file = fs::File::create(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut gz = GzWriter::new(io::BufWriter::with_capacity(1 << 20, file))?;
    let mut header = format!("%%MatrixMarket matrix coordinate {field} general\n");
    for c in comments {
        header.push('%');
        header.push_str(c);
        header.push('\n');
    }
    header.push_str(&format!("{nrow} {ncol} {nnz}\n"));
    gz.write_chunk(&Chunk::compress(header.as_bytes()))?;
    // pieces of whole columns, about PIECE entries each, formatted and compressed in parallel
    const PIECE: usize = 1 << 17;
    let mut bounds = vec![0usize];
    let mut j = 0;
    while j < ncol {
        let start = p[j];
        while j < ncol && p[j + 1] - start <= PIECE {
            j += 1;
        }
        if j == bounds[bounds.len() - 1] {
            j += 1; // one column larger than a piece
        }
        bounds.push(j);
    }
    let pieces: Vec<(usize, usize)> = bounds.windows(2).map(|w| (w[0], w[1])).collect();
    let format = |&(a, b): &(usize, usize)| {
        let mut out = Vec::with_capacity((p[b] - p[a]) * 24);
        for col in a..b {
            for idx in p[col]..p[col + 1] {
                push_uint(&mut out, i[idx] as u64 + 1);
                out.push(b' ');
                push_uint(&mut out, col as u64 + 1);
                out.push(b' ');
                format_value(idx, &mut out);
                out.push(b'\n');
            }
        }
        Chunk::compress(&out)
    };
    let threads = threads.max(1);
    for batch in pieces.chunks(threads * 2) {
        let chunks: Vec<Chunk> = if threads == 1 || batch.len() == 1 {
            batch.iter().map(format).collect()
        } else {
            std::thread::scope(|s| {
                let handles: Vec<_> = batch.iter().map(|piece| s.spawn(move || format(piece))).collect();
                handles.into_iter().map(|h| h.join().unwrap()).collect()
            })
        };
        for c in &chunks {
            gz.write_chunk(c)?;
        }
    }
    gz.finish()?.flush()?;
    Ok(())
}

/// How two Matrix Market files compare on their data lines: the banner, the size line and the
/// entries, in order; comment lines are skipped (an R reference's name R's script, a Rust run's
/// scrust's command).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataLines {
    /// The banners (`%%MatrixMarket ...`), the size lines.
    pub banner: (String, String),
    pub size: (String, String),
    /// Entry lines in each file, and the positions where the two differ (or where one has an
    /// entry the other lacks).
    pub entries: (usize, usize),
    pub differing: usize,
    /// The first differing entry, as `(position, line in a, line in b)`.
    pub first: Option<(usize, String, String)>,
}

impl DataLines {
    pub fn identical(&self) -> bool {
        self.banner.0 == self.banner.1 && self.size.0 == self.size.1 && self.entries.0 == self.entries.1 && self.differing == 0
    }
}

/// The banner, the size line and the entry lines of a Matrix Market file's text.
fn split_data_lines(text: &[u8]) -> (String, String, Vec<&[u8]>) {
    let mut lines = text.split(|&b| b == b'\n');
    let banner = String::from_utf8_lossy(lines.next().unwrap_or(b"")).into_owned();
    let mut size = String::new();
    for l in lines.by_ref() {
        if l.first() == Some(&b'%') {
            continue;
        }
        size = String::from_utf8_lossy(l).into_owned();
        break;
    }
    let entries: Vec<&[u8]> = lines.filter(|l| !l.is_empty() && l.first() != Some(&b'%')).collect();
    (banner, size, entries)
}

/// Compares two Matrix Market files (gzipped or not) on their data lines, byte for byte.
pub fn compare_data_lines(a: &Path, b: &Path) -> io::Result<DataLines> {
    let ta = read_maybe_gz(a)?;
    let tb = read_maybe_gz(b)?;
    let (ba, sa, ea) = split_data_lines(&ta);
    let (bb, sb, eb) = split_data_lines(&tb);
    let m = ea.len().min(eb.len());
    let mut differing = 0usize;
    let mut first = None;
    for k in 0..m {
        if ea[k] != eb[k] {
            differing += 1;
            if first.is_none() {
                first = Some((k + 1, String::from_utf8_lossy(ea[k]).into_owned(), String::from_utf8_lossy(eb[k]).into_owned()));
            }
        }
    }
    differing += ea.len().max(eb.len()) - m;
    Ok(DataLines { banner: (ba, bb), size: (sa, sb), entries: (ea.len(), eb.len()), differing, first })
}

/// A coordinate Matrix Market file's entries in file order, 0-based.
#[derive(Debug, Clone, PartialEq)]
pub struct Entries {
    pub field: String,
    pub nrow: usize,
    pub ncol: usize,
    pub rows: Vec<u32>,
    pub cols: Vec<u32>,
    pub vals: Vec<f64>,
}

/// Reads a `coordinate` `integer` or `real` `general` Matrix Market file (gzipped or not; `path`
/// with or without its `.gz`), keeping the entries in the file's order. Values are parsed
/// exactly (integers directly, others with a correctly rounded parser, which recovers every
/// `%.17g`).
pub fn read_entries(path: &Path) -> io::Result<Entries> {
    let path = match (path.exists(), path.parent(), path.file_name().and_then(|n| n.to_str())) {
        (true, _, _) => path.to_path_buf(),
        (false, Some(dir), Some(name)) => find(dir, name.trim_end_matches(".gz"))?,
        _ => return Err(io::Error::new(io::ErrorKind::NotFound, format!("{}: not found", path.display()))),
    };
    let bytes = read_maybe_gz(&path)?;
    let name = path.display().to_string();
    let text = std::str::from_utf8(&bytes).map_err(|_| invalid(format!("{name}: not text")))?;
    let mut lines = text.lines();
    let banner = lines.next().unwrap_or("");
    let words: Vec<String> = banner.split_whitespace().map(|w| w.to_ascii_lowercase()).collect();
    if words.len() < 5 || words[0] != "%%matrixmarket" || words[1] != "matrix" || words[2] != "coordinate"
        || !(words[3] == "integer" || words[3] == "real") || words[4] != "general"
    {
        return Err(invalid(format!("{name}: unsupported header {banner:?}")));
    }
    let size = loop {
        match lines.next() {
            None => return Err(invalid(format!("{name}: no size line"))),
            Some(l) if l.starts_with('%') => continue,
            Some(l) => break l,
        }
    };
    let dims: Vec<usize> = size.split_whitespace().filter_map(|t| t.parse().ok()).collect();
    if dims.len() != 3 {
        return Err(invalid(format!("{name}: bad size line {size:?}")));
    }
    let nnz = dims[2];
    let mut e = Entries { field: words[3].clone(), nrow: dims[0], ncol: dims[1], rows: Vec::with_capacity(nnz), cols: Vec::with_capacity(nnz), vals: Vec::with_capacity(nnz) };
    for line in lines {
        let mut t = line.split_whitespace();
        let Some(a) = t.next() else { continue };
        let (Some(b), Some(c)) = (t.next(), t.next()) else {
            return Err(invalid(format!("{name}: short entry line {line:?}")));
        };
        let r: usize = a.parse().map_err(|_| invalid(format!("{name}: bad row {a:?}")))?;
        let k: usize = b.parse().map_err(|_| invalid(format!("{name}: bad column {b:?}")))?;
        let v: f64 = c.parse().map_err(|_| invalid(format!("{name}: bad value {c:?}")))?;
        if r == 0 || r > e.nrow || k == 0 || k > e.ncol {
            return Err(invalid(format!("{name}: entry ({r}, {k}) out of range")));
        }
        e.rows.push((r - 1) as u32);
        e.cols.push((k - 1) as u32);
        e.vals.push(v);
    }
    if e.vals.len() != nnz {
        return Err(invalid(format!("{name}: expected {nnz} entries, read {}", e.vals.len())));
    }
    Ok(e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_keeps_order_and_bits() {
        let dir = crate::test_dir("mtx");
        let m = CscMatrix { nrow: 4, ncol: 3, p: vec![0, 2, 2, 5], i: vec![0, 3, 1, 2, 3], x: vec![1.0, 0.29032258064516131, 2.0, 1.0 / 3.0, 0.0] };
        let p = dir.join("g.mtx.gz");
        for t in [1, 4] {
            assert_eq!(write_gz(&p, &m, &["comment one".into()], t).unwrap(), "real");
            let e = read_entries(&p).unwrap();
            assert_eq!((e.nrow, e.ncol, e.field.as_str()), (4, 3, "real"));
            assert_eq!(e.rows, vec![0, 3, 1, 2, 3]);
            assert_eq!(e.cols, vec![0, 0, 2, 2, 2]);
            assert_eq!(e.vals.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), m.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            let text = String::from_utf8(crate::gz::gunzip(&fs::read(&p).unwrap()).unwrap()).unwrap();
            assert!(text.starts_with("%%MatrixMarket matrix coordinate real general\n%comment one\n4 3 5\n1 1 1\n4 1 0.29032258064516131\n"));
        }
        let ones = CscMatrix { nrow: 2, ncol: 2, p: vec![0, 1, 2], i: vec![1, 0], x: vec![1.0, 1.0] };
        assert_eq!(write_gz(&p, &ones, &[], 1).unwrap(), "integer");
        // a matrix spanning many pieces gives the same bytes at any thread count
        let n = 600usize;
        let mut big = CscMatrix::zeros(n, n);
        big.p.clear();
        big.p.push(0);
        for j in 0..n {
            for i in (j % 7..n).step_by(3) {
                big.i.push(i as u32);
                big.x.push((i * j) as f64 / 7.0);
            }
            big.p.push(big.i.len());
        }
        let (a, b) = (dir.join("a.mtx.gz"), dir.join("b.mtx.gz"));
        write_gz(&a, &big, &[], 1).unwrap();
        write_gz(&b, &big, &[], 6).unwrap();
        assert_eq!(fs::read(&a).unwrap(), fs::read(&b).unwrap());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn computed_values_and_data_lines() {
        let dir = crate::test_dir("mtx_with");
        // integer values computed from u32 storage, then their log1p from the same structure
        let (p, i, x) = (vec![0usize, 2, 2, 4], vec![0u32, 2, 1, 2], vec![1u32, 3, 2, 7]);
        let a = dir.join("a.mtx.gz");
        let field = write_gz_with(&a, 3, 3, &p, &i, &|k| x[k] as f64, &[" counts".into()], 2).unwrap();
        assert_eq!(field, "integer");
        let text = String::from_utf8(crate::gz::gunzip(&fs::read(&a).unwrap()).unwrap()).unwrap();
        assert_eq!(text, "%%MatrixMarket matrix coordinate integer general\n% counts\n3 3 4\n1 1 1\n3 1 3\n2 3 2\n3 3 7\n");
        let m = CscMatrix { nrow: 3, ncol: 3, p: p.clone(), i: i.clone(), x: x.iter().map(|&v| v as f64).collect() };
        let b = dir.join("b.mtx.gz");
        write_gz(&b, &m, &[" another comment".into(), " and another".into()], 1).unwrap();
        let d = compare_data_lines(&a, &b).unwrap();
        assert!(d.identical(), "{d:?}");
        assert_eq!(d.entries, (4, 4));
        let l = dir.join("l.mtx.gz");
        assert_eq!(write_gz_with(&l, 3, 3, &p, &i, &|k| (x[k] as f64).ln_1p(), &[], 1).unwrap(), "real");
        let d = compare_data_lines(&a, &l).unwrap();
        assert!(!d.identical());
        assert_eq!((d.differing, d.first.as_ref().map(|f| f.0)), (4, Some(1)));
        assert_ne!(d.banner.0, d.banner.1);
        // an empty matrix: every value integer-valued, as R's all(logical(0)) says
        assert_eq!(write_gz_with(&dir.join("e.mtx.gz"), 2, 2, &[0, 0, 0], &[], &|_| 0.5, &[], 1).unwrap(), "integer");
        fs::remove_dir_all(&dir).unwrap();
    }
}
