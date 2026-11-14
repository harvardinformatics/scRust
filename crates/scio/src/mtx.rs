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
    let mut integer = true;
    for &v in &m.x {
        if !v.is_finite() {
            return Err(invalid(format!("{}: a non-finite value cannot be written", path.display())));
        }
        integer &= v == v.trunc();
    }
    let field = if integer { "integer" } else { "real" };
    let file = fs::File::create(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut gz = GzWriter::new(io::BufWriter::with_capacity(1 << 20, file))?;
    let mut header = format!("%%MatrixMarket matrix coordinate {field} general\n");
    for c in comments {
        header.push('%');
        header.push_str(c);
        header.push('\n');
    }
    header.push_str(&format!("{} {} {}\n", m.nrow, m.ncol, m.nnz()));
    gz.write_chunk(&Chunk::compress(header.as_bytes()))?;
    // pieces of whole columns, about PIECE entries each, formatted and compressed in parallel
    const PIECE: usize = 1 << 17;
    let mut bounds = vec![0usize];
    let mut j = 0;
    while j < m.ncol {
        let start = m.p[j];
        while j < m.ncol && m.p[j + 1] - start <= PIECE {
            j += 1;
        }
        if j == bounds[bounds.len() - 1] {
            j += 1; // one column larger than a piece
        }
        bounds.push(j);
    }
    let pieces: Vec<(usize, usize)> = bounds.windows(2).map(|w| (w[0], w[1])).collect();
    let format = |&(a, b): &(usize, usize)| {
        let mut out = Vec::with_capacity((m.p[b] - m.p[a]) * 24);
        for col in a..b {
            for idx in m.p[col]..m.p[col + 1] {
                push_uint(&mut out, m.i[idx] as u64 + 1);
                out.push(b' ');
                push_uint(&mut out, col as u64 + 1);
                out.push(b' ');
                push_value(&mut out, m.x[idx]);
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
    Ok(field)
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
}
