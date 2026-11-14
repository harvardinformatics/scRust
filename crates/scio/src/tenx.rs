//! 10x Matrix Market directories (Cell Ranger v3 layout): `matrix.mtx(.gz)`,
//! `barcodes.tsv(.gz)` and `features.tsv(.gz)`. Reading covers the barcodes, the features and
//! the matrix's dimensions; [`crate::sparse::read_mtx`] reads the entries. [`write_columns`]
//! writes the columns a tool keeps as a new, gzipped directory (docs/PLAN.md A3).

use crate::deflate::{gzip, Chunk, GzWriter};
use crate::gz;
use crate::sparse::CscMatrix;
use crate::text::fmt17;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

fn invalid(m: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m)
}

/// The file `name` in `dir`, gzipped (`name.gz`) or not.
pub fn find(dir: &Path, name: &str) -> io::Result<PathBuf> {
    for candidate in [dir.join(format!("{name}.gz")), dir.join(name)] {
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    Err(io::Error::new(io::ErrorKind::NotFound, format!("{}: no {name}(.gz)", dir.display())))
}

/// The bytes of a file, decompressed if it is gzipped.
pub fn read_maybe_gz(path: &Path) -> io::Result<Vec<u8>> {
    let bytes = fs::read(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    if gz::is_gzip(&bytes) {
        gz::gunzip(&bytes).map_err(|e| invalid(format!("{}: {e}", path.display())))
    } else {
        Ok(bytes)
    }
}

/// The lines of a text file (gzipped or not).
pub fn read_lines(path: &Path) -> io::Result<Vec<String>> {
    let bytes = read_maybe_gz(path)?;
    let text = String::from_utf8(bytes).map_err(|_| invalid(format!("{}: not UTF-8", path.display())))?;
    Ok(text.lines().map(str::to_owned).collect())
}

/// The cell barcodes, in column order.
pub fn read_barcodes(dir: &Path) -> io::Result<Vec<String>> {
    read_lines(&find(dir, "barcodes.tsv")?)
}

/// One feature: id, name and type (columns 1 to 3 of `features.tsv`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub id: String,
    pub name: String,
    pub kind: String,
}

/// The features, in row order.
pub fn read_features(dir: &Path) -> io::Result<Vec<Feature>> {
    let lines = read_lines(&find(dir, "features.tsv")?)?;
    Ok(lines
        .iter()
        .map(|l| {
            let mut f = l.split('\t');
            let id = f.next().unwrap_or("").to_string();
            let name = f.next().unwrap_or(&id).to_string();
            let kind = f.next().unwrap_or("").to_string();
            Feature { id, name, kind }
        })
        .collect())
}

/// The header of `matrix.mtx`: field type ("integer" or "real") and dimensions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MtxHeader {
    pub field: String,
    pub nrow: usize,
    pub ncol: usize,
    pub nnz: usize,
}

/// Reads the banner and size line of `matrix.mtx(.gz)`.
pub fn read_mtx_header(dir: &Path) -> io::Result<MtxHeader> {
    let path = find(dir, "matrix.mtx")?;
    let bytes = read_maybe_gz(&path)?;
    let text = String::from_utf8_lossy(&bytes[..bytes.len().min(1 << 16)]);
    let mut lines = text.lines();
    let banner = lines.next().unwrap_or("");
    let words: Vec<&str> = banner.split_whitespace().collect();
    if words.len() < 5 || words[0] != "%%MatrixMarket" || words[2] != "coordinate" {
        return Err(invalid(format!("{}: unsupported header {banner:?}", path.display())));
    }
    let size = lines
        .find(|l| !l.starts_with('%'))
        .ok_or_else(|| invalid(format!("{}: no size line", path.display())))?;
    let dims: Vec<usize> = size.split_whitespace().filter_map(|w| w.parse().ok()).collect();
    if dims.len() != 3 {
        return Err(invalid(format!("{}: bad size line {size:?}", path.display())));
    }
    Ok(MtxHeader { field: words[3].to_string(), nrow: dims[0], ncol: dims[1], nnz: dims[2] })
}

/// What [`write_columns`] wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Written {
    pub nrow: usize,
    pub ncol: usize,
    pub nnz: usize,
    /// The header's field: "integer" when every value written is integer-valued, else "real".
    pub field: &'static str,
}

/// Entries per piece of `matrix.mtx.gz`: each piece is formatted and compressed on its own
/// ([`crate::deflate`]), on whichever thread is free.
pub const PIECE_ENTRIES: usize = 1 << 17;

/// Writes the columns `cols` of `m`, in that order, as a gzipped 10x (Cell Ranger v3) directory
/// `dir`, which is created if needed:
/// - `matrix.mtx.gz`: coordinate Matrix Market, 1-based. Each column's entries are written
///   exactly as stored, explicit zeros included, in column-major order. After the banner comes
///   the comment line Cell Ranger writes, `%metadata_json`, naming `software`;
/// - `barcodes.tsv.gz`: `barcodes`, one per column written;
/// - `features.tsv.gz`: `features`, one line per row, as given (usually the input's lines).
///
/// Values are written as C's `%.17g` writes them, integer-valued ones as integers, like the R
/// reference scripts' values: a correctly rounded reader recovers every double. The header says
/// "integer" when every value written is integer-valued, "real" otherwise. The bytes do not
/// depend on `threads`.
pub fn write_columns(dir: &Path, m: &CscMatrix, cols: &[usize], barcodes: &[String], features: &[String], software: &str, threads: usize) -> io::Result<Written> {
    write_columns_in(dir, m, cols, barcodes, features, software, threads, PIECE_ENTRIES)
}

#[allow(clippy::too_many_arguments)]
fn write_columns_in(dir: &Path, m: &CscMatrix, cols: &[usize], barcodes: &[String], features: &[String], software: &str, threads: usize, piece: usize) -> io::Result<Written> {
    if barcodes.len() != cols.len() {
        return Err(invalid(format!("{} barcodes for {} columns", barcodes.len(), cols.len())));
    }
    if features.len() != m.nrow {
        return Err(invalid(format!("{} features for {} rows", features.len(), m.nrow)));
    }
    if let Some(&j) = cols.iter().find(|&&j| j >= m.ncol) {
        return Err(invalid(format!("column {j} of a matrix with {} columns", m.ncol)));
    }
    // where each written column's entries start in the output, and the header's field
    let mut starts = Vec::with_capacity(cols.len() + 1);
    starts.push(0usize);
    let mut integer = true;
    for &j in cols {
        let (_, x) = m.col(j);
        for &v in x {
            if !v.is_finite() {
                return Err(invalid(format!("column {j}: a non-finite value cannot be written")));
            }
            integer &= v == v.trunc();
        }
        starts.push(starts.last().unwrap() + x.len());
    }
    let nnz = *starts.last().unwrap();
    let field = if integer { "integer" } else { "real" };
    fs::create_dir_all(dir)?;

    let path = dir.join("matrix.mtx.gz");
    let file = fs::File::create(&path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut gz = GzWriter::new(io::BufWriter::with_capacity(1 << 20, file))?;
    let software = software.replace('\\', "\\\\").replace('"', "\\\"");
    let header = format!(
        "%%MatrixMarket matrix coordinate {field} general\n%metadata_json: {{\"software_version\": \"{software}\", \"format_version\": 2}}\n{} {} {nnz}\n",
        m.nrow,
        cols.len()
    );
    gz.write_chunk(&Chunk::compress(header.as_bytes()))?;
    let n_pieces = nnz.div_ceil(piece.max(1));
    let threads = threads.max(1);
    let batch = threads * 2;
    let mut first = 0;
    while first < n_pieces {
        let last = (first + batch).min(n_pieces);
        let compress = |k: usize| Chunk::compress(&format_piece(m, cols, &starts, k * piece, ((k + 1) * piece).min(nnz)));
        let chunks: Vec<Chunk> = if threads == 1 {
            (first..last).map(compress).collect()
        } else {
            let next = AtomicUsize::new(first);
            let slots: Vec<Mutex<Option<Chunk>>> = (first..last).map(|_| Mutex::new(None)).collect();
            std::thread::scope(|s| {
                for _ in 0..threads.min(last - first) {
                    s.spawn(|| loop {
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        if k >= last {
                            break;
                        }
                        *slots[k - first].lock().unwrap() = Some(compress(k));
                    });
                }
            });
            slots.into_iter().map(|s| s.into_inner().unwrap().unwrap()).collect()
        };
        for c in &chunks {
            gz.write_chunk(c)?;
        }
        first = last;
    }
    gz.finish()?.flush()?;

    for (name, lines) in [("barcodes.tsv.gz", barcodes), ("features.tsv.gz", features)] {
        let mut text = Vec::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
        for l in lines {
            text.extend_from_slice(l.as_bytes());
            text.push(b'\n');
        }
        let path = dir.join(name);
        fs::write(&path, gzip(&text)).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    }
    Ok(Written { nrow: m.nrow, ncol: cols.len(), nnz, field })
}

/// The Matrix Market lines of output entries `a..b` (`starts`: where each written column's
/// entries begin).
fn format_piece(m: &CscMatrix, cols: &[usize], starts: &[usize], a: usize, b: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity((b - a) * 16);
    if a >= b {
        return out;
    }
    // the written column holding entry a (the last start at or before a)
    let mut c = starts.partition_point(|&s| s <= a) - 1;
    while c < cols.len() && starts[c] < b {
        let (lo, hi) = (a.max(starts[c]), b.min(starts[c + 1]));
        let j = cols[c];
        // output entries lo..hi of column c are the input column's entries from m.p[j] on;
        // starts[c] <= lo <= hi, while m.p[j] may be smaller than starts[c] when the columns are
        // written out of input order
        for idx in m.p[j] + (lo - starts[c])..m.p[j] + (hi - starts[c]) {
            push_uint(&mut out, m.i[idx] as u64 + 1);
            out.push(b' ');
            push_uint(&mut out, c as u64 + 1);
            out.push(b' ');
            push_value(&mut out, m.x[idx]);
            out.push(b'\n');
        }
        c += 1;
    }
    out
}

pub(crate) fn push_uint(out: &mut Vec<u8>, mut v: u64) {
    let mut buf = [0u8; 20];
    let mut k = buf.len();
    loop {
        k -= 1;
        buf[k] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[k..]);
}

/// `%.17g` of a finite value: integers below 1e15 in magnitude directly, anything else (and a
/// negative zero) through [`fmt17`].
pub(crate) fn push_value(out: &mut Vec<u8>, v: f64) {
    if v == v.trunc() && v.abs() < 1e15 && !(v == 0.0 && v.is_sign_negative()) {
        if v < 0.0 {
            out.push(b'-');
        }
        push_uint(out, v.abs() as u64);
    } else {
        out.extend_from_slice(fmt17(v).as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::read_mtx;

    /// 6 x 5: explicit zeros, non-integer values, an empty column.
    fn small() -> CscMatrix {
        CscMatrix {
            nrow: 6,
            ncol: 5,
            p: vec![0, 3, 3, 6, 8, 10],
            i: vec![0, 2, 5, 1, 2, 3, 0, 4, 2, 5],
            x: vec![1.0, 0.0, 7.0, 0.1, 2.5e-7, 123456.0, 0.0, 1e300, 3.0, -2.0],
        }
    }

    fn lines(prefix: &str, n: usize) -> Vec<String> {
        (0..n).map(|k| format!("{prefix}{k}\t{prefix}{k}\tGene Expression")).collect()
    }

    /// The columns `cols` of `m`, as a matrix of their own.
    fn subset(m: &CscMatrix, cols: &[usize]) -> CscMatrix {
        let mut s = CscMatrix::zeros(m.nrow, 0);
        for &j in cols {
            let (r, x) = m.col(j);
            s.i.extend_from_slice(r);
            s.x.extend_from_slice(x);
            s.p.push(s.x.len());
            s.ncol += 1;
        }
        s
    }

    #[test]
    fn writes_what_it_reads_back() {
        let m = small();
        let features = lines("g", 6);
        for cols in [vec![0, 1, 2, 3, 4], vec![0, 2, 4], vec![1], vec![]] {
            let barcodes: Vec<String> = cols.iter().map(|j| format!("cell{j}")).collect();
            // pieces of every size, so that columns (the empty one too) span piece boundaries
            for (threads, piece) in [(1, PIECE_ENTRIES), (1, 1), (3, 2), (2, 3), (4, 100)] {
                let dir = crate::test_dir("tenx_write").join(format!("c{}_t{threads}_p{piece}", cols.len()));
                let w = write_columns_in(&dir, &m, &cols, &barcodes, &features, "scio test", threads, piece).unwrap();
                let back = read_mtx(&dir).unwrap();
                let want = subset(&m, &cols);
                assert_eq!((back.nrow, back.ncol, &back.p, &back.i), (want.nrow, want.ncol, &want.p, &want.i));
                assert!(back.x.iter().zip(&want.x).all(|(a, b)| a.to_bits() == b.to_bits()), "{cols:?}");
                let field = if want.x.iter().all(|v| v.fract() == 0.0) { "integer" } else { "real" };
                assert_eq!((w.ncol, w.nnz, w.field), (cols.len(), want.nnz(), field));
                assert_eq!(read_barcodes(&dir).unwrap(), barcodes);
                assert_eq!(read_lines(&find(&dir, "features.tsv").unwrap()).unwrap(), features);
                assert_eq!(read_mtx_header(&dir).unwrap().nnz, want.nnz());
            }
        }
    }

    #[test]
    fn same_bytes_for_any_thread_count() {
        let m = small();
        let cols = vec![4, 0, 2, 3];
        let barcodes: Vec<String> = cols.iter().map(|j| format!("c{j}")).collect();
        let features = lines("f", 6);
        let mut first: Option<Vec<u8>> = None;
        for threads in [1, 2, 3, 8] {
            let dir = crate::test_dir("tenx_threads").join(format!("t{threads}"));
            write_columns_in(&dir, &m, &cols, &barcodes, &features, "scio test", threads, 2).unwrap();
            let bytes = fs::read(dir.join("matrix.mtx.gz")).unwrap();
            match &first {
                None => first = Some(bytes),
                Some(f) => assert_eq!(f, &bytes, "{threads} threads"),
            }
        }
    }

    #[test]
    fn integer_header_and_text() {
        let mut m = small();
        m.x = vec![1.0, 0.0, 7.0, 2.0, 2.0, 123456.0, 0.0, 9.0, 3.0, 4.0];
        let dir = crate::test_dir("tenx_int");
        let cols: Vec<usize> = (0..5).collect();
        let barcodes: Vec<String> = cols.iter().map(|j| format!("c{j}")).collect();
        let w = write_columns(&dir, &m, &cols, &barcodes, &lines("g", 6), "scrust test", 1).unwrap();
        assert_eq!(w.field, "integer");
        let text = String::from_utf8(gz::gunzip(&fs::read(dir.join("matrix.mtx.gz")).unwrap()).unwrap()).unwrap();
        let want = "%%MatrixMarket matrix coordinate integer general\n\
                    %metadata_json: {\"software_version\": \"scrust test\", \"format_version\": 2}\n\
                    6 5 10\n1 1 1\n3 1 0\n6 1 7\n2 3 2\n3 3 2\n4 3 123456\n1 4 0\n5 4 9\n3 5 3\n6 5 4\n";
        assert_eq!(text, want);
        let mut bad = small();
        bad.x[0] = f64::NAN;
        assert!(write_columns(&dir, &bad, &cols, &barcodes, &lines("g", 6), "x", 1).is_err());
    }

    /// A tier-T input written back whole: the same matrix (its explicit zeros included), the
    /// same barcodes and feature lines.
    #[test]
    fn rewrites_a_tier_t_input() {
        let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder/soupx_tenx_test/input");
        let m = read_mtx(&input).unwrap();
        let barcodes = read_barcodes(&input).unwrap();
        let features = read_lines(&find(&input, "features.tsv").unwrap()).unwrap();
        let cols: Vec<usize> = (0..m.ncol).collect();
        let dir = crate::test_dir("tenx_tier_t");
        let w = write_columns(&dir, &m, &cols, &barcodes, &features, "scio test", 2).unwrap();
        assert_eq!(w, Written { nrow: 2000, ncol: 300, nnz: 64934, field: "integer" });
        assert_eq!(read_mtx(&dir).unwrap(), m);
        assert_eq!(read_barcodes(&dir).unwrap(), barcodes);
        assert_eq!(read_lines(&find(&dir, "features.tsv").unwrap()).unwrap(), features);
        // the decompressed feature and barcode files are the input's, byte for byte
        for f in ["barcodes.tsv.gz", "features.tsv.gz"] {
            assert_eq!(read_maybe_gz(&dir.join(f)).unwrap(), read_maybe_gz(&input.join(f)).unwrap(), "{f}");
        }
    }
}
