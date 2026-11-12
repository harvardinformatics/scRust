//! `.npy` arrays (NumPy format 1.0 and 2.0), as `validation/R/lib/io.R` writes them: dtype `<f8`
//! or `<i4`, little-endian, matrices in Fortran (column-major) order with shape `(nrow, ncol)`,
//! vectors with shape `(n,)`, the header padded with spaces to a multiple of 64 bytes.

use std::fmt;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

/// The element type of an array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dtype {
    /// `<f8`: little-endian double.
    F8,
    /// `<i4`: little-endian 32-bit integer.
    I4,
}

impl Dtype {
    pub fn descr(self) -> &'static str {
        match self {
            Dtype::F8 => "<f8",
            Dtype::I4 => "<i4",
        }
    }

    fn size(self) -> usize {
        match self {
            Dtype::F8 => 8,
            Dtype::I4 => 4,
        }
    }
}

/// An array read from a `.npy` file. Matrices are kept in the file's order; [`Npy::col`] and
/// [`Npy::get`] address them by (row, column) either way.
#[derive(Debug, Clone)]
pub struct Npy {
    pub dtype: Dtype,
    pub shape: Vec<usize>,
    pub fortran_order: bool,
    data: Vec<u8>,
}

#[derive(Debug)]
pub enum NpyError {
    Io(io::Error),
    Format(String),
}

impl fmt::Display for NpyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NpyError::Io(e) => write!(f, "{e}"),
            NpyError::Format(m) => write!(f, "{m}"),
        }
    }
}

impl std::error::Error for NpyError {}

impl From<io::Error> for NpyError {
    fn from(e: io::Error) -> Self {
        NpyError::Io(e)
    }
}

fn bad(path: &Path, what: &str) -> NpyError {
    NpyError::Format(format!("{}: {what}", path.display()))
}

/// The value of `'key': value` in a .npy header dictionary.
fn header_field<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("'{key}':");
    let start = header.find(&pat)? + pat.len();
    Some(header[start..].trim_start())
}

impl Npy {
    pub fn read(path: &Path) -> Result<Npy, NpyError> {
        let bytes = fs::read(path).map_err(|e| NpyError::Format(format!("{}: {e}", path.display())))?;
        Npy::parse(&bytes).map_err(|m| bad(path, &m))
    }

    pub fn parse(bytes: &[u8]) -> Result<Npy, String> {
        if bytes.len() < 10 || &bytes[..6] != b"\x93NUMPY" {
            return Err("not a .npy file".into());
        }
        let (hlen, hstart) = match bytes[6] {
            1 => (u16::from_le_bytes([bytes[8], bytes[9]]) as usize, 10),
            2 | 3 => {
                if bytes.len() < 12 {
                    return Err("truncated header".into());
                }
                (u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize, 12)
            }
            v => return Err(format!("unsupported .npy version {v}")),
        };
        let header = std::str::from_utf8(bytes.get(hstart..hstart + hlen).ok_or("truncated header")?)
            .map_err(|_| "header is not text")?;
        let descr = header_field(header, "descr").ok_or("no descr")?;
        let dtype = if descr.starts_with("'<f8'") {
            Dtype::F8
        } else if descr.starts_with("'<i4'") {
            Dtype::I4
        } else {
            return Err(format!("unsupported dtype {}", descr.split(',').next().unwrap_or("")));
        };
        let fortran_order = header_field(header, "fortran_order").ok_or("no fortran_order")?.starts_with("True");
        let shape_txt = header_field(header, "shape").ok_or("no shape")?;
        let close = shape_txt.find(')').ok_or("bad shape")?;
        let shape: Vec<usize> = shape_txt[1..close]
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().map_err(|_| format!("bad shape {shape_txt}")))
            .collect::<Result<_, _>>()?;
        let n: usize = shape.iter().product();
        let data = bytes[hstart + hlen..].to_vec();
        if data.len() != n * dtype.size() {
            return Err(format!("{} data bytes for shape {shape:?} of {}", data.len(), dtype.descr()));
        }
        Ok(Npy { dtype, shape, fortran_order, data })
    }

    /// Rows: the first dimension (1 for a scalar).
    pub fn nrow(&self) -> usize {
        self.shape.first().copied().unwrap_or(1)
    }

    /// Columns: the second dimension (1 for a vector).
    pub fn ncol(&self) -> usize {
        if self.shape.len() > 1 {
            self.shape[1]
        } else {
            1
        }
    }

    fn index(&self, row: usize, col: usize) -> usize {
        if self.fortran_order || self.shape.len() < 2 {
            col * self.nrow() + row
        } else {
            row * self.ncol() + col
        }
    }

    fn f8(&self, i: usize) -> f64 {
        f64::from_le_bytes(self.data[8 * i..8 * i + 8].try_into().unwrap())
    }

    fn i4(&self, i: usize) -> i32 {
        i32::from_le_bytes(self.data[4 * i..4 * i + 4].try_into().unwrap())
    }

    /// Element (row, col) as a double (an `<i4` value converted exactly).
    pub fn get(&self, row: usize, col: usize) -> f64 {
        let i = self.index(row, col);
        match self.dtype {
            Dtype::F8 => self.f8(i),
            Dtype::I4 => self.i4(i) as f64,
        }
    }

    /// Element (row, col) of an `<i4` array.
    pub fn get_i32(&self, row: usize, col: usize) -> i32 {
        assert_eq!(self.dtype, Dtype::I4, "get_i32 on a {} array", self.dtype.descr());
        self.i4(self.index(row, col))
    }

    /// Rows `start..start + len` of column `col`, as doubles.
    pub fn col_f64(&self, col: usize, start: usize, len: usize) -> Vec<f64> {
        (start..start + len).map(|r| self.get(r, col)).collect()
    }

    /// Rows `start..start + len` of column `col` of an `<i4` array.
    pub fn col_i32(&self, col: usize, start: usize, len: usize) -> Vec<i32> {
        (start..start + len).map(|r| self.get_i32(r, col)).collect()
    }

    /// Every element in file order, as doubles.
    pub fn f64s(&self) -> Vec<f64> {
        let n: usize = self.shape.iter().product();
        (0..n)
            .map(|i| match self.dtype {
                Dtype::F8 => self.f8(i),
                Dtype::I4 => self.i4(i) as f64,
            })
            .collect()
    }

    /// Rows `start..start + len` as a column-major block (`len` x ncol), as doubles.
    pub fn block_f64(&self, start: usize, len: usize) -> Vec<f64> {
        let mut out = Vec::with_capacity(len * self.ncol());
        for c in 0..self.ncol() {
            out.extend(self.col_f64(c, start, len));
        }
        out
    }
}

/// Reads a `<f8` matrix (`nrow x ncol`) row by row: returns its elements in row-major order
/// (row `r` at `r * ncol..`), with `(nrow, ncol)`. A column-major file (R's `npy_write`) is
/// transposed while it streams in blocks of columns, so memory stays at about one copy of the
/// matrix.
pub fn read_f64_row_major(path: &Path) -> Result<(Vec<f64>, usize, usize), NpyError> {
    use std::io::Read;
    let mut file = fs::File::open(path).map_err(|e| bad(path, &e.to_string()))?;
    let mut pre = [0u8; 12];
    file.read_exact(&mut pre[..10]).map_err(|e| bad(path, &e.to_string()))?;
    if &pre[..6] != b"\x93NUMPY" {
        return Err(bad(path, "not a .npy file"));
    }
    let hlen = match pre[6] {
        1 => u16::from_le_bytes([pre[8], pre[9]]) as usize,
        2 | 3 => {
            file.read_exact(&mut pre[10..12]).map_err(|e| bad(path, &e.to_string()))?;
            u32::from_le_bytes(pre[8..12].try_into().unwrap()) as usize
        }
        v => return Err(bad(path, &format!("unsupported .npy version {v}"))),
    };
    let mut hbytes = vec![0u8; hlen];
    file.read_exact(&mut hbytes).map_err(|e| bad(path, &e.to_string()))?;
    let header = std::str::from_utf8(&hbytes).map_err(|_| bad(path, "header is not text"))?;
    let descr = header_field(header, "descr").ok_or_else(|| bad(path, "no descr"))?;
    if !descr.starts_with("'<f8'") {
        return Err(bad(path, "not a <f8 array"));
    }
    let fortran = header_field(header, "fortran_order").ok_or_else(|| bad(path, "no fortran_order"))?.starts_with("True");
    let shape_txt = header_field(header, "shape").ok_or_else(|| bad(path, "no shape"))?;
    let close = shape_txt.find(')').ok_or_else(|| bad(path, "bad shape"))?;
    let shape: Vec<usize> = shape_txt[1..close].split(',').map(str::trim).filter(|s| !s.is_empty()).map(|s| s.parse().map_err(|_| bad(path, "bad shape"))).collect::<Result<_, _>>()?;
    if shape.len() != 2 {
        return Err(bad(path, &format!("not a matrix: shape {shape:?}")));
    }
    let (nrow, ncol) = (shape[0], shape[1]);
    let mut out = vec![0.0f64; nrow * ncol];
    let mut reader = io::BufReader::with_capacity(1 << 22, file);
    if !fortran {
        let mut buf = vec![0u8; 8 * ncol.max(1)];
        for r in 0..nrow {
            reader.read_exact(&mut buf[..8 * ncol]).map_err(|e| bad(path, &e.to_string()))?;
            for (c, b) in buf[..8 * ncol].chunks_exact(8).enumerate() {
                out[r * ncol + c] = f64::from_le_bytes(b.try_into().unwrap());
            }
        }
    } else {
        // blocks of columns: each row's values of the block are written contiguously
        let block = (1usize << 21).checked_div(8 * nrow.max(1)).unwrap_or(1).clamp(1, ncol.max(1));
        let mut buf = vec![0u8; 8 * nrow * block];
        let mut c0 = 0;
        while c0 < ncol {
            let nb = block.min(ncol - c0);
            reader.read_exact(&mut buf[..8 * nrow * nb]).map_err(|e| bad(path, &e.to_string()))?;
            for r in 0..nrow {
                for k in 0..nb {
                    let i = 8 * (k * nrow + r);
                    out[r * ncol + c0 + k] = f64::from_le_bytes(buf[i..i + 8].try_into().unwrap());
                }
            }
            c0 += nb;
        }
    }
    let mut rest = [0u8; 1];
    if reader.read(&mut rest).map_err(|e| bad(path, &e.to_string()))? != 0 {
        return Err(bad(path, "more data than the shape holds"));
    }
    Ok((out, nrow, ncol))
}

fn header_bytes(dtype: Dtype, shape: &[usize]) -> Vec<u8> {
    let shape_txt = match shape {
        [n] => format!("({n},)"),
        dims => format!("({})", dims.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(", ")),
    };
    let fortran = if shape.len() > 1 { "True" } else { "False" };
    let mut header = format!("{{'descr': '{}', 'fortran_order': {fortran}, 'shape': {shape_txt}, }}", dtype.descr());
    let total = 10 + header.len() + 1;
    header.push_str(&" ".repeat((64 - total % 64) % 64));
    header.push('\n');
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend_from_slice(&(header.len() as u16).to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    out
}

/// Writes doubles in column-major order with the given shape, byte for byte as R's `npy_write`.
pub fn write_f64(path: &Path, data: &[f64], shape: &[usize]) -> io::Result<()> {
    assert_eq!(data.len(), shape.iter().product::<usize>());
    let mut out = header_bytes(Dtype::F8, shape);
    out.reserve(data.len() * 8);
    for v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    fs::File::create(path)?.write_all(&out)
}

/// Writes 32-bit integers in column-major order with the given shape, as R's `npy_write`.
pub fn write_i32(path: &Path, data: &[i32], shape: &[usize]) -> io::Result<()> {
    assert_eq!(data.len(), shape.iter().product::<usize>());
    let mut out = header_bytes(Dtype::I4, shape);
    out.reserve(data.len() * 4);
    for v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    fs::File::create(path)?.write_all(&out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_in_r_layout() {
        let dir = crate::test_dir("npy");
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("m.npy");
        write_f64(&p, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[3, 2]).unwrap();
        let bytes = fs::read(&p).unwrap();
        assert_eq!(bytes.len(), 128 + 48); // the header padded to 128 bytes, and 6 doubles
        let a = Npy::read(&p).unwrap();
        assert_eq!((a.nrow(), a.ncol(), a.fortran_order), (3, 2, true));
        assert_eq!(a.get(2, 1), 6.0);
        assert_eq!(a.col_f64(1, 1, 2), vec![5.0, 6.0]);
        let q = dir.join("v.npy");
        write_i32(&q, &[7, -8], &[2]).unwrap();
        let b = Npy::read(&q).unwrap();
        assert_eq!((b.nrow(), b.ncol(), b.get_i32(1, 0)), (2, 1, -8));
        // row-major streaming read of a column-major matrix
        let (rm, nr, nc) = read_f64_row_major(&p).unwrap();
        assert_eq!((nr, nc), (3, 2));
        assert_eq!(rm, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        assert!(read_f64_row_major(&q).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }
}
