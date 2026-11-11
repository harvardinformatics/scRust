//! 10x Matrix Market directories (Cell Ranger v3 layout): `matrix.mtx(.gz)`,
//! `barcodes.tsv(.gz)` and `features.tsv(.gz)`: the barcodes, the features and the matrix's
//! dimensions; the entries are read by [`crate::sparse::read_mtx`].

use crate::gz;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

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
