//! Reading an SNN graph from plain files (docs/PLAN.md C3, `scrust cluster --snn <dir>`): the
//! `snn.mtx.gz` that `scrust neighbors` writes (with its `knn_index.tsv`), or an R run's
//! `graphs/snn.mtx.gz` as M8 saved them (with the run's `metadata.tsv`).

use scio::{CscMatrix, Table};
use std::path::{Path, PathBuf};

/// An SNN graph with its cells.
#[derive(Debug, Clone)]
pub struct SnnInput {
    pub snn: CscMatrix,
    pub barcodes: Vec<String>,
    /// The Matrix Market file read, and where the barcodes came from.
    pub source: PathBuf,
    pub barcodes_source: PathBuf,
}

/// Builds the dgCMatrix of a Matrix Market file's entries, which must be in its stored order:
/// column by column, rows strictly increasing within each column (as Seurat stores a Graph and
/// as `scrust neighbors` and M8's exports write it). Entries out of that order are refused
/// rather than sorted: the network's neighbour lists follow the stored order.
pub fn csc_from_entries(e: &scio::mtx::Entries) -> Result<CscMatrix, String> {
    let n = e.vals.len();
    let mut m = CscMatrix { nrow: e.nrow, ncol: e.ncol, p: vec![0; e.ncol + 1], i: Vec::with_capacity(n), x: Vec::with_capacity(n) };
    let mut prev: Option<(u32, u32)> = None;
    for t in 0..n {
        let (r, c) = (e.rows[t], e.cols[t]);
        if let Some((pr, pc)) = prev {
            if c < pc || (c == pc && r <= pr) {
                return Err(format!("entry {} ({}, {}) is out of the stored order (columns in turn, rows increasing)", t + 1, r + 1, c + 1));
            }
        }
        prev = Some((r, c));
        m.i.push(r);
        m.x.push(e.vals[t]);
        m.p[c as usize + 1] += 1;
    }
    for j in 0..e.ncol {
        m.p[j + 1] += m.p[j];
    }
    Ok(m)
}

/// The barcodes of a knn table (`scrust neighbors`' `knn_index.tsv`) or of a run's
/// `metadata.tsv`: the `barcode` column (else the first).
fn table_barcodes(path: &Path) -> Result<Vec<String>, String> {
    let bytes = scio::tenx::read_maybe_gz(path).map_err(|e| e.to_string())?;
    let text = String::from_utf8(bytes).map_err(|_| format!("{}: not text", path.display()))?;
    let t = Table::parse(&text).ok_or_else(|| format!("{}: empty", path.display()))?;
    let c = t.col("barcode").unwrap_or(0);
    Ok(t.rows.iter().map(|r| r.get(c).cloned().unwrap_or_default()).collect())
}

/// Reads the SNN graph of `dir`: `snn.mtx(.gz)` in `dir`, or in `dir/graphs` or `dir/neighbors`
/// (an R run's or a chain run's directory). The cells' barcodes come from `barcodes.tsv(.gz)`
/// or `knn_index.tsv` beside the graph, else from the `metadata.tsv(.gz)` of the directory
/// above it (an R run's layout).
pub fn read_snn_dir(dir: &Path) -> Result<SnnInput, String> {
    let mut found = None;
    for d in [dir.to_path_buf(), dir.join("graphs"), dir.join("neighbors")] {
        if let Ok(p) = scio::tenx::find(&d, "snn.mtx") {
            found = Some((d, p));
            break;
        }
    }
    let (gdir, path) = found.ok_or_else(|| format!("{}: no snn.mtx(.gz) here or in graphs/ or neighbors/", dir.display()))?;
    let e = scio::mtx::read_entries(&path).map_err(|e| e.to_string())?;
    let snn = csc_from_entries(&e).map_err(|m| format!("{}: {m}", path.display()))?;
    let mut barcodes_source = None;
    for name in ["barcodes.tsv", "knn_index.tsv"] {
        if let Ok(p) = scio::tenx::find(&gdir, name) {
            barcodes_source = Some(p);
            break;
        }
    }
    if barcodes_source.is_none() {
        if let Some(parent) = gdir.parent() {
            barcodes_source = scio::tenx::find(parent, "metadata.tsv").ok();
        }
    }
    let bpath = barcodes_source.ok_or_else(|| format!("{}: no barcodes (barcodes.tsv or knn_index.tsv beside the graph, or ../metadata.tsv)", gdir.display()))?;
    let barcodes = if bpath.file_name().and_then(|f| f.to_str()).is_some_and(|f| f.starts_with("barcodes.tsv")) {
        scio::tenx::read_lines(&bpath).map_err(|e| e.to_string())?
    } else {
        table_barcodes(&bpath)?
    };
    if barcodes.len() != snn.ncol || snn.nrow != snn.ncol {
        return Err(format!("{}: a {} x {} graph, but {} barcodes in {}", path.display(), snn.nrow, snn.ncol, barcodes.len(), bpath.display()));
    }
    Ok(SnnInput { snn, barcodes, source: path, barcodes_source: bpath })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entries_must_be_in_stored_order() {
        let e = scio::mtx::Entries { field: "real".into(), nrow: 3, ncol: 3, rows: vec![0, 2, 1, 0], cols: vec![0, 0, 1, 2], vals: vec![1.0, 0.5, 1.0, 0.25] };
        let m = csc_from_entries(&e).unwrap();
        assert_eq!(m.p, vec![0, 2, 3, 4]);
        assert_eq!(m.i, vec![0, 2, 1, 0]);
        let bad = scio::mtx::Entries { rows: vec![2, 0, 1, 0], ..e.clone() };
        assert!(csc_from_entries(&bad).is_err());
        let bad = scio::mtx::Entries { cols: vec![0, 1, 0, 2], ..e };
        assert!(csc_from_entries(&bad).is_err());
    }
}
