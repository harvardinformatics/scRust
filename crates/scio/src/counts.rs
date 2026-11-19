//! A count matrix with its cells and features: what a 10x directory holds, and what one tool hands
//! to the next in memory (docs/PLAN.md A3, "In-memory hand-offs"; milestone M15).
//!
//! [`Counts::read`] reads a 10x directory as the tools read it: every stored entry kept, explicit
//! zeros included ([`crate::sparse::read_mtx`]), the barcodes, and the features file's lines as
//! they are. [`Counts::write`] writes it as a new gzipped 10x directory
//! ([`crate::tenx::write_columns`]), carrying the feature lines over unchanged. A tool that keeps
//! cells (DoubletFinder's singlets) selects them in memory with [`Counts::keep_cells`]: writing
//! that and reading it back gives the same `Counts`, so the next tool can take it from memory
//! instead of from the files.

use crate::sparse::{read_mtx, CscMatrix};
use crate::tenx::{self, Feature, Written};
use std::io;
use std::path::Path;

/// Genes x cells counts with their identifiers.
#[derive(Debug, Clone, PartialEq)]
pub struct Counts {
    /// Genes x cells, in R's `dgCMatrix` layout, every stored entry as read.
    pub matrix: CscMatrix,
    /// The cells' barcodes, one per column.
    pub barcodes: Vec<String>,
    /// The features file's lines (id, name, type, tab-separated), one per row, as read: a tool
    /// that writes count data carries them over unchanged.
    pub feature_lines: Vec<String>,
}

fn invalid(m: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m)
}

impl Counts {
    /// Checks that the identifiers fit the matrix.
    pub fn new(matrix: CscMatrix, barcodes: Vec<String>, feature_lines: Vec<String>) -> io::Result<Counts> {
        if matrix.ncol != barcodes.len() || matrix.nrow != feature_lines.len() {
            return Err(invalid(format!("the matrix is {} x {}, with {} features and {} barcodes", matrix.nrow, matrix.ncol, feature_lines.len(), barcodes.len())));
        }
        Ok(Counts { matrix, barcodes, feature_lines })
    }

    /// Reads a 10x Matrix Market directory (`matrix.mtx(.gz)`, `barcodes.tsv(.gz)`,
    /// `features.tsv(.gz)`).
    pub fn read(dir: &Path) -> io::Result<Counts> {
        let barcodes = tenx::read_barcodes(dir)?;
        let feature_lines = tenx::read_lines(&tenx::find(dir, "features.tsv")?)?;
        let matrix = read_mtx(dir)?;
        Counts::new(matrix, barcodes, feature_lines).map_err(|e| invalid(format!("{}: {e}", dir.display())))
    }

    pub fn n_cells(&self) -> usize {
        self.matrix.ncol
    }

    pub fn n_features(&self) -> usize {
        self.matrix.nrow
    }

    /// The features (id, name, type), as [`tenx::read_features`] reads them.
    pub fn features(&self) -> Vec<Feature> {
        self.feature_lines.iter().map(|l| Feature::parse(l)).collect()
    }

    /// Keeps the cells whose `keep` is true, in their order, each column exactly as stored. The
    /// matrix is compacted in place, so no second copy of the entries is made.
    pub fn keep_cells(mut self, keep: &[bool]) -> Counts {
        assert_eq!(keep.len(), self.matrix.ncol, "{} flags for {} cells", keep.len(), self.matrix.ncol);
        let m = &mut self.matrix;
        let mut w = 0usize;
        let mut p = Vec::with_capacity(keep.iter().filter(|&&k| k).count() + 1);
        p.push(0usize);
        for (j, &k) in keep.iter().enumerate() {
            let (a, b) = (m.p[j], m.p[j + 1]);
            if !k {
                continue;
            }
            // w <= a: a column moves towards the start, never over entries not yet moved
            m.i.copy_within(a..b, w);
            m.x.copy_within(a..b, w);
            w += b - a;
            p.push(w);
        }
        m.i.truncate(w);
        m.x.truncate(w);
        m.i.shrink_to_fit();
        m.x.shrink_to_fit();
        m.ncol = p.len() - 1;
        m.p = p;
        let mut k = keep.iter();
        self.barcodes.retain(|_| *k.next().unwrap());
        self
    }

    /// Writes the counts as a gzipped 10x directory `dir` (created if needed); `software` goes in
    /// the matrix's `%metadata_json` line. The bytes do not depend on `threads`.
    pub fn write(&self, dir: &Path, software: &str, threads: usize) -> io::Result<Written> {
        let cols: Vec<usize> = (0..self.matrix.ncol).collect();
        tenx::write_columns(dir, &self.matrix, &cols, &self.barcodes, &self.feature_lines, software, threads)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder/soupx_tenx_test/input")
    }

    #[test]
    fn keeping_cells_in_place_is_selecting_columns() {
        let c = Counts::read(&fixture()).unwrap();
        assert_eq!((c.n_features(), c.n_cells()), (2000, 300));
        let keep: Vec<bool> = (0..c.n_cells()).map(|j| j % 3 != 1 && j != 299).collect();
        let cols: Vec<usize> = (0..c.n_cells()).filter(|&j| keep[j]).collect();
        let want = c.matrix.select_cols(&cols);
        let bc: Vec<String> = cols.iter().map(|&j| c.barcodes[j].clone()).collect();
        let kept = c.clone().keep_cells(&keep);
        kept.matrix.validate().unwrap();
        assert_eq!(kept.matrix, want);
        assert_eq!(kept.barcodes, bc);
        assert_eq!(kept.feature_lines, c.feature_lines);
        // the explicit zeros stay
        assert_eq!(kept.matrix.x.iter().filter(|&&v| v == 0.0).count(), want.x.iter().filter(|&&v| v == 0.0).count());
        // none kept, all kept
        assert_eq!(c.clone().keep_cells(&vec![false; 300]).matrix, CscMatrix::zeros(2000, 0));
        assert_eq!(c.clone().keep_cells(&vec![true; 300]), c);
    }

    #[test]
    fn written_and_read_back_is_the_same() {
        let c = Counts::read(&fixture()).unwrap();
        let keep: Vec<bool> = (0..c.n_cells()).map(|j| j % 2 == 0).collect();
        let kept = c.keep_cells(&keep);
        let dir = crate::test_dir("counts_round_trip");
        let _ = std::fs::remove_dir_all(&dir);
        let w = kept.write(&dir, "scio test", 3).unwrap();
        assert_eq!((w.nrow, w.ncol, w.nnz, w.field), (2000, 150, kept.matrix.nnz(), "integer"));
        let back = Counts::read(&dir).unwrap();
        assert_eq!(back, kept);
        assert_eq!(back.features(), tenx::read_features(&dir).unwrap());
    }
}
