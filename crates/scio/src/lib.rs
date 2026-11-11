//! Sparse matrix and table I/O shared by every tool (docs/PLAN.md A3 and A4).
//!
//! Tools read and write plain files, never Seurat objects:
//!
//! - [`tenx`]: 10x Matrix Market directories (`matrix.mtx`, `barcodes.tsv`, `features.tsv`,
//!   optionally gzipped, through [`gz`], a standard-library gzip reader);
//! - [`sparse`]: the counts matrix in R's `dgCMatrix` layout, read from `matrix.mtx` with every
//!   stored entry kept, explicit zeros included (M4);
//! - [`tsv`]: per-cell and per-parameter tables (cluster labels in; `cells.tsv` and the sweep
//!   tables out);
//! - [`text`]: doubles as `%.17g` text, as the R reference scripts write them, and back;
//! - [`json`]: `summary.json` in the R scripts' layout, and a small parser;
//! - [`npy`]: `.npy` arrays, the format of the R scripts' dense intermediates and RNG
//!   checkpoints.

pub mod gz;
pub mod json;
pub mod npy;
pub mod sparse;
pub mod tenx;
pub mod text;
pub mod tsv;

pub use json::Json;
pub use npy::Npy;
pub use sparse::{CscMatrix, CsrMatrix};
pub use text::{fmt17, parse17};
pub use tsv::Table;

/// A scratch directory under the workspace's `target/` for unit tests.
#[cfg(test)]
pub(crate) fn test_dir(name: &str) -> std::path::PathBuf {
    let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/test-tmp/scio").join(name);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    #[test]
    fn reads_a_tier_t_input_directory() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder/soupx_tenx_test/input");
        let barcodes = super::tenx::read_barcodes(&dir).unwrap();
        let features = super::tenx::read_features(&dir).unwrap();
        let h = super::tenx::read_mtx_header(&dir).unwrap();
        assert_eq!((barcodes.len(), features.len()), (300, 2000));
        assert_eq!((h.nrow, h.ncol, h.nnz, h.field.as_str()), (2000, 300, 64934, "integer"));
        assert_eq!(barcodes[0], "AAACGAAAGGATTCCT-1");
    }
}
