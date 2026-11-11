//! Sparse matrix and table I/O shared by every tool (docs/PLAN.md A3 and A4).
//!
//! Tools read and write plain files, never Seurat objects. This crate will hold:
//!
//! - 10x Matrix Market directories (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, optionally
//!   gzipped), read and written; 10x HDF5 later;
//! - per-cell and per-parameter tables as TSV: cluster labels in; `cells.tsv` and sweep tables
//!   out;
//! - run metadata as `summary.json`: parameters, seed, versions, timings and the environment
//!   fingerprint (A5b);
//! - `.npy` arrays, the format the R export scripts use for dense intermediates and RNG
//!   checkpoints.

#[cfg(test)]
mod tests {
    /// Placeholder until the first reader lands: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "scio");
    }
}
