//! DoubletFinder 2.0.6 (commit 1b244d8), reimplemented without Seurat objects
//! (docs/PLAN.md Part B; milestones M3 to M7).
//!
//! It will reproduce the workflow's `doubletfinder.R` call sequence: `paramSweep` over the
//! pN x pK grid, `summarizeSweep` and `find.pK` (the bimodality coefficient), nExp with the
//! homotypic adjustment, and the final `doubletFinder` classification. It reads a counts matrix
//! and cluster labels, and writes pANN and the classification for every cell, the sweep tables,
//! and a `summary.json` that includes the column names R would write. It builds on the `rcompat`,
//! `scio`, `sctransform`, `pca`, `knn` and `kde` crates, and runs as `scrust doubletfinder`.

/// The original implementation this crate reproduces (docs/PLAN.md B1).
pub const REFERENCE: &str = "DoubletFinder 2.0.6 (commit 1b244d8)";

#[cfg(test)]
mod tests {
    /// Placeholder until M3: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "doubletfinder");
        assert!(super::REFERENCE.contains("1b244d8"));
    }
}
