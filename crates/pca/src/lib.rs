//! Truncated PCA, a port of irlba 2.4.1 as Seurat 5.5.1's `RunPCA` calls it
//! (docs/PLAN.md B3b; milestone M5).
//!
//! It will reproduce `set.seed(42)` followed by `irlba(t(x), nv = npcs)` with irlba's defaults:
//! the Lanczos bidiagonalisation loop with its convergence tests and restarts, every matrix
//! product through the same OpenBLAS routine R uses, and the small SVD through LAPACK `dgesdd`.
//! It returns the embeddings (`u %*% diag(d)`) and standard deviations, and consumes the same
//! random draws as R, so the stream continues correctly afterwards.

#[cfg(test)]
mod tests {
    /// Placeholder until M5: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "pca");
    }
}
