//! SCTransform v2, ported from sctransform 0.4.3 with glmGamPoi 1.22.0's fitting
//! (docs/PLAN.md B3a; milestones M4a to M4d).
//!
//! It will reproduce `sctransform::vst(vst.flavor = "v2")` as Seurat 5.5.1's `SCTransform` calls
//! it, through a typed API along the lines of `vst(&counts, &SctParams, &mut RRng) -> SctFit`,
//! where `SctParams` mirrors R's arguments. Each stage can also start from R's own intermediates,
//! so a difference shows up in the stage that caused it:
//!
//! - S2a setup: `log_umi`, the gene filter, geometric means, fitting cells and genes, and the
//!   500-gene bins;
//! - S2b fits: glmGamPoi's Newton fits, the overdispersion MLE through R's `nlminb` with the
//!   Cox-Reid adjustment, and `loc_median_fit`;
//! - S2c regularisation: outlier removal and `ksmooth` against the genes' geometric means;
//! - S2d residuals: clipped Pearson residuals, residual variance, the variable features and
//!   `scale.data`.
//!
//! Shared by DoubletFinder (7 calls per object) and the clustering crate.

#[cfg(test)]
mod tests {
    /// Placeholder until M4a: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "sctransform");
    }
}
