//! SCTransform v2, ported from sctransform 0.4.3 with glmGamPoi 1.22.0's fitting, as Seurat
//! 5.5.1's `SCTransform()` calls it with its defaults (docs/PLAN.md B3a; milestones M4a to M4d).
//!
//! The call chain is `SCTransform.Seurat` → `.StdAssay` → `.Assay` → `.default` (three of which
//! call `set.seed(1448145)`), then `sctransform::vst(vst.flavor = "v2")`, which fits
//! `glmGamPoi::glm_gp(~ 1, offset = log(umi))` per gene. The port follows that code, arithmetic
//! operation by arithmetic operation, so that its results are R's bit for bit. Each stage can
//! start from R's own intermediates (teacher forcing), so a difference shows up in the stage that
//! caused it:
//!
//! - [`setup`] (S2a): cell attributes, the gene filter, geometric means, the fitting cells and
//!   genes (density-weighted sampling), the bins;
//! - S2b fits: glmGamPoi's Newton fits, the overdispersion MLE through R's `nlminb` with the
//!   Cox-Reid adjustment, `loc_median_fit`, the theta rules;
//! - S2c regularisation: outlier removal and `ksmooth` against the genes' geometric means;
//! - S2d residuals: clipped Pearson residuals, residual variance, the variable features and
//!   `scale.data`.
//!
//! Shared by DoubletFinder (7 calls per object) and, later, the clustering crate.

pub mod matrix;
pub mod params;
pub mod setup;

pub use matrix::GeneRows;
pub use params::{MinVariance, SctParams};
pub use setup::{setup, CellAttr, Setup};

use std::fmt;

/// The original implementation this crate reproduces.
pub const REFERENCE: &str = "sctransform 0.4.3 (vst.flavor = \"v2\") with glmGamPoi 1.22.0, as Seurat 5.5.1's SCTransform calls it";

/// Why a run stops.
#[derive(Debug)]
pub enum SctError {
    /// Invalid input or parameters, or a condition under which R stops.
    Input(String),
}

impl fmt::Display for SctError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SctError::Input(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for SctError {}
