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
//! - [`fit`] (S2b): glmGamPoi's Newton fits, the overdispersion MLE through R's `nlminb` with
//!   the Cox-Reid adjustment, `loc_median_fit`, the theta rules;
//! - [`reg`] (S2c): the Poisson genes, outlier removal and `ksmooth` against the genes'
//!   geometric means;
//! - [`residuals`] (S2d): clipped Pearson residuals, residual variance, the variable features and
//!   `scale.data`;
//! - [`correct`] (K10, milestone M13): `correct()`'s corrected counts, made in S2d's pass when
//!   asked for, their `log1p`, `nCount_SCT` and `nFeature_SCT`;
//! - [`regress`] (milestone M15b): `vars.to.regress`, the second regression of `scale.data`'s
//!   clipped rows on per-cell covariates (`ScaleData`'s `RegressOutMatrix`, through R's LINPACK
//!   QR, `rcompat::linpack`), before they are centred.
//!
//! Shared by DoubletFinder (7 calls per object) and the clustering crate (M11).

pub mod correct;
pub mod fit;
pub mod glm;
pub mod matrix;
pub mod output;
pub mod params;
pub mod record;
pub mod reg;
pub mod regress;
pub mod residuals;
pub mod setup;
pub mod vst;

pub use fit::{fit, FitInput, Fits};
pub use matrix::GeneRows;
pub use params::{MinVariance, SctParams};
pub use reg::{regularize, Reg, RegInput};
pub use correct::Corrected;
pub use regress::Regression;
pub use residuals::{residuals, residuals_corrected, residuals_with, ResidualInput, Residuals};
pub use setup::{setup, CellAttr, Setup};
pub use vst::{finish, finish_full, finish_with, sctransform, start, Sct, Started};

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
