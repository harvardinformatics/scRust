//! The whole of `SCTransform()` as the port reproduces it, stage after stage on its own outputs
//! (S2a to S2d).

use crate::fit::{fit, FitInput, Fits};
use crate::matrix::GeneRows;
use crate::params::SctParams;
use crate::reg::{mean_cell_sum, regularize, Reg, RegInput};
use crate::residuals::{residuals, ResidualInput, Residuals};
use crate::setup::{setup, Setup};
use crate::SctError;
use rcompat::RRng;

/// Everything one `SCTransform()` call computes that the port reproduces.
#[derive(Debug, Clone)]
pub struct Sct {
    pub setup: Setup,
    pub rows: GeneRows,
    pub fits: Fits,
    pub mean_cell_sum: f64,
    pub reg: Reg,
    pub residuals: Residuals,
    /// `.Random.seed` when `vst` (and `SCTransform`) return: the fits' exit state
    pub rng_exit: Vec<i32>,
}

/// `SCTransform()` on a genes x cells counts matrix: `set.seed(seed_use)` (when set) on `rng`,
/// then `vst(vst.flavor = "v2")`, the variable features and `scale.data`. `rng` ends in R's state
/// at `SCTransform`'s exit. `keep_y` keeps the full residual matrix (small inputs only).
pub fn sctransform(counts: &scio::CscMatrix, params: &SctParams, rng: &mut RRng, keep_y: bool, threads: usize) -> Result<Sct, SctError> {
    params.check().map_err(SctError::Input)?;
    // SCTransform.Seurat, .StdAssay and .default each call set.seed(seed.use)
    if let Some(s) = params.seed_use {
        rng.set_seed(s);
    }
    let (setup, rows) = setup(counts, params, rng)?;
    let fits = fit(&rows, &FitInput::from_setup(&setup), params, rng, threads)?;
    let mcs = mean_cell_sum(&rows);
    let reg = regularize(
        &RegInput {
            model_pars: &fits.model_pars,
            genes_step1: &setup.genes_step1,
            log_gmean_step1: &setup.log_gmean_step1,
            log_gmean: &setup.log_gmean,
            amean: &setup.amean,
            var: &setup.var,
            mean_cell_sum: mcs,
        },
        params,
    )?;
    let residuals = residuals(&rows, &ResidualInput { model_pars_fit: &reg.model_pars_fit, log_umi: &setup.cell_attr.log_umi }, params, keep_y, threads)?;
    let rng_exit = rng.random_seed();
    Ok(Sct { setup, rows, fits, mean_cell_sum: mcs, reg, residuals, rng_exit })
}
