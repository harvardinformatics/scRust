//! The whole of `SCTransform()` as the port reproduces it, stage after stage on its own outputs
//! (S2a to S2d).
//!
//! [`sctransform`] runs it in one call. [`start`] and [`finish`] split it where the counts stop
//! being read: `start` makes `SCTransform`'s `set.seed` calls and stage S2a, which keeps the
//! rows of the genes `vst` keeps ([`GeneRows`]); `finish` runs S2b to S2d on those rows. A caller
//! that owns a large input (DoubletFinder's counts with their artificial doublets) can free it
//! in between. The results are the same either way.

use crate::correct::Corrected;
use crate::fit::{fit, FitInput, Fits};
use crate::matrix::GeneRows;
use crate::params::SctParams;
use crate::reg::{mean_cell_sum, regularize, Reg, RegInput};
use crate::regress::Regression;
use crate::residuals::{residuals_with, ResidualInput, Residuals};
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
    /// `correct()`'s corrected counts, when asked for ([`finish_with`]; M13)
    pub corrected: Option<Corrected>,
    /// `vars.to.regress`' regression of `scale.data`, when asked for ([`finish_full`]; M15b)
    pub regression: Option<Regression>,
    /// `.Random.seed` when `vst` starts, after `SCTransform`'s `set.seed` calls (R's checkpoint
    /// "vst_entry")
    pub rng_entry: Vec<i32>,
    /// `.Random.seed` when `vst` (and `SCTransform`) return: the fits' exit state
    pub rng_exit: Vec<i32>,
}

/// What [`start`] hands to [`finish`]: stage S2a, the kept genes' rows and the state at `vst`'s
/// entry.
#[derive(Debug, Clone)]
pub struct Started {
    pub setup: Setup,
    pub rows: GeneRows,
    pub rng_entry: Vec<i32>,
}

/// `SCTransform()` on a genes x cells counts matrix: `set.seed(seed_use)` (when set) on `rng`,
/// then `vst(vst.flavor = "v2")`, the variable features and `scale.data`. `rng` ends in R's state
/// at `SCTransform`'s exit. `keep_y` keeps the full residual matrix (small inputs only).
pub fn sctransform(counts: &scio::CscMatrix, params: &SctParams, rng: &mut RRng, keep_y: bool, threads: usize) -> Result<Sct, SctError> {
    let started = start(counts, params, rng)?;
    finish(started, params, rng, keep_y, threads)
}

/// The part of [`sctransform`] that reads the counts: `SCTransform.Seurat`, `.Assay` and
/// `.default` each call `set.seed(seed.use)`, then `vst`'s setup (S2a: cell attributes, the gene
/// filter, gene statistics, the fitting cells and genes), which keeps the rows of the genes it
/// keeps. The counts are not needed after it.
pub fn start(counts: &scio::CscMatrix, params: &SctParams, rng: &mut RRng) -> Result<Started, SctError> {
    params.check().map_err(SctError::Input)?;
    if let Some(s) = params.seed_use {
        rng.set_seed(s);
    }
    let rng_entry = rng.random_seed();
    let (setup, rows) = setup(counts, params, rng)?;
    Ok(Started { setup, rows, rng_entry })
}

/// The rest of [`sctransform`] after [`start`]: the fits (S2b), the regularisation (S2c), the
/// residuals, the variable features and `scale.data` (S2d).
pub fn finish(started: Started, params: &SctParams, rng: &mut RRng, keep_y: bool, threads: usize) -> Result<Sct, SctError> {
    finish_with(started, params, rng, keep_y, false, threads)
}

/// [`finish`], and with `correct` also the corrected counts of `correct()`, which `vst()` makes
/// from its residuals before clipping them (`do.correct.umi = TRUE`, Seurat's default; M13). They
/// change nothing else: `correct()` draws no random number and its result feeds no later stage.
pub fn finish_with(started: Started, params: &SctParams, rng: &mut RRng, keep_y: bool, correct: bool, threads: usize) -> Result<Sct, SctError> {
    finish_full(started, params, rng, keep_y, correct, None, threads)
}

/// [`finish_with`], and with `regress` SCTransform's second regression (`vars.to.regress`,
/// milestone M15b): `scale.data`'s clipped rows regressed on the covariates before they are
/// centred ([`crate::regress`]). It draws no random number and changes nothing but `scale.data`.
pub fn finish_full(started: Started, params: &SctParams, rng: &mut RRng, keep_y: bool, correct: bool, regress: Option<Regression>, threads: usize) -> Result<Sct, SctError> {
    let Started { setup, rows, rng_entry } = started;
    if let Some(r) = &regress {
        if r.qr.n != setup.cell_attr.log_umi.len() {
            return Err(SctError::Input(format!("vars.to.regress: covariates for {} cells, but the counts have {}", r.qr.n, setup.cell_attr.log_umi.len())));
        }
    }
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
    let input = ResidualInput { model_pars_fit: &reg.model_pars_fit, log_umi: &setup.cell_attr.log_umi };
    let (residuals, corrected) = residuals_with(&rows, &input, params, keep_y, threads, correct, regress.as_ref())?;
    let rng_exit = rng.random_seed();
    Ok(Sct { setup, rows, fits, mean_cell_sum: mcs, reg, residuals, corrected, regression: regress, rng_entry, rng_exit })
}
