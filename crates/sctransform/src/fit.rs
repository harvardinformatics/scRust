//! Stage S2b, sctransform 0.4.3's `get_model_pars` for `method = "glmGamPoi_offset"` and the
//! theta rules after it (`R/vst.R`):
//!
//! - for each bin of `bin_size` fitting genes, `future_lapply(..., future.seed = TRUE)` over one
//!   chunk (a sequential plan has one worker), which advances R's main random-number stream by
//!   exactly two draws: future 1.82.0's `make_rng_seeds` and future.apply 1.21.0's
//!   `future_xapply` each call `next_random_seed()` (`sample.int(1)`, one uniform) and restore
//!   that forwarded state on exit; the L'Ecuyer-CMRG kind switch in between is undone. Then
//!   `fit_glmGamPoi_offset`: [`crate::glm::glm_gp_offset`] on the bin's counts over the fitting
//!   cells with offset `log(10^log_umi)`, `theta = 1 / overdispersion`;
//! - `exclude_poisson`: genes whose moment estimate of theta, `amean^2 / (var - amean)`, is below
//!   1/1000 of the fitted one get `theta = Inf`. The rule goes through
//!   `rownames(model_pars[diff_theta < 1e-3, ])`, which is NULL when exactly one gene qualifies
//!   (the one-row result drops to a vector): that gene keeps its theta, as in R;
//! - thetas below 1e-7 are raised to 1e-7.

use crate::glm::{glm_gp_offset, GeneGlm};
use crate::matrix::{positions, GeneRows};
use crate::setup::{bins, Setup};
use crate::params::SctParams;
use crate::SctError;
use rcompat::libm::{log, pow};
use rcompat::RRng;

/// What S2b reads: indices into the kept genes (`genes`) and the vst input's cells.
#[derive(Debug, Clone)]
pub struct FitInput {
    /// fitting genes, positions in the kept genes, in draw order
    pub genes_step1: Vec<usize>,
    /// fitting cells, in draw order
    pub cells_step1: Vec<usize>,
    /// `log10(colSums(umi))` of every cell
    pub log_umi: Vec<f64>,
    /// `rowMeans(umi)` and `row_var(umi)` of the fitting genes
    pub amean_step1: Vec<f64>,
    pub var_step1: Vec<f64>,
}

impl FitInput {
    pub fn from_setup(s: &Setup) -> FitInput {
        FitInput {
            genes_step1: s.genes_step1.clone(),
            cells_step1: s.cells_step1.clone(),
            log_umi: s.cell_attr.log_umi.clone(),
            amean_step1: s.genes_step1.iter().map(|&k| s.amean[k]).collect(),
            var_step1: s.genes_step1.iter().map(|&k| s.var[k]).collect(),
        }
    }
}

/// The result of S2b; one entry per fitting gene, in draw order.
#[derive(Debug, Clone)]
pub struct Fits {
    pub glm: Vec<GeneGlm>,
    /// `get_model_pars`' result: theta (after the Inf rule), `(Intercept)`, `log_umi` (= log(10))
    pub model_pars_raw: Vec<[f64; 3]>,
    /// after the 1e-7 floor
    pub model_pars: Vec<[f64; 3]>,
    /// `.Random.seed` at each bin's `future_lapply` entry, and at the fits' exit
    pub rng_bins: Vec<Vec<i32>>,
    pub rng_exit: Vec<i32>,
}

/// Stage S2b. `rows` are the kept genes' counts; `rng` is R's state at `get_model_pars` entry.
pub fn fit(rows: &GeneRows, input: &FitInput, params: &SctParams, rng: &mut RRng, threads: usize) -> Result<Fits, SctError> {
    let n1 = input.genes_step1.len();
    let nc = input.cells_step1.len();
    let pos = positions(&input.cells_step1, rows.ncol);
    // log_umi <- log(10^log10_umi) over data_step1 = cell_attr[cells_step1, ]
    let offset: Vec<f64> = input.cells_step1.iter().map(|&c| log(pow(10.0, input.log_umi[c]))).collect();
    let mut glm = Vec::with_capacity(n1);
    let mut rng_bins = Vec::new();
    for b in bins(n1, params.bin_size) {
        rng_bins.push(rng.random_seed());
        // future.seed = TRUE: make_rng_seeds() and future_xapply() each call next_random_seed()
        rng.unif_index(1.0);
        rng.unif_index(1.0);
        let ng = b.len();
        let mut y = vec![0.0f64; ng * nc];
        for (k, gi) in b.clone().enumerate() {
            rows.dense_row_cells(input.genes_step1[gi], &pos, &mut y[k * nc..(k + 1) * nc]);
        }
        let res = glm_gp_offset(&y, ng, nc, &offset, threads).map_err(|e| SctError::Input(format!("bin {}: {e}", rng_bins.len())))?;
        glm.extend(res);
    }
    let rng_exit = rng.random_seed();
    let log10 = log(10.0);
    let mut model_pars_raw: Vec<[f64; 3]> = glm.iter().map(|g| [1.0 / g.disp_mle, g.beta2, log10]).collect();
    // exclude_poisson: theta = Inf where amean^2 / (var - amean) / theta < 1e-3
    let diff: Vec<f64> = (0..n1)
        .map(|k| {
            let a = input.amean_step1[k];
            let predicted = a * a / (input.var_step1[k] - a);
            predicted / model_pars_raw[k][0]
        })
        .collect();
    if diff.iter().any(|d| d.is_nan()) {
        return Err(SctError::Input("theta_mm / theta_mle is NaN for a fitting gene (R stops or misbehaves here)".into()));
    }
    let flagged: Vec<usize> = (0..n1).filter(|&k| diff[k] < 1e-3).collect();
    if flagged.len() >= 2 {
        for &k in &flagged {
            model_pars_raw[k][0] = f64::INFINITY;
        }
    }
    let mut model_pars = model_pars_raw.clone();
    if model_pars.iter().any(|p| p[0] < 1e-7) {
        for p in model_pars.iter_mut() {
            // pmax(theta, 1e-7)
            if !p[0].is_nan() && 1e-7 > p[0] {
                p[0] = 1e-7;
            }
        }
    }
    Ok(Fits { glm, model_pars_raw, model_pars, rng_bins, rng_exit })
}
