//! Stage S2c, sctransform 0.4.3's `reg_model_pars` as `vst(vst.flavor = "v2")` calls it
//! (`theta_regularization = "od_factor"`, `exclude_poisson = TRUE`, no batch variable, neither
//! slope nor intercept fixed; `R/vst.R` and `R/utils.R`):
//!
//! 1. Poisson genes: kept genes whose variance does not exceed their mean, then those with mean
//!    below 1e-3 (`union()`, so in that order); their parameters come from the offset model,
//!    `theta = Inf`, `(Intercept) = log(amean) - log(mean(colSums(umi)))`, `log_umi = log(10)`;
//! 2. the fitting genes' dispersion parameter `log10(1 + 10^log_gmean / theta)`, next to their
//!    intercept and slope; outliers in any of the three against the log geometric means
//!    (`is_outlier`: robust z-scores within bins of width `range * bw.SJ / 2`, on two shifted
//!    grids, the smaller one above 10), and genes with an infinite theta, are left out, then the
//!    Poisson genes;
//! 3. each parameter is smoothed against the remaining genes' log geometric means with a normal
//!    kernel (`ksmooth`, bandwidth `bw.SJ * bw_adjust`) at every kept gene's log geometric mean,
//!    clamped to their range;
//! 4. the dispersion parameter back to theta, `10^log_gmean / (10^dispersion_par - 1)`, and the
//!    Poisson genes' parameters put in.
//!
//! In R, `poisson_genes_step1` and `overdispersed_genes_step1` only feed messages, so they are not
//! computed here.

use crate::matrix::GeneRows;
use crate::params::SctParams;
use crate::SctError;
use rcompat::bandwidth::bw_sj;
use rcompat::libm::{log, log10};
use rcompat::quantile::{mad, median};
use rcompat::smooth::{bincode, ksmooth, Kernel};
use rcompat::{order, r_pow, seq_by};

/// What S2c reads.
#[derive(Debug, Clone)]
pub struct RegInput<'a> {
    /// `model_pars` after the 1e-7 floor (theta, `(Intercept)`, `log_umi`), per fitting gene
    pub model_pars: &'a [[f64; 3]],
    /// fitting genes, positions in the kept genes, in draw order
    pub genes_step1: &'a [usize],
    pub log_gmean_step1: &'a [f64],
    /// per kept gene: `log10(row_gmean)`, `rowMeans`, `row_var`
    pub log_gmean: &'a [f64],
    pub amean: &'a [f64],
    pub var: &'a [f64],
    /// `mean(colSums(umi))` of the kept genes' counts ([`mean_cell_sum`])
    pub mean_cell_sum: f64,
}

/// The result of S2c.
#[derive(Debug, Clone)]
pub struct Reg {
    /// per kept gene: theta, `(Intercept)`, `log_umi` (`model_pars_fit`)
    pub model_pars_fit: Vec<[f64; 3]>,
    /// per fitting gene: an outlier (or an infinite theta)
    pub outliers: Vec<bool>,
    /// `all_poisson_genes`, positions in the kept genes, in R's order
    pub poisson_genes: Vec<usize>,
    /// the fitting genes whose parameters are smoothed, in order
    pub smoothing_genes: Vec<usize>,
    /// `bw.SJ(log_gmean of the smoothing genes) * bw_adjust`
    pub bw: f64,
}

/// `mean(colSums(umi))` over the kept genes: Matrix's `colSums` (each column's stored values in
/// row order, in double), then base R's `mean` (long double).
pub fn mean_cell_sum(rows: &GeneRows) -> f64 {
    let mut cs = vec![0.0f64; rows.ncol];
    for k in 0..rows.nrow() {
        let (cols, vals) = rows.row(k);
        for (&j, &v) in cols.iter().zip(vals) {
            cs[j as usize] += v;
        }
    }
    rcompat::mean(&cs, false)
}

fn r_min(x: &[f64]) -> f64 {
    let mut s = f64::INFINITY;
    for &v in x {
        if v.is_nan() {
            return v;
        }
        if v < s {
            s = v;
        }
    }
    s
}

fn r_max(x: &[f64]) -> f64 {
    let mut s = f64::NEG_INFINITY;
    for &v in x {
        if v.is_nan() {
            return v;
        }
        if v > s {
            s = v;
        }
    }
    s
}

/// `robust_scale_binned(y, x, breaks)`: within each bin of `cut(x, breaks)`, `(y - median(y)) /
/// (mad(y) + .Machine$double.eps)`. (`aggregate` returns the bins in order and each bin's values
/// in their order, which `score[order(bins)] <- unlist(...)` puts back in place.)
fn robust_scale_binned(y: &[f64], x: &[f64], breaks: &[f64]) -> Result<Vec<f64>, SctError> {
    if breaks.windows(2).any(|w| w[0] == w[1]) {
        return Err(SctError::Input("is_outlier: 'breaks' are not unique".into()));
    }
    let codes = bincode(x, breaks, true, false);
    let nbins = breaks.len() - 1;
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); nbins];
    for (i, c) in codes.iter().enumerate() {
        match c {
            Some(b) => members[b - 1].push(i),
            // R's aggregate() drops the NA bin and the scores no longer line up
            None => return Err(SctError::Input("is_outlier: a value falls outside the bins".into())),
        }
    }
    let mut score = vec![0.0f64; x.len()];
    for m in members.iter().filter(|m| !m.is_empty()) {
        let v: Vec<f64> = m.iter().map(|&i| y[i]).collect();
        let med = median(&v);
        let denom = mad(&v) + f64::EPSILON;
        for &i in m {
            score[i] = (y[i] - med) / denom;
        }
    }
    Ok(score)
}

/// sctransform's `is_outlier(y, x, th = 10)`.
pub fn is_outlier(y: &[f64], x: &[f64]) -> Result<Vec<bool>, SctError> {
    let bw_x = bw_sj(x).map_err(|e| SctError::Input(format!("is_outlier: bw.SJ: {e}")))?;
    is_outlier_bw(y, x, bw_x)
}

/// `is_outlier(y, x, th = 10)`, with `bw.SJ(x)` given (it depends on `x` only).
fn is_outlier_bw(y: &[f64], x: &[f64], bw_x: f64) -> Result<Vec<bool>, SctError> {
    let (mn, mx) = (r_min(x), r_max(x));
    let bin_width = (mx - mn) * bw_x / 2.0;
    let eps = f64::EPSILON * 10.0;
    let seq = |from: f64| seq_by(from, mx + bin_width, bin_width).map_err(|e| SctError::Input(format!("is_outlier: seq(): {e}")));
    let breaks1 = seq(mn - eps)?;
    let breaks2 = seq(mn - eps - bin_width / 2.0)?;
    let s1 = robust_scale_binned(y, x, &breaks1)?;
    let s2 = robust_scale_binned(y, x, &breaks2)?;
    Ok(s1
        .iter()
        .zip(&s2)
        .map(|(&a, &b)| {
            // pmin(abs(score1), abs(score2)) > th
            let (a, b) = (a.abs(), b.abs());
            let m = if b < a { b } else { a };
            m > 10.0
        })
        .collect())
}

/// Stage S2c.
pub fn regularize(input: &RegInput, params: &SctParams) -> Result<Reg, SctError> {
    let ng = input.log_gmean.len();
    let n1 = input.genes_step1.len();
    assert!(input.amean.len() == ng && input.var.len() == ng && input.model_pars.len() == n1 && input.log_gmean_step1.len() == n1);
    // 1. all_poisson_genes <- union(genes[var - amean <= 0], genes[amean < 1e-3])
    let mut is_poisson = vec![false; ng];
    let mut poisson_genes = Vec::new();
    for g in 0..ng {
        let od = input.var[g] - input.amean[g];
        if od.is_nan() {
            return Err(SctError::Input("reg_model_pars: NaN overdispersion factor".into()));
        }
        if od <= 0.0 {
            is_poisson[g] = true;
            poisson_genes.push(g);
        }
    }
    for g in 0..ng {
        if input.amean[g] < 1e-3 && !is_poisson[g] {
            is_poisson[g] = true;
            poisson_genes.push(g);
        }
    }
    // 2. dispersion_par <- log10(1 + 10^genes_log_gmean_step1 / theta); the parameters' columns
    let cols: [Vec<f64>; 3] = [
        (0..n1).map(|k| log10(1.0 + r_pow(10.0, input.log_gmean_step1[k]) / input.model_pars[k][0])).collect(),
        input.model_pars.iter().map(|p| p[1]).collect(),
        input.model_pars.iter().map(|p| p[2]).collect(),
    ];
    // outliers in any column, and theta = Inf
    let bw_outlier = bw_sj(input.log_gmean_step1).map_err(|e| SctError::Input(format!("is_outlier: bw.SJ: {e}")))?;
    let mut outliers = vec![false; n1];
    for col in &cols {
        for (o, flag) in outliers.iter_mut().zip(is_outlier_bw(col, input.log_gmean_step1, bw_outlier)?) {
            *o |= flag;
        }
    }
    for (k, o) in outliers.iter_mut().enumerate() {
        *o |= !input.model_pars[k][0].is_finite();
    }
    // the smoothing genes: not outliers, then not Poisson genes (setdiff keeps the order)
    let keep: Vec<usize> = (0..n1).filter(|&k| !outliers[k] && !is_poisson[input.genes_step1[k]]).collect();
    if keep.len() < 2 {
        return Err(SctError::Input(format!("reg_model_pars: {} genes left to smooth", keep.len())));
    }
    let smoothing_genes: Vec<usize> = keep.iter().map(|&k| input.genes_step1[k]).collect();
    let lg: Vec<f64> = keep.iter().map(|&k| input.log_gmean_step1[k]).collect();
    // 3. bw <- bw.SJ(genes_log_gmean_step1) * bw_adjust
    let bw = bw_sj(&lg).map_err(|e| SctError::Input(format!("reg_model_pars: bw.SJ: {e}")))? * params.bw_adjust;
    // x_points <- pmin(pmax(genes_log_gmean, min(lg)), max(lg))
    let (mn, mx) = (r_min(&lg), r_max(&lg));
    let x_points: Vec<f64> = input
        .log_gmean
        .iter()
        .map(|&v| {
            let v = if mn > v { mn } else { v };
            if mx < v {
                mx
            } else {
                v
            }
        })
        .collect();
    let o = order(&x_points, false);
    let mut fit = vec![[f64::NAN; 3]; ng];
    for (c, col) in cols.iter().enumerate() {
        let y: Vec<f64> = keep.iter().map(|&k| col[k]).collect();
        let (_, yp) = ksmooth(&lg, &y, Kernel::Normal, bw, &x_points);
        for (i, &g) in o.iter().enumerate() {
            fit[g][c] = yp[i];
        }
    }
    // 4. Poisson genes: dispersion_par 0; theta <- 10^log_gmean / (10^dispersion_par - 1)
    for &g in &poisson_genes {
        fit[g][0] = 0.0;
    }
    let log_mcs = log(input.mean_cell_sum);
    let log10_ = log(10.0);
    let mut model_pars_fit: Vec<[f64; 3]> = (0..ng)
        .map(|g| [r_pow(10.0, input.log_gmean[g]) / (r_pow(10.0, fit[g][0]) - 1.0), fit[g][1], fit[g][2]])
        .collect();
    for &g in &poisson_genes {
        model_pars_fit[g] = [f64::INFINITY, log(input.amean[g]) - log_mcs, log10_];
    }
    Ok(Reg { model_pars_fit, outliers, poisson_genes, smoothing_genes, bw })
}
