//! `density()` with a Gaussian kernel and no weights (`density.default` in R 4.5.3's
//! `src/library/stats/R/density.R`, its binning `BinDist` in `massdist.c`), and the bandwidth rule
//! `bw.nrd` (`bandwidths.R`). sctransform calls `density(x, bw = "nrd", adjust = 1)` to weight its
//! gene sampling.
//!
//! The steps are R's, in R's order: linear binning onto `2n` points over `[lo, up]` (the upper
//! half stays zero), the kernel evaluated at R's 4.4+ coordinates (`old.coords = FALSE`), the
//! circular convolution through `fft(fft(y) * Conj(fft(kords)), inverse = TRUE)` with R's own FFT
//! and C99 complex multiplication, and linear interpolation of the result onto the output grid.

use crate::approx::{ApproxFun, Method};
use crate::arith::{r_pow, seq_length_out};
use crate::fft::{fft, Complex};
use crate::normal::dnorm;
use crate::quantile::quantile7;
use crate::summary::var;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum DensityError {
    /// fewer than 2 finite points to choose a bandwidth from
    TooFewPoints,
    /// "non-finite 'bw'" or "'bw' is not positive."
    BadBandwidth(f64),
    /// "'x' contains missing values"
    MissingValues,
}

impl fmt::Display for DensityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DensityError::TooFewPoints => f.write_str("need at least 2 points to select a bandwidth automatically"),
            DensityError::BadBandwidth(b) => write!(f, "bandwidth {b} is not positive and finite"),
            DensityError::MissingValues => f.write_str("'x' contains missing values"),
        }
    }
}

impl std::error::Error for DensityError {}

/// `bw.nrd(x)`: `1.06 * min(sqrt(var(x)), IQR(x) / 1.34) * length(x)^(-1/5)`.
pub fn bw_nrd(x: &[f64]) -> Result<f64, DensityError> {
    if x.len() < 2 {
        return Err(DensityError::TooFewPoints);
    }
    let r = quantile7(x, &[0.25, 0.75]);
    let h = (r[1] - r[0]) / 1.34;
    let s = var(x, false).sqrt();
    // min(): NaN-propagating, else the smaller (R's do_summary keeps the first of equal values)
    let m = if s.is_nan() || h.is_nan() { f64::NAN } else if h < s { h } else { s };
    Ok(1.06 * m * r_pow(x.len() as f64, -1.0 / 5.0))
}

/// The result of [`density_nrd`]: `x`, `y` and `bw` of R's `density` object.
#[derive(Debug, Clone)]
pub struct Density {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub bw: f64,
    pub n: usize,
}

/// `BinDist(x, w, lo, hi, n)`: linear binning onto `n` grid points, returned in a vector of `2n`
/// whose upper half is zero.
pub fn bin_dist(x: &[f64], w: &[f64], xlo: f64, xhi: f64, n: usize) -> Vec<f64> {
    let mut y = vec![0.0; 2 * n];
    let ixmin: i64 = 0;
    let ixmax: i64 = n as i64 - 2;
    let xdelta = (xhi - xlo) / (n as i32 - 1) as f64;
    for (&xi, &wi) in x.iter().zip(w) {
        if xi.is_finite() {
            let xpos = (xi - xlo) / xdelta;
            if xpos > i32::MAX as f64 || xpos < i32::MIN as f64 {
                continue;
            }
            let ix = xpos.floor() as i32 as i64;
            let fx = xpos - ix as f64;
            if ixmin <= ix && ix <= ixmax {
                y[ix as usize] += (1.0 - fx) * wi;
                y[ix as usize + 1] += fx * wi;
            } else if ix == -1 {
                y[0] += fx * wi;
            } else if ix == ixmax + 1 {
                y[ix as usize] += (1.0 - fx) * wi;
            }
        }
    }
    y
}

/// `density(x, bw = "nrd", adjust, n = n_user)` with the Gaussian kernel, `cut = 3`, `ext = 4` and
/// no weights. Every element of `x` must be finite (sctransform passes log geometric means).
pub fn density_nrd(x: &[f64], adjust: f64, n_user: usize) -> Result<Density, DensityError> {
    if x.iter().any(|v| v.is_nan()) {
        return Err(DensityError::MissingValues);
    }
    let big_n = x.len();
    let xs: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    let nx = xs.len();
    // weights <- rep.int(1/nx, nx); totMass <- nx/N
    let w = vec![1.0 / nx as f64; nx];
    let tot_mass = nx as f64 / big_n as f64;
    let mut n = n_user.max(512);
    if n > 512 {
        n = r_pow(2.0, (n as f64).log2().ceil()) as usize;
    }
    if nx < 2 {
        return Err(DensityError::TooFewPoints);
    }
    let bw = adjust * bw_nrd(&xs)?;
    if !bw.is_finite() || bw <= 0.0 {
        return Err(DensityError::BadBandwidth(bw));
    }
    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for &v in &xs {
        if v < xmin {
            xmin = v;
        }
        if v > xmax {
            xmax = v;
        }
    }
    let cut = 3.0;
    let ext = 4.0;
    let from = xmin - cut * bw;
    let to = xmax + cut * bw;
    let lo = from - ext * bw;
    let up = to + ext * bw;
    let y: Vec<f64> = bin_dist(&xs, &w, lo, up, n).into_iter().map(|v| v * tot_mass).collect();
    // kords <- seq.int(0, (2L*n-1)/(n-1) * (up-lo), length.out = 2L * n)
    let nf = n as f64;
    let span = (2.0 * nf - 1.0) / (nf - 1.0) * (up - lo);
    let mut kords = seq_length_out(0.0, span, 2.0 * nf).expect("finite grid");
    // kords[(n + 2):(2 * n)] <- -kords[n:2]
    for t in 0..n - 1 {
        kords[n + 1 + t] = -kords[n - 1 - t];
    }
    let kern: Vec<f64> = kords.iter().map(|&k| dnorm(k, 0.0, bw, false)).collect();
    let to_c = |v: &[f64]| v.iter().map(|&re| Complex::new(re, 0.0)).collect::<Vec<_>>();
    let fy = fft(&to_c(&y), false).expect("fft");
    let fk = fft(&to_c(&kern), false).expect("fft");
    let prod: Vec<Complex> = fy.iter().zip(&fk).map(|(&a, &b)| a * Complex::new(b.re, -b.im)).collect();
    let conv = fft(&prod, true).expect("fft");
    let len_y = y.len() as f64;
    let dens: Vec<f64> = conv[..n]
        .iter()
        .map(|c| {
            let v = c.re / len_y;
            // pmax.int(0, v)
            if v.is_nan() {
                v
            } else if v > 0.0 {
                v
            } else {
                0.0
            }
        })
        .collect();
    let xords = seq_length_out(lo, up, nf).expect("finite grid");
    let xout = seq_length_out(from, to, n_user as f64).expect("finite grid");
    let f = ApproxFun::new(&xords, &dens, Method::Linear).expect("approx");
    let yout = f.eval_all(&xout);
    Ok(Density { x: xout, y: yout, bw, n: big_n })
}
