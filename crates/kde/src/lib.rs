//! Kernel density estimation and moments, as DoubletFinder uses them to choose pK
//! (docs/PLAN.md B3; milestone M3), ported so that every result is bit-identical to R's:
//!
//! - [`bkde`]: `KernSmooth::bkde(x, kernel = "normal")` with its defaults (KernSmooth 2.23-27):
//!   the normal-reference bandwidth, a 401-point grid from `seq(length.out = )`, Fortran linear
//!   binning ([`linbin`]), kernel weights truncated at 4 bandwidths, and the convolution through
//!   R's `fft` ([`rcompat::fft`]);
//! - [`density_function`]: R's `approxfun()` of that estimate (linear, `rule = 1`, `ties = mean`);
//! - [`skewness`], [`kurtosis`] and [`bimodality_coefficient`]: DoubletFinder 2.0.6's own
//!   functions (taken from the `modes` package), with R's long-double `sum` and `mean` and R's
//!   `^` (`x^2` is `x*x`; other powers are glibc's `pow`).

use rcompat::approx::{ApproxError, ApproxFun, Method};
use rcompat::arith::{r_pow, seq_length_out};
use rcompat::fft::{fft, fft_real, Complex};
use rcompat::{dnorm, libm, mean, sum, var};
use std::fmt;

/// `bkde`'s `gridsize` default.
pub const GRIDSIZE: usize = 401;

/// A binned kernel density estimate: the grid and the density on it.
#[derive(Debug, Clone, PartialEq)]
pub struct Bkde {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Why `bkde` stops in R.
#[derive(Debug, Clone, PartialEq)]
pub enum BkdeError {
    /// Fewer than two values: `var(x)` is NA, so the grid's ends are not finite
    /// ("'from' must be a finite number" in `seq`).
    NotFinite,
    /// Every value is the same: the bandwidth is 0, the grid spacing over the bandwidth is
    /// 0/0, and `if (L == 0)` stops with "missing value where TRUE/FALSE needed". DoubletFinder
    /// stops there when every pANN of a (pN, pK) pair is identical.
    ZeroBandwidth,
    /// `fft` could not factor the length (impossible for these sizes).
    Fft,
}

impl fmt::Display for BkdeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            BkdeError::NotFinite => "bkde: 'from' must be a finite number (fewer than two values)",
            BkdeError::ZeroBandwidth => {
                "bkde: missing value where TRUE/FALSE needed (all values identical: bandwidth 0)"
            }
            BkdeError::Fft => "bkde: fft factorization error",
        })
    }
}

impl std::error::Error for BkdeError {}

/// Fortran's `INT()` of a double as gfortran compiles it on x86-64 (`cvttsd2si`): truncation
/// toward zero, and the "integer indefinite" value `i32::MIN` for NaN and out-of-range values.
fn fortran_int(v: f64) -> i32 {
    if v.is_nan() || v >= 2147483648.0 || v <= -2147483649.0 {
        i32::MIN
    } else {
        v as i32
    }
}

/// KernSmooth's Fortran `linbin(X, n, a, b, M, trun, gcnts)`: linear binning of `x` onto `m`
/// equally spaced points from `a` to `b`. With `truncate`, points outside the grid are dropped;
/// otherwise their weight goes to the end points.
pub fn linbin(x: &[f64], a: f64, b: f64, m: usize, truncate: bool) -> Vec<f64> {
    let mut gcnts = vec![0.0f64; m];
    let mi = m as i32;
    let delta = (b - a) / (m as i32 - 1) as f64;
    for &xi in x {
        let lxi = ((xi - a) / delta) + 1.0;
        // find integer part of lxi
        let li = fortran_int(lxi);
        let rem = lxi - li as f64;
        if li >= 1 && li < mi {
            gcnts[(li - 1) as usize] += 1.0 - rem;
            gcnts[li as usize] += rem;
        }
        if li < 1 && !truncate {
            gcnts[0] += 1.0;
        }
        if li >= mi && !truncate {
            gcnts[m - 1] += 1.0;
        }
    }
    gcnts
}

/// `min(x)` and `max(x)` of a numeric vector without NA.
fn range(x: &[f64]) -> (f64, f64) {
    x.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| (lo.min(v), hi.max(v)))
}

/// `KernSmooth::bkde(x, kernel = "normal")` with every other argument at its default
/// (`canonical = FALSE`, bandwidth `del0 * (243/(35 n))^(1/5) * sqrt(var(x))`, `gridsize = 401`,
/// range `[min(x) - 4h, max(x) + 4h]`, `truncate = TRUE`).
pub fn bkde(x: &[f64]) -> Result<Bkde, BkdeError> {
    let n = x.len();
    let m = GRIDSIZE;
    let nf = n as f64;
    // canonical scaling factor of the normal kernel: (1/(4*pi))^(1/10)
    let del0 = r_pow(1.0 / (4.0 * std::f64::consts::PI), 1.0 / 10.0);
    // the default bandwidth: del0 * (243/(35*n))^(1/5) * sqrt(var(x))
    let h = del0 * r_pow(243.0 / (35.0 * nf), 1.0 / 5.0) * var(x, false).sqrt();
    let tau = 4.0;
    let (lo, hi) = range(x);
    let a = lo - tau * h;
    let b = hi + tau * h;
    // set up grid points and bin the data
    let gpoints = seq_length_out(a, b, m as f64).map_err(|_| BkdeError::NotFinite)?;
    let gcounts = linbin(x, gpoints[0], gpoints[m - 1], m, true);
    // kernel weights
    let delta = (b - a) / (h * (m as i32 - 1) as f64);
    let l = (tau / delta).floor();
    if l.is_nan() {
        return Err(BkdeError::ZeroBandwidth); // `if (L == 0)` on NA
    }
    let l = l.min(m as f64);
    // (L == 0 only warns: "Binning grid too coarse for current (small) bandwidth")
    let li = l as usize;
    let kappa: Vec<f64> = (0..=li).map(|i| dnorm(i as f64 * delta, 0.0, 1.0, false) / (nf * h)).collect();
    // combine weights and counts: P = 2^ceiling(log(M + L + 1)/log(2)) >= 2L + 1
    let p = r_pow(2.0, (libm::log(m as f64 + l + 1.0) / libm::log(2.0)).ceil());
    let pu = p as usize;
    let mut kap = Vec::with_capacity(pu);
    kap.extend_from_slice(&kappa);
    kap.resize(kap.len() + (pu - 2 * li - 1), 0.0);
    kap.extend(kappa[1..].iter().rev());
    // should have total weight one
    let tot = sum(&kap, false) * (b - a) / (m as i32 - 1) as f64 * nf;
    let mut gc = gcounts;
    gc.resize(pu, 0.0);
    let kap_scaled: Vec<f64> = kap.iter().map(|&k| k / tot).collect();
    let kf = fft_real(&kap_scaled, false).map_err(|_| BkdeError::Fft)?;
    let gf = fft_real(&gc, false).map_err(|_| BkdeError::Fft)?;
    let prod: Vec<Complex> = kf.iter().zip(&gf).map(|(&k, &g)| k * g).collect();
    let conv = fft(&prod, true).map_err(|_| BkdeError::Fft)?;
    let y = conv[..m].iter().map(|c| c.re / p).collect();
    Ok(Bkde { x: gpoints, y })
}

/// `approxfun(bkde(...))`: linear interpolation of the estimate, NA outside the grid.
pub fn density_function(est: &Bkde) -> Result<ApproxFun, ApproxError> {
    ApproxFun::new(&est.x, &est.y, Method::Linear)
}

/// `sum((x - m)^p)` with R's `^` and long-double `sum`.
fn central_power_sum(x: &[f64], m: f64, p: f64) -> f64 {
    let v: Vec<f64> = x.iter().map(|&xi| r_pow(xi - m, p)).collect();
    sum(&v, false)
}

/// DoubletFinder's `skewness(x)`:
/// `S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)`, then
/// `S*(sqrt(n*(n-1)))/(n-2)`.
pub fn skewness(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let m = mean(x, false);
    let s = (1.0 / n) * central_power_sum(x, m, 3.0) / r_pow((1.0 / n) * central_power_sum(x, m, 2.0), 1.5);
    s * (n * (n - 1.0)).sqrt() / (n - 2.0)
}

/// DoubletFinder's `kurtosis(x)`:
/// `K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3`, then
/// `((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3`.
pub fn kurtosis(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let m = mean(x, false);
    let k = (1.0 / n) * central_power_sum(x, m, 4.0) / r_pow((1.0 / n) * central_power_sum(x, m, 2.0), 2.0) - 3.0;
    ((n - 1.0) * ((n + 1.0) * k - 3.0 * (n - 1.0)) / ((n - 2.0) * (n - 3.0))) + 3.0
}

/// DoubletFinder's `bimodality_coefficient(x)`:
/// `((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))` with G the skewness and K the kurtosis.
pub fn bimodality_coefficient(x: &[f64]) -> f64 {
    let g = skewness(x);
    let k = kurtosis(x);
    let n = x.len() as f64;
    (r_pow(g, 2.0) + 1.0) / (k + ((3.0 * r_pow(n - 1.0, 2.0)) / ((n - 2.0) * (n - 3.0))))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_values_stop_like_r() {
        assert_eq!(bkde(&[0.25; 50]), Err(BkdeError::ZeroBandwidth));
        assert_eq!(bkde(&[0.25]), Err(BkdeError::NotFinite));
    }

    #[test]
    fn density_integrates_to_about_one() {
        let x: Vec<f64> = (0..1000).map(|i| ((i * 7919) % 1000) as f64 / 1000.0).collect();
        let est = bkde(&x).unwrap();
        assert_eq!(est.x.len(), GRIDSIZE);
        let dx = est.x[1] - est.x[0];
        let total: f64 = est.y.iter().sum::<f64>() * dx;
        assert!((total - 1.0).abs() < 1e-3, "{total}");
    }

    #[test]
    fn fortran_truncation() {
        assert_eq!(fortran_int(2.9), 2);
        assert_eq!(fortran_int(-2.9), -2);
        assert_eq!(fortran_int(f64::NAN), i32::MIN);
        assert_eq!(fortran_int(3e9), i32::MIN);
    }

    #[test]
    fn moments_of_a_known_sample() {
        let x = [1.0, 2.0, 2.0, 3.0, 9.0];
        // the definitions evaluated by hand: m2 = 8.24, m3 = 31.248, m4 = 204.8672
        let g = skewness(&x);
        let k = kurtosis(&x);
        assert!((g - 1.969363).abs() < 1e-5, "{g}");
        assert!((k - -4.930813).abs() < 1e-5, "{k}");
        assert!(bimodality_coefficient(&x).is_finite());
    }
}
