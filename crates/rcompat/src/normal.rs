//! Normal quantiles and deviates: `qnorm5` (R 4.5.3's `src/nmath/qnorm.c`: Wichura's AS 241,
//! with Maechler's asymptotic formulas for the extreme tail) and `rnorm` by inversion
//! (`src/nmath/snorm.c`, `src/nmath/rnorm.c`).
//!
//! Expressions keep R's evaluation order and constants; Rust does not contract them into fused
//! multiply-adds, and `log`, `log1p` and `expm1` are glibc's.

// R's coefficients are kept exactly as written in qnorm.c; they parse to the same doubles.
#![allow(clippy::excessive_precision)]

use crate::libm;
use crate::rng::RRng;

/// 2^27: `norm_rand` combines two uniforms into one with more precision.
const BIG: f64 = 134217728.0;

/// `qnorm(p, mu, sigma, lower.tail, log.p)`.
pub fn qnorm(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    if p.is_nan() || mu.is_nan() || sigma.is_nan() {
        return p + mu + sigma;
    }
    // R_Q_P01_boundaries(p, -Inf, Inf)
    if log_p {
        if p > 0.0 {
            return f64::NAN;
        }
        if p == 0.0 {
            return if lower_tail { f64::INFINITY } else { f64::NEG_INFINITY };
        }
        if p == f64::NEG_INFINITY {
            return if lower_tail { f64::NEG_INFINITY } else { f64::INFINITY };
        }
    } else {
        if !(0.0..=1.0).contains(&p) {
            return f64::NAN;
        }
        if p == 0.0 {
            return if lower_tail { f64::NEG_INFINITY } else { f64::INFINITY };
        }
        if p == 1.0 {
            return if lower_tail { f64::INFINITY } else { f64::NEG_INFINITY };
        }
    }
    if sigma < 0.0 {
        return f64::NAN;
    }
    if sigma == 0.0 {
        return mu;
    }

    // R_DT_qIv(p): the lower-tail probability
    let p_ = if log_p {
        if lower_tail {
            libm::exp(p)
        } else {
            -libm::expm1(p)
        }
    } else if lower_tail {
        p
    } else {
        0.5 - p + 0.5
    };
    let q = p_ - 0.5;
    let mut val;

    if q.abs() <= 0.425 {
        let r = 0.180625 - q * q;
        val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r
            + 67265.770927008700853) * r + 45921.953931549871457) * r
            + 13731.693765509461125) * r + 1971.5909503065514427) * r
            + 133.14166789178437745) * r + 3.387132872796366608)
            / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r
                + 39307.89580009271061) * r + 21213.794301586595867) * r
                + 5394.1960214247511077) * r + 687.1870074920579083) * r
                + 42.313330701600911252) * r + 1.);
    } else {
        // r = sqrt(-log(min(p~, 1 - p~)))
        let lp = if log_p && ((lower_tail && q <= 0.0) || (!lower_tail && q > 0.0)) {
            p
        } else {
            libm::log(if q > 0.0 {
                // R_DT_CIv(p): 1 - p~
                if log_p {
                    if lower_tail {
                        -libm::expm1(p)
                    } else {
                        libm::exp(p)
                    }
                } else if lower_tail {
                    0.5 - p + 0.5
                } else {
                    p
                }
            } else {
                p_
            })
        };
        let mut r = (-lp).sqrt();
        if r <= 5. {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r
                + 0.24178072517745061177) * r + 1.27045825245236838258) * r
                + 3.64784832476320460504) * r + 5.7694972214606914055) * r
                + 4.6303378461565452959) * r + 1.42343711074968357734)
                / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r
                    + 0.0151986665636164571966) * r + 0.14810397642748007459) * r
                    + 0.68976733498510000455) * r + 1.6763848301838038494) * r
                    + 2.05319162663775882187) * r + 1.);
        } else if r <= 27. {
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r
                + 0.0012426609473880784386) * r + 0.026532189526576123093) * r
                + 0.29656057182850489123) * r + 1.7848265399172913358) * r
                + 5.4637849111641143699) * r + 6.6579046435011037772)
                / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r
                    + 1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r
                    + 0.0148753612908506148525) * r + 0.13692988092273580531) * r
                    + 0.59983220655588793769) * r + 1.);
        } else if r >= 6.4e8 {
            // p is extremely close to 0 or 1: the "0-th order" asymptotic formula
            val = r * std::f64::consts::SQRT_2;
        } else {
            // r in (27, 6.4e8): Maechler (2022), asymptotic tail formulas
            // M_2PI (6.283185307179586476925286766559) is the same double as TAU
            const M_2PI: f64 = std::f64::consts::TAU;
            let s2 = -(lp * 2.0); // -ldexp(lp, 1)
            let mut x2 = s2 - libm::log(M_2PI * s2);
            if r < 36000. {
                x2 = s2 - libm::log(M_2PI * x2) - 2. / (2. + x2);
                if r < 840. {
                    x2 = s2 - libm::log(M_2PI * x2)
                        + 2. * libm::log1p(-(1. - 1. / (4. + x2)) / (2. + x2));
                    if r < 109. {
                        x2 = s2 - libm::log(M_2PI * x2)
                            + 2. * libm::log1p(-(1. - (1. - 5. / (6. + x2)) / (4. + x2)) / (2. + x2));
                        if r < 55. {
                            x2 = s2 - libm::log(M_2PI * x2)
                                + 2. * libm::log1p(
                                    -(1. - (1. - (5. - 9. / (8. + x2)) / (6. + x2)) / (4. + x2))
                                        / (2. + x2),
                                );
                        }
                    }
                }
            }
            val = x2.sqrt();
        }
        if q < 0.0 {
            val = -val;
        }
    }
    mu + sigma * val
}

/// `norm_rand()` with `normal.kind = "Inversion"`: two uniforms per deviate,
/// `u = (int)(2^27 * u1) + u2`, then `qnorm(u / 2^27)`.
pub fn norm_rand(rng: &mut RRng) -> f64 {
    let u1 = rng.unif_rand();
    let u = (BIG * u1) as i32 as f64 + rng.unif_rand();
    qnorm(u / BIG, 0.0, 1.0, true, false)
}

/// `rnorm(1, mu, sigma)`: NaN if mu is NaN or sigma is not finite or negative; `mu` without a draw
/// if sigma is 0 or mu infinite; else `mu + sigma * norm_rand()`.
pub fn rnorm_one(rng: &mut RRng, mu: f64, sigma: f64) -> f64 {
    if mu.is_nan() || !sigma.is_finite() || sigma < 0.0 {
        return f64::NAN;
    }
    if sigma == 0.0 || !mu.is_finite() {
        mu
    } else {
        mu + sigma * norm_rand(rng)
    }
}

/// `rnorm(n, mean, sd)` with scalar mean and sd.
pub fn rnorm(rng: &mut RRng, n: usize, mean: f64, sd: f64) -> Vec<f64> {
    (0..n).map(|_| rnorm_one(rng, mean, sd)).collect()
}

/// `M_1_SQRT_2PI` from Rmath.h: 1/sqrt(2 pi).
const M_1_SQRT_2PI: f64 = 0.398942280401432677939946059934;
/// `M_LN_SQRT_2PI` from Rmath.h: log(sqrt(2 pi)).
const M_LN_SQRT_2PI: f64 = 0.918938533204672741780329736406;

/// `dnorm(x, mu, sigma, log)` (`dnorm4` in `src/nmath/dnorm.c`): `M_1_SQRT_2PI * exp(-x*x/2) /
/// sigma` for |x| < 5 standard deviations; beyond, x is split into a part with 16 fractional bits
/// and a remainder, so that `x*x` loses no accuracy (PR#15620).
///
/// As R's `dnorm()` calls it (`math3_1` in `src/main/arithmetic.c`): `NA` in any argument gives
/// `NA`, any other NaN gives `NaN`, before `dnorm4` runs.
pub fn dnorm(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    use crate::summary::{is_na, NA_REAL};
    if is_na(x) || is_na(mu) || is_na(sigma) {
        return NA_REAL;
    }
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return f64::from_bits(0x7ff8_0000_0000_0000); // R_NaN
    }
    let d0 = if give_log { f64::NEG_INFINITY } else { 0.0 };
    if sigma < 0.0 {
        return f64::NAN;
    }
    if !sigma.is_finite() {
        return d0;
    }
    if !x.is_finite() && mu == x {
        return f64::NAN; // x - mu is NaN
    }
    if sigma == 0.0 {
        return if x == mu { f64::INFINITY } else { d0 };
    }
    let x = (x - mu) / sigma;
    if !x.is_finite() {
        return d0;
    }
    let x = x.abs();
    if x >= 2.0 * f64::MAX.sqrt() {
        return d0;
    }
    if give_log {
        return -(M_LN_SQRT_2PI + 0.5 * x * x + libm::log(sigma));
    }
    if x < 5.0 {
        return M_1_SQRT_2PI * libm::exp(-0.5 * x * x) / sigma;
    }
    // underflow boundary: sqrt(-2 log(2) (DBL_MIN_EXP + 1 - DBL_MANT_DIG)), about 38.59
    if x > (-2.0 * std::f64::consts::LN_2 * (-1021.0 + 1.0 - 53.0)).sqrt() {
        return 0.0;
    }
    // x1 = ldexp(nearbyint(ldexp(x, 16)), -16): exact scalings by 2^16
    let x1 = (x * 65536.0).round_ties_even() / 65536.0;
    let x2 = x - x1;
    M_1_SQRT_2PI / sigma * (libm::exp(-0.5 * x1 * x1) * libm::exp((-0.5 * x2 - x1) * x2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_quantiles() {
        assert_eq!(qnorm(0.5, 0.0, 1.0, true, false), 0.0);
        assert!((qnorm(0.975, 0.0, 1.0, true, false) - 1.959963984540054).abs() < 1e-15);
        assert_eq!(qnorm(0.0, 0.0, 1.0, true, false), f64::NEG_INFINITY);
        assert!(qnorm(1.5, 0.0, 1.0, true, false).is_nan());
        assert_eq!(qnorm(0.3, 2.0, 0.0, true, false), 2.0);
    }

    #[test]
    fn two_uniforms_per_normal() {
        let mut a = RRng::new(42);
        let mut b = a.clone();
        rnorm(&mut a, 5, 0.0, 1.0);
        for _ in 0..10 {
            b.unif_rand();
        }
        assert_eq!(a.unif_rand(), b.unif_rand());
    }
}

