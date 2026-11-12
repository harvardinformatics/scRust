//! R 4.5.3's negative binomial log density in the mean parametrisation, `dnbinom(x, size, mu =
//! mu, log = TRUE)` (`src/nmath/dnbinom.c`), and what it calls: `dbinom_raw` (log form),
//! `stirlerr`, `bd0`, and `lgamma1p` with `logcf` and `log1pmx` (`pgamma.c`). glmGamPoi's
//! `optimize()` fallback for a Newton fit that fails evaluates it.
//!
//! Inside R's math library `lgamma` is the C library's (glibc's) and `lgammafn` R's own; `log1p`
//! and `ldexp` are glibc's. Only the log densities are ported (`give_log = TRUE`), and not the
//! Poisson limit `size = Inf`, which needs `dpois_raw` and `ebd0`.

use crate::gamma::{lgamma_glibc, lgammafn};
use crate::libm::{ldexp, log, log1p};

const M_LN_2PI: f64 = 1.837877066409345483560659472811;
const M_LN_SQRT_2PI: f64 = 0.918938533204672741780329736406;

// Generated from R 4.5.3 src/nmath/stirlerr.c and pgamma.c (the decimal literals as written there).
const S0: f64 = 0.083333333333333333333;
const S1: f64 = 0.00277777777777777777778;
const S2: f64 = 0.00079365079365079365079365;
const S3: f64 = 0.000595238095238095238095238;
const S4: f64 = 0.0008417508417508417508417508;
const S5: f64 = 0.0019175269175269175269175262;
const S6: f64 = 0.0064102564102564102564102561;
const S7: f64 = 0.029550653594771241830065352;
const S8: f64 = 0.17964437236883057316493850;
const S9: f64 = 1.3924322169059011164274315;
const S10: f64 = 13.402864044168391994478957;
const S11: f64 = 156.84828462600201730636509;
const S12: f64 = 2193.1033333333333333333333;
const S13: f64 = 36108.771253724989357173269;
const S14: f64 = 691472.26885131306710839498;
const S15: f64 = 15238221.539407416192283370;
const S16: f64 = 382900751.39141414141414141;
const SFERR_HALVES: [f64; 31] = [
    0.0,
    0.1534264097200273452913848,
    0.0810614667953272582196702,
    0.0548141210519176538961390,
    0.0413406959554092940938221,
    0.03316287351993628748511048,
    0.02767792568499833914878929,
    0.02374616365629749597132920,
    0.02079067210376509311152277,
    0.01848845053267318523077934,
    0.01664469118982119216319487,
    0.01513497322191737887351255,
    0.01387612882307074799874573,
    0.01281046524292022692424986,
    0.01189670994589177009505572,
    0.01110455975820691732662991,
    0.010411265261972096497478567,
    0.009799416126158803298389475,
    0.009255462182712732917728637,
    0.008768700134139385462952823,
    0.008330563433362871256469318,
    0.007934114564314020547248100,
    0.007573675487951840794972024,
    0.007244554301320383179543912,
    0.006942840107209529865664152,
    0.006665247032707682442354394,
    0.006408994188004207068439631,
    0.006171712263039457647532867,
    0.005951370112758847735624416,
    0.005746216513010115682023589,
    0.005554733551962801371038690,
];
const LGAMMA1P_COEFFS: [f64; 40] = [
    0.3224670334241132182362075833230126e-0,
    0.6735230105319809513324605383715000e-1,
    0.2058080842778454787900092413529198e-1,
    0.7385551028673985266273097291406834e-2,
    0.2890510330741523285752988298486755e-2,
    0.1192753911703260977113935692828109e-2,
    0.5096695247430424223356548135815582e-3,
    0.2231547584535793797614188036013401e-3,
    0.9945751278180853371459589003190170e-4,
    0.4492623673813314170020750240635786e-4,
    0.2050721277567069155316650397830591e-4,
    0.9439488275268395903987425104415055e-5,
    0.4374866789907487804181793223952411e-5,
    0.2039215753801366236781900709670839e-5,
    0.9551412130407419832857179772951265e-6,
    0.4492469198764566043294290331193655e-6,
    0.2120718480555466586923135901077628e-6,
    0.1004322482396809960872083050053344e-6,
    0.4769810169363980565760193417246730e-7,
    0.2271109460894316491031998116062124e-7,
    0.1083865921489695409107491757968159e-7,
    0.5183475041970046655121248647057669e-8,
    0.2483674543802478317185008663991718e-8,
    0.1192140140586091207442548202774640e-8,
    0.5731367241678862013330194857961011e-9,
    0.2759522885124233145178149692816341e-9,
    0.1330476437424448948149715720858008e-9,
    0.6422964563838100022082448087644648e-10,
    0.3104424774732227276239215783404066e-10,
    0.1502138408075414217093301048780668e-10,
    0.7275974480239079662504549924814047e-11,
    0.3527742476575915083615072228655483e-11,
    0.1711991790559617908601084114443031e-11,
    0.8315385841420284819798357793954418e-12,
    0.4042200525289440065536008957032895e-12,
    0.1966475631096616490411045679010286e-12,
    0.9573630387838555763782200936508615e-13,
    0.4664076026428374224576492565974577e-13,
    0.2273736960065972320633279596737272e-13,
    0.1109139947083452201658320007192334e-13,
];

/// `stirlerr(n) = log(n!) - log(sqrt(2*pi*n) * (n/e)^n)`, for `n >= 0`.
pub fn stirlerr(n: f64) -> f64 {
    if n <= 23.5 {
        let nn = n + n;
        if n <= 15.0 && nn == (nn as i32) as f64 {
            return SFERR_HALVES[nn as i32 as usize];
        }
        if n <= 5.25 {
            if n >= 1.0 {
                // "MM2"
                let l_n = log(n);
                return lgamma_glibc(n) + n * (1.0 - l_n) + ldexp(l_n - M_LN_2PI, -1);
            }
            return lgamma1p(n) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI;
        }
        // 5.25 < n <= 23.5
        let nn = n * n;
        if n > 12.8 {
            return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - S6 / nn) / nn) / nn) / nn) / nn) / nn) / n;
        }
        if n > 12.3 {
            return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - S7 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
        }
        if n > 8.9 {
            return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - S8 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
        }
        if n > 7.3 {
            return (S0
                - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - S10 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn)
                    / nn)
                / n;
        }
        if n > 6.6 {
            return (S0
                - (S1
                    - (S2
                        - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - (S10 - (S11 - S12 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn)
                            / nn)
                        / nn)
                    / nn)
                / n;
        }
        if n > 6.1 {
            return (S0
                - (S1
                    - (S2
                        - (S3
                            - (S4
                                - (S5 - (S6 - (S7 - (S8 - (S9 - (S10 - (S11 - (S12 - (S13 - S14 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / n;
        }
        // 6.1 >= n > 5.25
        return (S0
            - (S1
                - (S2
                    - (S3
                        - (S4
                            - (S5
                                - (S6
                                    - (S7 - (S8 - (S9 - (S10 - (S11 - (S12 - (S13 - (S14 - (S15 - S16 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn)
                                        / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / nn)
            / n;
    }
    // n > 23.5
    let nn = n * n;
    if n > 15.7e6 {
        return S0 / n;
    }
    if n > 6180.0 {
        return (S0 - S1 / nn) / n;
    }
    if n > 205.0 {
        return (S0 - (S1 - S2 / nn) / nn) / n;
    }
    if n > 86.0 {
        return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
    }
    if n > 27.0 {
        return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
    }
    // 23.5 < n <= 27
    (S0 - (S1 - (S2 - (S3 - (S4 - S5 / nn) / nn) / nn) / nn) / nn) / n
}

/// `bd0(x, np) = x log(x/np) + np - x`, the deviance term, computed stably.
pub fn bd0(x: f64, np: f64) -> f64 {
    if !x.is_finite() || !np.is_finite() || np == 0.0 {
        return f64::NAN; // ML_WARN_return_NAN
    }
    if (x - np).abs() < 0.1 * (x + np) {
        let d = x - np;
        let mut v = d / (x + np);
        if d != 0.0 && v == 0.0 {
            // v has underflown to 0 (x + np = Inf)
            let x_ = ldexp(x, -2);
            let n_ = ldexp(np, -2);
            v = (x_ - n_) / (x_ + n_);
        }
        let mut s = ldexp(d, -1) * v;
        if ldexp(s, 1).abs() < f64::MIN_POSITIVE {
            return ldexp(s, 1);
        }
        let mut ej = x * v;
        v *= v;
        for j in 1..1000 {
            // Taylor series
            ej *= v;
            let s_ = s;
            s += ej / ((j << 1) + 1) as f64;
            if s == s_ {
                return ldexp(s, 1);
            }
        }
        // (R warns that the series did not converge and goes on below)
    }
    let lg_x_n = if (x / np).is_finite() { log(x / np) } else { log(x) - log(np) };
    if x > np {
        x * (lg_x_n - 1.0) + np
    } else {
        x * lg_x_n + np - x
    }
}

const SCALEFACTOR: f64 = 4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0;

/// The continued fraction `sum_{k>=0} x^k / (i + k d)` (`logcf`).
fn logcf(x: f64, i: f64, d: f64, eps: f64) -> f64 {
    let mut c1 = 2.0 * d;
    let mut c2 = i + d;
    let mut c4 = c2 + d;
    let mut a1 = c2;
    let mut b1 = i * (c2 - i * x);
    let mut b2 = d * d * x;
    let mut a2 = c4 * c2 - b2;
    b2 = c4 * b1 - i * b2;
    while (a2 * b1 - a1 * b2).abs() > (eps * b1 * b2).abs() {
        let mut c3 = c2 * c2 * x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;
        c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;
        if b2.abs() > SCALEFACTOR {
            a1 /= SCALEFACTOR;
            b1 /= SCALEFACTOR;
            a2 /= SCALEFACTOR;
            b2 /= SCALEFACTOR;
        } else if b2.abs() < 1.0 / SCALEFACTOR {
            a1 *= SCALEFACTOR;
            b1 *= SCALEFACTOR;
            a2 *= SCALEFACTOR;
            b2 *= SCALEFACTOR;
        }
    }
    a2 / b2
}

/// `log1pmx(x) = log(1 + x) - x`, accurate for small `x`.
pub fn log1pmx(x: f64) -> f64 {
    const MIN_LOG1_VALUE: f64 = -0.79149064;
    if x > 1.0 || x < MIN_LOG1_VALUE {
        return log1p(x) - x;
    }
    // expand in [x/(2+x)]^2 =: y
    let r = x / (2.0 + x);
    let y = r * r;
    if x.abs() < 1e-2 {
        const TWO: f64 = 2.0;
        r * ((((TWO / 9.0 * y + TWO / 7.0) * y + TWO / 5.0) * y + TWO / 3.0) * y - x)
    } else {
        r * (2.0 * y * logcf(y, 3.0, 2.0, 1e-14) - x)
    }
}

/// `lgamma1p(a) = log(gamma(a + 1))`, accurate also for small `a`.
pub fn lgamma1p(a: f64) -> f64 {
    if a.abs() >= 0.5 {
        return lgammafn(a + 1.0);
    }
    const EULERS_CONST: f64 = 0.5772156649015328606065120900824024;
    const N: usize = 40;
    const C: f64 = 0.2273736845824652515226821577978691e-12; // zeta(N+2)-1
    const TOL_LOGCF: f64 = 1e-14;
    let mut lgam = C * logcf(-a / 2.0, (N + 2) as f64, 1.0, TOL_LOGCF);
    for i in (0..N).rev() {
        lgam = LGAMMA1P_COEFFS[i] - a * lgam;
    }
    (a * lgam - EULERS_CONST) * a - log1pmx(a)
}

/// `dbinom_raw(x, n, p, q, log = TRUE)`.
pub fn dbinom_raw_log(x: f64, n: f64, p: f64, q: f64) -> f64 {
    if p == 0.0 {
        return if x == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if q == 0.0 {
        return if x == n { 0.0 } else { f64::NEG_INFINITY };
    }
    if x == 0.0 {
        if n == 0.0 {
            return 0.0;
        }
        return if p > q { n * log(q) } else { n * log1p(-p) };
    }
    if x == n {
        return if p > q { n * log1p(-q) } else { n * log(p) };
    }
    if x < 0.0 || x > n {
        return f64::NEG_INFINITY;
    }
    let mut n = n;
    if !n.is_finite() {
        if x.is_finite() {
            return f64::NEG_INFINITY;
        }
        n = f64::MAX;
    }
    let lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
    let lf = M_LN_2PI + log(x) + log1p(-x / n);
    lc - 0.5 * lf
}

/// R's `R_nonint(x)`.
fn r_nonint(x: f64) -> bool {
    let ax = x.abs();
    // fmax2(1., fabs(x))
    let m = if ax.is_nan() { 1.0 + ax } else if 1.0 < ax { ax } else { 1.0 };
    (x - x.round_ties_even()).abs() > 1e-9 * m
}

/// `dnbinom(x, size, mu = mu, log = TRUE)` for one value. Panics for `size = Inf` (the Poisson
/// limit, not ported).
pub fn dnbinom_mu_log(x: f64, size: f64, mu: f64) -> f64 {
    if x.is_nan() || size.is_nan() || mu.is_nan() {
        return x + size + mu;
    }
    if mu < 0.0 || size < 0.0 {
        return f64::NAN; // with a warning
    }
    if r_nonint(x) {
        return f64::NEG_INFINITY; // with a warning: non-integer x
    }
    if x < 0.0 || !x.is_finite() {
        return f64::NEG_INFINITY;
    }
    if x == 0.0 && size == 0.0 {
        return 0.0;
    }
    let x = x.round_ties_even(); // R_forceint: nearbyint
    assert!(size.is_finite(), "dnbinom_mu_log: size = Inf (the Poisson limit) is not ported");
    if x == 0.0 {
        return size * if size < mu { log(size / (size + mu)) } else { log1p(-mu / (size + mu)) };
    }
    if x < 1e-10 * size {
        // MM's formula instead of dbinom_raw()
        let p = if size < mu { log(size / (1.0 + size / mu)) } else { log(mu / (1.0 + mu / size)) };
        let xx2s = if x < f64::MAX.sqrt() { ldexp(x * (x - 1.0), -1) / size } else { x * (ldexp(x, -1) / size) };
        return x * p - mu - lgamma1p(x) + log1p(xx2s);
    }
    let p = if x < size { log1p(-x / (size + x)) } else { log(size / (size + x)) };
    let ans = dbinom_raw_log(size, x + size, size / (size + mu), mu / (size + mu));
    p + ans
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_values() {
        // log P(X = 0) for NB(size 1, mu 1) is log(1/2)
        assert_eq!(dnbinom_mu_log(0.0, 1.0, 1.0), log1p(-0.5));
        assert_eq!(stirlerr(1.0), 0.0810614667953272582196702);
        assert_eq!(dnbinom_mu_log(0.5, 1.0, 1.0), f64::NEG_INFINITY);
    }
}
