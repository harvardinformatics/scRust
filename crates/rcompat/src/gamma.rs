//! The gamma function family, as glmGamPoi's C++ calls it (milestone M4):
//!
//! - R's own `lgammafn`, `gammafn` (for `|x| <= 10`), `digamma` and `trigamma` (R 4.5.3's
//!   `src/nmath/lgamma.c`, `gamma.c`, `lgammacor.c`, `chebyshev.c` and `polygamma.c`'s `dpsifn`),
//!   ported statement by statement: Rcpp sugar's vectorised `lgamma`, `digamma` and `trigamma`
//!   and the `Rf_digamma`/`Rf_trigamma` calls go there;
//! - glibc's `lgamma` ([`lgamma_glibc`]): a plain `lgamma(double)` in C++ is the C library's,
//!   which glmGamPoi's log-likelihood calls on its count table. It is called through FFI
//!   (`lgamma_r`, the re-entrant form of the same function).

use crate::libm::{cos, exp, log, pow, sin};
use std::os::raw::c_int;

extern "C" {
    fn lgamma_r(x: f64, signp: *mut c_int) -> f64;
}

/// glibc's `lgamma(x)` (through `lgamma_r`, which computes the same value without touching the
/// global `signgam`).
pub fn lgamma_glibc(x: f64) -> f64 {
    let mut s: c_int = 0;
    // SAFETY: lgamma_r writes one int through the pointer, which is valid.
    unsafe { lgamma_r(x, &mut s) }
}

const M_LN_SQRT_2PI: f64 = 0.918938533204672741780329736406;
const M_LN_SQRT_PID2: f64 = 0.225791352644727432363097614947;
const M_LOG10_2: f64 = 0.301029995663981195213738894724;
const M_PI: f64 = std::f64::consts::PI;

/// `chebyshev_eval(x, a, n)`.
fn chebyshev_eval(x: f64, a: &[f64], n: usize) -> f64 {
    if !(1..=1000).contains(&n) || !(-1.1..=1.1).contains(&x) {
        return f64::NAN;
    }
    let twox = x * 2.0;
    let (mut b0, mut b1, mut b2) = (0.0, 0.0, 0.0);
    for i in 1..=n {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    (b0 - b2) * 0.5
}

const ALGMCS: [f64; 15] = [
    0.1666389480451863247205729650822e+0,
    -0.1384948176067563840732986059135e-4,
    0.9810825646924729426157171547487e-8,
    -0.1809129475572494194263306266719e-10,
    0.6221098041892605227126015543416e-13,
    -0.3399615005417721944303330599666e-15,
    0.2683181998482698748957538846666e-17,
    -0.2868042435334643284144622399999e-19,
    0.3962837061046434803679306666666e-21,
    -0.6831888753985766870111999999999e-23,
    0.1429227355942498147573333333333e-24,
    -0.3547598158101070547199999999999e-26,
    0.1025680058010470912000000000000e-27,
    -0.3401102254316748799999999999999e-29,
    0.1276642195630062933333333333333e-30,
];

/// `lgammacor(x)` for `x >= 10`.
fn lgammacor(x: f64) -> f64 {
    const NALGM: usize = 5;
    const XBIG: f64 = 94906265.62425156;
    if x < 10.0 {
        f64::NAN
    } else if x < XBIG {
        let tmp = 10.0 / x;
        chebyshev_eval(tmp * tmp * 2.0 - 1.0, &ALGMCS, NALGM) / x
    } else {
        1.0 / (x * 12.0)
    }
}

const GAMCS: [f64; 42] = [
    0.8571195590989331421920062399942e-2,
    0.4415381324841006757191315771652e-2,
    0.5685043681599363378632664588789e-1,
    -0.4219835396418560501012500186624e-2,
    0.1326808181212460220584006796352e-2,
    -0.1893024529798880432523947023886e-3,
    0.3606925327441245256578082217225e-4,
    -0.6056761904460864218485548290365e-5,
    0.1055829546302283344731823509093e-5,
    -0.1811967365542384048291855891166e-6,
    0.3117724964715322277790254593169e-7,
    -0.5354219639019687140874081024347e-8,
    0.9193275519859588946887786825940e-9,
    -0.1577941280288339761767423273953e-9,
    0.2707980622934954543266540433089e-10,
    -0.4646818653825730144081661058933e-11,
    0.7973350192007419656460767175359e-12,
    -0.1368078209830916025799499172309e-12,
    0.2347319486563800657233471771688e-13,
    -0.4027432614949066932766570534699e-14,
    0.6910051747372100912138336975257e-15,
    -0.1185584500221992907052387126192e-15,
    0.2034148542496373955201026051932e-16,
    -0.3490054341717405849274012949108e-17,
    0.5987993856485305567135051066026e-18,
    -0.1027378057872228074490069778431e-18,
    0.1762702816060529824942759660748e-19,
    -0.3024320653735306260958772112042e-20,
    0.5188914660218397839717833550506e-21,
    -0.8902770842456576692449251601066e-22,
    0.1527474068493342602274596891306e-22,
    -0.2620731256187362900257328332799e-23,
    0.4496464047830538670331046570666e-24,
    -0.7714712731336877911703901525333e-25,
    0.1323635453126044036486572714666e-25,
    -0.2270999412942928816702313813333e-26,
    0.3896418998003991449320816639999e-27,
    -0.6685198115125953327792127999999e-28,
    0.1146998663140024384347613866666e-28,
    -0.1967938586345134677295103999999e-29,
    0.3376448816585338090334890666666e-30,
    -0.5793070335782135784625493333333e-31,
];

/// R's `gammafn(x)` for `|x| <= 10` (the only range `lgammafn` uses it in); `None` outside.
pub fn gammafn_small(x: f64) -> Option<f64> {
    const NGAM: usize = 22;
    const XSML: f64 = 2.2474362225598545e-308;
    if x.is_nan() {
        return Some(x);
    }
    if x == 0.0 || (x < 0.0 && x == x.round()) {
        return Some(f64::NAN);
    }
    if x.abs() > 10.0 {
        return None;
    }
    let mut n = x as i32;
    if x < 0.0 {
        n -= 1;
    }
    let y = x - n as f64;
    n -= 1;
    let mut value = chebyshev_eval(y * 2.0 - 1.0, &GAMCS, NGAM) + 0.9375;
    if n == 0 {
        return Some(value);
    }
    if n < 0 {
        if y < XSML {
            return Some(if x > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY });
        }
        let n = -n;
        for i in 0..n {
            value /= x + i as f64;
        }
        Some(value)
    } else {
        for i in 1..=n {
            value *= y + i as f64;
        }
        Some(value)
    }
}

/// `sinpi(x)` (R's `src/nmath/cospi.c`).
fn sinpi(x: f64) -> f64 {
    if x.is_nan() || !x.is_finite() {
        return f64::NAN;
    }
    let mut x = x % 2.0; // fmod(x, 2.): sign of x
    if x <= -1.0 {
        x += 2.0;
    } else if x > 1.0 {
        x -= 2.0;
    }
    if x == 0.0 || x == 1.0 {
        return 0.0;
    }
    if x == 0.5 {
        return 1.0;
    }
    if x == -0.5 {
        return -1.0;
    }
    sin(M_PI * x)
}

/// R's `lgammafn(x)`: `log|Gamma(x)|`.
pub fn lgammafn(x: f64) -> f64 {
    const XMAX: f64 = 2.5327372760800758e+305;
    if x.is_nan() {
        return x;
    }
    if x <= 0.0 && x == x.trunc() {
        return f64::INFINITY;
    }
    let y = x.abs();
    if y < 1e-306 {
        return -log(y);
    }
    if y <= 10.0 {
        return log(gammafn_small(x).unwrap().abs());
    }
    if y > XMAX {
        return f64::INFINITY;
    }
    if x > 0.0 {
        if x > 1e17 {
            return x * (log(x) - 1.0);
        } else if x > 4934720.0 {
            return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x;
        } else {
            return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
        }
    }
    let sinpiy = sinpi(y).abs();
    if sinpiy == 0.0 {
        return f64::NAN;
    }
    M_LN_SQRT_PID2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y)
}

/// `d_n_cot(x, n)`: the n-th derivative of `cot` at `x` (polygamma.c), n in 0..=5.
fn d_n_cot(x: f64, n: i32) -> f64 {
    use crate::arith::r_pow_di;
    match n {
        0 => cos(x) / sin(x),
        1 => -1.0 / r_pow_di(sin(x), 2),
        2 => 2.0 * cos(x) / r_pow_di(sin(x), 3),
        3 => {
            let sin2 = r_pow_di(sin(x), 2);
            -2.0 * (3.0 - 2.0 * sin2) / r_pow_di(sin2, 2)
        }
        4 => {
            let co = cos(x);
            8.0 * co * (r_pow_di(co, 2) + 2.0) / r_pow_di(sin(x), 5)
        }
        5 => {
            let co2 = r_pow_di(cos(x), 2);
            -8.0 * (2.0 * r_pow_di(co2, 2) + 11.0 * co2 + 2.0) / r_pow_di(sin(x), 6)
        }
        _ => f64::NAN,
    }
}

const BVALUES: [f64; 22] = [
    1.00000000000000000e+00,
    -5.00000000000000000e-01,
    1.66666666666666667e-01,
    -3.33333333333333333e-02,
    2.38095238095238095e-02,
    -3.33333333333333333e-02,
    7.57575757575757576e-02,
    -2.53113553113553114e-01,
    1.16666666666666667e+00,
    -7.09215686274509804e+00,
    5.49711779448621554e+01,
    -5.29124242424242424e+02,
    6.19212318840579710e+03,
    -8.65802531135531136e+04,
    1.42551716666666667e+06,
    -2.72982310678160920e+07,
    6.01580873900642368e+08,
    -1.51163157670921569e+10,
    4.29614643061166667e+11,
    -1.37116552050883328e+13,
    4.88332318973593167e+14,
    -1.92965793419400681e+16,
];

const N_MAX: usize = 100;

/// `dpsifn(x, n, kode = 1, m = 1, ...)`: `(-1)^(n+1) / gamma(n+1) * psi(n, x)`, or the error
/// code `ierr`.
fn dpsifn1(x: f64, n: i32) -> Result<f64, i32> {
    let (kode, m) = (1, 1usize);
    let mut ans = [0.0f64; 1];
    if n < 0 {
        return Err(1);
    }
    if x <= 0.0 {
        if x == x.round() {
            return Ok(if n % 2 != 0 { f64::INFINITY } else { f64::NAN });
        }
        let mut a = dpsifn1(1.0 - x, n)?;
        if n > 5 {
            return Err(4);
        }
        let xp = x * M_PI;
        let (mut t1, mut t2, mut s) = (1.0f64, 1.0f64, 1.0f64);
        let mut k = 0i32;
        let mut j = k - n;
        while j < m as i32 {
            t1 *= M_PI;
            if k >= 2 {
                t2 *= k as f64;
            }
            if j >= 0 {
                a = s * (a + t1 / t2 * d_n_cot(xp, k));
            }
            k += 1;
            j += 1;
            s = -s;
        }
        return Ok(a);
    }
    let xln = log(x);
    if kode == 1 && m == 1 {
        let lrg = 1.0 / (2.0 * f64::EPSILON);
        if n == 0 && x * xln > lrg {
            return Ok(-xln);
        } else if n >= 1 && x > n as f64 * lrg {
            return Ok(exp(-(n as f64) * xln) / n as f64);
        }
    }
    let mut nx: i32 = 1021; // imin2(-i1mach(15), i1mach(16))
    let r1m5 = M_LOG10_2;
    let r1m4 = f64::EPSILON * 0.5;
    let wdtol = if r1m4 > 0.5e-18 { r1m4 } else { 0.5e-18 };
    let elim = 2.302 * (nx as f64 * r1m5 - 3.0);
    let rln = {
        let v = r1m5 * 53.0;
        if v < 18.06 { v } else { 18.06 }
    };
    let mut fln = (if rln > 3.0 { rln } else { 3.0 }) - 3.0;
    let yint = 3.50 + 0.40 * fln;
    let slope = 0.21 + fln * (0.0006038 * fln + 0.008677);
    let mut mm = m;
    let mut trm = [0.0f64; 23];
    let mut trmr = vec![0.0f64; N_MAX + 1];
    let (mut nn, mut fn_);
    let (mut t, mut t1, mut t2, mut tk, mut s);
    let (mut xdmy, mut xdmln, mut xinc);
    let xmin;
    loop {
        nn = n + mm as i32 - 1;
        fn_ = nn as f64;
        t = (fn_ + 1.0) * xln;
        if t.abs() > elim {
            if t <= 0.0 {
                return Err(2);
            }
        } else {
            if x < wdtol {
                // R_pow_di(x, -n-1)
                return Ok(crate::arith::r_pow_di(x, -n - 1));
            }
            let mut xm = yint + slope * fn_;
            let mx = xm as i32 + 1;
            let xmin_try = mx as f64;
            if n != 0 {
                xm = -2.302 * rln - (if 0.0 < xln { 0.0 } else { xln });
                let arg = { let v = xm / n as f64; if 0.0 < v { 0.0 } else { v } };
                let eps = exp(arg);
                xm = if arg.abs() < 1.0e-3 { -arg } else { 1.0 - eps };
                fln = x * xm / eps;
                xm = xmin_try - x;
                if xm > 7.0 && fln < 15.0 {
                    // series branch below
                    return Ok(series(x, n, xln, fln));
                }
            }
            xdmy = x;
            xdmln = xln;
            xinc = 0.0;
            if x < xmin_try {
                nx = x as i32;
                xinc = xmin_try - nx as f64;
                xdmy = x + xinc;
                xdmln = log(xdmy);
            }
            t = fn_ * xdmln;
            t1 = xdmln + xdmln;
            t2 = t + xdmln;
            tk = t.abs().max_r(t1.abs().max_r(t2.abs()));
            if tk <= elim {
                xmin = xmin_try;
                break;
            }
        }
        // underflow: nz++, mm--, ans[mm] = 0
        mm -= 1;
        ans[mm] = 0.0;
        if mm == 0 {
            return Ok(ans[0]);
        }
    }
    let _ = xmin;
    // L10:
    let tss = exp(-t);
    let tt = 0.5 / xdmy;
    t1 = tt;
    let tst = wdtol * tt;
    if nn != 0 {
        t1 = tt + 1.0 / fn_;
    }
    let rxsq = 1.0 / (xdmy * xdmy);
    let ta = 0.5 * rxsq;
    t = (fn_ + 1.0) * ta;
    s = t * BVALUES[2];
    if s.abs() >= tst {
        tk = 2.0;
        for k in 4..=22 {
            t = t * ((tk + fn_ + 1.0) / (tk + 1.0)) * ((tk + fn_) / (tk + 2.0)) * rxsq;
            trm[k] = t * BVALUES[k - 1];
            if trm[k].abs() < tst {
                break;
            }
            s += trm[k];
            tk += 2.0;
        }
    }
    s = (s + t1) * tss;
    if xinc != 0.0 {
        nx = xinc as i32;
        let np = nn + 1;
        if nx as usize > N_MAX {
            return Err(3);
        }
        if nn == 0 {
            // L20
            for i in 1..=nx {
                s += 1.0 / (x + (nx - i) as f64);
            }
            return Ok(s - xdmln);
        }
        let mut xm = xinc - 1.0;
        let mut fx = x + xm;
        for i in 1..=nx as usize {
            trmr[i] = pow(fx, -(np as f64));
            s += trmr[i];
            xm -= 1.0;
            fx = x + xm;
        }
    }
    ans[mm - 1] = s;
    if fn_ == 0.0 {
        // L30
        return Ok(s - xdmln);
    }
    // m == 1: no lower derivatives to generate
    Ok(ans[0])
}

/// The series branch of `dpsifn` (large `n` relative to `x`), for `m = 1`.
fn series(x: f64, n: i32, xln: f64, fln: f64) -> f64 {
    let nn = fln as i32 + 1;
    let np = n + 1;
    let t1 = (n + 1) as f64 * xln;
    let t = exp(-t1);
    let mut s = t;
    let mut den = x;
    for _ in 1..=nn {
        den += 1.0;
        s += pow(den, -(np as f64));
    }
    s
}

trait MaxR {
    fn max_r(self, other: f64) -> f64;
}

impl MaxR for f64 {
    /// R's `fmax2(x, y)`: NaN if either is NaN, else the larger (`x < y ? y : x`).
    fn max_r(self, other: f64) -> f64 {
        if self.is_nan() || other.is_nan() {
            return self + other;
        }
        if self < other { other } else { self }
    }
}

/// R's `digamma(x)`.
pub fn digamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    match dpsifn1(x, 0) {
        Ok(a) => -a,
        Err(_) => f64::NAN,
    }
}

/// R's `trigamma(x)`.
pub fn trigamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    dpsifn1(x, 1).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn familiar_values() {
        assert_eq!(lgammafn(1.0), 0.0);
        assert_eq!(lgammafn(2.0), 0.0);
        assert!((lgammafn(0.5) - 0.5723649429247001).abs() < 1e-15);
        assert!((digamma(1.0) + 0.5772156649015329).abs() < 1e-15);
        assert!((trigamma(1.0) - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-15);
        assert!((lgamma_glibc(10.5) - lgammafn(10.5)).abs() < 1e-13);
    }
}
