//! R arithmetic: `round()` (`src/nmath/fround.c`), `^` (`R_POW`, `R_pow` and `R_pow_di` in
//! `src/main/arithmetic.c`), and `seq()` (`seq.default` in `src/library/base/R/seq.R`), from
//! R 4.5.3.

use crate::libm;
use crate::summary::{is_na, NA_REAL};
use std::fmt;

/// R's `NaN` (`R_NaN`), as R's constant and `math2` write it.
pub const R_NAN: f64 = f64::from_bits(0x7ff8_0000_0000_0000);

/// `math2`'s NA handling, which R applies before `round(x, digits)` reaches `fround`: `NA` if
/// either argument is `NA`, else `NaN` if either is NaN.
fn math2_nan(a: f64, b: f64) -> Option<f64> {
    if is_na(a) || is_na(b) {
        Some(NA_REAL)
    } else if a.is_nan() || b.is_nan() {
        Some(R_NAN)
    } else {
        None
    }
}

/// `round(x)`: `nearbyint`, i.e. half to even (`round(0.5) == 0`, `round(2.5) == 2`).
pub fn round(x: f64) -> f64 {
    math2_nan(x, 0.0).unwrap_or_else(|| x.round_ties_even())
}

/// `round(x, digits)`: R 4's algorithm, which picks the closer of the two candidates
/// `floor(x * 10^d) / 10^d` and `ceil(x * 10^d) / 10^d` (the even one on a tie).
pub fn round_digits(x: f64, digits: f64) -> f64 {
    const DBL_DIG: f64 = 15.0;
    const MAX10E: i32 = 308; // DBL_MAX_10_EXP
    const MAX_DIGITS: f64 = 308.0 + 15.0;
    // M_LOG10_2 (0.301029995663981195213738894724) is the same double as LOG10_2
    const M_LOG10_2: f64 = std::f64::consts::LOG10_2;
    if let Some(nan) = math2_nan(x, digits) {
        return nan;
    }
    if !x.is_finite() {
        return x;
    }
    if digits > MAX_DIGITS || x == 0.0 {
        return x;
    } else if digits < -(MAX10E as f64) {
        return 0.0;
    } else if digits == 0.0 {
        return round(x);
    }
    let dig = (digits + 0.5).floor() as i32;
    let (sgn, x) = if x < 0.0 { (-1.0, -x) } else { (1.0, x) };
    let l10x = M_LOG10_2 * (0.5 + libm::logb(x));
    if l10x + dig as f64 > DBL_DIG {
        return sgn * x;
    }
    let (i10, xd, xu);
    if dig <= MAX10E {
        let pow10 = r_pow_di(10., dig);
        let x10 = x * pow10;
        i10 = x10.floor();
        xd = i10 / pow10;
        xu = x10.ceil() / pow10;
    } else {
        let e10 = dig - MAX10E;
        let p10 = r_pow_di(10., e10);
        let pow10 = r_pow_di(10., MAX10E);
        let x10 = (x * pow10) * p10;
        i10 = x10.floor();
        xd = i10 / pow10 / p10;
        xu = x10.ceil() / pow10 / p10;
    }
    let du = xu - x;
    let dd = x - xd;
    sgn * if du < dd || (i10 % 2.0 == 1.0 && du == dd) { xu } else { xd }
}

/// `x ^ y` for doubles: `x * x` when y is 2, R's own rules for 1, 0, zero and non-finite
/// operands, and glibc's `pow` for everything else.
pub fn r_pow(x: f64, y: f64) -> f64 {
    if y == 2.0 {
        return x * x;
    }
    if x == 1. || y == 0. {
        return 1.;
    }
    if x == 0. {
        return if y > 0. {
            0.
        } else if y < 0. {
            f64::INFINITY
        } else {
            y // NA or NaN
        };
    }
    if x.is_finite() && y.is_finite() {
        return libm::pow(x, y);
    }
    if x.is_nan() || y.is_nan() {
        // `return(x + y)`: R's build evaluates this as SSE y + x, which returns its first NaN
        // operand quietened, so y's NaN when both are NaN (e.g. NaN^NA is NA)
        let quiet = |v: f64| f64::from_bits(v.to_bits() | (1 << 51));
        return if y.is_nan() { quiet(y) } else { quiet(x) };
    }
    if !x.is_finite() {
        if x > 0. {
            return if y < 0. { 0. } else { f64::INFINITY };
        } else if y.is_finite() && y == y.floor() {
            // (-Inf) ^ n
            return if y < 0. {
                0.
            } else if myfmod(y, 2.) != 0. {
                x
            } else {
                -x
            };
        }
    }
    if !y.is_finite() && x >= 0. {
        return if y > 0. {
            if x >= 1. { f64::INFINITY } else { 0. }
        } else if x < 1. {
            f64::INFINITY
        } else {
            0.
        };
    }
    f64::NAN // (-Inf)^(+-Inf or non-integer), (negative)^(+-Inf)
}

/// `myfmod(x1, x2)` from arithmetic.c, for the finite integer-valued `x1` and `x2 = 2` that
/// `R_pow` passes it: equivalent to `fmod` there.
fn myfmod(x1: f64, x2: f64) -> f64 {
    let q = x1 / x2;
    let tmp = x1 - q.floor() * x2;
    tmp - (tmp / x2).floor() * x2
}

/// `R_pow_di(x, n)`: x to an integer power by repeated squaring.
pub fn r_pow_di(mut x: f64, mut n: i32) -> f64 {
    let mut xn = 1.0;
    if x.is_nan() {
        return x;
    }
    if n != 0 {
        if !x.is_finite() {
            return r_pow(x, n as f64);
        }
        let is_neg = n < 0;
        if is_neg {
            n = -n;
        }
        loop {
            if n & 1 != 0 {
                xn *= x;
            }
            n >>= 1;
            if n != 0 {
                x *= x;
            } else {
                break;
            }
        }
        if is_neg {
            xn = 1. / xn;
        }
    }
    xn
}

/// The conditions under which `seq()` stops with an error.
#[derive(Debug, Clone, PartialEq)]
pub enum SeqError {
    /// "'from' must be a finite number" (or 'to')
    NotFinite,
    /// "invalid '(to - from)/by'"
    InvalidBy,
    /// "wrong sign in 'by' argument"
    WrongSign,
    /// "'by' argument is much too small"
    ByTooSmall,
    /// "'length.out' must be a non-negative number"
    InvalidLength,
}

impl fmt::Display for SeqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            SeqError::NotFinite => "'from' and 'to' must be finite numbers",
            SeqError::InvalidBy => "invalid '(to - from)/by'",
            SeqError::WrongSign => "wrong sign in 'by' argument",
            SeqError::ByTooSmall => "'by' argument is much too small",
            SeqError::InvalidLength => "'length.out' must be a non-negative number",
        })
    }
}

impl std::error::Error for SeqError {}

/// `seq(from, to, by)` for doubles: `from + (0:n) * by` with `n = as.integer((to - from)/by +
/// 1e-10)`, then `pmin(x, to)` (or `pmax` for negative `by`) against overshoot. A span below 100
/// machine epsilons of the endpoints' magnitude gives just `from`.
pub fn seq_by(from: f64, to: f64, by: f64) -> Result<Vec<f64>, SeqError> {
    if !from.is_finite() || !to.is_finite() {
        return Err(SeqError::NotFinite);
    }
    let del = to - from;
    if del == 0.0 && to == 0.0 {
        return Ok(vec![to]);
    }
    let finite_del = del.is_finite();
    let n = if finite_del { del / by } else { to / by - from / by };
    if !n.is_finite() {
        if !by.is_nan() && by == 0.0 && del == 0.0 {
            return Ok(vec![from]);
        }
        return Err(SeqError::InvalidBy);
    }
    if n < 0.0 {
        return Err(SeqError::WrongSign);
    }
    if n > i32::MAX as f64 {
        return Err(SeqError::ByTooSmall);
    }
    if finite_del && del.abs() / to.abs().max(from.abs()) < 100.0 * f64::EPSILON {
        return Ok(vec![from]);
    }
    let n = (n + 1e-10) as i32; // as.integer truncates
    let x = (0..=n).map(|i| {
        if finite_del {
            from + i as f64 * by
        } else {
            (from / 4.0 + i as f64 * (by / 4.0)) * 4.0
        }
    });
    // pmin(x, to) replaces x[i] only when to < x[i]; pmax only when to > x[i]
    Ok(if by > 0.0 {
        x.map(|v| if to < v { to } else { v }).collect()
    } else {
        x.map(|v| if to > v { to } else { v }).collect()
    })
}

/// `seq(from, to, length.out = n)` for doubles: `c(from, from + (1:(n-2)) * ((to - from)/(n-1)),
/// to)`, so both ends are exact; `rep(from, n)` if `from == to`; `c(from, to)[1:n]` for n <= 2.
/// A non-integer `length_out` is rounded up, as R does.
pub fn seq_length_out(from: f64, to: f64, length_out: f64) -> Result<Vec<f64>, SeqError> {
    if !from.is_finite() || !to.is_finite() {
        return Err(SeqError::NotFinite);
    }
    if !length_out.is_finite() || length_out < 0.0 {
        return Err(SeqError::InvalidLength);
    }
    let len = length_out.ceil();
    if len == 0.0 {
        return Ok(Vec::new());
    }
    let len_u = len as usize;
    if len_u <= 2 {
        return Ok([from, to][..len_u].to_vec());
    }
    if from == to {
        return Ok(vec![from; len_u]);
    }
    let n1 = len - 1.0;
    let del = to - from;
    let mut out = Vec::with_capacity(len_u);
    if del.is_finite() {
        let by = del / n1;
        out.push(from);
        out.extend((1..len_u - 1).map(|i| from + i as f64 * by));
        out.push(to);
    } else {
        let (from4, to4) = (from / 4.0, to / 4.0);
        let by = (to4 - from4) / n1;
        out.push(from4);
        out.extend((1..len_u - 1).map(|i| from4 + i as f64 * by));
        out.push(to4);
        for v in out.iter_mut() {
            *v *= 4.0;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_half_even() {
        assert_eq!(round(0.5), 0.0);
        assert_eq!(round(1.5), 2.0);
        assert_eq!(round(2.5), 2.0);
        assert_eq!(round(-2.5), -2.0);
        // DoubletFinder: n = 1122, pN = 0.2 gives 280 artificial doublets
        assert_eq!(round(1122.0 / (1.0 - 0.2) - 1122.0), 280.0);
        assert_eq!(round_digits(2.675, 2.0), 2.67);
        assert_eq!(round_digits(0.125, 2.0), 0.12);
    }

    #[test]
    fn powers() {
        assert_eq!(r_pow(3.0, 2.0), 9.0);
        assert_eq!(r_pow(-0.0, -1.0), f64::INFINITY);
        assert!(r_pow(-1.0, f64::INFINITY).is_nan());
        assert!(r_pow(f64::NEG_INFINITY, 0.5).is_nan());
        assert_eq!(r_pow(f64::NEG_INFINITY, 3.0), f64::NEG_INFINITY);
        assert_eq!(r_pow(f64::NAN, 0.0), 1.0);
        assert_eq!(r_pow_di(10.0, -2), 1.0 / 100.0);
    }

    #[test]
    fn sequences() {
        let pk = seq_by(0.01, 0.3, 0.01).unwrap();
        assert_eq!(pk.len(), 30);
        assert_eq!(*pk.last().unwrap(), 0.3);
        assert_eq!(seq_by(1.0, 2.0, 3.0).unwrap(), vec![1.0]);
        assert_eq!(seq_by(1.0, 1.0, 0.0).unwrap(), vec![1.0]);
        assert_eq!(seq_by(1.0, 0.0, 1.0), Err(SeqError::WrongSign));
        let g = seq_length_out(-1.0, 2.0, 401.0).unwrap();
        assert_eq!((g.len(), g[0], g[400]), (401, -1.0, 2.0));
        assert_eq!(seq_length_out(3.0, 3.0, 4.0).unwrap(), vec![3.0; 4]);
        assert_eq!(seq_length_out(0.0, 1.0, 1.0).unwrap(), vec![0.0]);
    }
}
