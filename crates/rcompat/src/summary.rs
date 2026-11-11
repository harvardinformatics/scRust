//! `sum()`, `mean()`, `var()` and `sd()` of double vectors, with R 4.5.3's long-double
//! arithmetic (`rsum` and `real_mean` in `src/main/summary.c`; `cov_na_1` and `cov_complete1`
//! in `src/library/stats/src/cov.c`), emulated exactly by [`F80`].

use crate::ldouble::F80;

/// R's `NA_real_`: a NaN whose low word is 1954.
pub const NA_REAL: f64 = f64::from_bits(0x7ff0_0000_0000_07a2);

/// `ISNA(x)`: true for R's `NA` (and not for other NaNs).
pub fn is_na(x: f64) -> bool {
    x.is_nan() && (x.to_bits() & 0xffff_ffff) == 1954
}

fn ld(x: f64) -> F80 {
    F80::from_f64(x)
}

fn ld_n(n: usize) -> F80 {
    F80::from_u64(n as u64)
}

/// `sum(x, na.rm)`: accumulates in long double in index order; a total beyond the double range
/// becomes +-Inf, anything else is rounded to double.
pub fn sum(x: &[f64], na_rm: bool) -> f64 {
    let mut s = F80::ZERO;
    for &v in x {
        if !na_rm || !v.is_nan() {
            s = s.add_f64(v);
        }
    }
    if s > ld(f64::MAX) {
        f64::INFINITY
    } else if s < ld(-f64::MAX) {
        f64::NEG_INFINITY
    } else {
        s.to_f64()
    }
}

/// `mean(x)` (`mean.default` without trimming; with `na_rm`, NA and NaN are dropped first).
///
/// `real_mean`: the long-double sum divided by n, then a second pass adds the mean of the
/// residuals, `t += x[i] - s` in long double. If the sum overflows the double range, the mean is
/// recomputed from `x[i]/n` (each quotient in double) and corrected with `(x[i] - s)/n`.
pub fn mean(x: &[f64], na_rm: bool) -> f64 {
    if na_rm && x.iter().any(|v| v.is_nan()) {
        let kept: Vec<f64> = x.iter().copied().filter(|v| !v.is_nan()).collect();
        return real_mean(&kept);
    }
    real_mean(x)
}

fn real_mean(x: &[f64]) -> f64 {
    let n = x.len();
    let nl = ld_n(n);
    let mut s = F80::ZERO;
    for &v in x {
        s = s.add_f64(v);
    }
    let finite_s = s.to_f64().is_finite();
    if finite_s {
        s = s / nl;
    } else {
        // each x[i]/n in double (which quietens a signalling NaN), added in long double
        s = F80::ZERO;
        for &v in x {
            s = s.add_f64(v / n as f64);
        }
    }
    if finite_s && s.to_f64().is_finite() {
        let mut t = F80::ZERO;
        for &v in x {
            t = t + (ld(v) - s);
        }
        s = s + t / nl;
    } else if s.to_f64().is_finite() {
        let mut t = F80::ZERO;
        for &v in x {
            t = t + (ld(v) - s) / nl;
        }
        s = s + t;
    }
    s.to_f64()
}

/// `var(x, na.rm)` of a vector: NA if fewer than two observations, or if any is NA/NaN and
/// `na_rm` is false.
///
/// Two passes, all in long double: the mean as `cov.c`'s MEAN macro computes it (the sum over n,
/// plus the mean residual; no overflow fallback), rounded to double; then the sum of
/// `(x[i] - m) * (x[i] - m)`, where `m` is held in a long-double variable, so each difference and
/// each product is rounded to 64 bits; finally `/ (n - 1)` in long double.
pub fn var(x: &[f64], na_rm: bool) -> f64 {
    let kept: Vec<f64>;
    let x = if x.iter().any(|v| v.is_nan()) {
        if !na_rm {
            return NA_REAL;
        }
        kept = x.iter().copied().filter(|v| !v.is_nan()).collect();
        &kept[..]
    } else {
        x
    };
    let n = x.len();
    if n <= 1 {
        return NA_REAL;
    }
    let nl = ld_n(n);
    let mut sum = F80::ZERO;
    for &v in x {
        sum = sum.add_f64(v);
    }
    let mut tmp = sum / nl;
    if tmp.to_f64().is_finite() {
        sum = F80::ZERO;
        for &v in x {
            sum = sum + (ld(v) - tmp);
        }
        tmp = tmp + sum / nl;
    }
    let xxm = ld(tmp.to_f64());
    let mut sum = F80::ZERO;
    for &v in x {
        let d = ld(v) - xxm;
        sum = sum + d * d;
    }
    (sum / ld_n(n - 1)).to_f64()
}

/// `sd(x, na.rm)`: `sqrt(var(x, na.rm))`. R's `sqrt` (`math1`) gives back an incoming NaN
/// unchanged, so `NA` stays exactly `NA`.
pub fn sd(x: &[f64], na_rm: bool) -> f64 {
    let v = var(x, na_rm);
    if v.is_nan() {
        v
    } else {
        v.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_double_accumulation() {
        // 1 + 10 * 1e-16: double accumulation would stay at 1
        let mut x = vec![1.0];
        x.extend(std::iter::repeat_n(1e-16, 10));
        assert_eq!(sum(&x, false), 1.000000000000001);
        // an intermediate beyond the double range is fine in long double
        assert_eq!(sum(&[1e308, 1e308, -1e308], false), 1e308);
        assert_eq!(sum(&[f64::MAX, f64::MAX], false), f64::INFINITY);
        assert_eq!(mean(&[1e308, 1e308], false), 1e308);
        assert!(sum(&[], false) == 0.0 && mean(&[], false).is_nan());
    }

    #[test]
    fn na_handling() {
        assert!(is_na(var(&[1.0, NA_REAL], false)));
        assert!(is_na(var(&[1.0], false)));
        assert_eq!(var(&[1.0, f64::NAN, 3.0], true), 2.0);
        assert!(is_na(sum(&[1.0, NA_REAL], false)));
        assert_eq!(sum(&[1.0, NA_REAL], true), 1.0);
        assert_eq!(sd(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0], false), 2.138089935299395);
    }
}
