//! Order statistics: `quantile(type = 7)`, `IQR`, `median` and `mad` (R 4.5.3's
//! `src/library/stats/R/quantile.R`, `median.R`, `mad.R` and `IQR`).
//!
//! R selects order statistics with a partial sort (`sort(x, partial = ...)`); the selected values
//! are the ones a full sort puts at those positions, so a full sort gives the same results. (Where
//! `-0` and `+0` both occur, which of them a partial sort selects is not reproduced; it only
//! changes the sign of a zero.) Inputs must not contain NaN: R stops (`quantile`) or returns `NA`
//! (`median`), and the callers here never pass any.

use crate::summary::mean;

fn sorted(x: &[f64]) -> Vec<f64> {
    assert!(x.iter().all(|v| !v.is_nan()), "order statistics of NaN");
    let mut s = x.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    s
}

/// `quantile(x, probs, type = 7, names = FALSE)`.
///
/// `index = 1 + (n - 1) * p`, `lo = floor(index)`, `hi = ceiling(index)`; the result is `x[lo]`,
/// replaced by `(1 - h) * x[lo] + h * x[hi]` with `h = index - lo` where `index > lo` and
/// `x[hi] != x[lo]`.
pub fn quantile7(x: &[f64], probs: &[f64]) -> Vec<f64> {
    let s = sorted(x);
    let n = s.len();
    probs
        .iter()
        .map(|&p| {
            // probs <- pmax(0, pmin(1, probs))
            let p = if 1.0 < p { 1.0 } else { p };
            let p = if 0.0 > p { 0.0 } else { p };
            let index = 1.0 + (n.max(1) - 1) as f64 * p;
            let lo = index.floor();
            let hi = index.ceil();
            if n == 0 {
                return f64::NAN;
            }
            let qs = s[lo as usize - 1];
            let xhi = s[hi as usize - 1];
            if index > lo && xhi != qs {
                let h = index - lo;
                (1.0 - h) * qs + h * xhi
            } else {
                qs
            }
        })
        .collect()
}

/// `IQR(x)`: `diff(quantile(x, c(0.25, 0.75)))`.
pub fn iqr(x: &[f64]) -> f64 {
    let q = quantile7(x, &[0.25, 0.75]);
    q[1] - q[0]
}

/// `median(x)`: the middle order statistic, or the long-double `mean()` of the two middle ones.
/// (Selected in linear time; the values are those of a full sort.)
pub fn median(x: &[f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return crate::summary::NA_REAL;
    }
    assert!(x.iter().all(|v| !v.is_nan()), "order statistics of NaN");
    let mut s = x.to_vec();
    let half = n.div_ceil(2); // (n + 1) %/% 2, 1-based
    let cmp = |a: &f64, b: &f64| a.partial_cmp(b).unwrap();
    let (_, &mut lo, upper) = s.select_nth_unstable_by(half - 1, cmp);
    if n % 2 == 1 {
        lo
    } else {
        let hi = upper.iter().copied().fold(f64::INFINITY, |m, v| if v < m { v } else { m });
        mean(&[lo, hi], false)
    }
}

/// `mad(x)`: `1.4826 * median(abs(x - median(x)))`.
pub fn mad(x: &[f64]) -> f64 {
    let center = median(x);
    let dev: Vec<f64> = x.iter().map(|&v| (v - center).abs()).collect();
    1.4826 * median(&dev)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_cases() {
        assert_eq!(quantile7(&[3.0, 1.0, 2.0, 4.0], &[0.25, 0.75]), vec![1.75, 3.25]);
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[3.0, 1.0, 2.0, 10.0]), 2.5);
        assert_eq!(mad(&[1.0, 2.0, 3.0, 4.0, 100.0]), 1.4826);
        assert_eq!(iqr(&[1.0, 2.0, 3.0, 4.0, 5.0]), 2.0);
    }
}
