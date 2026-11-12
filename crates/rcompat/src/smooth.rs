//! R 4.5.3's kernel regression smoother `ksmooth()` (`src/library/stats/R/ksmooth.R`, C routine
//! `BDRksmooth` in `src/library/stats/src/ksmooth.c`) and the bin codes of `cut()`
//! (`.bincode`, `bincode` in `src/main/util.c`).

use crate::libm::exp;
use crate::order::order;
use crate::summary::NA_REAL;

/// `ksmooth()`'s kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    Box,
    Normal,
}

/// `ksmooth(x, y, kernel, bandwidth, x.points = x_points)`: returns `sort(x.points)` and the
/// smoothed values there (NA where no point is within the kernel's reach). `x.points` must hold
/// no NaN (R's `sort` would drop it).
pub fn ksmooth(x: &[f64], y: &[f64], kernel: Kernel, bandwidth: f64, x_points: &[f64]) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(x.len(), y.len());
    assert!(x_points.iter().all(|v| !v.is_nan()), "ksmooth: NaN in x.points");
    // sort(x.points): the values in increasing order (stable; equal values are equal)
    let xp: Vec<f64> = order(x_points, false).into_iter().map(|i| x_points[i]).collect();
    let ord = order(x, false);
    let xs: Vec<f64> = ord.iter().map(|&i| x[i]).collect();
    let ys: Vec<f64> = ord.iter().map(|&i| y[i]).collect();
    let yp = bdr_ksmooth(&xs, &ys, &xp, kernel, bandwidth);
    (xp, yp)
}

fn dokern(x: f64, kernel: Kernel) -> f64 {
    match kernel {
        Kernel::Box => 1.0,
        Kernel::Normal => exp(-0.5 * x * x),
    }
}

/// `BDRksmooth`: `x` sorted increasingly, `xp` too.
fn bdr_ksmooth(x: &[f64], y: &[f64], xp: &[f64], kernel: Kernel, bw: f64) -> Vec<f64> {
    let n = x.len();
    let mut imin = 0usize;
    // the bandwidth is in units of half the inter-quartile range
    let (bw, cutoff) = match kernel {
        Kernel::Box => {
            let b = bw * 0.5;
            (b, b)
        }
        Kernel::Normal => {
            let b = bw * 0.3706506;
            (b, 4.0 * b)
        }
    };
    let mut yp = vec![0.0; xp.len()];
    if xp.is_empty() {
        return yp;
    }
    while imin < n && x[imin] < xp[0] - cutoff {
        imin += 1;
    }
    for (j, &x0) in xp.iter().enumerate() {
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        let mut i = imin;
        while i < n {
            if x[i] < x0 - cutoff {
                imin = i;
            } else {
                if x[i] > x0 + cutoff {
                    break;
                }
                let w = dokern((x[i] - x0).abs() / bw, kernel);
                num += w * y[i];
                den += w;
            }
            i += 1;
        }
        yp[j] = if den > 0.0 { num / den } else { NA_REAL };
    }
    yp
}

/// `.bincode(x, breaks, right, include.lowest)`, as `cut()` computes its codes: 1-based interval
/// numbers, `None` for NA (outside the breaks, or NaN). `breaks` must be sorted increasingly.
pub fn bincode(x: &[f64], breaks: &[f64], right: bool, include_lowest: bool) -> Vec<Option<usize>> {
    let nb = breaks.len();
    assert!(nb >= 2 && breaks.windows(2).all(|w| w[0] <= w[1]), "bincode: 'breaks' is not sorted");
    let nb1 = nb - 1;
    let lft = !right;
    x.iter()
        .map(|&xi| {
            if xi.is_nan() {
                return None;
            }
            let (mut lo, mut hi) = (0usize, nb1);
            if xi < breaks[lo] || breaks[hi] < xi || (xi == breaks[if lft { hi } else { lo }] && !include_lowest) {
                return None;
            }
            while hi - lo >= 2 {
                let new = (hi + lo) / 2;
                if xi > breaks[new] || (lft && xi == breaks[new]) {
                    lo = new;
                } else {
                    hi = new;
                }
            }
            Some(lo + 1)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_like_cut() {
        let b = [0.0, 1.0, 2.0, 3.0];
        assert_eq!(bincode(&[0.0, 0.5, 1.0, 2.5, 3.0, 3.5], &b, true, false), vec![None, Some(1), Some(1), Some(3), Some(3), None]);
        assert_eq!(bincode(&[0.0, 3.0], &b, true, true), vec![Some(1), Some(3)]);
        assert_eq!(bincode(&[0.0, 1.0, 3.0], &b, false, false), vec![Some(1), Some(2), None]);
    }

    #[test]
    fn constant_is_smoothed_to_itself_roughly() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let (_, yp) = ksmooth(&x, &[2.0; 4], Kernel::Normal, 1.0, &[2.5, 1.0]);
        assert!(yp.iter().all(|v| (v - 2.0).abs() < 1e-12));
    }
}
