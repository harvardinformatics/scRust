//! Linear and constant interpolation as R's `approx()` and `approxfun()` compute it: the R-level
//! preparation of `regularize.values` (`src/library/stats/R/approx.R`) and `approx1` and
//! `R_approxfun` in `src/library/stats/src/approx.c`, R 4.5.3.

use crate::order::order;
use crate::summary::{mean, NA_REAL};

/// `method = "linear"` (1) or `"constant"` (2, with `f`: the weight of the right-hand value).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method {
    Linear,
    Constant { f: f64 },
}

/// Why `approxfun()` stops.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ApproxError {
    /// "need at least two non-NA values to interpolate" (linear) or "zero non-NA points".
    TooFewPoints,
    /// "approx(): invalid f value" (constant method, f outside [0, 1]).
    InvalidF,
    /// x and y differ in length.
    Lengths,
}

impl std::fmt::Display for ApproxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ApproxError::TooFewPoints => "need at least two non-NA values to interpolate",
            ApproxError::InvalidF => "approx(): invalid f value",
            ApproxError::Lengths => "x and y lengths differ",
        })
    }
}

impl std::error::Error for ApproxError {}

/// The interpolating function `approxfun(x, y, method, yleft, yright, rule, f, ties = mean,
/// na.rm = TRUE)` returns.
#[derive(Debug, Clone)]
pub struct ApproxFun {
    x: Vec<f64>,
    y: Vec<f64>,
    method: Method,
    yleft: f64,
    yright: f64,
}

/// `regularize.values(x, y, ties = mean, na.rm = TRUE)`: drop pairs with an NA (or NaN) in x or
/// y, sort by x if x is unsorted (stable `order`), and collapse tied x values to the mean of
/// their y values (R's long-double `mean`), keeping the sorted unique x.
pub fn regularize_values(x: &[f64], y: &[f64]) -> Result<(Vec<f64>, Vec<f64>), ApproxError> {
    if x.len() != y.len() {
        return Err(ApproxError::Lengths);
    }
    let (mut xs, mut ys): (Vec<f64>, Vec<f64>) = x
        .iter()
        .zip(y)
        .filter(|(a, b)| !a.is_nan() && !b.is_nan())
        .map(|(&a, &b)| (a, b))
        .unzip();
    if xs.windows(2).any(|w| w[1] < w[0]) {
        let o = order(&xs, false);
        xs = o.iter().map(|&i| xs[i]).collect();
        ys = o.iter().map(|&i| ys[i]).collect();
    }
    if xs.windows(2).any(|w| w[1] == w[0]) {
        // sorted, so ties are adjacent; groups in order of first appearance
        let (mut ux, mut uy) = (Vec::new(), Vec::new());
        let mut start = 0;
        while start < xs.len() {
            let mut end = start + 1;
            while end < xs.len() && xs[end] == xs[start] {
                end += 1;
            }
            ux.push(xs[start]);
            uy.push(mean(&ys[start..end], false));
            start = end;
        }
        xs = ux;
        ys = uy;
    }
    Ok((xs, ys))
}

impl ApproxFun {
    /// `approxfun(x, y, method, rule = 1)`: NA outside the range of x on both sides.
    pub fn new(x: &[f64], y: &[f64], method: Method) -> Result<ApproxFun, ApproxError> {
        ApproxFun::with_rule(x, y, method, [1, 1])
    }

    /// `approxfun(x, y, method, rule = c(left, right))`: rule 1 gives NA outside the range of x,
    /// rule 2 the value at the closest end.
    pub fn with_rule(x: &[f64], y: &[f64], method: Method, rule: [u8; 2]) -> Result<ApproxFun, ApproxError> {
        let (x, y) = regularize_values(x, y)?;
        let nx = x.len();
        if nx <= 1 && (method == Method::Linear || nx == 0) {
            return Err(ApproxError::TooFewPoints);
        }
        if let Method::Constant { f } = method {
            if !f.is_finite() || !(0.0..=1.0).contains(&f) {
                return Err(ApproxError::InvalidF);
            }
        }
        let yleft = if rule[0] == 1 { NA_REAL } else { y[0] };
        let yright = if rule[1] == 1 { NA_REAL } else { y[nx - 1] };
        Ok(ApproxFun { x, y, method, yleft, yright })
    }

    /// The regularized abscissae and ordinates.
    pub fn points(&self) -> (&[f64], &[f64]) {
        (&self.x, &self.y)
    }

    /// The function at `v` (`R_approxfun`: NaN and NA pass through).
    pub fn eval(&self, v: f64) -> f64 {
        if v.is_nan() {
            return v;
        }
        approx1(v, &self.x, &self.y, self.method, self.yleft, self.yright)
    }

    /// The function at each of `v`.
    pub fn eval_all(&self, v: &[f64]) -> Vec<f64> {
        v.iter().map(|&t| self.eval(t)).collect()
    }
}

/// `approx1(v, x, y, n, Meth)`: bisection for the interval, then linear interpolation
/// `y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i]))`, or the constant method's
/// `y[i] * f1 + y[j] * f2` (each term skipped when its weight is 0).
fn approx1(v: f64, x: &[f64], y: &[f64], method: Method, ylow: f64, yhigh: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN; // R_NaN
    }
    let (mut i, mut j) = (0usize, n - 1);
    // handle out-of-domain points
    if v < x[i] {
        return ylow;
    }
    if v > x[j] {
        return yhigh;
    }
    // find the correct interval by bisection
    while i + 1 < j {
        let ij = (i + j) / 2;
        if v < x[ij] {
            j = ij;
        } else {
            i = ij;
        }
    }
    // interpolation
    if v == x[j] {
        return y[j];
    }
    if v == x[i] {
        return y[i];
    }
    match method {
        Method::Linear => y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i])),
        Method::Constant { f } => {
            let (f1, f2) = (1.0 - f, f);
            (if f1 != 0.0 { y[i] * f1 } else { 0.0 }) + (if f2 != 0.0 { y[j] * f2 } else { 0.0 })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::summary::is_na;

    #[test]
    fn linear_and_constant() {
        let fun = ApproxFun::new(&[1.0, 2.0, 4.0], &[10.0, 20.0, 0.0], Method::Linear).unwrap();
        assert_eq!(fun.eval_all(&[1.0, 1.5, 3.0, 4.0]), vec![10.0, 15.0, 10.0, 0.0]);
        assert!(is_na(fun.eval(0.5)) && is_na(fun.eval(4.5)));
        let c = ApproxFun::new(&[1.0, 2.0], &[10.0, 20.0], Method::Constant { f: 0.0 }).unwrap();
        assert_eq!(c.eval(1.5), 10.0);
        let r2 = ApproxFun::with_rule(&[1.0, 2.0], &[10.0, 20.0], Method::Linear, [2, 2]).unwrap();
        assert_eq!((r2.eval(0.0), r2.eval(3.0)), (10.0, 20.0));
    }

    #[test]
    fn regularizes_like_r() {
        // unsorted with a tie: sorted, and the tie collapsed to the mean
        let fun = ApproxFun::new(&[3.0, 1.0, 3.0, f64::NAN], &[1.0, 0.0, 2.0, 5.0], Method::Linear).unwrap();
        assert_eq!(fun.points(), (&[1.0, 3.0][..], &[0.0, 1.5][..]));
        assert_eq!(ApproxFun::new(&[1.0], &[1.0], Method::Linear).unwrap_err(), ApproxError::TooFewPoints);
    }
}
