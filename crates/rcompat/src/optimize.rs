//! R 4.5.3's one-dimensional optimisation, `optimize()` (`src/library/stats/R/nlm.R`, the C
//! routine `Brent_fmin` and its callback `fcn1` in `src/library/stats/src/optimize.c`).

/// `.Machine$double.eps^0.25`, `optimize()`'s default tolerance (2^-13).
pub const OPTIMIZE_TOL: f64 = 0.0001220703125;

/// What `fcn1` makes of the function's value: -Inf becomes `-DBL_MAX`, NaN and +Inf `DBL_MAX`
/// (each with a warning in R).
fn fcn1_value(v: f64) -> f64 {
    if v.is_finite() {
        v
    } else if v == f64::NEG_INFINITY {
        -f64::MAX
    } else {
        f64::MAX
    }
}

/// `Brent_fmin(ax, bx, f, tol)`: Brent's minimiser without derivatives on `[ax, bx]`. `f` gets
/// each evaluation point in R's order.
pub fn brent_fmin(ax: f64, bx: f64, f: &mut dyn FnMut(f64) -> f64, tol: f64) -> f64 {
    // c is the squared inverse of the golden ratio
    let c = (3.0 - 5.0f64.sqrt()) * 0.5;
    let mut eps = f64::EPSILON;
    eps = eps.sqrt();
    let mut a = ax;
    let mut b = bx;
    let mut v = a + c * (b - a);
    let mut w = v;
    let mut x = v;
    let mut d = 0.0f64;
    let mut e = 0.0f64;
    let mut fx = f(x);
    let mut fv = fx;
    let mut fw = fx;
    let tol3 = tol / 3.0;
    loop {
        let xm = (a + b) * 0.5;
        let tol1 = eps * x.abs() + tol3;
        let t2 = tol1 * 2.0;
        // stopping criterion
        if (x - xm).abs() <= t2 - (b - a) * 0.5 {
            break;
        }
        let mut p = 0.0f64;
        let mut q = 0.0f64;
        let mut r = 0.0f64;
        if e.abs() > tol1 {
            // fit a parabola
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.0;
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
        }
        let mut u;
        if p.abs() >= (q * 0.5 * r).abs() || p <= q * (a - x) || p >= q * (b - x) {
            // a golden-section step
            e = if x < xm { b - x } else { a - x };
            d = c * e;
        } else {
            // a parabolic-interpolation step
            d = p / q;
            u = x + d;
            // f must not be evaluated too close to ax or bx
            if u - a < t2 || b - u < t2 {
                d = tol1;
                if x >= xm {
                    d = -d;
                }
            }
        }
        // f must not be evaluated too close to x
        u = if d.abs() >= tol1 {
            x + d
        } else if d > 0.0 {
            x + tol1
        } else {
            x - tol1
        };
        let fu = f(u);
        // update a, b, v, w, and x
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    x
}

/// `optimize(f, lower = lower, upper = upper, maximum = TRUE, tol = tol)$maximum`: Brent's
/// minimiser of `-f`, each value passed through `fcn1`.
pub fn optimize_maximum(f: &mut dyn FnMut(f64) -> f64, lower: f64, upper: f64, tol: f64) -> f64 {
    assert!(lower.is_finite() && upper.is_finite() && lower < upper && tol.is_finite() && tol > 0.0, "optimize: invalid interval or tolerance");
    brent_fmin(lower, upper, &mut |x| fcn1_value(-f(x)), tol)
}

/// `optimize(f, lower = lower, upper = upper, tol = tol)$minimum`.
pub fn optimize_minimum(f: &mut dyn FnMut(f64) -> f64, lower: f64, upper: f64, tol: f64) -> f64 {
    assert!(lower.is_finite() && upper.is_finite() && lower < upper && tol.is_finite() && tol > 0.0, "optimize: invalid interval or tolerance");
    brent_fmin(lower, upper, &mut |x| fcn1_value(f(x)), tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_a_parabola_minimum() {
        let m = optimize_minimum(&mut |x| (x - 1.0) * (x - 1.0), -3.0, 5.0, OPTIMIZE_TOL);
        assert!((m - 1.0).abs() < 1e-4);
        assert_eq!(OPTIMIZE_TOL, (2.0f64).powi(-13));
    }
}
