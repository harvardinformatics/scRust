//! R 4.5.3's root finder, `uniroot()` (`src/library/stats/R/nlm.R`), over the C routine
//! `R_zeroin2` (`src/library/stats/src/zeroin.c`) and its callback `fcn2` (`optimize.c`):
//! Brent's method from a bracketing interval whose function values are already known.

/// `uniroot()`'s result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Uniroot {
    pub root: f64,
    /// iterations used (`maxiter` when it did not converge)
    pub iter: i32,
    /// whether it converged
    pub converged: bool,
    /// the estimated precision, `|c - b|`
    pub estim_prec: f64,
}

/// What `fcn2` makes of the function's value: -Inf becomes `-DBL_MAX`, NaN and +Inf `DBL_MAX`.
fn fcn2_value(v: f64) -> f64 {
    if v.is_finite() {
        v
    } else if v == f64::NEG_INFINITY {
        -f64::MAX
    } else {
        f64::MAX
    }
}

/// `R_zeroin2(ax, bx, fa, fb, f, &tol, &maxit)`: returns the root estimate, the iterations used
/// (-1 when `maxit + 1` iterations were not enough) and the final `|c - b|`.
pub fn zeroin2(ax: f64, bx: f64, fa: f64, fb: f64, f: &mut dyn FnMut(f64) -> f64, tol: f64, maxit: i32) -> (f64, i32, f64) {
    const EPSILON: f64 = f64::EPSILON;
    let (mut a, mut b) = (ax, bx);
    let (mut fa, mut fb) = (fa, fb);
    let mut c = a;
    let mut fc = fa;
    let mut left = maxit + 1;
    // a root at an endpoint
    if fa == 0.0 {
        return (a, 0, 0.0);
    }
    if fb == 0.0 {
        return (b, 0, 0.0);
    }
    while left > 0 {
        left -= 1;
        let prev_step = b - a;
        if fc.abs() < fb.abs() {
            // swap data for b to be the best approximation
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol_act = 2.0 * EPSILON * b.abs() + tol / 2.0;
        let mut new_step = (c - b) / 2.0;
        if new_step.abs() <= tol_act || fb == 0.0 {
            // *Maxit -= maxit: `left` is C's loop counter after its post-decrement
            return (b, maxit - left, (c - b).abs());
        }
        // decide if the interpolation can be tried
        if prev_step.abs() >= tol_act && fa.abs() > fb.abs() {
            let cb = c - b;
            let (mut p, mut q);
            if a == c {
                // linear interpolation
                let t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                // inverse quadratic interpolation
                q = fa / fc;
                let t1 = fb / fc;
                let t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            if p < (0.75 * cb * q - (tol_act * q).abs() / 2.0) && p < (prev_step * q / 2.0).abs() {
                new_step = p / q;
            }
        }
        if new_step.abs() < tol_act {
            // not less than the tolerance
            new_step = if new_step > 0.0 { tol_act } else { -tol_act };
        }
        a = b;
        fa = fb;
        b += new_step;
        fb = f(b);
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            // c gets the sign opposite to b's
            c = a;
            fc = fa;
        }
    }
    // failed
    (b, -1, (c - b).abs())
}

/// `uniroot(f, c(lower, upper), tol = tol, maxiter = maxiter)` with `extendInt = "no"`: `f` is
/// evaluated at `lower` and `upper` first (as the defaults of `f.lower` and `f.upper`), then by
/// [`zeroin2`] (through `fcn2`). Errors where R stops.
pub fn uniroot(f: &mut dyn FnMut(f64) -> f64, lower: f64, upper: f64, tol: f64, maxiter: i32) -> Result<Uniroot, String> {
    if !(lower < upper) {
        return Err("lower < upper  is not fulfilled".into());
    }
    let f_lower = f(lower);
    let f_upper = f(upper);
    if f_lower.is_nan() {
        return Err("f.lower = f(lower) is NA".into());
    }
    if f_upper.is_nan() {
        return Err("f.upper = f(upper) is NA".into());
    }
    // truncate(): pmax.int(pmin(x, .Machine$double.xmax), -.Machine$double.xmax)
    let truncate = |x: f64| x.min(f64::MAX).max(-f64::MAX);
    let (f_low, f_upp) = (truncate(f_lower), truncate(f_upper));
    let sign = |x: f64| if x > 0.0 { 1.0 } else if x < 0.0 { -1.0 } else { 0.0 };
    if !(sign(f_lower) * sign(f_upper) <= 0.0) {
        return Err("f() values at end points not of opposite sign".into());
    }
    if !(tol.is_finite() && tol > 0.0) || maxiter <= 0 {
        return Err("invalid 'tol' or 'maxiter'".into());
    }
    let (root, iter, prec) = zeroin2(lower, upper, f_low, f_upp, &mut |x| fcn2_value(f(x)), tol, maxiter);
    Ok(if iter < 0 {
        Uniroot { root, iter: maxiter, converged: false, estim_prec: prec }
    } else {
        Uniroot { root, iter, converged: true, estim_prec: prec }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_sqrt2() {
        let r = uniroot(&mut |x| x * x - 2.0, 0.0, 2.0, 1e-12, 1000).unwrap();
        assert!((r.root - 2f64.sqrt()).abs() < 1e-11 && r.converged);
    }
}
