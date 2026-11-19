//! R's LINPACK QR decomposition, as `lm()` makes it and `qr.resid()` applies it (milestone M15b:
//! SCTransform's `vars.to.regress`, docs/PLAN.md C3).
//!
//! Seurat 5.5.1's `RegressOutMatrix()` regresses every row of `scale.data` on the same covariates:
//! it takes `lm(GENE ~ <vars>, data, qr = TRUE)$qr` once, then `qr.resid(qr, y = <row>)` per row.
//! In R 4.5.3 that is:
//!
//! - **the QR** ([`lm_qr`]): `lm.fit()` calls `.Call(C_Cdqrls, x, y, tol = 1e-07, FALSE)`
//!   (`src/library/stats/src/lm.c`), which stops on a non-finite value of `x` ("NA/NaN/Inf in
//!   'x'") and calls `dqrls` (`src/appl/dqrls.f`), which calls [`dqrdc2`] (`src/appl/dqrdc2.f`):
//!   R's modification of LINPACK's `dqrdc`, Householder reflections with a limited column
//!   pivoting that moves each column whose norm falls below `tol` times its original norm to the
//!   right-hand edge, and computes the rank. `x` is the model matrix: the intercept (1.0), then the
//!   covariates in the formula's order;
//! - **the residuals** ([`QrResid`]): `qr.resid()` (`src/library/base/R/qr.R`) returns `y` when the
//!   rank is 0, else calls `dqrrsd` (`src/appl/dqrutl.f`), which calls `dqrsl(..., job = 10)`
//!   (`src/appl/dqrsl.f`): `qty = Q' y` in place, then `rsd` = `qty` with its first `rank`
//!   entries zeroed, multiplied by `Q`.
//!
//! `dqrdc2` and `dqrsl` are Fortran in `libR.so`, compiled without FMA (`-march=nocona`; no
//! `vfmadd` in their code), and call BLAS `dnrm2`, `dscal`, `ddot` and `daxpy`, which `libR.so`
//! takes from `libblas.so.3`, the conda env's OpenBLAS 0.3.34 (`ldd`, `nm -D`). Its kernels differ
//! by CPU (Haswell's `ddot` and `daxpy` fuse multiply and add), so the port calls the same
//! routines with the same lengths ([`crate::blas`]): the results are R's bit for bit under the
//! same kernel. OpenBLAS's level-1 kernels load unaligned and peel nothing for alignment, so a
//! result depends on the vector's length and values, not on where it lies in memory: the port
//! passes each Householder vector as a contiguous copy with its leading element replaced, where
//! `dqrsl` swaps `qraux(j)` into `x(j,j)` for the call and back. `dcopy` copies values unchanged
//! and is done in Rust.

use crate::blas::{daxpy, ddot, dnrm2, dscal};

/// `lm.fit()`'s default tolerance (`tol = 1e-07`).
pub const LM_TOL: f64 = 1e-7;

/// R's `qr` object as `lm(..., qr = TRUE)$qr` holds it.
#[derive(Debug, Clone, PartialEq)]
pub struct Qr {
    /// Rows (cells) and columns (the intercept and the covariates).
    pub n: usize,
    pub p: usize,
    /// `qr$qr`: the compact decomposition, `n x p`, column-major.
    pub qr: Vec<f64>,
    /// `qr$qraux`.
    pub qraux: Vec<f64>,
    /// `qr$pivot`: 1-based column numbers of `x`, in their order after pivoting.
    pub pivot: Vec<usize>,
    /// `qr$rank`.
    pub rank: usize,
    /// `qr$tol`.
    pub tol: f64,
}

/// LINPACK `dqrdc2` as R 4.5.3 has it (`src/appl/dqrdc2.f`, "modernized to F77, 2024-05"), on
/// `x` (`n x p`, column-major, overwritten with the decomposition) with `ldx = n`. Returns
/// `qraux`, `jpvt` (1-based; it starts as `1..p`, as `Cdqrls` sets it) and the rank `k`.
pub fn dqrdc2(x: &mut [f64], n: usize, p: usize, tol: f64) -> (Vec<f64>, Vec<usize>, usize) {
    assert_eq!(x.len(), n * p, "dqrdc2: x is not n x p");
    let mut qraux = vec![0.0f64; p];
    // work(p, 2): the columns' norms as updated, and their original norms (1 for a zero column)
    let mut work1 = vec![0.0f64; p];
    let mut work2 = vec![0.0f64; p];
    let mut jpvt: Vec<usize> = (1..=p).collect();
    let at = |i: usize, j: usize| i + j * n;
    // compute the norms of the columns of x
    if n > 0 {
        for j in 0..p {
            qraux[j] = dnrm2(&x[at(0, j)..at(0, j) + n]);
            work1[j] = qraux[j];
            work2[j] = qraux[j];
            if work2[j] == 0.0 {
                work2[j] = 1.0;
            }
        }
    }
    // perform the householder reduction of x; `l` and `k` are Fortran's (1-based)
    let lup = n.min(p);
    let mut k = p + 1;
    for l in 1..=lup {
        let c = l - 1; // column l, 0-based
        // cycle the columns from l to p left-to-right until one with non-negligible norm is
        // located; the check for l .ge. k avoids infinite cycling
        loop {
            if l >= k || qraux[c] >= work2[c] * tol {
                break;
            }
            for i in 0..n {
                let t = x[at(i, c)];
                for j in (l + 1)..=p {
                    x[at(i, j - 2)] = x[at(i, j - 1)];
                }
                x[at(i, p - 1)] = t;
            }
            let i = jpvt[c];
            let t = qraux[c];
            let tt = work1[c];
            let ttt = work2[c];
            for j in (l + 1)..=p {
                jpvt[j - 2] = jpvt[j - 1];
                qraux[j - 2] = qraux[j - 1];
                work1[j - 2] = work1[j - 1];
                work2[j - 2] = work2[j - 1];
            }
            jpvt[p - 1] = i;
            qraux[p - 1] = t;
            work1[p - 1] = tt;
            work2[p - 1] = ttt;
            k -= 1;
        }
        if l != n {
            // compute the householder transformation for column l: rows l..n
            let (head, tail) = x.split_at_mut(at(0, c) + n);
            let col_l = &mut head[at(c, c)..];
            let mut nrmxl = dnrm2(col_l);
            if nrmxl != 0.0 {
                // Fortran's sign(a, b): |a| with b's sign
                if col_l[0] != 0.0 {
                    nrmxl = nrmxl.abs().copysign(col_l[0]);
                }
                dscal(1.0 / nrmxl, col_l);
                col_l[0] = 1.0 + col_l[0];
                // apply the transformation to the remaining columns, updating the norms
                let xll = col_l[0];
                for j in (l + 1)..=p {
                    let cj = j - 1;
                    let off = at(0, cj) - (at(0, c) + n);
                    let col_j = &mut tail[off + c..off + n];
                    // t = -ddot(n-l+1, x(l,l), 1, x(l,j), 1) / x(l,l): divsd, then the sign flip
                    let t = -(ddot(col_l, col_j) / xll);
                    daxpy(t, col_l, col_j);
                    if qraux[cj] != 0.0 {
                        // tt = 1 - (|x(l,j)| / qraux(j))**2 (y**2 is y*y); max(tt, 0) as maxsd
                        let q = col_j[0].abs() / qraux[cj];
                        let tt = 1.0 - q * q;
                        let t = if tt > 0.0 { tt } else { 0.0 };
                        // re-compute the norm when it has dropped a lot (relative to the original)
                        if t.abs() >= 1e-6 {
                            qraux[cj] *= t.sqrt();
                        } else {
                            qraux[cj] = dnrm2(&col_j[1..]);
                            work1[cj] = qraux[cj];
                        }
                    }
                }
                // save the transformation
                qraux[c] = col_l[0];
                col_l[0] = -nrmxl;
            }
        }
    }
    let k = (k - 1).min(n);
    (qraux, jpvt, k)
}

/// The QR that `lm.fit(x, y)` returns (`Cdqrls` with `tol = 1e-07`): `x` is the model matrix,
/// `n x p`, column-major. R's checks first: `n == 0` stops ("0 (non-NA) cases"), and a non-finite
/// value of `x` stops ("NA/NaN/Inf in 'x'").
pub fn lm_qr(mut x: Vec<f64>, n: usize, p: usize) -> Result<Qr, String> {
    if x.len() != n * p {
        return Err(format!("lm: the model matrix holds {} values, not {n} x {p}", x.len()));
    }
    if n == 0 {
        return Err("lm.fit: 0 (non-NA) cases".into());
    }
    if p == 0 {
        return Err("lm.fit: a model matrix without columns (the null model) is not supported".into());
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err("lm.fit: NA/NaN/Inf in 'x'".into());
    }
    let (qraux, pivot, rank) = dqrdc2(&mut x, n, p, LM_TOL);
    Ok(Qr { n, p, qr: x, qraux, pivot, rank, tol: LM_TOL })
}

/// `model.matrix(GENE ~ v1 + v2 + ...)` for numeric covariates: the intercept (1.0), then each
/// covariate's values as they are, `n x (1 + q)`, column-major.
pub fn model_matrix(covariates: &[&[f64]], n: usize) -> Vec<f64> {
    let mut x = vec![1.0f64; n * (1 + covariates.len())];
    for (j, c) in covariates.iter().enumerate() {
        assert_eq!(c.len(), n, "model_matrix: a covariate does not have one value per row");
        x[(j + 1) * n..(j + 2) * n].copy_from_slice(c);
    }
    x
}

/// `qr.resid(qr, y)` for many `y` with one QR: the Householder vectors `dqrsl` uses, prepared
/// once. It only reads, so threads may share it.
#[derive(Debug, Clone)]
pub struct QrResid<'a> {
    qr: &'a Qr,
    /// `ju = min(rank, n - 1)` vectors: column `j`'s rows `j..n` with `qraux(j)` in place of
    /// `x(j,j)`, as `dqrsl` passes them to `ddot` and `daxpy` (empty where `qraux(j) == 0`,
    /// which `dqrsl` skips)
    v: Vec<Vec<f64>>,
}

impl<'a> QrResid<'a> {
    pub fn new(qr: &'a Qr) -> QrResid<'a> {
        let n = qr.n;
        let ju = if n == 0 { 0 } else { qr.rank.min(n - 1) };
        let v = (0..ju)
            .map(|j| {
                if qr.qraux[j] == 0.0 {
                    return Vec::new();
                }
                let mut v = qr.qr[j + j * n..(j + 1) * n].to_vec();
                v[0] = qr.qraux[j];
                v
            })
            .collect();
        QrResid { qr, v }
    }

    /// `qr.resid(qr, y)` into `rsd`; `qty` is the work vector `dqrsl` computes `Q' y` in (`y`'s
    /// own copy, in `dqrrsd`). All three hold `n` values.
    pub fn resid(&self, y: &[f64], qty: &mut [f64], rsd: &mut [f64]) {
        let (n, k) = (self.qr.n, self.qr.rank);
        assert!(y.len() == n && qty.len() == n && rsd.len() == n, "qr.resid: 'qr' and 'y' must have the same number of rows");
        if k == 0 {
            // qr.resid: if (k == 0) return(y)
            rsd.copy_from_slice(y);
            return;
        }
        // dqrsl, job = 10: qty, then rsd (cqty and cr; not qy, b or xb)
        let ju = k.min(n - 1);
        if ju == 0 {
            // special action when n = 1
            qty[0] = y[0];
            rsd[0] = 0.0;
            return;
        }
        // dcopy(n, y, 1, qty, 1)
        qty.copy_from_slice(y);
        // compute trans(q) * y
        for j in 0..ju {
            if self.qr.qraux[j] == 0.0 {
                continue;
            }
            let v = &self.v[j];
            let t = -(ddot(v, &qty[j..]) / v[0]);
            daxpy(t, v, &mut qty[j..]);
        }
        // set up to compute rsd: dcopy(n-k, qty(k+1), 1, rsd(k+1), 1), then rsd(1..k) = 0
        if k < n {
            rsd[k..].copy_from_slice(&qty[k..]);
        }
        for r in rsd.iter_mut().take(k) {
            *r = 0.0;
        }
        // compute rsd
        for j in (0..ju).rev() {
            if self.qr.qraux[j] == 0.0 {
                continue;
            }
            let v = &self.v[j];
            let t = -(ddot(v, &rsd[j..]) / v[0]);
            daxpy(t, v, &mut rsd[j..]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The residuals are orthogonal to the model's columns, and y minus them lies in its span.
    #[test]
    fn residuals_of_an_exact_fit_and_of_noise() {
        let n = 50;
        let v: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin() * 3.0).collect();
        let x = model_matrix(&[&v], n);
        let qr = lm_qr(x, n, 2).unwrap();
        assert_eq!((qr.rank, qr.pivot.clone()), (2, vec![1, 2]));
        let r = QrResid::new(&qr);
        // y = 2 + 5 v: residuals ~ 0
        let y: Vec<f64> = v.iter().map(|a| 2.0 + 5.0 * a).collect();
        let (mut qty, mut rsd) = (vec![0.0; n], vec![0.0; n]);
        r.resid(&y, &mut qty, &mut rsd);
        assert!(rsd.iter().all(|e| e.abs() < 1e-12), "{rsd:?}");
        // noise: residuals orthogonal to 1 and v
        let y: Vec<f64> = (0..n).map(|i| ((i * 7919) % 101) as f64 / 10.0).collect();
        r.resid(&y, &mut qty, &mut rsd);
        let s1: f64 = rsd.iter().sum();
        let s2: f64 = rsd.iter().zip(&v).map(|(a, b)| a * b).sum();
        assert!(s1.abs() < 1e-10 && s2.abs() < 1e-10, "{s1} {s2}");
    }

    /// A zero column is moved to the edge and leaves rank 1: the regression on the intercept
    /// alone, so the residuals are the centred values.
    #[test]
    fn a_zero_column_is_pivoted_out() {
        let n = 9;
        let zero = vec![0.0; n];
        let qr = lm_qr(model_matrix(&[&zero], n), n, 2).unwrap();
        assert_eq!((qr.rank, qr.pivot.clone()), (1, vec![1, 2]));
        let r = QrResid::new(&qr);
        let y: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let (mut qty, mut rsd) = (vec![0.0; n], vec![0.0; n]);
        r.resid(&y, &mut qty, &mut rsd);
        for (i, e) in rsd.iter().enumerate() {
            assert!((e - (i as f64 - 4.0)).abs() < 1e-12);
        }
        // the zero column first: it moves behind the intercept
        let x: Vec<f64> = zero.iter().copied().chain(std::iter::repeat(1.0).take(n)).collect();
        let qr = lm_qr(x, n, 2).unwrap();
        assert_eq!((qr.rank, qr.pivot.clone()), (1, vec![2, 1]));
    }

    #[test]
    fn r_stops_on_non_finite_values_and_no_rows() {
        assert!(lm_qr(vec![1.0, f64::NAN], 1, 2).unwrap_err().contains("NA/NaN/Inf"));
        assert!(lm_qr(vec![1.0, f64::INFINITY], 1, 2).unwrap_err().contains("NA/NaN/Inf"));
        assert!(lm_qr(Vec::new(), 0, 2).unwrap_err().contains("0 (non-NA) cases"));
    }
}
