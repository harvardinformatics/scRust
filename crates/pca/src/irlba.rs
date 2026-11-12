//! irlba 2.4.1's `irlba()` (`R/irlba.R`, `R/utility.R`), as `RunPCA` calls it: a dense real
//! matrix, `nu = nv`, and the defaults `maxit = 1000`, `work = nv + 16`, `reorth = TRUE`,
//! `tol = svtol = 1e-7`, no `v`, `scale`, `center` or `shift`, `smallest = FALSE`,
//! `right_only = FALSE`. irlba 2.4.1 runs only its R code (`fastpath` is deprecated).
//!
//! The augmented implicitly restarted Lanczos bidiagonalisation is ported statement by
//! statement: the Lanczos steps with full reorthogonalisation of both bases, the
//! near-invariant-subspace branches that draw `rnorm()` vectors from R's stream, `convtests`,
//! the "clustered spectrum" test (whose `rel_change_k` is always 0, as `lastsv` is updated just
//! before it), the working-dimension expansions, and the restarts with the Ritz vectors.
//!
//! Every product is the BLAS call R makes, with the same routine, shape, transposition flags
//! and leading dimensions:
//!
//! - `A %*% v` and `t(w) %*% A`: `dgemm` through `direct_dgemm_c` ([`Products`]);
//! - `norm2(x)`: `dgemm("T", "N", 1, 1, length(x), ...)` (`gemm(x, x, TRUE)`), then
//!   `sqrt(abs(.))`;
//! - `orthog(Y, X)` for a matrix `X`: `dgemm("T", "N", 1, ncol(X), nrow(X), Y, X)` for the
//!   coefficients, then R's `X %*% doty`, which is `dgemv("N")` unless NaN or Inf may be
//!   present; when R has dropped `X` to a plain vector (`W[, 1:j]` with `j = 1`), `cross()`
//!   falls back to `t(Y) %*% X` (`dgemv("N")` with one row), then `X %*% doty`;
//! - the restarts and the final vectors: R's `%*%` (`dgemm("N", "N")`, or `dgemv` for one
//!   column);
//! - `svd(B)`: LAPACK `dgesdd` with R's workspace query ([`rcompat::blas::la_svd`]).
//!
//! Elementwise arithmetic (`F - S * V[, j]`, `F / R`, ...) is IEEE double, as in R.

use crate::products::Products;
use crate::PcaError;
use rcompat::blas;
use rcompat::RRng;

/// irlba's arguments beyond `nv` (`irlba()`'s defaults).
#[derive(Debug, Clone, PartialEq)]
pub struct IrlbaParams {
    /// `maxit`: maximum number of restarts (default 1000)
    pub maxit: usize,
    /// `work`: working subspace dimension (default `nv + 16`)
    pub work: Option<usize>,
    /// `tol` (default 1e-7)
    pub tol: f64,
    /// `svtol` (default `tol`)
    pub svtol: Option<f64>,
    /// `reorth`: full reorthogonalisation of both bases (default `TRUE`)
    pub reorth: bool,
}

impl Default for IrlbaParams {
    fn default() -> Self {
        IrlbaParams { maxit: 1000, work: None, tol: 1e-7, svtol: None, reorth: true }
    }
}

/// `irlba()`'s value: `d`, `u`, `v`, `iter`, `mprod`, and what else happened.
#[derive(Debug, Clone)]
pub struct Irlba {
    /// the `nv` largest singular values, decreasing
    pub d: Vec<f64>,
    /// left singular vectors, `m x nv`, column-major
    pub u: Vec<f64>,
    /// right singular vectors, `n x nv`, column-major
    pub v: Vec<f64>,
    /// `iter` and `mprod` as irlba returns them
    pub iter: usize,
    pub mprod: usize,
    /// false when `maxit` was reached (R warns "did not converge")
    pub converged: bool,
    /// the working dimension at the end (more than at the start after an expansion)
    pub work: usize,
    /// vectors drawn with `rnorm()` in the near-invariant-subspace branches (beyond the start
    /// vector)
    pub invariant_draws: usize,
    /// the tiny-problem path (`min(m, n) < 6`: R's `svd()` instead of the iteration)
    pub tiny: bool,
}

/// A small dense matrix, column-major, for irlba's `B`.
#[derive(Debug, Clone)]
struct Mat {
    nrow: usize,
    ncol: usize,
    x: Vec<f64>,
}

impl Mat {
    fn at(&self, i: usize, j: usize) -> f64 {
        self.x[i + j * self.nrow]
    }

    /// `rbind(cbind(B, 0), c(rep(0, ncol(B) - 1), S, R))`
    fn grow(&self, s: f64, r: f64) -> Mat {
        let (nr, nc) = (self.nrow + 1, self.ncol + 1);
        let mut x = vec![0.0; nr * nc];
        for j in 0..self.ncol {
            for i in 0..self.nrow {
                x[i + j * nr] = self.at(i, j);
            }
        }
        x[self.nrow + (self.ncol - 1) * nr] = s;
        x[self.nrow + self.ncol * nr] = r;
        Mat { nrow: nr, ncol: nc, x }
    }

    /// `rbind(B, c(rep(0, j - 1), S))`, `j = ncol(B)`
    fn last_row(&self, s: f64) -> Mat {
        let nr = self.nrow + 1;
        let mut x = vec![0.0; nr * self.ncol];
        for j in 0..self.ncol {
            for i in 0..self.nrow {
                x[i + j * nr] = self.at(i, j);
            }
        }
        x[self.nrow + (self.ncol - 1) * nr] = s;
        Mat { nrow: nr, ncol: self.ncol, x }
    }
}

/// `svd(B)` as irlba uses it: `d`, `u` and `v = t(vt)`, all `bsz x bsz`.
struct BSvd {
    d: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    bsz: usize,
}

impl BSvd {
    fn of(b: &Mat) -> Result<BSvd, String> {
        let s = blas::la_svd(&b.x, b.nrow, b.ncol)?;
        // B is square here (bsz x bsz); v = t(vt)
        let (np, p) = (s.np(), s.p);
        let mut v = vec![0.0; p * np];
        for i in 0..np {
            for j in 0..p {
                v[j + i * p] = s.vt[i + j * np];
            }
        }
        Ok(BSvd { d: s.d, u: s.u, v, bsz: b.nrow })
    }

    /// The first `k` columns of `u` (or `v`): `bsz x k`, column-major.
    fn cols(x: &[f64], bsz: usize, k: usize) -> &[f64] {
        &x[..bsz * k]
    }
}

/// `norm2(x)`: `sqrt(abs(drop(gemm(x, x, TRUE))))`.
fn norm2(x: &[f64]) -> f64 {
    let n = x.len();
    let mut c = [0.0];
    blas::dgemm(b'T', b'N', 1, 1, n, 1.0, x, n, x, n, 0.0, &mut c, 1);
    c[0].abs().sqrt()
}

/// `Y <- orthog(Y, X)` for a `p x q` matrix `X` (an R matrix, `q >= 1`):
/// `doty <- Conj(t(gemm(Y, X, TRUE)))`, then `Y - X %*% doty`.
fn orthog_mat(y: &mut [f64], x: &[f64], p: usize, q: usize) {
    debug_assert!(y.len() == p && x.len() >= p * q);
    let mut c = vec![0.0; q];
    blas::dgemm(b'T', b'N', 1, q, p, 1.0, y, p, &x[..p * q], p, 0.0, &mut c, 1);
    let z = blas::matprod(&x[..p * q], p, q, &c, q, 1);
    for (yi, zi) in y.iter_mut().zip(&z) {
        *yi -= zi;
    }
}

/// `Y <- orthog(Y, X)` when `X` is a plain vector (R drops `W[, 1:1]`): `cross(Y, X)` is then
/// `Conj(t(Y)) %*% X` (a `1 x p` matrix times a vector), and `X %*% doty` a column vector times
/// a `1 x 1` matrix.
fn orthog_vec(y: &mut [f64], x: &[f64]) {
    let p = y.len();
    let c = blas::matprod(y, 1, p, x, p, 1);
    let z = blas::matprod(x, p, 1, &c, 1, 1);
    for (yi, zi) in y.iter_mut().zip(&z) {
        *yi -= zi;
    }
}

/// `W[, j] <- orthog(W[, j], W[, 1:(j - 1)])` with 0-based column `col` of a `p`-row matrix
/// `w`: against columns `0..col`, a plain vector when there is one (`drop = TRUE`).
fn orthog_against_previous(w: &mut [f64], p: usize, col: usize, drop: bool) {
    let (left, right) = w.split_at_mut(col * p);
    let y = &mut right[..p];
    if col == 1 && drop {
        orthog_vec(y, &left[..p]);
    } else {
        orthog_mat(y, left, p, col);
    }
}

/// `rnorm(len)`: irlba's `random`.
fn random(rng: &mut RRng, len: usize) -> Vec<f64> {
    rcompat::rnorm(rng, len, 0.0, 1.0)
}

/// R's `x < y` as a logical: `None` (NA) when either is NaN.
fn lt(x: f64, y: f64) -> Option<bool> {
    if x.is_nan() || y.is_nan() {
        None
    } else {
        Some(x < y)
    }
}

/// R's `&` on logicals.
fn and(a: Option<bool>, b: Option<bool>) -> Option<bool> {
    match (a, b) {
        (Some(false), _) | (_, Some(false)) => Some(false),
        (Some(true), Some(true)) => Some(true),
        _ => None,
    }
}

/// `convtests(Bsz, tol, k_org, Bsvd, residuals, k, Smax, lastsv, svtol, maxritz, S)`:
/// (converged, k).
#[allow(clippy::too_many_arguments)]
fn convtests(bsz: usize, tol: f64, k_org: usize, d: &[f64], residuals: &[f64], k: usize, smax: f64, lastsv: Option<&[f64]>, svtol: f64, s: f64) -> (bool, usize) {
    // subspace_converged <- residuals[1:k_org] < tol * Smax
    let thr = tol * smax;
    // delta_converged <- (abs(Bsvd$d[1:k_org] - lastsv[1:k_org]) / Bsvd$d[1:k_org]) < svtol,
    // with lastsv <- 0 * Bsvd$d the first time
    let mut len_res = 0usize;
    let mut na = false;
    for i in 0..k_org {
        let sub = lt(residuals[i], thr);
        let last = lastsv.map_or(0.0 * d[i], |l| l[i]);
        let delta = lt((d[i] - last).abs() / d[i], svtol);
        match and(sub, delta) {
            Some(true) => len_res += 1,
            Some(false) => {}
            None => na = true,
        }
    }
    // len_res <- sum(...); if (is.na(len_res)) len_res <- 0
    if na {
        len_res = 0;
    }
    if len_res >= k_org || s == 0.0 {
        return (true, k);
    }
    let k = k.max(k_org + len_res);
    let k = (k as i64).min((bsz as i64 - 3).max(1)) as usize;
    (false, k)
}

/// `irlba(A, nv)` with [`IrlbaParams`], `A` being `m x n` (`prod`); draws from `rng` exactly
/// what R draws.
pub fn irlba(prod: &Products, nv: usize, params: &IrlbaParams, rng: &mut RRng) -> Result<Irlba, PcaError> {
    let err = |m: &str| Err(PcaError::Input(format!("irlba: {m}")));
    let (m, n) = (prod.m, prod.n);
    let a = prod.a();
    // oknum(A): okatomic is TRUE for a single element, else every value must be finite
    if a.len() != 1 && a.iter().any(|v| !v.is_finite()) {
        return err("Data matrix A must contain finite, non-missing numeric (real or complex) values (atomic matrix- or Matrix-classed).");
    }
    if m.min(n) < 1 {
        return err("A is not matrix-like.");
    }
    let nu = nv;
    let mut k = nu.max(nv);
    if k == 0 {
        return err("nv and/or nu must be finite positive integer values");
    }
    if k > (m - 1).min(n - 1) {
        return err("max(nu, nv) must be strictly less than min(nrow(A), ncol(A))");
    }
    let mut work = params.work.unwrap_or(nv + 16);
    if work <= 1 {
        return err("work must be greater than 1");
    }
    let tol = params.tol;
    if tol.is_nan() || tol < 0.0 || !tol.is_finite() {
        return err("tol must be a non-negative real number");
    }
    if params.maxit == 0 {
        return err("maxit must be positive");
    }
    let svtol = params.svtol.unwrap_or(tol);
    let mut reorth = params.reorth;
    if work <= k {
        work = k + 1;
    }
    if work >= n.min(m) {
        work = (n.min(m) as i64 - 2).max(1) as usize;
        if work <= k {
            work = k + 1;
        }
    }
    let k_org = k;
    let w_dim = work;

    // Tiny problem: R's svd() instead of the iteration
    if m.min(n) < 6 {
        let s = blas::la_svd(a, m, n).map_err(|e| PcaError::Input(format!("irlba (svd): {e}")))?;
        let np = s.np();
        let d = s.d[..k].to_vec();
        let u = s.u[..m * nu].to_vec();
        let mut v = vec![0.0; n * nv];
        for i in 0..nv {
            for j in 0..n {
                v[j + i * n] = s.vt[i + j * np];
            }
        }
        return Ok(Irlba { d, u, v, iter: 0, mprod: 0, converged: true, work, invariant_draws: 0, tiny: true });
    }

    let mut w = vec![0.0; m * w_dim];
    let mut f = vec![0.0; n];
    let mut v = vec![0.0; n * work];
    v[..n].copy_from_slice(&random(rng, n));

    let eps = f64::EPSILON;
    let eps2 = rcompat::r_pow(eps, 4.0 / 5.0);
    let eps23 = rcompat::r_pow(eps, 2.0 / 3.0);
    let sqrteps = eps.sqrt();
    let mut b: Option<Mat> = None;
    let mut iter = 1usize;
    let mut mprod = 0usize;
    let mut smax = 1.0f64;
    let mut smin = f64::NAN;
    let mut lastsv: Option<Vec<f64>> = None;
    let mut clustered_count = 0usize;
    let mut s: f64;
    let mut invariant_draws = 0usize;
    let mut bsvd: Option<BSvd> = None;
    let mut converged = false;

    while iter <= params.maxit {
        // j and the column indices below are 1-based, as in R
        let mut j = 1usize;
        if iter == 1 {
            let nrm = norm2(&v[..n]);
            for x in &mut v[..n] {
                *x /= nrm;
            }
        } else {
            j = k + 1;
        }
        // W[, j_w] <- drop(mult(A, V[, j]))
        {
            let vj = &v[(j - 1) * n..j * n];
            prod.av(vj, &mut w[(j - 1) * m..j * m]);
        }
        mprod += 1;
        if iter != 1 && w_dim > 1 && reorth {
            // W[, j] <- orthog(W[, j, drop=FALSE], W[, 1:(j - 1), drop=FALSE])
            orthog_against_previous(&mut w, m, j - 1, false);
        }
        s = norm2(&w[(j - 1) * m..j * m]);
        if s.is_nan() || (s < eps2 && j == 1) {
            return err("starting vector near the null space");
        }
        if s.is_nan() || s < eps2 {
            // W[, j_w] <- random(nrow(W)); W[, j] <- orthog(W[, j], W[, 1:(j - 1)])
            w[(j - 1) * m..j * m].copy_from_slice(&random(rng, m));
            invariant_draws += 1;
            if w_dim > 1 {
                orthog_against_previous(&mut w, m, j - 1, true);
            }
            let nrm = norm2(&w[(j - 1) * m..j * m]);
            for x in &mut w[(j - 1) * m..j * m] {
                *x /= nrm;
            }
            s = 0.0;
        } else {
            for x in &mut w[(j - 1) * m..j * m] {
                *x /= s;
            }
        }

        // Lanczos process
        while j <= work {
            let mut invs = 0;
            // F <- t(mult(W[, j_w], A, tx=TRUE))
            prod.atw(&w[(j - 1) * m..j * m], &mut f);
            mprod += 1;
            // F <- drop(F - S * V[, j])
            for (fi, vi) in f.iter_mut().zip(&v[(j - 1) * n..j * n]) {
                *fi -= s * vi;
            }
            // F <- orthog(F, V[, 1:j, drop=FALSE])
            orthog_mat(&mut f, &v[..j * n], n, j);
            if j < work {
                let mut r = norm2(&f);
                if r.is_nan() || r < eps2 {
                    invs += 1;
                    // F <- matrix(random(dim(V)[1]), dim(V)[1], 1); F <- orthog(F, V[, 1:j, drop=FALSE])
                    f = random(rng, n);
                    invariant_draws += 1;
                    orthog_mat(&mut f, &v[..j * n], n, j);
                    let nrm = norm2(&f);
                    for (vi, fi) in v[j * n..(j + 1) * n].iter_mut().zip(&f) {
                        *vi = fi / nrm;
                    }
                    r = 0.0;
                } else {
                    for (vi, fi) in v[j * n..(j + 1) * n].iter_mut().zip(&f) {
                        *vi = fi / r;
                    }
                }
                b = Some(match b.take() {
                    None => Mat { nrow: 1, ncol: 2, x: vec![s, r] },
                    Some(bm) => bm.grow(s, r),
                });
                let w_old = w[(j - 1) * m..j * m].to_vec();
                // W[, jp1_w] <- drop(mult(A, V[, j + 1]))
                prod.av(&v[j * n..(j + 1) * n], &mut w[j * m..(j + 1) * m]);
                mprod += 1;
                // W[, jp1_w] <- W[, jp1_w] - R * w_old
                for (wi, oi) in w[j * m..(j + 1) * m].iter_mut().zip(&w_old) {
                    *wi -= r * oi;
                }
                if reorth && w_dim > 1 {
                    // W[, j + 1] <- orthog(W[, j + 1], W[, 1:j]): a plain vector when j = 1
                    orthog_against_previous(&mut w, m, j, true);
                }
                s = norm2(&w[j * m..(j + 1) * m]);
                if s.is_nan() || s < eps2 {
                    invs += 1;
                    w[j * m..(j + 1) * m].copy_from_slice(&random(rng, m));
                    invariant_draws += 1;
                    if w_dim > 1 {
                        orthog_against_previous(&mut w, m, j, true);
                    }
                    let nrm = norm2(&w[j * m..(j + 1) * m]);
                    for x in &mut w[j * m..(j + 1) * m] {
                        *x /= nrm;
                    }
                    s = 0.0;
                } else {
                    for x in &mut w[j * m..(j + 1) * m] {
                        *x /= s;
                    }
                }
            } else {
                // B <- rbind(B, c(rep(0, j - 1), S))
                let bm = b.take().ok_or_else(|| PcaError::Input("irlba: B is empty".into()))?;
                debug_assert_eq!(bm.ncol, j);
                b = Some(bm.last_row(s));
            }
            if invs > 3 {
                // Near-invariant subspace in Lanczos bidiagonalization (irlba resets invs at
                // each step, so this never happens; ported as written)
                let work_new = (work + 5).min(m.min(n) - 1);
                if work_new > work {
                    v.resize(n * work_new, 0.0);
                    if w_dim > 1 {
                        w.resize(m * work_new, 0.0);
                    }
                    work = work_new;
                }
            }
            j += 1;
        }
        let bm = b.as_ref().ok_or_else(|| PcaError::Input("irlba: B is empty".into()))?;
        let bsz = bm.nrow;
        let r_f = norm2(&f);
        for x in &mut f {
            *x /= r_f;
        }
        // Bsvd <- tryCatch(svd(B), error=function(e) FALSE); if(identical(Bsvd, FALSE)) break
        let sv = match BSvd::of(bm) {
            Ok(sv) => sv,
            Err(e) => {
                // R then stops at ct$converged or Bsvd$d
                return Err(PcaError::Input(format!("irlba: svd(B) failed ({e})")));
            }
        };
        if iter == 1 {
            smax = sv.d[0];
            smin = sv.d[bsz - 1];
        } else {
            smax = smax.max(sv.d[0]);
            smin = smin.min(sv.d[bsz - 1]);
        }
        smax = eps23.max(smax);
        if !reorth && smin / smax < sqrteps {
            // "The matrix is ill-conditioned. Basis will be reorthogonalized."
            reorth = true;
        }
        // R <- R_F * Bsvd$u[Bsz, , drop=FALSE]
        let rr: Vec<f64> = (0..bsz).map(|c| r_f * sv.u[(bsz - 1) + c * bsz]).collect();
        let abs_r: Vec<f64> = rr.iter().map(|x| x.abs()).collect();
        let (ct_converged, ct_k) = convtests(bsz, tol, k_org, &sv.d, &abs_r, k, smax, lastsv.as_deref(), svtol, s);
        lastsv = Some(sv.d.clone());
        if !ct_converged && k_org < bsz {
            // gap_abs <- min(abs(diff(Bsvd$d)))
            let gap_abs = sv.d.windows(2).map(|p| (p[1] - p[0]).abs()).fold(f64::INFINITY, f64::min);
            let denom = sv.d[k_org - 1].abs().max(eps * smax);
            let ext_gap = gap_abs / denom;
            // lastsv is Bsvd$d already, so the change is 0
            let rel_change_k = (sv.d[k_org - 1] - sv.d[k_org - 1]).abs() / denom;
            match and(lt(ext_gap, sqrteps), lt(rel_change_k, svtol)) {
                Some(true) => clustered_count += 1,
                Some(false) => clustered_count = 0,
                None => return err("missing value where TRUE/FALSE needed"),
            }
            if clustered_count > 1 {
                let work_new = (work + 5).min(m.min(n) - 1);
                if work_new > work {
                    v.resize(n * work_new, 0.0);
                    if w_dim > 1 {
                        w.resize(m * work_new, 0.0);
                    }
                    work = work_new;
                    clustered_count = 0;
                }
            }
        }
        k = ct_k;
        converged = ct_converged;
        if ct_converged || iter >= params.maxit {
            bsvd = Some(sv);
            break;
        }
        // V[, 1:(k + 1)] <- cbind(V[, 1:Bsz, drop=FALSE] %*% Bsvd$v[, 1:k], F)
        let vk = blas::matprod(&v[..n * bsz], n, bsz, BSvd::cols(&sv.v, bsz, k), bsz, k);
        v[..n * k].copy_from_slice(&vk);
        v[n * k..n * (k + 1)].copy_from_slice(&f);
        // B <- cbind(diag(Bsvd$d[1:k], nrow=k), R[1:k])
        let mut bx = vec![0.0; k * (k + 1)];
        for i in 0..k {
            bx[i + i * k] = sv.d[i];
            bx[i + k * k] = rr[i];
        }
        b = Some(Mat { nrow: k, ncol: k + 1, x: bx });
        // W[, 1:k] <- W[, 1:Bsz, drop=FALSE] %*% Bsvd$u[, 1:k]
        let wk = blas::matprod(&w[..m * bsz], m, bsz, BSvd::cols(&sv.u, bsz, k), bsz, k);
        w[..m * k].copy_from_slice(&wk);
        bsvd = Some(sv);
        iter += 1;
    }
    let sv = bsvd.ok_or_else(|| PcaError::Input("irlba: no iteration".into()))?;
    let bsz = sv.bsz;
    let d = sv.d[..k_org].to_vec();
    // U <- W[, 1:Bsz, drop=FALSE] %*% Bsvd$u[, 1:k_org, drop=FALSE]; V likewise
    let u = blas::matprod(&w[..m * bsz], m, bsz, BSvd::cols(&sv.u, bsz, k_org), bsz, k_org);
    let vv = blas::matprod(&v[..n * bsz], n, bsz, BSvd::cols(&sv.v, bsz, k_org), bsz, k_org);
    Ok(Irlba { d, u: u[..m * nu].to_vec(), v: vv[..n * nv].to_vec(), iter, mprod, converged, work, invariant_draws, tiny: false })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convtests_like_r() {
        let d = [3.0, 2.0, 1.0, 0.5];
        // first iteration: lastsv NULL, nothing converges; k is capped at max(1, Bsz - 3)
        assert_eq!(convtests(26, 1e-7, 2, &d, &[0.0, 0.0, 0.0, 0.0], 2, 3.0, None, 1e-7, 1.0), (false, 2));
        assert_eq!(convtests(4, 1e-7, 2, &d, &[0.0, 0.0, 0.0, 0.0], 2, 3.0, None, 1e-7, 1.0), (false, 1));
        // both criteria met for the first value: k grows by one, capped at max(1, Bsz - 3)
        let r = [0.0, 1.0, 0.0, 0.0];
        assert_eq!(convtests(26, 1e-7, 10, &[1.0; 26], &[0.0; 26], 10, 1.0, Some(&[1.0; 26]), 1e-7, 1.0), (true, 10));
        assert_eq!(convtests(26, 1e-7, 2, &d, &r, 2, 3.0, Some(&d), 1e-7, 1.0), (false, 3));
        assert_eq!(convtests(4, 1e-7, 2, &d, &r, 2, 3.0, Some(&d), 1e-7, 1.0), (false, 1));
        // S == 0 converges
        assert_eq!(convtests(4, 1e-7, 2, &d, &[1.0; 4], 2, 3.0, None, 1e-7, 0.0), (true, 2));
        // NA (0/0) anywhere makes len_res 0
        let dz = [0.0, 0.0, 0.0];
        assert_eq!(convtests(26, 1e-7, 2, &dz, &[0.0; 3], 5, 3.0, Some(&dz), 1e-7, 1.0), (false, 5));
    }

    #[test]
    fn b_grows_like_r() {
        let b = Mat { nrow: 1, ncol: 2, x: vec![5.0, 6.0] };
        let b = b.grow(7.0, 8.0);
        assert_eq!((b.nrow, b.ncol), (2, 3));
        assert_eq!(b.x, vec![5.0, 0.0, 6.0, 7.0, 0.0, 8.0]);
        let b = b.last_row(9.0);
        assert_eq!((b.nrow, b.ncol), (3, 3));
        assert_eq!(b.x, vec![5.0, 0.0, 0.0, 6.0, 7.0, 0.0, 0.0, 8.0, 9.0]);
    }
}
