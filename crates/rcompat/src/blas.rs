//! The OpenBLAS routines R and its packages call, through FFI (docs/PLAN.md A5, "BLAS and LAPACK").
//!
//! OpenBLAS chooses its kernel per CPU at run time (`OPENBLAS_CORETYPE` pins it, `env.sh`), and
//! kernels round differently, so a port reproduces R bit for bit only by calling the same library
//! with the same routine, the same arguments and the same matrix shapes. The bindings are the
//! Fortran entry points (`dgemm_` etc.) with 32-bit integers, as R's `F77_CALL`s pass them, and the
//! hidden string-length arguments R adds (`FCONE`).
//!
//! OpenBLAS must run single-threaded (`OPENBLAS_NUM_THREADS=1`, set by `env.sh`); several Rust
//! threads may call it at once.

use std::os::raw::{c_char, c_int};

mod sys {
    use std::os::raw::{c_char, c_int};
    extern "C" {
        pub fn dgemm_(
            transa: *const c_char,
            transb: *const c_char,
            m: *const c_int,
            n: *const c_int,
            k: *const c_int,
            alpha: *const f64,
            a: *const f64,
            lda: *const c_int,
            b: *const f64,
            ldb: *const c_int,
            beta: *const f64,
            c: *mut f64,
            ldc: *const c_int,
            la: usize,
            lb: usize,
        );
        pub fn dgemv_(
            trans: *const c_char,
            m: *const c_int,
            n: *const c_int,
            alpha: *const f64,
            a: *const f64,
            lda: *const c_int,
            x: *const f64,
            incx: *const c_int,
            beta: *const f64,
            y: *mut f64,
            incy: *const c_int,
            lt: usize,
        );
        pub fn ddot_(n: *const c_int, x: *const f64, incx: *const c_int, y: *const f64, incy: *const c_int) -> f64;
        pub fn dnrm2_(n: *const c_int, x: *const f64, incx: *const c_int) -> f64;
    }
}

fn int(n: usize, what: &str) -> c_int {
    c_int::try_from(n).unwrap_or_else(|_| panic!("BLAS: {what} = {n} exceeds the 32-bit range"))
}

/// Number of elements a column-major `rows x cols` operand with leading dimension `ld` spans.
fn span(rows: usize, cols: usize, ld: usize) -> usize {
    if rows == 0 || cols == 0 {
        0
    } else {
        ld * (cols - 1) + rows
    }
}

/// `dgemm`: `C := alpha * op(A) * op(B) + beta * C`, column-major, `op(X) = X` for `b'N'` and
/// `X^T` for `b'T'`. `C` is `m x n`, `op(A)` is `m x k`.
#[allow(clippy::too_many_arguments)]
pub fn dgemm(
    transa: u8,
    transb: u8,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) {
    let (ar, ac) = if transa == b'N' { (m, k) } else { (k, m) };
    let (br, bc) = if transb == b'N' { (k, n) } else { (n, k) };
    assert!(transa == b'N' || transa == b'T', "dgemm: transa must be N or T");
    assert!(transb == b'N' || transb == b'T', "dgemm: transb must be N or T");
    assert!(lda >= ar.max(1) && ldb >= br.max(1) && ldc >= m.max(1), "dgemm: leading dimension too small");
    assert!(a.len() >= span(ar, ac, lda) && b.len() >= span(br, bc, ldb) && c.len() >= span(m, n, ldc), "dgemm: operand too short");
    let (ta, tb) = (transa as c_char, transb as c_char);
    let (m, n, k) = (int(m, "m"), int(n, "n"), int(k, "k"));
    let (lda, ldb, ldc) = (int(lda, "lda"), int(ldb, "ldb"), int(ldc, "ldc"));
    // SAFETY: the shapes and leading dimensions were checked against the slices above.
    unsafe {
        sys::dgemm_(&ta, &tb, &m, &n, &k, &alpha, a.as_ptr(), &lda, b.as_ptr(), &ldb, &beta, c.as_mut_ptr(), &ldc, 1, 1);
    }
}

/// `dgemv`: `y := alpha * op(A) * x + beta * y`, `A` is `m x n` column-major.
#[allow(clippy::too_many_arguments)]
pub fn dgemv(trans: u8, m: usize, n: usize, alpha: f64, a: &[f64], lda: usize, x: &[f64], incx: usize, beta: f64, y: &mut [f64], incy: usize) {
    assert!(trans == b'N' || trans == b'T', "dgemv: trans must be N or T");
    let (lx, ly) = if trans == b'N' { (n, m) } else { (m, n) };
    assert!(lda >= m.max(1) && incx >= 1 && incy >= 1, "dgemv: bad leading dimension or increment");
    assert!(a.len() >= span(m, n, lda), "dgemv: A too short");
    assert!(lx == 0 || x.len() > (lx - 1) * incx, "dgemv: x too short");
    assert!(ly == 0 || y.len() > (ly - 1) * incy, "dgemv: y too short");
    let t = trans as c_char;
    let (m, n, lda, incx, incy) = (int(m, "m"), int(n, "n"), int(lda, "lda"), int(incx, "incx"), int(incy, "incy"));
    // SAFETY: checked above.
    unsafe { sys::dgemv_(&t, &m, &n, &alpha, a.as_ptr(), &lda, x.as_ptr(), &incx, &beta, y.as_mut_ptr(), &incy, 1) }
}

/// `ddot(n, x, 1, y, 1)`.
pub fn ddot(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "ddot: lengths differ");
    let (n, one) = (int(x.len(), "n"), 1 as c_int);
    // SAFETY: both slices hold n elements.
    unsafe { sys::ddot_(&n, x.as_ptr(), &one, y.as_ptr(), &one) }
}

/// `dnrm2(n, x, 1)`.
pub fn dnrm2(x: &[f64]) -> f64 {
    let (n, one) = (int(x.len(), "n"), 1 as c_int);
    // SAFETY: the slice holds n elements.
    unsafe { sys::dnrm2_(&n, x.as_ptr(), &one) }
}

/// R's imprecise test for NaN or infinite values before a matrix product (`mayHaveNaNOrInf` in
/// `src/main/array.c`): an odd first element alone, then the sums of pairs.
pub fn may_have_nan_or_inf(x: &[f64]) -> bool {
    let n = x.len();
    if n & 1 != 0 && !x[0].is_finite() {
        return true;
    }
    let mut i = n & 1;
    while i < n {
        if !(x[i] + x[i + 1]).is_finite() {
            return true;
        }
        i += 2;
    }
    false
}

/// `tcrossprod(x, y)` = `x %*% t(y)` as R 4.5.3 computes it with its default `matprod` option
/// (`tcrossprod` in `src/main/array.c`): `x` is `nrx x ncx` and `y` is `nry x ncx`, column-major;
/// the result is `nrx x nry`. If either operand may hold NaN or Inf, R's own triple loop
/// (`simple_tcrossprod`); otherwise `dgemv` when one of them has a single row, else
/// `dgemm("N", "T")` with exactly R's arguments.
pub fn tcrossprod(x: &[f64], nrx: usize, ncx: usize, y: &[f64], nry: usize, ncy: usize) -> Vec<f64> {
    let mut z = vec![0.0; nrx * nry];
    tcrossprod_into(x, nrx, ncx, y, nry, ncy, &mut z);
    z
}

/// [`tcrossprod`] into `z` (`nrx * nry` elements), so that a caller can reuse one buffer for many
/// products. What `z` held before does not matter: R's result is uninitialised memory too, which
/// `dgemm` with `beta = 0` never reads (BLAS: "when BETA is supplied as zero then C need not be
/// set on input"); `z` is cleared first for the other routes.
pub fn tcrossprod_into(x: &[f64], nrx: usize, ncx: usize, y: &[f64], nry: usize, ncy: usize, z: &mut [f64]) {
    assert_eq!(ncx, ncy, "tcrossprod: non-conformable arguments");
    assert_eq!(x.len(), nrx * ncx, "tcrossprod: x has the wrong length");
    assert_eq!(y.len(), nry * ncy, "tcrossprod: y has the wrong length");
    assert_eq!(z.len(), nrx * nry, "tcrossprod: z has the wrong length");
    if nrx == 0 || ncx == 0 || nry == 0 || ncy == 0 {
        // zero-extent operations should return zeroes
        z.fill(0.0);
        return;
    }
    if may_have_nan_or_inf(x) || may_have_nan_or_inf(y) {
        for i in 0..nrx {
            for k in 0..nry {
                let mut sum = 0.0;
                for j in 0..ncx {
                    sum += x[i + j * nrx] * y[k + j * nry];
                }
                z[i + k * nrx] = sum;
            }
        }
        return;
    }
    if nry == 1 {
        z.fill(0.0);
        dgemv(b'N', nrx, ncx, 1.0, x, nrx, y, 1, 0.0, z, 1);
    } else if nrx == 1 {
        z.fill(0.0);
        dgemv(b'N', nry, ncy, 1.0, y, nry, x, 1, 0.0, z, 1);
    } else {
        dgemm(b'N', b'T', nrx, nry, ncx, 1.0, x, nrx, y, nry, 0.0, z, nrx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_products_are_exact() {
        // 2 x 2 times its transpose: integers, exact in any kernel
        let x = [1.0, 3.0, 2.0, 4.0]; // [[1, 2], [3, 4]]
        let z = tcrossprod(&x, 2, 2, &x, 2, 2);
        assert_eq!(z, vec![5.0, 11.0, 11.0, 25.0]);
        assert_eq!(ddot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]), 32.0);
        assert_eq!(dnrm2(&[3.0, 4.0]), 5.0);
        let mut y = [0.0];
        dgemv(b'T', 3, 1, 1.0, &[0.5, 0.25, 0.125], 3, &[1.0, 1.0, 1.0], 1, 0.0, &mut y, 1);
        assert_eq!(y[0], 0.875);
    }

    #[test]
    fn products_into_a_used_buffer_ignore_its_contents() {
        // dgemm (beta = 0) must not read C: a NaN left in the buffer would propagate
        let (nrx, nry, k) = (37, 53, 2);
        let x: Vec<f64> = (0..nrx * k).map(|i| (i as f64 * 0.37).sin() - 1.5).collect();
        let y: Vec<f64> = (0..nry * k).map(|i| if i < nry { 1.0 } else { 2.0 + (i as f64 * 0.11).cos() }).collect();
        let fresh = tcrossprod(&x, nrx, k, &y, nry, k);
        let mut z = vec![f64::NAN; nrx * nry];
        tcrossprod_into(&x, nrx, k, &y, nry, k, &mut z);
        assert!(z.iter().zip(&fresh).all(|(a, b)| a.to_bits() == b.to_bits()));
        // the dgemv routes (one row on either side) and the zero-extent case
        let fresh = tcrossprod(&x[..k], 1, k, &y, nry, k);
        let mut z = vec![f64::NAN; nry];
        tcrossprod_into(&x[..k], 1, k, &y, nry, k, &mut z);
        assert!(z.iter().zip(&fresh).all(|(a, b)| a.to_bits() == b.to_bits()));
        let fresh = tcrossprod(&x, nrx, k, &y[..k], 1, k);
        let mut z = vec![f64::NAN; nrx];
        tcrossprod_into(&x, nrx, k, &y[..k], 1, k, &mut z);
        assert!(z.iter().zip(&fresh).all(|(a, b)| a.to_bits() == b.to_bits()));
        let mut z = vec![f64::NAN; nrx * nry];
        tcrossprod_into(&[], nrx, 0, &[], nry, 0, &mut z);
        assert!(z.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn nan_falls_back_to_r_loop() {
        let x = [f64::NAN, 1.0];
        let z = tcrossprod(&x, 2, 1, &[2.0, 3.0], 2, 1);
        assert!(z[0].is_nan() && z[2].is_nan());
        assert_eq!((z[1], z[3]), (2.0, 3.0));
        assert!(may_have_nan_or_inf(&[f64::INFINITY]));
        assert!(!may_have_nan_or_inf(&[1.0, 2.0, 3.0]));
    }
}
