//! irlba's two large products, `A %*% v` and `t(w) %*% A`, as R makes them, and the only way to
//! split them across threads that keeps R's bits.
//!
//! irlba 2.4.1's `mult()` calls `direct_dgemm_c` (its `src/irlb.c`) for a dense matrix, which
//! calls BLAS `dgemm` directly, with no NaN check:
//!
//! - `A %*% v`: `dgemm("N", "N", m, 1, n, 1, A, m, v, n, 0, w, m)`;
//! - `t(w) %*% A`: `dgemm("T", "N", 1, n, m, 1, w, m, A, m, 0, f, 1)`.
//!
//! OpenBLAS 0.3.34 does not forward these one-row or one-column calls to `gemv` on x86-64
//! (`GEMM_GEMV_FORWARD` is set for ARM64, POWER and RISC-V only). The SkylakeX kernels (also used
//! for Cooperlake and SapphireRapids) send a call with `M N K <= 10^6` (and, for `"T", "N"`,
//! `M N <= 1200` and `K >= 32`) to small-matrix kernels (`dgemm_small_kernel_permit_skylakex.c`);
//! every other call, and every Haswell call, goes through the level-3 driver
//! (`driver/level3/level3.c`). That driver blocks K by `GEMM_Q` and `GEMM_UNROLL_M` only, so each
//! element of the result is the same chain of K blocks however M and N are cut; within a block
//! the micro-kernel's tile decides the summation order:
//!
//! - Haswell (`dgemm_kernel_4x8_haswell.S`): the 4-row tiles of the `N = 1` product keep four
//!   partial sums (`KERNEL4x1`, k modulo 4), the 2- and 1-row tails one. Every `M = 1` tile
//!   (`1x12` ... `1x1`) keeps one sum per column.
//! - SkylakeX (`dgemm_kernel_16x2_skylakex.c`): every `N = 1` tile keeps one sum per row, and so
//!   does every `M = 1` tile the driver reaches. The kernel's 12-column tile (`KERNEL_k1m1n12`)
//!   keeps four (k modulo 4), but the driver only hands it 12 columns in its SkylakeX branch
//!   (`6 * GEMM_UNROLL_N`), which a `DYNAMIC_ARCH` build compiles out: there the driver is built
//!   once for the base target and only `kernel/` per core, so it passes `3 * GEMM_UNROLL_N = 6`
//!   columns at a time. Its column blocks (`GEMM_R`, about 43,000 at run time) are far wider than
//!   any input here.
//!
//! So a row of `A %*% v` keeps its tile if the pieces start at multiples of 4 rows, the last
//! piece holding the call's own tail; the columns of `t(w) %*% A` may be cut anywhere, and are
//! cut at multiples of 12 anyway, in case a build reaches the 12-column tile. On the SkylakeX
//! family the call and every piece must also stay on the level-3 path: a piece small enough for
//! the small-matrix kernels rounds differently. Under those rules each piece, a `dgemm` call with
//! the same K and leading dimensions on a block of `A`, gives exactly the elements the single
//! call gives, whatever the number of pieces; this module's tests check it under both reference
//! kernels, and check that cuts breaking the rules do change bits. On any other kernel the
//! products are not split.

use rcompat::blas;

/// Pieces of `A %*% v` start at multiples of this many rows (Haswell's 4-row tile divides it).
pub const ROW_ALIGN: usize = 16;
/// Pieces of `t(w) %*% A` start at multiples of this many columns (SkylakeX's 12-column tile, a
/// margin: see the module documentation).
pub const COL_ALIGN: usize = 12;
/// Wider products are not split by column: below SkylakeX's column blocking (`GEMM_R`).
pub const MAX_SPLIT_COLS: usize = 16_384;
/// Below this many elements of `A`, a product is not split (for speed only).
pub const MIN_SPLIT_ELEMENTS: usize = 1 << 20;

/// The OpenBLAS kernels whose tiles this module was checked against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// `Haswell` (`dgemm_kernel_4x8_haswell.S`, no small-matrix kernels)
    Haswell,
    /// `SkylakeX`, `Cooperlake`, `SapphireRapids` (`dgemm_kernel_16x2_skylakex.c`)
    SkylakeX,
    /// any other kernel: products are never split
    Other,
}

impl Family {
    pub fn of(core: &str) -> Family {
        match core {
            "Haswell" => Family::Haswell,
            "SkylakeX" | "Cooperlake" | "SapphireRapids" => Family::SkylakeX,
            _ => Family::Other,
        }
    }

    /// Whether OpenBLAS sends this `dgemm` call to a small-matrix kernel instead of the level-3
    /// driver (`GEMM_SMALL_MATRIX_PERMIT` in `interface/gemm.c`).
    pub fn small_kernel(self, transa: bool, m: usize, n: usize, k: usize) -> bool {
        match self {
            // kernel/generic/gemm_small_matrix_permit.c: never
            Family::Haswell | Family::Other => false,
            // kernel/x86_64/dgemm_small_kernel_permit_skylakex.c
            Family::SkylakeX => {
                let mnk = m as f64 * n as f64 * k as f64;
                if mnk > 100.0 * 100.0 * 100.0 {
                    return false;
                }
                !(transa && (m * n > 1200 || k < 32))
            }
        }
    }
}

/// Cuts `0..len` into at most `pieces` ranges whose starts are multiples of `align`, each
/// accepted by `ok`; fewer pieces when needed, one piece (`0..len`) if nothing else works.
fn plan(len: usize, pieces: usize, align: usize, ok: impl Fn(usize) -> bool) -> Vec<(usize, usize)> {
    let mut t = pieces.min(len / align);
    while t > 1 {
        let mut bounds = vec![0];
        for i in 1..t {
            let b = (i * len / t) / align * align;
            if b > *bounds.last().unwrap() {
                bounds.push(b);
            }
        }
        bounds.push(len);
        let chunks: Vec<(usize, usize)> = bounds.windows(2).map(|w| (w[0], w[1])).collect();
        if chunks.len() > 1 && chunks.iter().all(|&(a, b)| b > a && ok(b - a)) {
            return chunks;
        }
        t -= 1;
    }
    vec![(0, len)]
}

/// Pieces of `A %*% v` (row ranges) for `A` of `m x n` on `threads` threads.
pub fn plan_rows(m: usize, n: usize, threads: usize, family: Family, min_elements: usize) -> Vec<(usize, usize)> {
    let one = vec![(0, m)];
    if threads <= 1 || family == Family::Other || m * n < min_elements || family.small_kernel(false, m, 1, n) {
        return one;
    }
    plan(m, threads, ROW_ALIGN, |rows| !family.small_kernel(false, rows, 1, n))
}

/// Pieces of `t(w) %*% A` (column ranges) for `A` of `m x n` on `threads` threads.
pub fn plan_cols(m: usize, n: usize, threads: usize, family: Family, min_elements: usize) -> Vec<(usize, usize)> {
    let one = vec![(0, n)];
    if threads <= 1 || family == Family::Other || m * n < min_elements || n > MAX_SPLIT_COLS || family.small_kernel(true, 1, n, m) {
        return one;
    }
    plan(n, threads, COL_ALIGN, |cols| !family.small_kernel(true, 1, cols, m))
}

/// The matrix `A` (`m x n`, column-major) of irlba's products, and how they are split.
pub struct Products<'a> {
    a: &'a [f64],
    pub m: usize,
    pub n: usize,
    rows: Vec<(usize, usize)>,
    cols: Vec<(usize, usize)>,
}

impl<'a> Products<'a> {
    /// Products of `a` on `threads` threads, split as the OpenBLAS kernel in use allows.
    pub fn new(a: &'a [f64], m: usize, n: usize, threads: usize) -> Products<'a> {
        let family = Family::of(&blas::corename());
        let rows = plan_rows(m, n, threads, family, MIN_SPLIT_ELEMENTS);
        let cols = plan_cols(m, n, threads, family, MIN_SPLIT_ELEMENTS);
        Products::with_pieces(a, m, n, rows, cols)
    }

    /// Products with the given pieces (row ranges of `A %*% v`, column ranges of `t(w) %*% A`).
    pub fn with_pieces(a: &'a [f64], m: usize, n: usize, rows: Vec<(usize, usize)>, cols: Vec<(usize, usize)>) -> Products<'a> {
        assert_eq!(a.len(), m * n, "Products: A has the wrong length");
        for (pieces, len) in [(&rows, m), (&cols, n)] {
            assert!(!pieces.is_empty() && pieces[0].0 == 0 && pieces.last().unwrap().1 == len, "Products: pieces must cover the range");
            assert!(pieces.windows(2).all(|w| w[0].1 == w[1].0) && pieces.iter().all(|&(a, b)| b > a), "Products: pieces must be contiguous");
        }
        Products { a, m, n, rows, cols }
    }

    /// `A` itself.
    pub fn a(&self) -> &'a [f64] {
        self.a
    }

    /// The row pieces of `A %*% v` and the column pieces of `t(w) %*% A`.
    pub fn pieces(&self) -> (&[(usize, usize)], &[(usize, usize)]) {
        (&self.rows, &self.cols)
    }

    /// `out = A %*% v` (`v`: n values, `out`: m), as `mult(A, v)`.
    pub fn av(&self, v: &[f64], out: &mut [f64]) {
        let (m, n, a) = (self.m, self.n, self.a);
        assert!(v.len() == n && out.len() == m, "Products::av: wrong lengths");
        // dgemm("N", "N", rows, 1, n, 1, A[r0, ], m, v, n, 0, out[r0], m) for each piece
        run(&self.rows, out, |(r0, r1), o| blas::dgemm(b'N', b'N', r1 - r0, 1, n, 1.0, &a[r0..], m, v, n, 0.0, o, m));
    }

    /// `out = t(w) %*% A` (`w`: m values, `out`: n), as `t(mult(w, A, tx = TRUE))`.
    pub fn atw(&self, w: &[f64], out: &mut [f64]) {
        let (m, a) = (self.m, self.a);
        assert!(w.len() == m && out.len() == self.n, "Products::atw: wrong lengths");
        // dgemm("T", "N", 1, cols, m, 1, w, m, A[, c0], m, 0, out[c0], 1) for each piece
        run(&self.cols, out, |(c0, c1), o| blas::dgemm(b'T', b'N', 1, c1 - c0, m, 1.0, w, m, &a[c0 * m..], m, 0.0, o, 1));
    }
}

/// Runs `f` on each piece of `out` (contiguous ranges), one thread per piece.
fn run(pieces: &[(usize, usize)], out: &mut [f64], f: impl Fn((usize, usize), &mut [f64]) + Sync) {
    if pieces.len() == 1 {
        f(pieces[0], out);
        return;
    }
    let mut parts = Vec::with_capacity(pieces.len());
    let mut rest = out;
    for &(a, b) in pieces {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(b - a);
        parts.push(((a, b), head));
        rest = tail;
    }
    std::thread::scope(|s| {
        let f = &f;
        let mut parts = parts.into_iter();
        let first = parts.next().unwrap();
        for (piece, o) in parts {
            s.spawn(move || f(piece, o));
        }
        f(first.0, first.1);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cheap, reproducible values of both signs and many magnitudes.
    fn values(len: usize, seed: u64) -> Vec<f64> {
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let u = (x >> 11) as f64 / (1u64 << 53) as f64;
                (u - 0.5) * (1.0 + (x & 0xff) as f64)
            })
            .collect()
    }

    fn bits(x: &[f64]) -> Vec<u64> {
        x.iter().map(|v| v.to_bits()).collect()
    }

    /// Every plan, at every thread count, gives the single call's bits; checked on whichever
    /// OpenBLAS kernel runs the tests (Haswell under env.sh; with `OPENBLAS_CORETYPE` unset on a
    /// Sapphire Rapids node, SapphireRapids), at shapes that cover every tile and tail, and
    /// sizes on both sides of SkylakeX's small-matrix limit.
    #[test]
    fn pieces_give_the_single_call_bits() {
        let family = Family::of(&blas::corename());
        let shapes = [(37, 29), (64, 12), (1003, 517), (4099, 251), (2050, 611), (333, 3000), (4000, 1201), (5003, 777)];
        let mut split_rows = 0;
        let mut split_cols = 0;
        for (s, &(m, n)) in shapes.iter().enumerate() {
            let a = values(m * n, 1 + s as u64);
            let v = values(n, 100 + s as u64);
            let w = values(m, 200 + s as u64);
            let one = Products::with_pieces(&a, m, n, vec![(0, m)], vec![(0, n)]);
            let (mut av1, mut atw1) = (vec![0.0; m], vec![0.0; n]);
            one.av(&v, &mut av1);
            one.atw(&w, &mut atw1);
            for t in 2..=13 {
                // the plans for this kernel, with no minimum size, so that small shapes split too
                let rows = plan_rows(m, n, t, family, 0);
                let cols = plan_cols(m, n, t, family, 0);
                split_rows += (rows.len() > 1) as usize;
                split_cols += (cols.len() > 1) as usize;
                let p = Products::with_pieces(&a, m, n, rows.clone(), cols.clone());
                let (mut av, mut atw) = (vec![f64::NAN; m], vec![f64::NAN; n]);
                p.av(&v, &mut av);
                p.atw(&w, &mut atw);
                assert_eq!(bits(&av), bits(&av1), "A v, {m} x {n}, pieces {rows:?}");
                assert_eq!(bits(&atw), bits(&atw1), "t(w) A, {m} x {n}, pieces {cols:?}");
            }
        }
        if family != Family::Other {
            assert!(split_rows > 20 && split_cols > 20, "too few split plans tested: {split_rows}, {split_cols}");
        }
    }

    /// The rules are needed, and no more than they say: on Haswell a cut that moves rows of
    /// `A %*% v` out of a 4-row tile changes their bits; on the SkylakeX family a piece that
    /// takes the small-matrix kernel changes bits, and any other cut changes none.
    #[test]
    fn rule_breaking_pieces_change_bits() {
        let family = Family::of(&blas::corename());
        eprintln!("OpenBLAS core {} ({family:?})", blas::corename());
        let (m, n) = (4000, 1201);
        let a = values(m * n, 7);
        let (v, w) = (values(n, 8), values(m, 9));
        let one = Products::with_pieces(&a, m, n, vec![(0, m)], vec![(0, n)]);
        let (mut av, mut atw) = (vec![0.0; m], vec![0.0; n]);
        one.av(&v, &mut av);
        one.atw(&w, &mut atw);
        let differ = |x: &[f64], y: &[f64]| x.iter().zip(y).filter(|(a, b)| a.to_bits() != b.to_bits()).count();
        let rows_cut = |cut: usize| {
            let mut c = vec![0.0; m];
            Products::with_pieces(&a, m, n, vec![(0, cut), (cut, m)], vec![(0, n)]).av(&v, &mut c);
            differ(&av, &c)
        };
        let cols_cut = |cut: usize| {
            let mut c = vec![0.0; n];
            Products::with_pieces(&a, m, n, vec![(0, m)], vec![(0, cut), (cut, n)]).atw(&w, &mut c);
            differ(&atw, &c)
        };
        match family {
            Family::Haswell => {
                // rows 1000-1001 and the last two move from 4-row tiles to 2-row tails
                assert_eq!(rows_cut(1002), 4);
                assert_eq!(rows_cut(1008), 0);
                // every one-row tile keeps one sum per column
                assert_eq!(cols_cut(606) + cols_cut(1196) + cols_cut(240), 0);
            }
            Family::SkylakeX => {
                // level-3 path on both sides: nothing changes, aligned or not
                assert_eq!(rows_cut(1002) + cols_cut(606) + cols_cut(607), 0);
                // a piece small enough for the small-matrix kernel changes
                assert!(!family.small_kernel(false, 1002, 1, n) && family.small_kernel(false, 200, 1, n));
                assert!(rows_cut(200) > 0 && rows_cut(3800) > 0);
                assert!(family.small_kernel(true, 1, 240, m) && family.small_kernel(true, 1, 5, m));
                assert!(cols_cut(240) > 0 && cols_cut(1196) > 0);
            }
            Family::Other => {}
        }
    }

    #[test]
    fn plans_follow_the_rules() {
        // Haswell: any size splits, at multiples of 16 rows and 12 columns
        let r = plan_rows(1000, 50, 4, Family::Haswell, 0);
        assert_eq!(r, vec![(0, 240), (240, 496), (496, 736), (736, 1000)]);
        let c = plan_cols(1000, 50, 4, Family::Haswell, 0);
        assert_eq!(c, vec![(0, 12), (12, 24), (24, 36), (36, 50)]);
        // SkylakeX: no piece may take the small-matrix kernel
        assert_eq!(plan_rows(1000, 50, 4, Family::SkylakeX, 0), vec![(0, 1000)]);
        let r = plan_rows(4000, 3000, 32, Family::SkylakeX, 0);
        assert!(r.len() > 1 && r.iter().all(|&(a, b)| (b - a) * 3000 > 1_000_000 && a % ROW_ALIGN == 0));
        let c = plan_cols(4144, 3000, 32, Family::SkylakeX, 0);
        assert!(c.len() > 1 && c.iter().all(|&(a, b)| (b - a) * 4144 > 1_000_000 && a % COL_ALIGN == 0));
        assert_eq!(plan_cols(400, 850, 8, Family::SkylakeX, 0), vec![(0, 850)]);
        // other kernels, one thread or small products: one call
        assert_eq!(plan_rows(100_000, 3000, 32, Family::Other, 0), vec![(0, 100_000)]);
        assert_eq!(plan_cols(100_000, 3000, 1, Family::Haswell, 0), vec![(0, 3000)]);
        assert_eq!(plan_rows(400, 850, 8, Family::Haswell, MIN_SPLIT_ELEMENTS), vec![(0, 400)]);
        assert!(Family::SkylakeX.small_kernel(true, 1, 850, 400));
        assert!(!Family::SkylakeX.small_kernel(true, 1, 1201, 400));
        assert!(!Family::SkylakeX.small_kernel(true, 1, 850, 31));
        assert!(Family::SkylakeX.small_kernel(false, 400, 1, 850));
        assert!(!Family::Haswell.small_kernel(false, 4, 1, 4));
    }
}
