//! Exact nearest neighbours, memory-bounded and parallel over cells (docs/PLAN.md B3;
//! milestone M3).
//!
//! DoubletFinder ranks neighbours with `order(fields::rdist(pca)[, i])`: the Euclidean distance
//! from cell i to every cell, then R's stable `order()`, so ties go to the lower index. This
//! crate gives the same ranking without the N x N matrix:
//!
//! - **Distances** use the arithmetic of `fields::rdist` (fields 17.3, Fortran `rdist`): for a
//!   pair of cells, `s = (x[1] - y[1])^2`, then `s = (x[c] - y[c])^2 + s` for c = 2..dim in
//!   column order, then `sqrt(s)`, with `^2` as a multiplication and no fused multiply-add. The
//!   distances are compared after the square root, as R compares them (distinct sums can share
//!   a square root).
//! - **Selection**: per query cell, the distances to all N cells go into an O(N) buffer per
//!   thread, and the first L positions are selected by (distance, index), which is exactly
//!   R's `order()` followed by truncation to L. Position 1 is not assumed to be the cell itself:
//!   a cell with a bitwise-identical twin of lower index comes second.
//! - **Parallelism** is over query cells; each result depends only on its cell, so the output is
//!   the same for any number of threads.
//!
//! Memory is the points (N x dim doubles), two O(N) buffers per thread and whatever the caller
//! keeps per query.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Points in R's matrix layout: `n` rows (cells) by `dim` columns (coordinates), column-major,
/// so coordinate c of point i is `data[c * n + i]` (as `cell.embeddings` is stored).
#[derive(Debug, Clone, Copy)]
pub struct Points<'a> {
    n: usize,
    dim: usize,
    data: &'a [f64],
}

impl<'a> Points<'a> {
    pub fn new(data: &'a [f64], n: usize, dim: usize) -> Points<'a> {
        assert_eq!(data.len(), n * dim, "points: {} values for {n} x {dim}", data.len());
        assert!(dim >= 1 && n >= 1);
        assert!(n <= u32::MAX as usize);
        Points { n, dim, data }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    fn col(&self, c: usize) -> &'a [f64] {
        &self.data[c * self.n..(c + 1) * self.n]
    }

    /// Coordinate c of point i.
    pub fn get(&self, i: usize, c: usize) -> f64 {
        self.data[c * self.n + i]
    }
}

/// How the first L positions are handed to the caller.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    /// In R's `order()` order: by distance, ties to the lower index.
    Sorted,
    /// Position 1 (the smallest by (distance, index)) first, then the other L - 1 positions in
    /// no particular order. Enough to count neighbours at positions 2..L, without a sort.
    FirstThenSet,
}

/// Distances from point `q` to every point, as `fields::rdist` computes them.
#[allow(clippy::assign_op_pattern)] // rdist's `k = (x - y)**2 + k`, kept as written
pub fn distances(points: &Points, q: usize, out: &mut [f64]) {
    let n = points.n;
    debug_assert_eq!(out.len(), n);
    // blocks keep the running sums in L1; the inner loops vectorize without changing any
    // element's operations (no reduction across elements, no FMA)
    const BLOCK: usize = 256;
    let mut s = [0.0f64; BLOCK];
    let mut start = 0;
    while start < n {
        let end = (start + BLOCK).min(n);
        let len = end - start;
        let s = &mut s[..len];
        let x0 = &points.col(0)[start..end];
        let q0 = points.get(q, 0);
        for (si, &xi) in s.iter_mut().zip(x0) {
            let d = xi - q0;
            *si = d * d;
        }
        for c in 1..points.dim {
            let xc = &points.col(c)[start..end];
            let qc = points.get(q, c);
            for (si, &xi) in s.iter_mut().zip(xc) {
                let d = xi - qc;
                *si = d * d + *si;
            }
        }
        for (o, &si) in out[start..end].iter_mut().zip(s.iter()) {
            *o = si.sqrt();
        }
        start = end;
    }
}

/// A distance as an unsigned key in R's `order()` order: non-negative doubles order like their
/// bit patterns; NaN (which R puts last) is largest.
#[inline]
fn key(d: f64) -> u64 {
    if d.is_nan() {
        u64::MAX
    } else if d == 0.0 {
        0
    } else {
        d.to_bits()
    }
}

/// Per-thread buffers.
struct Scratch {
    dist: Vec<f64>,
    keys: Vec<u64>,
    sel: Vec<u32>,
}

/// The first `l` positions (0-based indices) of `order(d)` for the distances in `scratch.dist`,
/// in the requested `order`, left in `scratch.sel`.
fn select(scratch: &mut Scratch, l: usize, order: Order) {
    let Scratch { dist, keys, sel } = scratch;
    let n = dist.len();
    sel.clear();
    if l == 0 {
        return;
    }
    if l >= n {
        sel.extend(0..n as u32);
    } else {
        // threshold: the l-th smallest distance key
        keys.clear();
        keys.extend(dist.iter().map(|&d| key(d)));
        let (_, &mut t, _) = keys.select_nth_unstable(l - 1);
        // every index below the threshold, then threshold ties in index order
        let mut ties = 0usize;
        for (i, &d) in dist.iter().enumerate() {
            let k = key(d);
            if k < t {
                sel.push(i as u32);
            } else if k == t {
                ties += 1;
            }
        }
        let need = l - sel.len();
        debug_assert!(need >= 1 && need <= ties);
        let mut taken = 0;
        for (i, &d) in dist.iter().enumerate() {
            if taken == need {
                break;
            }
            if key(d) == t {
                sel.push(i as u32);
                taken += 1;
            }
        }
    }
    match order {
        Order::Sorted => sel.sort_unstable_by_key(|&i| (key(dist[i as usize]), i)),
        Order::FirstThenSet => {
            let first = (0..sel.len()).min_by_key(|&j| (key(dist[sel[j] as usize]), sel[j])).unwrap();
            sel.swap(0, first);
        }
    }
}

/// For each query point `q` in `0..n_query`, the first `l` positions of R's
/// `order(rdist(points)[, q + 1])` (0-based indices of points), passed to `f(q, positions)`;
/// returns `f`'s results in query order. `l` is capped at the number of points (R's positions
/// beyond N are NA). Uses `threads` worker threads (at least one); the result does not depend
/// on it.
pub fn neighbours<R, F>(points: &Points, n_query: usize, l: usize, order: Order, threads: usize, f: F) -> Vec<R>
where
    R: Send,
    F: Fn(usize, &[u32]) -> R + Sync,
{
    assert!(n_query <= points.n);
    let l = l.min(points.n);
    let threads = threads.max(1).min(n_query.max(1));
    let chunk = 16usize;
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<(usize, Vec<R>)>> = Mutex::new(Vec::new());
    let work = || {
        let mut scratch = Scratch { dist: vec![0.0; points.n], keys: Vec::with_capacity(points.n), sel: Vec::with_capacity(l) };
        loop {
            let start = next.fetch_add(chunk, Ordering::Relaxed);
            if start >= n_query {
                break;
            }
            let end = (start + chunk).min(n_query);
            let mut out = Vec::with_capacity(end - start);
            for q in start..end {
                distances(points, q, &mut scratch.dist);
                select(&mut scratch, l, order);
                out.push(f(q, &scratch.sel));
            }
            results.lock().unwrap().push((start, out));
        }
    };
    if threads == 1 {
        work();
    } else {
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(work);
            }
        });
    }
    let mut parts = results.into_inner().unwrap();
    parts.sort_unstable_by_key(|(start, _)| *start);
    let mut all = Vec::with_capacity(n_query);
    for (_, part) in parts {
        all.extend(part);
    }
    debug_assert_eq!(all.len(), n_query);
    all
}

/// The digest of one cell's order positions that DoubletFinder's instrumented R runs dump as
/// stage S4 (`dfi_digest_matrix` in `validation/R/lib/doubletfinder_hooks.R`), for 1-based
/// positions `pos[1..L]`: `h1 = sum_r r * pos[r]` and `h2 = (sum_r (pos[r] * w[r] mod P)) mod
/// P` with `P = 2^31 - 1` and `w[r] = 16807^r mod P`. Exact in integers (and in R's doubles).
#[derive(Debug, Clone)]
pub struct Digest {
    weights: Vec<u64>,
}

/// `2^31 - 1`.
pub const DIGEST_P: u64 = 2_147_483_647;

impl Digest {
    pub fn new(l: usize) -> Digest {
        let mut weights = Vec::with_capacity(l);
        let mut v = 1u64;
        for _ in 0..l {
            v = (v * 16807) % DIGEST_P;
            weights.push(v);
        }
        Digest { weights }
    }

    /// `(h1, h2)` of 0-based positions (converted to R's 1-based ones), at most L of them.
    pub fn of(&self, positions0: &[u32]) -> (u64, u64) {
        let (mut h1, mut h2) = (0u64, 0u64);
        for (r, &p0) in positions0.iter().enumerate() {
            let pos = p0 as u64 + 1;
            h1 += (r as u64 + 1) * pos;
            h2 += (pos * self.weights[r]) % DIGEST_P;
        }
        (h1, h2 % DIGEST_P)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// R's `order(d)[1:l]` by a stable sort of the full distance vector.
    fn brute(points: &Points, q: usize, l: usize) -> Vec<u32> {
        let mut d = vec![0.0; points.n()];
        distances(points, q, &mut d);
        let o = rcompat::order(&d, false);
        o[..l.min(o.len())].iter().map(|&i| i as u32).collect()
    }

    /// Points with many exact ties: coordinates on a coarse grid, duplicated cells, and a
    /// cluster of identical points.
    fn tied_points(n: usize, dim: usize, seed: u64) -> Vec<f64> {
        let mut g = Lcg(seed);
        let mut rows: Vec<Vec<f64>> = (0..n).map(|_| (0..dim).map(|_| (g.next() * 4.0).floor() - 1.5).collect()).collect();
        for i in (5..n).step_by(7) {
            rows[i] = rows[i / 3].clone(); // bitwise twins, some of lower index
        }
        let mut data = vec![0.0; n * dim];
        for (i, r) in rows.iter().enumerate() {
            for (c, &v) in r.iter().enumerate() {
                data[c * n + i] = v;
            }
        }
        data
    }

    #[test]
    fn rdist_arithmetic() {
        // 3 points in 2 dimensions; d(0, 1) = sqrt(3^2 + 4^2)
        let data = [0.0, 3.0, 1.0, 0.0, 4.0, 1.0];
        let p = Points::new(&data, 3, 2);
        let mut d = [0.0; 3];
        distances(&p, 0, &mut d);
        assert_eq!(d, [0.0, 5.0, 2f64.sqrt()]);
    }

    #[test]
    fn selection_equals_rs_order_with_ties() {
        for (n, dim, l) in [(60, 2, 10), (200, 3, 37), (500, 10, 120), (64, 1, 64), (50, 2, 80)] {
            let data = tied_points(n, dim, n as u64);
            let p = Points::new(&data, n, dim);
            let nq = n / 2;
            let got = neighbours(&p, nq, l, Order::Sorted, 1, |_, pos| pos.to_vec());
            for (q, g) in got.iter().enumerate() {
                assert_eq!(g, &brute(&p, q, l), "n={n} q={q}");
            }
            let sets = neighbours(&p, nq, l, Order::FirstThenSet, 3, |_, pos| pos.to_vec());
            for (q, s) in sets.iter().enumerate() {
                let want = brute(&p, q, l);
                assert_eq!(s[0], want[0], "first position, n={n} q={q}");
                let (mut a, mut b) = (s.clone(), want.clone());
                a.sort_unstable();
                b.sort_unstable();
                assert_eq!(a, b, "set, n={n} q={q}");
            }
        }
    }

    #[test]
    fn a_lower_twin_comes_first() {
        // points 0 and 2 identical: for query 2, position 1 is point 0
        let data = [1.0, 5.0, 1.0, 2.0, 6.0, 2.0];
        let p = Points::new(&data, 3, 2);
        let got = neighbours(&p, 3, 3, Order::Sorted, 1, |_, pos| pos.to_vec());
        assert_eq!(got[2], vec![0, 2, 1]);
        assert_eq!(got[0], vec![0, 2, 1]);
    }

    #[test]
    fn independent_of_thread_count() {
        let n = 3000;
        let data = tied_points(n, 10, 7);
        let p = Points::new(&data, n, 10);
        let dg = Digest::new(400);
        let one = neighbours(&p, 2000, 400, Order::Sorted, 1, |_, pos| dg.of(pos));
        for t in [2, 3, 8, 13] {
            assert_eq!(neighbours(&p, 2000, 400, Order::Sorted, t, |_, pos| dg.of(pos)), one, "{t} threads");
        }
    }

    #[test]
    fn digest_definition() {
        let d = Digest::new(3);
        // positions 3, 1, 2 (1-based): h1 = 1*3 + 2*1 + 3*2; w = 16807, 16807^2 mod P, ...
        let w1 = 16807u64;
        let w2 = w1 * w1 % DIGEST_P;
        let w3 = w2 * w1 % DIGEST_P;
        let h2 = ((3 * w1) % DIGEST_P + (w2) % DIGEST_P + (2 * w3) % DIGEST_P) % DIGEST_P;
        assert_eq!(d.of(&[2, 0, 1]), (11, h2));
    }
}
