//! The exact search `bluster::neighborPurity` makes (docs/PLAN.md C3e; milestone M12), with
//! knncolle's arithmetic, and a certificate that KMKNN, the algorithm BiocNeighbors uses, finds
//! the same.
//!
//! *The arithmetic is ported from knncolle's `EuclideanDistance` (Aaron Lun, MIT; the version
//! vendored by assorthead 1.4.0), and the k-th distance and the range query from BiocNeighbors
//! 2.4.0's `generic_find_knn` and `generic_find_all` (Aaron Lun, GPL-3). NOTICE lists what was
//! ported and the changes.*
//!
//! The workflow calls `neighborPurity(emb[, 1:30], clusters, k = 50, BNPARAM = NULL)`, which
//! builds BiocNeighbors' default index, `KmknnParam()`, and makes two searches:
//!
//! 1. `findDistance(idx, k)`: each cell's distance to its k-th nearest other cell. knncolle
//!    fills a queue of the k + 1 nearest cells (the cell itself among them), drops the cell
//!    itself (or, when more than k + 1 duplicates sit at distance 0, the farthest), and returns
//!    the last. Whatever the ties, that is the (k + 1)-th smallest raw distance to all cells,
//!    the cell itself included at 0, which is the k-th smallest to the other cells. BiocNeighbors
//!    caps k at n - 1 (with a warning), and returns 0 when no neighbour is left.
//! 2. `median()` of those distances is the radius, one for all cells.
//! 3. `findNeighbors(idx, threshold = radius)`: each cell's other cells whose raw distance is
//!    `<=` `radius * radius`.
//!
//! **Arithmetic.** knncolle's `EuclideanDistance::raw` sums `delta * delta` over the dims in
//! order, in double, from 0 (`output += delta * delta`); `normalize` takes `sqrt`, `denormalize`
//! squares. Searches compare raw values. Here each pair's raw distance is computed the same way
//! (the order of a subtraction's operands changes nothing: `(a - b)^2` and `(b - a)^2` are equal),
//! which is the sum `fields::rdist` makes before its square root ([`crate::distances`]).
//!
//! **Exhaustive, not KMKNN.** KMKNN visits clusters of cells around k-means centres and skips
//! those that the triangle inequality puts out of reach. It computes every visited pair's raw
//! distance as above and compares it as above, so it returns what an exhaustive search returns,
//! except where its floating-point bounds skip a cell they should not. Searching exhaustively
//! makes the result independent of KMKNN's k-means (`std::mt19937_64` seeded with 6523,
//! k-means++, Hartigan-Wong), and [`Certificate`] reports every pair close enough to a threshold
//! for KMKNN's bounds to treat it differently.
//!
//! # The certificate
//!
//! KMKNN (`knncolle_kmknn.hpp`, as BiocNeighbors 2.4.0 builds it: `KNNCOLLE_KMKNN_USE_UPPER`
//! undefined, so only the lower bound is used) keeps, for each cell `j` of a cluster with centre
//! `c`, `r_j = sqrt(raw(x_j, c))`, computed as above; for a query `t` it computes
//! `D = sqrt(raw(t, c))` and skips every cell of the cluster with `r_j < D - S` (`lower_bd`,
//! by `std::lower_bound` on the sorted `r_j`, or the whole cluster when its largest `r_j` is below
//! it), where `S` is the search's current distance threshold: the radius in step 3, and in step 1
//! the square root of the largest raw distance in its queue once the queue is full.
//!
//! In exact arithmetic `D <= d + r_j` (`d = |t - x_j|`), so no cell within the threshold is
//! skipped. In floating point, with `u = 2^-53` and `m` dims:
//!
//! - a raw distance has a relative error of at most `g = gamma(m + 2) = (m + 2)u / (1 - (m + 2)u)`
//!   (one rounding for the subtraction, two for the square, `m - 1` for the sum; the terms are
//!   non-negative, so the bound holds for the sum), and a distance, after `sqrt`, at most
//!   `e = g/2 + g^2/2 + u + ug`;
//! - so a skip (`r_j < fl(D - S)`) implies `S < d (1 + e) + kappa r_j`, with
//!   `kappa = 2e + u + eu`, after using `D <= d + r_j`;
//! - KMKNN's centres are means of cells (or cells): `kmeans::compute_centroids` sums and divides,
//!   and Hartigan-Wong's last centres come from it, so each lies in the cells' bounding box,
//!   widened per dim by its rounding, `(gamma(n + 2) + 2u) max|x|`. Then `r_j <= B`, the widened
//!   box's diagonal.
//!
//! **Step 3** (the range query, threshold `T`, the radius): a cell is in the result when its raw
//! distance is `<= T2 = fl(T * T)`. KMKNN can only drop cells, and only a cell with
//! `d > (T - kappa B) / (1 + e)`, so with raw distance above `L = ((T - kappa B)/(1 + e))^2 (1 - g)`.
//! Flagged: every pair (i, j), `j != i`, with `L < raw <= T2` and `raw > 0` (a duplicate of the
//! query, raw 0, has `r_j = D` bit for bit and is never skipped).
//!
//! **Step 1** (the k-th distance): KMKNN's result can differ only if, when it reaches the cluster
//! of a cell `j` among the true k + 1 nearest (raw `<= v`, the true (k + 1)-th smallest), its
//! threshold `theta` makes it skip `j`. `theta` is the largest raw distance in its full queue, so
//! `theta >= v`; if `theta = v` the queue already holds k + 1 values `<= v` and the result is `v`
//! whatever is skipped later. So `theta > v` is the raw distance of some other cell `j'`, and the
//! skip needs `sqrt(theta) (1 - u) <= S < d (1 + e) + kappa B` with `d <= sqrt(v / (1 - g))`,
//! that is `theta < H = ((sqrt(v / (1 - g)) (1 + e) + kappa B) / (1 - u))^2`. Flagged: every cell
//! `i` with another cell `j'` at a raw distance in `(v_i, H_i)`.
//!
//! If nothing is flagged, every result equals KMKNN's, whatever its k-means did. The bounds are
//! computed with a little slack (each constant enlarged by 1%, each threshold moved by 16 ulps
//! outwards), so the flagged set can only be larger than the exact one.
//!
//! **The median pair.** When the radius is a cell's own k-th distance (always, with an odd
//! number of cells: the median is the middle value), that cell and its k-th neighbour lie on the
//! radius by construction: their raw distance `v` gives `radius = sqrt(v)`, and the pair is within
//! the radius exactly when `fl(radius * radius) >= v`, which the last bits of the embeddings
//! decide. Such a pair is flagged whenever it is within the radius, both ways, and labelled
//! [`Kind::MedianPair`]; every other flag is [`Kind::Other`]. M12 found nothing but median pairs on
//! the 34 reference objects, and R's KMKNN kept every one of them (docs/validation/clustering.md);
//! the user chose to keep the exhaustive result and report the flagged pairs (M12's decision
//! point).
//!
//! **Memory and threads.** No n x n matrix: each query's distances are computed block by block
//! and consumed at once (a queue of k + 1 values in step 1, the cells within the radius in
//! step 3). Queries are processed in tiles of 8 that share each block of points, in parallel over
//! tiles; each result is placed by its query, so the output does not depend on the thread count.

use crate::Points;
use std::collections::BinaryHeap;
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Queries per tile.
const TQ: usize = 8;
/// Points per block.
const BP: usize = 256;
/// `u = 2^-53`, the unit roundoff of double precision.
const U: f64 = f64::EPSILON / 2.0;

/// The raw distances of a tile of queries to one block of points.
struct Tile {
    /// Query coordinates, dim-major: coordinate c of the tile's query t is `q[c * TQ + t]`.
    q: Vec<f64>,
    nt: usize,
    /// `acc[t * BP + j]`: the raw distance of query t to point `start + j`.
    acc: Vec<f64>,
}

impl Tile {
    fn new(dim: usize) -> Tile {
        Tile { q: vec![0.0; dim * TQ], nt: 0, acc: vec![0.0; TQ * BP] }
    }

    fn load(&mut self, points: &Points, queries: Range<usize>) {
        self.nt = queries.len();
        debug_assert!(self.nt <= TQ);
        for (t, qi) in queries.enumerate() {
            for c in 0..points.dim() {
                self.q[c * TQ + t] = points.get(qi, c);
            }
        }
    }

    /// The tile's raw distances to points `start..end` (at most `BP` of them).
    fn compute(&mut self, points: &Points, start: usize, end: usize) {
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx512f") {
                // SAFETY: the CPU supports AVX-512F, checked just above.
                unsafe { tile_avx512(&self.q, self.nt, &mut self.acc, points, start, end) };
                return;
            }
            if std::arch::is_x86_feature_detected!("avx2") {
                // SAFETY: the CPU supports AVX2, checked just above.
                unsafe { tile_avx2(&self.q, self.nt, &mut self.acc, points, start, end) };
                return;
            }
        }
        tile_body(&self.q, self.nt, &mut self.acc, points, start, end);
    }

    fn row(&self, t: usize, len: usize) -> &[f64] {
        &self.acc[t * BP..t * BP + len]
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn tile_avx512(q: &[f64], nt: usize, acc: &mut [f64], points: &Points, start: usize, end: usize) {
    tile_body(q, nt, acc, points, start, end);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn tile_avx2(q: &[f64], nt: usize, acc: &mut [f64], points: &Points, start: usize, end: usize) {
    tile_body(q, nt, acc, points, start, end);
}

/// knncolle's raw distance for every (query, point) pair of the tile and block: `delta * delta`
/// for the first dim, then `+ delta * delta` for each further dim, in order. Each element sees
/// the same IEEE operations whatever the vector width (subtract, multiply, add; Rust never fuses
/// a multiply and an add), so every code path gives the same bits.
#[inline(always)]
fn tile_body(q: &[f64], nt: usize, acc: &mut [f64], points: &Points, start: usize, end: usize) {
    let len = end - start;
    let col0 = &points.col(0)[start..end];
    for t in 0..nt {
        let qc = q[t];
        let a = &mut acc[t * BP..t * BP + len];
        for (ai, &xi) in a.iter_mut().zip(col0) {
            let d = qc - xi;
            *ai = d * d;
        }
    }
    for c in 1..points.dim() {
        let col = &points.col(c)[start..end];
        for t in 0..nt {
            let qc = q[c * TQ + t];
            let a = &mut acc[t * BP..t * BP + len];
            for (ai, &xi) in a.iter_mut().zip(col) {
                let d = qc - xi;
                *ai += d * d;
            }
        }
    }
}

/// knncolle's raw distance between points `i` and `j`, one pair at a time (the definition the
/// tiled kernel is tested against).
pub fn raw_distance(points: &Points, i: usize, j: usize) -> f64 {
    let mut out = 0.0;
    for c in 0..points.dim() {
        let delta = points.get(i, c) - points.get(j, c);
        out += delta * delta;
    }
    out
}

/// Runs `f` on every tile of queries `0..n`, in parallel over tiles, with one scratch per worker
/// (made by `init`); `f` appends one result per query. The results come back in query order.
fn run_tiles<S, R, I, F>(n: usize, threads: usize, init: I, f: F) -> Vec<R>
where
    R: Send,
    I: Fn() -> S + Sync,
    F: Fn(&mut S, Range<usize>, &mut Vec<R>) + Sync,
{
    let ntiles = n.div_ceil(TQ);
    let threads = threads.max(1).min(ntiles.max(1));
    let next = AtomicUsize::new(0);
    let parts: Mutex<Vec<(usize, Vec<R>)>> = Mutex::new(Vec::new());
    let work = || {
        let mut scratch = init();
        loop {
            let t = next.fetch_add(1, Ordering::Relaxed);
            if t >= ntiles {
                break;
            }
            let qs = t * TQ..((t + 1) * TQ).min(n);
            let mut out = Vec::with_capacity(qs.len());
            f(&mut scratch, qs, &mut out);
            parts.lock().unwrap().push((t, out));
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
    let mut parts = parts.into_inner().unwrap();
    parts.sort_unstable_by_key(|p| p.0);
    let mut all = Vec::with_capacity(n);
    for (_, p) in parts {
        all.extend(p);
    }
    debug_assert_eq!(all.len(), n);
    all
}

/// The smallest of a few raw distances (none is NaN).
#[inline(always)]
fn chunk_min(chunk: &[f64]) -> f64 {
    chunk.iter().fold(f64::INFINITY, |a, &v| if v < a { v } else { a })
}

/// The ordering key of a raw distance (non-negative, not NaN): its bits.
#[inline]
fn key(v: f64) -> u64 {
    v.to_bits()
}

/// The unit in the last place of a positive finite double.
pub fn ulp(x: f64) -> f64 {
    let a = x.abs();
    a.next_up() - a
}

/// The constants of the certificate (module documentation): how far KMKNN's floating-point
/// bounds can be off.
#[derive(Debug, Clone, PartialEq)]
pub struct Margin {
    /// Dims summed per distance.
    pub dims: usize,
    /// Relative error bound of a raw distance, `gamma(m + 2)`, enlarged by 1%.
    pub gamma: f64,
    /// Relative error bound of a distance after `sqrt`, enlarged by 1%.
    pub eps: f64,
    /// `kappa = 2e + u + eu`, enlarged by 1%.
    pub kappa: f64,
    /// `B`: a bound on the distance of any cell to any k-means centre (the diagonal of the
    /// cells' bounding box, each side widened by the centres' rounding).
    pub bound: f64,
}

impl Margin {
    /// The margin of the cells `points`.
    pub fn of(points: &Points) -> Margin {
        let m = points.dim();
        let n = points.n();
        let gamma_k = |k: f64| k * U / (1.0 - k * U);
        let g = gamma_k(m as f64 + 2.0);
        let eps = g / 2.0 + g * g / 2.0 + U + U * g;
        let kappa = 2.0 * eps + U + eps * U;
        // the centres: means computed by summing (n - 1 additions) and dividing
        let gn = gamma_k(n as f64 + 2.0) + 2.0 * U;
        let mut s = 0.0;
        for c in 0..m {
            let col = points.col(c);
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for &v in col {
                lo = lo.min(v);
                hi = hi.max(v);
            }
            let w = (hi - lo) + gn * lo.abs().max(hi.abs());
            s += w * w;
        }
        let bound = s.sqrt() * (1.0 + 1e-12);
        Margin { dims: m, gamma: g * 1.01, eps: eps * 1.01, kappa: kappa * 1.01, bound }
    }

    /// Step 3: pairs with a raw distance above this (and at most `radius * radius`) are flagged.
    /// Negative when the margin exceeds the radius (then every pair within it is flagged).
    pub fn radius_low(&self, radius: f64) -> f64 {
        let x = (radius - self.kappa * self.bound) / (1.0 + self.eps);
        if x > 0.0 {
            x * x * (1.0 - self.gamma) * (1.0 - 16.0 * U)
        } else {
            -1.0
        }
    }

    /// Step 1: a cell with raw k-th distance `v` is flagged when another cell's raw distance lies
    /// in `(v, kth_high(v))`.
    pub fn kth_high(&self, v: f64) -> f64 {
        let x = ((v / (1.0 - self.gamma)).sqrt() * (1.0 + self.eps) + self.kappa * self.bound) / (1.0 - U);
        x * x * (1.0 + 16.0 * U)
    }
}

/// Which of KMKNN's searches a flagged pair concerns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// `findDistance`: the pair's raw distance lies just above the cell's k-th.
    KthDistance,
    /// `findNeighbors`: the pair lies just within the radius.
    Radius,
}

impl Step {
    pub fn name(self) -> &'static str {
        match self {
            Step::KthDistance => "kth_distance",
            Step::Radius => "radius",
        }
    }
}

/// What kind of pair a flag is (module documentation, "The median pair").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A cell whose k-th distance is the radius, and its k-th neighbour: on the radius by
    /// construction (with an odd number of cells, the median cell's pair).
    MedianPair,
    /// Anything else.
    Other,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::MedianPair => "median_pair",
            Kind::Other => "other",
        }
    }
}

/// A pair the certificate flags (module documentation).
#[derive(Debug, Clone, PartialEq)]
pub struct Flag {
    pub step: Step,
    pub kind: Kind,
    /// The query cell and the other cell, 0-based.
    pub i: u32,
    pub j: u32,
    /// The pair's raw distance.
    pub raw: f64,
    /// The threshold it is close to: cell i's raw k-th distance (step 1), or `radius * radius`
    /// (step 3).
    pub threshold: f64,
    /// The end of the flagged band: `kth_high(v_i)` (step 1) or `radius_low(radius)` (step 3).
    pub band: f64,
}

impl Flag {
    /// How far the pair lies from its threshold, in units in the last place of the threshold
    /// (positive above it).
    pub fn ulps(&self) -> f64 {
        let u = ulp(self.threshold.max(f64::MIN_POSITIVE));
        (self.raw - self.threshold) / u
    }
}

/// What the certificate saw: the flagged pairs, and how close the pairs it did not flag came.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Certificate {
    pub flags: Vec<Flag>,
    /// Step 3: the gap between `radius * radius` and the largest raw distance of an included
    /// pair (other than duplicates of the query, at 0); `None` when no pair is included.
    pub radius_gap_below: Option<f64>,
    /// Step 3: the gap to the smallest raw distance above `radius * radius` among the pairs seen
    /// (those up to `2^-20` above it, relatively).
    pub radius_gap_above: Option<f64>,
    /// `radius * radius - radius_low`: the width of step 3's band.
    pub radius_band: f64,
    /// Step 1: the smallest ratio, over the cells, of the gap from the raw k-th distance to the
    /// next larger raw distance, to the width of the cell's band (`kth_high(v) - v`), among the
    /// cells whose next larger value lies within `2^-20` of `v`, relatively. Below 1 means
    /// flagged.
    pub kth_gap_ratio: Option<f64>,
    /// The widest and narrowest step-1 band, `kth_high(v) - v`, over the cells.
    pub kth_band_max: f64,
    pub kth_band_min: f64,
}

impl Certificate {
    pub fn n_flagged(&self) -> usize {
        self.flags.len()
    }

    pub fn n_flagged_step(&self, step: Step) -> usize {
        self.flags.iter().filter(|f| f.step == step).count()
    }

    pub fn n_flagged_kind(&self, kind: Kind) -> usize {
        self.flags.iter().filter(|f| f.kind == kind).count()
    }
}

/// The first two searches: each cell's k-th distance and the radius.
#[derive(Debug, Clone, PartialEq)]
pub struct KthDistances {
    pub n: usize,
    /// The `k` asked for.
    pub k_requested: usize,
    /// `k` after BiocNeighbors' cap at `n - 1` ("'k' capped at the number of observations minus
    /// 1").
    pub k: usize,
    /// Each cell's raw k-th distance (0 when `k = 0`).
    pub raw_k: Vec<f64>,
    /// `findDistance`: `sqrt` of `raw_k`.
    pub dist_k: Vec<f64>,
    /// `median(dist_k)` (R's).
    pub radius: f64,
    /// `radius * radius`, knncolle's `denormalize`.
    pub radius_raw: f64,
    pub margin: Margin,
}

impl KthDistances {
    /// True when BiocNeighbors would warn that it capped `k`.
    pub fn capped(&self) -> bool {
        self.k < self.k_requested
    }
}

/// Steps 1 and 2: every cell's distance to its k-th nearest other cell, and their median.
/// Panics on non-finite coordinates (BiocNeighbors refuses NA).
pub fn kth_distances(points: &Points, k: usize, threads: usize) -> KthDistances {
    let n = points.n();
    assert!(points_are_finite(points), "neighborPurity: the coordinates must be finite");
    let k_eff = if k < n { k } else { n - 1 };
    let m = k_eff + 1; // the queue holds the cell itself too
    let raw_k: Vec<f64> = run_tiles(
        n,
        threads,
        || (Tile::new(points.dim()), (0..TQ).map(|_| BinaryHeap::<u64>::with_capacity(m + 1)).collect::<Vec<_>>()),
        |(tile, heaps), qs, out| {
            tile.load(points, qs.clone());
            let nt = qs.len();
            for h in heaps.iter_mut() {
                h.clear();
            }
            let mut start = 0;
            while start < n {
                let end = (start + BP).min(n);
                tile.compute(points, start, end);
                for (t, h) in heaps.iter_mut().enumerate().take(nt) {
                    let row = tile.row(t, end - start);
                    let mut rest = row;
                    // fill the queue
                    while h.len() < m && !rest.is_empty() {
                        h.push(key(rest[0]));
                        rest = &rest[1..];
                    }
                    if rest.is_empty() {
                        continue;
                    }
                    let mut thr = f64::from_bits(*h.peek().unwrap());
                    for chunk in rest.chunks(8) {
                        // most chunks hold nothing below the threshold
                        if chunk_min(chunk) >= thr {
                            continue;
                        }
                        for &v in chunk {
                            if v < thr {
                                *h.peek_mut().unwrap() = key(v);
                                thr = f64::from_bits(*h.peek().unwrap());
                            }
                        }
                    }
                }
                start = end;
            }
            for h in heaps.iter().take(nt) {
                out.push(f64::from_bits(*h.peek().expect("n >= 1")));
            }
        },
    );
    // the cell itself, at 0, is among the k + 1: its k-th is the (k + 1)-th smallest of all
    let dist_k: Vec<f64> = raw_k.iter().map(|&v| v.sqrt()).collect();
    let radius = rcompat::quantile::median(&dist_k);
    let radius_raw = radius * radius;
    KthDistances { n, k_requested: k, k: k_eff, raw_k, dist_k, radius, radius_raw, margin: Margin::of(points) }
}

fn points_are_finite(points: &Points) -> bool {
    (0..points.dim()).all(|c| points.col(c).iter().all(|v| v.is_finite()))
}

/// Per-query scratch of step 3.
struct Within {
    tile: Tile,
    /// Each query's neighbours: (raw distance, index).
    lists: Vec<Vec<(f64, u32)>>,
    sorted: Vec<u32>,
}

/// The relative width of the watch band above a threshold within which the certificate reports
/// how close the pairs came: `2^-20`.
const WATCH: f64 = 1.0 / 1048576.0;

impl KthDistances {
    /// Step 3, and the certificate: for each cell `i`, `f(i, neighbours)` with the other cells
    /// whose raw distance is `<= radius * radius` (0-based; in index order, or with `sorted` by
    /// raw distance and then index, as BiocNeighbors sorts them by raw distance). Returns `f`'s
    /// results in cell order, and the certificate of both searches.
    pub fn within<R, F>(&self, points: &Points, sorted: bool, threads: usize, f: F) -> (Vec<R>, Certificate)
    where
        R: Send,
        F: Fn(usize, &[u32]) -> R + Sync,
    {
        let n = self.n;
        assert_eq!(points.n(), n);
        let t2 = self.radius_raw;
        let low = self.margin.radius_low(self.radius);
        let watch_hi = t2 * (1.0 + WATCH);
        // per cell: the end of its step-1 band, and what it must see
        let hi: Vec<f64> = self.raw_k.iter().map(|&v| self.margin.kth_high(v)).collect();
        let limit: Vec<f64> = self.raw_k.iter().zip(&hi).map(|(&v, &h)| h.max(t2).max(watch_hi).max(v * (1.0 + WATCH))).collect();
        let (mut kth_band_max, mut kth_band_min) = (0.0f64, f64::INFINITY);
        for (&v, &h) in self.raw_k.iter().zip(&hi) {
            kth_band_max = kth_band_max.max(h - v);
            kth_band_min = kth_band_min.min(h - v);
        }
        let parts: Vec<(R, Vec<Flag>, Option<f64>, Option<f64>, Option<f64>)> = run_tiles(
            n,
            threads,
            || Within { tile: Tile::new(points.dim()), lists: vec![Vec::new(); TQ], sorted: Vec::new() },
            |w, qs, out| {
                w.tile.load(points, qs.clone());
                let nt = qs.len();
                let q0 = qs.start;
                for l in w.lists.iter_mut() {
                    l.clear();
                }
                // per query: step-1 tracking (the smallest raw distance above v), watch gaps
                let mut next_above: [f64; TQ] = [f64::INFINITY; TQ];
                let mut flags: Vec<Flag> = Vec::new();
                let (mut gap_below, mut gap_above): (Option<f64>, Option<f64>) = (None, None);
                let mut start = 0;
                while start < n {
                    let end = (start + BP).min(n);
                    w.tile.compute(points, start, end);
                    for t in 0..nt {
                        let i = q0 + t;
                        let lim = limit[i];
                        let v = self.raw_k[i];
                        let row = w.tile.row(t, end - start);
                        for (c, chunk) in row.chunks(8).enumerate() {
                            if chunk_min(chunk) > lim {
                                continue;
                            }
                            for (o, &raw) in chunk.iter().enumerate() {
                                if raw > lim {
                                    continue;
                                }
                                let j = start + c * 8 + o;
                                if j == i {
                                    continue;
                                }
                                if raw <= t2 {
                                    w.lists[t].push((raw, j as u32));
                                    if raw > 0.0 {
                                        let g = t2 - raw;
                                        gap_below = Some(gap_below.map_or(g, |x: f64| x.min(g)));
                                        if raw > low {
                                            // the median pair: the raw distance is the k-th of a cell whose
                                            // k-th distance is the radius
                                            let on = |c: usize| raw.to_bits() == self.raw_k[c].to_bits() && self.dist_k[c].to_bits() == self.radius.to_bits();
                                            let kind = if on(i) || on(j) { Kind::MedianPair } else { Kind::Other };
                                            flags.push(Flag { step: Step::Radius, kind, i: i as u32, j: j as u32, raw, threshold: t2, band: low });
                                        }
                                    }
                                } else if raw <= watch_hi {
                                    let g = raw - t2;
                                    gap_above = Some(gap_above.map_or(g, |x: f64| x.min(g)));
                                }
                                if raw > v {
                                    if raw < next_above[t] {
                                        next_above[t] = raw;
                                    }
                                    if raw < hi[i] {
                                        flags.push(Flag { step: Step::KthDistance, kind: Kind::Other, i: i as u32, j: j as u32, raw, threshold: v, band: hi[i] });
                                    }
                                }
                            }
                        }
                    }
                    start = end;
                }
                let mut kth_ratio: Option<f64> = None;
                for t in 0..nt {
                    let i = q0 + t;
                    let v = self.raw_k[i];
                    if next_above[t].is_finite() && next_above[t] <= v * (1.0 + WATCH) {
                        let r = (next_above[t] - v) / (hi[i] - v);
                        kth_ratio = Some(kth_ratio.map_or(r, |x: f64| x.min(r)));
                    }
                }
                for t in 0..nt {
                    let list = &mut w.lists[t];
                    w.sorted.clear();
                    if sorted {
                        list.sort_unstable_by(|a, b| key(a.0).cmp(&key(b.0)).then(a.1.cmp(&b.1)));
                    }
                    w.sorted.extend(list.iter().map(|p| p.1));
                    let r = f(q0 + t, &w.sorted);
                    // the tile's flags and gaps travel with its first query's result
                    if t == 0 {
                        out.push((r, std::mem::take(&mut flags), gap_below, gap_above, kth_ratio));
                    } else {
                        out.push((r, Vec::new(), None, None, None));
                    }
                }
            },
        );
        let mut results = Vec::with_capacity(n);
        let mut cert = Certificate { radius_band: if low < 0.0 { t2 } else { t2 - low }, kth_band_max, kth_band_min, ..Certificate::default() };
        for (r, fl, gb, ga, kr) in parts {
            results.push(r);
            cert.flags.extend(fl);
            if let Some(g) = gb {
                cert.radius_gap_below = Some(cert.radius_gap_below.map_or(g, |x| x.min(g)));
            }
            if let Some(g) = ga {
                cert.radius_gap_above = Some(cert.radius_gap_above.map_or(g, |x| x.min(g)));
            }
            if let Some(g) = kr {
                cert.kth_gap_ratio = Some(cert.kth_gap_ratio.map_or(g, |x| x.min(g)));
            }
        }
        // flags in (step, i, j) order: step 1 first, as the searches run
        cert.flags.sort_by(|a, b| (a.step == Step::Radius, a.i, a.j).cmp(&(b.step == Step::Radius, b.i, b.j)));
        (results, cert)
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

    fn random_points(n: usize, dim: usize, seed: u64, scale: f64) -> Vec<f64> {
        let mut g = Lcg(seed);
        (0..n * dim).map(|_| (g.next() * 2.0 - 1.0) * scale).collect()
    }

    /// The k-th smallest raw distance to the other cells, by a full sort.
    fn brute_kth(p: &Points, i: usize, k: usize) -> f64 {
        let mut d: Vec<f64> = (0..p.n()).filter(|&j| j != i).map(|j| raw_distance(p, i, j)).collect();
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if k == 0 {
            0.0
        } else {
            d[k - 1]
        }
    }

    #[test]
    fn every_code_path_gives_knncolle_s_bits() {
        let (n, dim) = (700usize, 30usize);
        let data = random_points(n, dim, 3, 40.0);
        let p = Points::new(&data, n, dim);
        let mut a = Tile::new(dim);
        let mut b = Tile::new(dim);
        for (q0, start) in [(0usize, 0usize), (8, 256), (693, 512)] {
            let qs = q0..(q0 + TQ).min(n);
            a.load(&p, qs.clone());
            b.load(&p, qs.clone());
            let end = (start + BP).min(n);
            a.compute(&p, start, end);
            tile_body(&b.q, b.nt, &mut b.acc, &p, start, end);
            for t in 0..qs.len() {
                for j in 0..end - start {
                    let want = raw_distance(&p, q0 + t, start + j);
                    assert_eq!(a.row(t, end - start)[j].to_bits(), want.to_bits(), "dispatched path, q {} j {}", q0 + t, start + j);
                    assert_eq!(b.row(t, end - start)[j].to_bits(), want.to_bits(), "scalar path, q {} j {}", q0 + t, start + j);
                }
            }
            #[cfg(target_arch = "x86_64")]
            if std::arch::is_x86_feature_detected!("avx2") {
                let mut c = Tile::new(dim);
                c.load(&p, qs.clone());
                unsafe { tile_avx2(&c.q, c.nt, &mut c.acc, &p, start, end) };
                for t in 0..qs.len() {
                    let (x, y) = (c.row(t, end - start), a.row(t, end - start));
                    assert!(x.iter().zip(y).all(|(u, v)| u.to_bits() == v.to_bits()), "AVX2 path, query {}", q0 + t);
                }
            }
        }
        // the sum is the one fields::rdist makes before its square root
        let mut d = vec![0.0; n];
        crate::distances(&p, 5, &mut d);
        for (j, dj) in d.iter().enumerate() {
            assert_eq!(dj.to_bits(), raw_distance(&p, 5, j).sqrt().to_bits());
        }
    }

    #[test]
    fn kth_distance_and_radius_by_brute_force() {
        for (n, dim, k, seed) in [(300usize, 30usize, 50usize, 1u64), (257, 5, 10, 2), (90, 3, 50, 3), (51, 2, 50, 4), (40, 2, 50, 5), (13, 4, 0, 6)] {
            let data = random_points(n, dim, seed, 10.0);
            let p = Points::new(&data, n, dim);
            for threads in [1, 3] {
                let s = kth_distances(&p, k, threads);
                assert_eq!(s.k, k.min(n - 1));
                assert_eq!(s.capped(), k >= n);
                for i in 0..n {
                    let want = brute_kth(&p, i, s.k);
                    assert_eq!(s.raw_k[i].to_bits(), want.to_bits(), "n {n} i {i}");
                    assert_eq!(s.dist_k[i].to_bits(), want.sqrt().to_bits());
                }
                assert_eq!(s.radius.to_bits(), rcompat::quantile::median(&s.dist_k).to_bits());
                assert_eq!(s.radius_raw, s.radius * s.radius);
                let (lists, _) = s.within(&p, true, threads, |_, nb| nb.to_vec());
                for i in 0..n {
                    let mut want: Vec<(f64, u32)> = (0..n).filter(|&j| j != i).map(|j| (raw_distance(&p, i, j), j as u32)).filter(|&(r, _)| r <= s.radius_raw).collect();
                    want.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
                    let want: Vec<u32> = want.into_iter().map(|w| w.1).collect();
                    assert_eq!(lists[i], want, "n {n} i {i}");
                }
                let (unsorted, _) = s.within(&p, false, threads, |_, nb| nb.to_vec());
                for i in 0..n {
                    let mut a = lists[i].clone();
                    a.sort_unstable();
                    assert_eq!(unsorted[i], a, "index order");
                }
            }
        }
    }

    #[test]
    fn duplicates_ties_and_the_self_exclusion() {
        // 30 cells, 12 of them at one point: their k-th distance is 0 up to k = 11
        let (n, dim) = (30usize, 2usize);
        let mut data = random_points(n, dim, 9, 5.0);
        for i in 0..12 {
            data[i] = 1.25;
            data[n + i] = -0.5;
        }
        let p = Points::new(&data, n, dim);
        for k in [1usize, 5, 11, 12, 20] {
            let s = kth_distances(&p, k, 2);
            for i in 0..n {
                assert_eq!(s.raw_k[i].to_bits(), brute_kth(&p, i, k).to_bits(), "k {k} i {i}");
            }
            if k <= 11 {
                assert!(s.raw_k[..12].iter().all(|&v| v == 0.0));
            }
            let (lists, cert) = s.within(&p, true, 1, |_, nb| nb.to_vec());
            // a duplicate at raw 0 is never flagged
            assert!(cert.flags.iter().all(|f| f.raw > 0.0));
            for i in 0..12 {
                assert!(lists[i].iter().take(11).all(|&j| j < 12), "the duplicates come first, by index");
            }
        }
    }

    #[test]
    fn independent_of_the_thread_count() {
        let (n, dim) = (2000usize, 30usize);
        let data = random_points(n, dim, 11, 20.0);
        let p = Points::new(&data, n, dim);
        let one = kth_distances(&p, 50, 1);
        let (l1, c1) = one.within(&p, true, 1, |i, nb| (i, nb.to_vec()));
        for t in [2, 5, 16] {
            let s = kth_distances(&p, 50, t);
            assert_eq!(s, one, "{t} threads");
            let (l, c) = s.within(&p, true, t, |i, nb| (i, nb.to_vec()));
            assert_eq!(l, l1);
            assert_eq!(c, c1);
        }
    }

    #[test]
    fn the_median_cell_s_kth_neighbour_sits_on_the_radius() {
        // odd n: the radius is one cell's k-th distance, so that pair lies exactly at it, and is
        // flagged whenever squaring the square root keeps it within the radius
        let (n, dim) = (301usize, 10usize);
        let data = random_points(n, dim, 21, 3.0);
        let p = Points::new(&data, n, dim);
        let s = kth_distances(&p, 50, 2);
        let m = (0..n).find(|&i| s.dist_k[i] == s.radius).unwrap();
        let (_, cert) = s.within(&p, false, 2, |_, _| ());
        let included = s.raw_k[m] <= s.radius_raw;
        let flagged = cert.flags.iter().any(|f| f.step == Step::Radius && f.kind == Kind::MedianPair && f.i as usize == m && f.raw == s.raw_k[m]);
        assert_eq!(flagged, included);
        // and the other way round, labelled the same
        let back = cert.flags.iter().any(|f| f.step == Step::Radius && f.kind == Kind::MedianPair && f.j as usize == m && f.raw == s.raw_k[m]);
        assert_eq!(back, included);
        for f in &cert.flags {
            match f.step {
                Step::Radius => assert!(f.raw <= s.radius_raw && f.raw > s.margin.radius_low(s.radius)),
                Step::KthDistance => assert!(f.raw > s.raw_k[f.i as usize] && f.raw < s.margin.kth_high(s.raw_k[f.i as usize])),
            }
        }
    }

    #[test]
    fn a_pair_just_above_the_kth_distance_is_flagged() {
        // cell 0 at the origin; cells 1..=3 at distance 1 along axis 0, cell 4 one ulp further:
        // with k = 3, cell 0's k-th raw distance is 1 and cell 4 sits just above it
        let x0 = [0.0, 1.0, -1.0, 0.0, 1.0f64.next_up(), 50.0, 60.0, 70.0];
        let x1 = [0.0, 0.0, 0.0, 1.0, 0.0, 50.0, 60.0, 70.0];
        let n = x0.len();
        let data: Vec<f64> = x0.iter().chain(x1.iter()).copied().collect();
        let p = Points::new(&data, n, 2);
        let s = kth_distances(&p, 3, 1);
        assert_eq!(s.raw_k[0], 1.0);
        let (_, cert) = s.within(&p, false, 1, |_, _| ());
        assert!(cert.flags.iter().any(|f| f.step == Step::KthDistance && f.kind == Kind::Other && f.i == 0 && f.j == 4), "{:?}", cert.flags);
    }

    #[test]
    fn the_margin_bounds_the_rounding() {
        // the bound on a raw distance's relative error holds on random data, against the exact
        // value in long double (F80): the sum of squares of exactly computed differences
        let (n, dim) = (200usize, 30usize);
        let data = random_points(n, dim, 33, 100.0);
        let p = Points::new(&data, n, dim);
        let mg = Margin::of(&p);
        for i in 0..20 {
            for j in 0..n {
                let raw = raw_distance(&p, i, j);
                let mut exact = rcompat::F80::ZERO;
                for c in 0..dim {
                    let d = rcompat::F80::from_f64(p.get(i, c)) - rcompat::F80::from_f64(p.get(j, c));
                    exact = exact + d * d;
                }
                let e = exact.to_f64();
                assert!((raw - e).abs() <= mg.gamma * e + 1e-300, "pair {i} {j}");
            }
        }
        // the band of step 3 widens with the data's extent
        assert!(mg.radius_low(10.0) < 100.0 && mg.radius_low(10.0) > 99.9);
        assert!(mg.kth_high(100.0) > 100.0 && mg.kth_high(100.0) < 100.1);
        assert_eq!(mg.radius_low(1e-300), -1.0);
    }
}
