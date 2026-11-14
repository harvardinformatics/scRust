//! Annoy 1.17.3, as RcppAnnoy 0.0.23 instantiates it for Seurat (docs/PLAN.md C3a):
//! `AnnoyIndex<int32_t, float, D, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy>` with
//! `D` = `Euclidean`, `Angular` (Seurat's `annoy.metric = "cosine"`) or `Manhattan`. Ported
//! statement by statement from `annoylib.h`; the comments name the C++ functions.
//!
//! *Ported from Annoy 1.17.3's `annoylib.h` (Copyright (c) 2013 Spotify AB, Apache License 2.0,
//! `third_party/licenses/Apache-2.0.txt`). This file is a modified version: translated into Rust,
//! for the three metrics Seurat uses; NOTICE lists the changes.*
//!
//! **Arithmetic.** Everything is IEEE single precision, as compiled for R (`-O2
//! -march=nocona`): no FMA, no AVX, every reduction a sequential loop (GCC vectorizes only the
//! element-wise products and differences, which does not change any value). Two places mix in
//! double precision, as written: `_split_imbalance` and the Angular distance
//! `2.0 - 2.0 * pq / sqrt(ppqq)`, whose `sqrt` is the float one (`sqrtss` in `RcppAnnoy.so`) and
//! whose quotient and difference are double, rounded to float on return.
//!
//! **Memory layout.** The nodes live in one array of 4-byte words, exactly as Annoy holds them
//! and `save()` writes them, so the file is the array itself ([`AnnoyIndex::to_bytes`]):
//!
//! | metric | words | `n_descendants` | offset `a` | `children` | `v` (f floats) |
//! |---|---|---|---|---|---|
//! | Euclidean, Manhattan (`Minkowski::Node`) | 4 + f | 0 | 1 | 2, 3 | 4.. |
//! | Angular (`Angular::Node`) | 3 + f | 0 | (`norm` at 1, in a union with `children[0]`) | 1, 2 | 3.. |
//!
//! A leaf stores up to K = f + 2 item ids from `children` on, into the space of `v`. Memory that
//! Annoy never writes is zero (`_reallocate_nodes` clears what it adds), so the offset `a` of an
//! item or a leaf is 0.
//!
//! **Node ids.** Items take 0..n-1. Then, tree by tree, each node takes the next id when it is
//! created: a leaf when it is reached, a split node after both its children; the smaller side is
//! built first (side 0 on equal sizes). A subset of one item is that item, not a new node. After
//! the last tree come copies of the roots, in tree order.

use crate::kiss::Kiss64Random;
use crate::stl::{pair_less, partial_sort, PriorityQueue};
use std::fmt;

/// The distance Annoy builds and searches with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Metric {
    /// `Annoy::Euclidean` (`AnnoyEuclidean`; Seurat's default).
    Euclidean,
    /// `Annoy::Angular` (`AnnoyAngular`; Seurat's `annoy.metric = "cosine"`).
    Angular,
    /// `Annoy::Manhattan` (`AnnoyManhattan`).
    Manhattan,
}

impl Metric {
    /// The name Seurat's `annoy.metric` uses: euclidean, cosine or manhattan.
    pub fn seurat_name(self) -> &'static str {
        match self {
            Metric::Euclidean => "euclidean",
            Metric::Angular => "cosine",
            Metric::Manhattan => "manhattan",
        }
    }

    /// Annoy's own name (`D::name()`): euclidean, angular or manhattan.
    pub fn annoy_name(self) -> &'static str {
        match self {
            Metric::Euclidean => "euclidean",
            Metric::Angular => "angular",
            Metric::Manhattan => "manhattan",
        }
    }

    /// Seurat's names, and Annoy's `angular` for `cosine`. Seurat's fourth metric, hamming
    /// (`AnnoyHamming`, on 64-bit integers), is another code path and is refused.
    pub fn parse(s: &str) -> Result<Metric, String> {
        match s {
            "euclidean" => Ok(Metric::Euclidean),
            "cosine" | "angular" => Ok(Metric::Angular),
            "manhattan" => Ok(Metric::Manhattan),
            "hamming" => Err("annoy.metric = \"hamming\" (Annoy's Hamming index on integers) is not ported".into()),
            _ => Err(format!("unknown metric '{s}': expected euclidean, cosine or manhattan")),
        }
    }

    fn layout(self, f: usize) -> Layout {
        match self {
            Metric::Euclidean | Metric::Manhattan => Layout { words: 4 + f, children: 2, v: 4 },
            Metric::Angular => Layout { words: 3 + f, children: 1, v: 3 },
        }
    }
}

impl fmt::Display for Metric {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.seurat_name())
    }
}

/// Word offsets of a node's fields.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Layout {
    /// Words per node: `_s / 4`.
    words: usize,
    /// `children` (and a leaf's item list).
    children: usize,
    /// `v`.
    v: usize,
}

/// What happened while building, for the reports: how often the rare branches ran.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BuildStats {
    pub trees: usize,
    /// Split nodes created (not counting the root copies).
    pub split_nodes: usize,
    /// Leaf nodes created (subsets of 2..K items, and roots of small indexes).
    pub leaves: usize,
    /// Splits accepted at the first, second and third attempt (`_split_imbalance < 0.95`).
    pub accepted_at_attempt: [usize; 3],
    /// Splits whose three attempts all stayed at or above 0.95, kept because the last one was at
    /// most 0.99.
    pub kept_imbalanced: usize,
    /// Splits that fell back to random sides (the last attempt above 0.99): the node keeps its
    /// last offset with a zero normal.
    pub random_fallbacks: usize,
    /// Rounds of random sides in those fallbacks (one `flip()` per point each).
    pub random_rounds: usize,
    /// Points whose margin was exactly 0 during a split, each costing one `flip()`.
    pub zero_margins: usize,
    /// Calls to `kiss()`.
    pub draws: u64,
    /// The deepest split node below a root (the root is depth 0).
    pub max_depth: usize,
}

/// An Annoy forest: [`AnnoyIndex::add_item`] for each item, then [`AnnoyIndex::build`] once.
#[derive(Debug, Clone)]
pub struct AnnoyIndex {
    metric: Metric,
    f: usize,
    lay: Layout,
    /// `_K`: the most item ids a leaf holds, f + 2.
    k_leaf: usize,
    n_items: usize,
    nodes: Vec<u32>,
    roots: Vec<i32>,
    seed: u64,
    built: bool,
    stats: BuildStats,
}

/// `&[u32]` words as `&[f32]`.
#[inline]
fn as_f32(w: &[u32]) -> &[f32] {
    // SAFETY: u32 and f32 have the same size and alignment, and every bit pattern is a valid
    // f32, so the reinterpretation is sound (it is what Annoy's union-based nodes do).
    unsafe { std::slice::from_raw_parts(w.as_ptr() as *const f32, w.len()) }
}

#[inline]
fn as_f32_mut(w: &mut [u32]) -> &mut [f32] {
    // SAFETY: as in `as_f32`; the borrow is exclusive.
    unsafe { std::slice::from_raw_parts_mut(w.as_mut_ptr() as *mut f32, w.len()) }
}

/// `dot(x, y, f)`: `s = 0; s += x[z] * y[z]` in order.
#[inline]
fn dot(x: &[f32], y: &[f32]) -> f32 {
    let mut s = 0.0f32;
    for (a, b) in x.iter().zip(y) {
        s += a * b;
    }
    s
}

/// `euclidean_distance(x, y, f)`: the squared distance, `d += tmp * tmp` in order.
#[inline]
fn euclidean_distance(x: &[f32], y: &[f32]) -> f32 {
    let mut d = 0.0f32;
    for (a, b) in x.iter().zip(y) {
        let t = a - b;
        d += t * t;
    }
    d
}

/// `manhattan_distance(x, y, f)`: `d += fabs(x[i] - y[i])` in order.
#[inline]
fn manhattan_distance(x: &[f32], y: &[f32]) -> f32 {
    let mut d = 0.0f32;
    for (a, b) in x.iter().zip(y) {
        d += (a - b).abs();
    }
    d
}

/// `get_norm(v, f)`: `sqrt(dot(v, v, f))`, in float.
#[inline]
fn get_norm(v: &[f32]) -> f32 {
    dot(v, v).sqrt()
}

/// `Base::normalize`: divides by the norm when it is positive.
#[inline]
fn normalize(v: &mut [f32]) {
    let norm = get_norm(v);
    if norm > 0.0 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
}

/// `Angular::distance(x, y, f)`: `2 - 2 cos`, with each side's stored `norm` (its squared
/// length) where it is nonzero. The quotient and the difference are double, as written
/// (`2.0` is a double literal); `sqrt` is the float one.
#[inline]
fn angular_distance(xv: &[f32], xnorm: f32, yv: &[f32], ynorm: f32) -> f32 {
    let pp = if xnorm != 0.0 { xnorm } else { dot(xv, xv) };
    let qq = if ynorm != 0.0 { ynorm } else { dot(yv, yv) };
    let pq = dot(xv, yv);
    let ppqq = pp * qq;
    if ppqq > 0.0 {
        (2.0f64 - 2.0f64 * pq as f64 / ppqq.sqrt() as f64) as f32
    } else {
        2.0
    }
}

/// `std::min(a, b)`: `(b < a) ? b : a`.
#[inline]
fn std_min(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}

/// `std::max(a, b)`: `(a < b) ? b : a`.
#[inline]
fn std_max(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}

/// `pq_distance(distance, margin, child_nr)` (Minkowski and Angular): the child's priority.
#[inline]
fn pq_distance(distance: f32, margin: f32, child_nr: usize) -> f32 {
    let margin = if child_nr == 0 { -margin } else { margin };
    std_min(distance, margin)
}

/// `_split_imbalance`: the larger side's share, mixing float and double as written.
#[inline]
pub(crate) fn split_imbalance(left: usize, right: usize) -> f64 {
    let ls = (left as f32) as f64;
    let rs = (right as f32) as f64;
    let f = (ls / (ls + rs + 1e-9)) as f32;
    std_max(f, 1.0f32 - f) as f64
}

/// A two-means centroid (`p`, `q` in `two_means`), with the `norm` that Angular keeps.
#[derive(Debug, Clone)]
struct Centroid {
    v: Vec<f32>,
    norm: f32,
}

/// Per-thread buffers for [`AnnoyIndex::get_nns_by_vector`].
pub struct SearchScratch {
    q: PriorityQueue<(f32, i32), fn(&(f32, i32), &(f32, i32)) -> bool>,
    nns: Vec<i32>,
    dist: Vec<(f32, i32)>,
}

impl SearchScratch {
    pub fn new() -> SearchScratch {
        SearchScratch { q: PriorityQueue::new(pair_less as fn(&(f32, i32), &(f32, i32)) -> bool), nns: Vec::new(), dist: Vec::new() }
    }
}

impl Default for SearchScratch {
    fn default() -> Self {
        SearchScratch::new()
    }
}

/// Why an index file cannot be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoadError(pub String);

impl fmt::Display for LoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for LoadError {}

impl AnnoyIndex {
    /// `AnnoyIndex(f)`: an empty index of `f`-dimensional items, seeded with
    /// [`Kiss64Random::DEFAULT_SEED`].
    pub fn new(metric: Metric, f: usize) -> AnnoyIndex {
        assert!(f >= 1, "Annoy: f must be at least 1");
        let lay = metric.layout(f);
        AnnoyIndex {
            metric,
            f,
            lay,
            k_leaf: lay.words - lay.children,
            n_items: 0,
            nodes: Vec::new(),
            roots: Vec::new(),
            seed: Kiss64Random::DEFAULT_SEED,
            built: false,
            stats: BuildStats::default(),
        }
    }

    pub fn metric(&self) -> Metric {
        self.metric
    }

    /// `get_f()`.
    pub fn f(&self) -> usize {
        self.f
    }

    /// `_s`: bytes per node, 16 + 4f (Euclidean, Manhattan) or 12 + 4f (Angular).
    pub fn node_bytes(&self) -> usize {
        4 * self.lay.words
    }

    /// `_K`: the most item ids a leaf holds.
    pub fn leaf_capacity(&self) -> usize {
        self.k_leaf
    }

    /// `get_n_items()`.
    pub fn n_items(&self) -> usize {
        self.n_items
    }

    /// `_n_nodes`: items, tree nodes and, once built, the root copies.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len() / self.lay.words
    }

    /// `_roots`, in tree order.
    pub fn roots(&self) -> &[i32] {
        &self.roots
    }

    /// `get_n_trees()`.
    pub fn n_trees(&self) -> usize {
        self.roots.len()
    }

    pub fn stats(&self) -> &BuildStats {
        &self.stats
    }

    /// `set_seed(seed)`: the seed `build` starts its generator from. Seurat never calls it.
    pub fn set_seed(&mut self, seed: u64) {
        self.seed = seed;
    }

    #[inline]
    fn base(&self, i: usize) -> usize {
        i * self.lay.words
    }

    #[inline]
    fn n_desc(&self, i: usize) -> i32 {
        self.nodes[self.base(i)] as i32
    }

    /// Item or node `i`'s vector `v`.
    #[inline]
    pub fn v(&self, i: usize) -> &[f32] {
        let b = self.base(i) + self.lay.v;
        as_f32(&self.nodes[b..b + self.f])
    }

    /// Angular's stored `norm` (word 1), 0 for the other metrics.
    #[inline]
    fn norm(&self, i: usize) -> f32 {
        match self.metric {
            Metric::Angular => f32::from_bits(self.nodes[self.base(i) + 1]),
            _ => 0.0,
        }
    }

    #[inline]
    fn child(&self, i: usize, side: usize) -> i32 {
        self.nodes[self.base(i) + self.lay.children + side] as i32
    }

    /// `get_item(item, v)`: the item's floats (`getItemsVector`).
    pub fn get_item(&self, item: usize) -> &[f32] {
        self.v(item)
    }

    /// `add_item(item, w)`: `w` (already rounded to float, as RcppAnnoy's `addItem` copies R's
    /// doubles into a `std::vector<float>`) becomes node `item`.
    pub fn add_item(&mut self, item: usize, w: &[f32]) -> Result<(), String> {
        if self.built {
            return Err("You can't add an item to a built index".into());
        }
        if w.len() != self.f {
            return Err(format!("item {item}: {} values for f = {}", w.len(), self.f));
        }
        if item > i32::MAX as usize - 1 {
            return Err(format!("item {item}: beyond int32_t"));
        }
        // _allocate_size(item + 1): new nodes are zero
        let need = (item + 1) * self.lay.words;
        if self.nodes.len() < need {
            self.nodes.resize(need, 0);
        }
        let b = self.base(item);
        let lay = self.lay;
        let node = &mut self.nodes[b..b + lay.words];
        // D::zero_value(n) does nothing for these metrics; children[0] = children[1] = 0
        node[lay.children] = 0;
        node[lay.children + 1] = 0;
        node[0] = 1; // n_descendants
        as_f32_mut(&mut node[lay.v..]).copy_from_slice(w);
        // D::init_node: Angular stores the squared norm in the union with children[0]
        if self.metric == Metric::Angular {
            let n2 = dot(w, w);
            node[1] = n2.to_bits();
        }
        if item >= self.n_items {
            self.n_items = item + 1;
        }
        Ok(())
    }

    /// `build(q)`: `q` trees (`q = -1`: until there are twice as many nodes as items), then the
    /// root copies. Single-threaded, as RcppAnnoy builds: one generator runs through every tree.
    pub fn build(&mut self, q: i32) -> Result<(), String> {
        if self.built {
            return Err("You can't build a built index".into());
        }
        if q < -1 {
            return Err(format!("build({q}): the number of trees must be -1 or at least 0"));
        }
        // D::preprocess does nothing for these metrics. _n_nodes = _n_items: n_nodes() already is.
        debug_assert_eq!(self.n_nodes(), self.n_items);
        // The trees are built recursively; a large stack keeps pathological splits (each peels
        // off a few points) from overflowing, whatever thread the caller is on.
        let mut stats = BuildStats::default();
        let mut thread_roots: Vec<i32> = Vec::new();
        std::thread::scope(|s| {
            std::thread::Builder::new()
                .name("annoy-build".into())
                .stack_size(256 << 20)
                .spawn_scoped(s, || self.thread_build(q, &mut thread_roots, &mut stats))
                .expect("annoy: cannot start the build thread")
                .join()
                .expect("annoy: the build thread panicked");
        });
        self.roots.extend_from_slice(&thread_roots);
        // copy the roots into the last segment of the array, in tree order
        for t in 0..self.roots.len() {
            let r = self.roots[t] as usize;
            let (b, w) = (self.base(r), self.lay.words);
            self.nodes.extend_from_within(b..b + w);
        }
        stats.trees = self.roots.len();
        self.stats = stats;
        self.built = true;
        Ok(())
    }

    /// `thread_build(q, 0, policy)`: the generator is `Random(_seed + thread_idx)`.
    fn thread_build(&mut self, q: i32, thread_roots: &mut Vec<i32>, stats: &mut BuildStats) {
        let mut rng = Kiss64Random::new(self.seed.wrapping_add(0));
        loop {
            if q == -1 {
                if self.n_nodes() >= 2 * self.n_items {
                    break;
                }
            } else if thread_roots.len() >= q as usize {
                break;
            }
            let indices: Vec<i32> = (0..self.n_items).filter(|&i| self.n_desc(i) >= 1).map(|i| i as i32).collect();
            let root = self.make_tree(&indices, true, &mut rng, stats, 0);
            thread_roots.push(root);
        }
    }

    /// Appends a zeroed node and returns its id (`_allocate_size(_n_nodes + 1); _n_nodes++`).
    fn alloc_node(&mut self) -> usize {
        let id = self.n_nodes();
        self.nodes.resize(self.nodes.len() + self.lay.words, 0);
        id
    }

    /// `_make_tree(indices, is_root, random, policy)`.
    fn make_tree(&mut self, indices: &[i32], is_root: bool, rng: &mut Kiss64Random, stats: &mut BuildStats, depth: usize) -> i32 {
        if indices.len() == 1 && !is_root {
            return indices[0];
        }
        let k = self.k_leaf;
        if indices.len() <= k && (!is_root || self.n_items <= k || indices.len() == 1) {
            let item = self.alloc_node();
            let b = self.base(item);
            self.nodes[b] = if is_root { self.n_items as u32 } else { indices.len() as u32 };
            let c = b + self.lay.children;
            for (t, &j) in indices.iter().enumerate() {
                self.nodes[c + t] = j as u32;
            }
            stats.leaves += 1;
            return item as i32;
        }

        stats.max_depth = stats.max_depth.max(depth);
        // Node* m = (Node*)alloca(_s): every word of it is written below
        let mut m = vec![0u32; self.lay.words];
        let mut sides: [Vec<i32>; 2] = [Vec::with_capacity(indices.len()), Vec::with_capacity(indices.len())];
        let mut p = Centroid { v: vec![0.0; self.f], norm: 0.0 };
        let mut q = Centroid { v: vec![0.0; self.f], norm: 0.0 };
        let mut accepted = None;
        for attempt in 0..3 {
            sides[0].clear();
            sides[1].clear();
            self.create_split(indices, rng, &mut m, &mut p, &mut q, stats);
            for &j in indices {
                let side = self.side(&m, j as usize, rng, stats);
                sides[side].push(j);
            }
            if split_imbalance(sides[0].len(), sides[1].len()) < 0.95 {
                accepted = Some(attempt);
                break;
            }
        }
        match accepted {
            Some(a) => stats.accepted_at_attempt[a] += 1,
            None => {
                if split_imbalance(sides[0].len(), sides[1].len()) > 0.99 {
                    stats.random_fallbacks += 1;
                } else {
                    stats.kept_imbalanced += 1;
                }
            }
        }

        // If we didn't find a hyperplane, just randomize sides as a last option
        while split_imbalance(sides[0].len(), sides[1].len()) > 0.99 {
            stats.random_rounds += 1;
            sides[0].clear();
            sides[1].clear();
            // the normal becomes 0; the offset of the last attempt stays
            for w in &mut m[self.lay.v..] {
                *w = 0.0f32.to_bits();
            }
            for &j in indices {
                stats.draws += 1;
                sides[rng.flip()].push(j);
            }
        }

        let flip = (sides[0].len() > sides[1].len()) as usize;
        m[0] = if is_root { self.n_items as u32 } else { indices.len() as u32 };
        for side in 0..2 {
            // the smaller side first (for cache locality, as Annoy says)
            let s = side ^ flip;
            let child = self.make_tree(&sides[s], false, rng, stats, depth + 1);
            m[self.lay.children + s] = child as u32;
        }
        let item = self.alloc_node();
        let b = self.base(item);
        self.nodes[b..b + self.lay.words].copy_from_slice(&m);
        stats.split_nodes += 1;
        item as i32
    }

    /// `D::side(m, y, f, random)`: by the margin's sign, a `flip()` when it is exactly 0.
    #[inline]
    fn side(&self, m: &[u32], j: usize, rng: &mut Kiss64Random, stats: &mut BuildStats) -> usize {
        let margin = self.margin_of(m, self.v(j));
        if margin != 0.0 {
            (margin > 0.0) as usize
        } else {
            stats.zero_margins += 1;
            stats.draws += 1;
            rng.flip()
        }
    }

    /// `D::margin(n, y, f)`: `a + dot(v, y)` (Minkowski) or `dot(v, y)` (Angular), for a node
    /// held in `m`.
    #[inline]
    fn margin_of(&self, m: &[u32], y: &[f32]) -> f32 {
        let v = as_f32(&m[self.lay.v..self.lay.v + self.f]);
        match self.metric {
            Metric::Angular => dot(v, y),
            _ => f32::from_bits(m[1]) + dot(v, y),
        }
    }

    /// `D::create_split(nodes, f, s, random, m)`.
    fn create_split(&self, indices: &[i32], rng: &mut Kiss64Random, m: &mut [u32], p: &mut Centroid, q: &mut Centroid, stats: &mut BuildStats) {
        let f = self.f;
        let cosine = self.metric == Metric::Angular;
        self.two_means(indices, rng, cosine, p, q, stats);
        let v_off = self.lay.v;
        {
            let v = as_f32_mut(&mut m[v_off..v_off + f]);
            for z in 0..f {
                v[z] = p.v[z] - q.v[z];
            }
            normalize(v);
        }
        if !cosine {
            // n->a = 0.0; n->a += -n->v[z] * (p->v[z] + q->v[z]) / 2
            let v = as_f32(&m[v_off..v_off + f]);
            let mut a = 0.0f32;
            for z in 0..f {
                a += -v[z] * (p.v[z] + q.v[z]) / 2.0f32;
            }
            m[1] = a.to_bits();
        }
    }

    /// `D::init_node` for a centroid: Angular keeps the squared norm.
    #[inline]
    fn init_centroid(&self, c: &mut Centroid) {
        if self.metric == Metric::Angular {
            c.norm = dot(&c.v, &c.v);
        }
    }

    /// `Distance::distance(centroid, nodes[k], f)`.
    #[inline]
    fn centroid_distance(&self, c: &Centroid, k: usize) -> f32 {
        match self.metric {
            Metric::Euclidean => euclidean_distance(&c.v, self.v(k)),
            Metric::Manhattan => manhattan_distance(&c.v, self.v(k)),
            Metric::Angular => angular_distance(&c.v, c.norm, self.v(k), self.norm(k)),
        }
    }

    /// `two_means(nodes, f, random, cosine, p, q)`: two distinct random items, then 200 random
    /// items, each moving the nearer centroid (distances weighted by the centroid's count;
    /// neither moves on a tie).
    fn two_means(&self, indices: &[i32], rng: &mut Kiss64Random, cosine: bool, p: &mut Centroid, q: &mut Centroid, stats: &mut BuildStats) {
        const ITERATION_STEPS: usize = 200;
        let count = indices.len();
        let i = rng.index(count);
        let mut j = rng.index(count - 1);
        j += (j >= i) as usize; // ensure that i != j
        stats.draws += 2;
        p.v.copy_from_slice(self.v(indices[i] as usize));
        q.v.copy_from_slice(self.v(indices[j] as usize));
        if cosine {
            normalize(&mut p.v);
            normalize(&mut q.v);
        }
        self.init_centroid(p);
        self.init_centroid(q);
        let (mut ic, mut jc) = (1i32, 1i32);
        for _ in 0..ITERATION_STEPS {
            let k = indices[rng.index(count)] as usize;
            stats.draws += 1;
            let di = ic as f32 * self.centroid_distance(p, k);
            let dj = jc as f32 * self.centroid_distance(q, k);
            let norm = if cosine { get_norm(self.v(k)) } else { 1.0f32 };
            if !(norm > 0.0) {
                continue;
            }
            let x = self.v(k);
            if di < dj {
                let (a, b) = (ic as f32, (ic + 1) as f32);
                for z in 0..self.f {
                    p.v[z] = (p.v[z] * a + x[z] / norm) / b;
                }
                self.init_centroid(p);
                ic += 1;
            } else if dj < di {
                let (a, b) = (jc as f32, (jc + 1) as f32);
                for z in 0..self.f {
                    q.v[z] = (q.v[z] * a + x[z] / norm) / b;
                }
                self.init_centroid(q);
                jc += 1;
            }
        }
    }

    /// `D::distance(v_node, item, f)` for a query whose `norm` (Angular) is `qnorm`.
    #[inline]
    fn query_distance(&self, w: &[f32], qnorm: f32, j: usize) -> f32 {
        match self.metric {
            Metric::Euclidean => euclidean_distance(w, self.v(j)),
            Metric::Manhattan => manhattan_distance(w, self.v(j)),
            Metric::Angular => angular_distance(w, qnorm, self.v(j), self.norm(j)),
        }
    }

    /// `D::normalized_distance(d)`: `sqrt(max(d, 0))` (Euclidean, Angular), `max(d, 0)`
    /// (Manhattan).
    #[inline]
    pub fn normalized_distance(&self, d: f32) -> f32 {
        match self.metric {
            Metric::Manhattan => std_max(d, 0.0),
            _ => std_max(d, 0.0).sqrt(),
        }
    }

    /// `D::margin(nd, v, f)` for node `i` of the array.
    #[inline]
    fn node_margin(&self, i: usize, w: &[f32]) -> f32 {
        match self.metric {
            Metric::Angular => dot(self.v(i), w),
            _ => f32::from_bits(self.nodes[self.base(i) + 1]) + dot(self.v(i), w),
        }
    }

    /// `get_nns_by_vector(w, n, search_k, result, distances)` (`_get_all_nns`): the `n` nearest
    /// items to `w` found by visiting nodes until `search_k` candidates are collected
    /// (`search_k = -1`: `n` times the number of trees), appended to `ids` and `dists`.
    pub fn get_nns_by_vector(&self, w: &[f32], n: usize, search_k: i32, s: &mut SearchScratch, ids: &mut Vec<i32>, dists: &mut Vec<f32>) {
        assert_eq!(w.len(), self.f, "query of {} values for f = {}", w.len(), self.f);
        // v_node: D::zero_value (nothing), the query, D::init_node (Angular's norm)
        let qnorm = if self.metric == Metric::Angular { dot(w, w) } else { 0.0 };
        let mut search_k = search_k;
        if search_k == -1 {
            // n * _roots.size() is size_t, assigned to an int
            search_k = n.wrapping_mul(self.roots.len()) as i32;
        }
        // nns.size() < (size_t)search_k: a negative int converts to a huge size_t
        let limit = search_k as i64 as u64;
        s.q.clear();
        for &r in &self.roots {
            s.q.push((f32::INFINITY, r));
        }
        s.nns.clear();
        while (s.nns.len() as u64) < limit && !s.q.is_empty() {
            let (d, i) = s.q.top();
            s.q.pop();
            let iu = i as usize;
            let nd = self.n_desc(iu);
            if nd == 1 && iu < self.n_items {
                s.nns.push(i);
            } else if nd <= self.k_leaf as i32 {
                if nd > 0 {
                    let c = self.base(iu) + self.lay.children;
                    s.nns.extend(self.nodes[c..c + nd as usize].iter().map(|&x| x as i32));
                }
            } else {
                let margin = self.node_margin(iu, w);
                s.q.push((pq_distance(d, margin, 1), self.child(iu, 1)));
                s.q.push((pq_distance(d, margin, 0), self.child(iu, 0)));
            }
        }
        // distances for all candidates, each once: sorted by id
        s.nns.sort_unstable();
        s.dist.clear();
        let mut last = -1i32;
        for &j in &s.nns {
            if j == last {
                continue;
            }
            last = j;
            if self.n_desc(j as usize) == 1 {
                s.dist.push((self.query_distance(w, qnorm, j as usize), j));
            }
        }
        let p = n.min(s.dist.len());
        if p == 0 {
            return;
        }
        partial_sort(&mut s.dist, p, &pair_less);
        for &(d, j) in &s.dist[..p] {
            dists.push(self.normalized_distance(d));
            ids.push(j);
        }
    }

    /// The index file `save()` writes: the node array, little-endian.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.nodes.len() * 4);
        for w in &self.nodes {
            out.extend_from_slice(&w.to_le_bytes());
        }
        out
    }

    /// Writes [`to_bytes`](Self::to_bytes) to `path`.
    pub fn save(&self, path: &std::path::Path) -> std::io::Result<()> {
        if !self.built {
            return Err(std::io::Error::other("You can't save an index that hasn't been built"));
        }
        use std::io::Write;
        let mut out = std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(path)?);
        for w in &self.nodes {
            out.write_all(&w.to_le_bytes())?;
        }
        out.flush()
    }

    /// `load(filename)`: an index read from Annoy's file. The roots are found as Annoy finds
    /// them, by scanning back from the end of the file for the nodes with the most descendants
    /// (and dropping the last root itself when it precedes its copy); the item count is their
    /// `n_descendants`.
    pub fn from_bytes(metric: Metric, f: usize, bytes: &[u8]) -> Result<AnnoyIndex, LoadError> {
        let mut idx = AnnoyIndex::new(metric, f);
        let s = idx.node_bytes();
        if bytes.is_empty() {
            return Err(LoadError("Size of file is zero".into()));
        }
        if bytes.len() % s != 0 {
            return Err(LoadError(format!("Index size {} is not a multiple of the node size {s} ({metric}, f = {f})", bytes.len())));
        }
        idx.nodes = bytes.chunks_exact(4).map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect();
        let n_nodes = idx.n_nodes();
        let mut m: i32 = -1;
        let mut roots = Vec::new();
        for i in (0..n_nodes).rev() {
            let k = idx.n_desc(i);
            if m == -1 || k == m {
                roots.push(i as i32);
                m = k;
            } else {
                break;
            }
        }
        if roots.len() > 1 && idx.child(roots[0] as usize, 0) == idx.child(*roots.last().unwrap() as usize, 0) {
            roots.pop();
        }
        idx.roots = roots;
        idx.n_items = m.max(0) as usize;
        idx.built = true;
        Ok(idx)
    }

    /// The fields of node `i` as a readable summary (for reports of differences).
    pub fn describe_node(&self, i: usize) -> String {
        let nd = self.n_desc(i);
        if i < self.n_items {
            format!("item {i}")
        } else if nd as usize <= self.k_leaf {
            format!("node {i}: leaf of {nd}")
        } else {
            let a = if self.metric == Metric::Angular { String::new() } else { format!(", a = {:e}", f32::from_bits(self.nodes[self.base(i) + 1])) };
            format!("node {i}: split of {nd}, children {} and {}{a}", self.child(i, 0), self.child(i, 1))
        }
    }

    /// Which field of a node a byte offset within it falls in.
    pub fn field_at(&self, byte_in_node: usize) -> &'static str {
        let w = byte_in_node / 4;
        match (self.metric, w) {
            (_, 0) => "n_descendants",
            (Metric::Angular, 1) => "children[0] / norm",
            (Metric::Angular, 2) => "children[1]",
            (_, 1) => "a",
            (_, 2) => "children[0]",
            (_, 3) => "children[1]",
            _ => "v (or a leaf's items)",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layouts_match_annoys_structs() {
        // Minkowski::Node<int32_t, float>: n_descendants, a, children[2], v -> offsetof(v) = 16
        let e = AnnoyIndex::new(Metric::Euclidean, 30);
        assert_eq!((e.node_bytes(), e.leaf_capacity()), (136, 32));
        // Angular::Node: n_descendants, union { children[2]; norm }, v -> offsetof(v) = 12
        let a = AnnoyIndex::new(Metric::Angular, 30);
        assert_eq!((a.node_bytes(), a.leaf_capacity()), (132, 32));
        let m = AnnoyIndex::new(Metric::Manhattan, 3);
        assert_eq!((m.node_bytes(), m.leaf_capacity()), (28, 5));
    }

    #[test]
    fn split_imbalance_mixes_float_and_double() {
        assert_eq!(split_imbalance(1, 1), 0.5f32 as f64);
        // 19 of 20 on one side: 19 / (20 + 1e-9) rounds to the float just below 0.95, so a 95%
        // split passes the `< 0.95` test; 96 of 100 does not
        assert_eq!(split_imbalance(19, 1), 0.95f32 as f64);
        assert!(split_imbalance(19, 1) < 0.95);
        assert!(split_imbalance(96, 4) > 0.95);
        assert!(split_imbalance(0, 33) > 0.99);
        assert!(split_imbalance(1, 32) < 0.99);
    }

    #[test]
    fn angular_distance_is_double_then_float() {
        let x = [1.0f32, 2.0, 3.0];
        let y = [3.0f32, -1.0, 0.5];
        let (nx, ny) = (dot(&x, &x), dot(&y, &y));
        let d = angular_distance(&x, nx, &y, ny);
        let want = (2.0f64 - 2.0 * dot(&x, &y) as f64 / ((nx * ny).sqrt() as f64)) as f32;
        assert_eq!(d.to_bits(), want.to_bits());
        // zero vectors: cos is taken as 0
        assert_eq!(angular_distance(&[0.0; 3], 0.0, &y, ny), 2.0);
    }

    #[test]
    fn small_index_has_one_leaf_per_tree() {
        let mut a = AnnoyIndex::new(Metric::Euclidean, 2);
        for i in 0..3 {
            a.add_item(i, &[i as f32, 1.0]).unwrap();
        }
        a.build(2).unwrap();
        // 3 items, 2 root leaves, 2 root copies
        assert_eq!(a.n_nodes(), 7);
        assert_eq!(a.roots(), &[3, 4]);
        let mut s = SearchScratch::new();
        let (mut ids, mut d) = (Vec::new(), Vec::new());
        a.get_nns_by_vector(&[2.0, 1.0], 3, -1, &mut s, &mut ids, &mut d);
        assert_eq!(ids, vec![2, 1, 0]);
        assert_eq!(d, vec![0.0, 1.0, 2.0]);
        // Annoy's load() on this file: every node from the first root leaf on has 3
        // descendants, and the "hacky fix" drops only the last one found, so it finds 3 roots
        let b = AnnoyIndex::from_bytes(Metric::Euclidean, 2, &a.to_bytes()).unwrap();
        assert_eq!(b.n_items(), 3);
        assert_eq!(b.roots(), &[6, 5, 4]);
    }
}
