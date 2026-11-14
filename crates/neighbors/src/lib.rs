//! Seurat 5.5.1's `FindNeighbors` without Seurat (docs/PLAN.md A7b, C1 step 4, C3, C3a, C3b;
//! milestone M9): the Annoy index RcppAnnoy 0.0.23 builds (Annoy 1.17.3), its search, the
//! `SCT_nn` graph and `ComputeSNN`'s shared-nearest-neighbour graph.
//!
//! The chain, as `FindNeighbors(obj, dims = 1:30)` runs it with its defaults:
//!
//! 1. **Items** ([`items`]): the embedding columns `dims`, each double rounded to float (to
//!    nearest even, as C++'s conversion in RcppAnnoy's `addItem` does), added in cell order.
//!    With fewer cells than `k.param`, k becomes the number of cells minus 1.
//! 2. **The forest** ([`AnnoyIndex::build`]): `n.trees` random-projection trees, one
//!    [`Kiss64Random`] stream through all of them from the default seed (Seurat never seeds it).
//! 3. **The search** ([`search`]): for each cell, its own float vector as the query, `k`
//!    neighbours with `search.k = -1`, so at least k × n.trees candidates are collected. Queries
//!    are independent and run in parallel; each result is placed by its cell.
//! 4. **The graphs** ([`graph`]): `SCT_nn`, the 0/1 kNN matrix, and `SCT_snn`, the Jaccard
//!    values s/(k + (k − s)) of shared neighbours, pruned below `prune.SNN`.
//!
//! None of it uses BLAS: every result is the same on any machine and for any thread count,
//! given the same embeddings. The standard library is the only dependency besides the
//! workspace's own crates.

pub mod annoy;
pub mod graph;
pub mod kiss;
pub mod params;
pub mod stl;
pub mod validate;

pub use annoy::{AnnoyIndex, BuildStats, Metric, SearchScratch};
pub use kiss::Kiss64Random;
pub use params::{parse_dims, parse_number, NeighborParams};
use scio::CscMatrix;
use std::fmt;
use std::sync::Mutex;
use std::time::Instant;

/// What this crate reproduces, for run records.
pub const REFERENCE: &str = "Seurat 5.5.1 FindNeighbors: RcppAnnoy 0.0.23 (Annoy 1.17.3, Kiss64Random), AnnoySearch, ComputeSNN";

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NeighborsError {
    /// The parameters do not fit the input (Seurat would stop too).
    Params(String),
    /// Annoy found fewer than k distinct neighbours for a cell; Seurat's `AnnoySearch` cannot
    /// store such a result in its n × k matrix (R stops, or recycles the list when its length
    /// divides k), so the port refuses it.
    TooFewNeighbours { cell: usize, found: usize, k: usize },
}

impl fmt::Display for NeighborsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NeighborsError::Params(m) => f.write_str(m),
            NeighborsError::TooFewNeighbours { cell, found, k } => {
                write!(f, "Annoy found only {found} of k = {k} neighbours for cell {} (Seurat's AnnoySearch would fail or recycle them); try more trees or a larger search.k", cell + 1)
            }
        }
    }
}

impl std::error::Error for NeighborsError {}

/// Wall-clock seconds of each step.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Timings {
    pub build: f64,
    pub search: f64,
    pub graphs: f64,
}

/// FindNeighbors' result.
#[derive(Debug, Clone)]
pub struct Neighbors {
    /// Cells.
    pub n: usize,
    /// k after Seurat's adjustment (`k.param`, or the number of cells minus 1).
    pub k: usize,
    /// The forest, whose [`AnnoyIndex::to_bytes`] is `annoy.ann`.
    pub index: AnnoyIndex,
    /// Each cell's k neighbours, 0-based, nearest first: `nn_idx[i * k + r]`.
    pub nn_idx: Vec<i32>,
    /// Their distances as Annoy returns them (`normalized_distance`, in float).
    pub nn_dist: Vec<f32>,
    /// `SCT_nn`.
    pub nn: CscMatrix,
    /// `SCT_snn`.
    pub snn: CscMatrix,
    pub timings: Timings,
}

/// The items Annoy holds: columns `dims` (1-based, in that order) of the embeddings, `n` rows
/// by `ncol` columns in column-major order (as R stores `cell.embeddings`), each double
/// rounded to float. Row-major: item i is `out[i * f..(i + 1) * f]`.
pub fn items(emb: &[f64], n: usize, ncol: usize, dims: &[usize]) -> Vec<f32> {
    assert_eq!(emb.len(), n * ncol, "embeddings: {} values for {n} x {ncol}", emb.len());
    let f = dims.len();
    let mut out = vec![0.0f32; n * f];
    for (z, &d) in dims.iter().enumerate() {
        assert!(d >= 1 && d <= ncol, "dimension {d} of {ncol}");
        let col = &emb[(d - 1) * n..d * n];
        for (i, &x) in col.iter().enumerate() {
            out[i * f + z] = x as f32;
        }
    }
    out
}

/// `AnnoyBuildIndex(data, metric, n.trees)`: a new index, `addItem(i - 1, data[i, ])` for every
/// row in order, `build(n.trees)`.
pub fn build_index(items: &[f32], n: usize, f: usize, metric: Metric, n_trees: i32) -> Result<AnnoyIndex, String> {
    assert_eq!(items.len(), n * f);
    let mut index = AnnoyIndex::new(metric, f);
    for i in 0..n {
        index.add_item(i, &items[i * f..(i + 1) * f])?;
    }
    index.build(n_trees)?;
    Ok(index)
}

/// `AnnoySearch(index, query, k, search.k)`: for each query row (row-major, f floats each),
/// `getNNsByVectorList(query[x, ], k, search.k, TRUE)`. Returns the ids (0-based) and distances,
/// `n_query × k` row-major, and how many neighbours each query found (at most k; unfilled slots
/// hold -1 and NaN). Runs on `threads` threads; the result does not depend on it.
pub fn search(index: &AnnoyIndex, queries: &[f32], k: usize, search_k: i32, threads: usize) -> (Vec<i32>, Vec<f32>, Vec<usize>) {
    let f = index.f();
    assert_eq!(queries.len() % f, 0);
    let nq = queries.len() / f;
    let mut ids = vec![-1i32; nq * k];
    let mut dists = vec![f32::NAN; nq * k];
    let mut counts = vec![0usize; nq];
    if nq == 0 || k == 0 {
        return (ids, dists, counts);
    }
    const CHUNK: usize = 64;
    // disjoint pieces of the outputs, handed out to the threads
    type Piece<'a> = (usize, &'a mut [i32], &'a mut [f32], &'a mut [usize]);
    let pieces: Vec<Piece> = ids
        .chunks_mut(CHUNK * k)
        .zip(dists.chunks_mut(CHUNK * k))
        .zip(counts.chunks_mut(CHUNK))
        .enumerate()
        .map(|(c, ((a, b), d))| (c * CHUNK, a, b, d))
        .rev()
        .collect();
    let work = Mutex::new(pieces);
    let run = || {
        let mut scratch = SearchScratch::new();
        let (mut ri, mut rd) = (Vec::with_capacity(k), Vec::with_capacity(k));
        loop {
            let piece = work.lock().unwrap().pop();
            let Some((start, a, b, d)) = piece else { break };
            for t in 0..d.len() {
                let q = start + t;
                ri.clear();
                rd.clear();
                index.get_nns_by_vector(&queries[q * f..(q + 1) * f], k, search_k, &mut scratch, &mut ri, &mut rd);
                d[t] = ri.len();
                a[t * k..t * k + ri.len()].copy_from_slice(&ri);
                b[t * k..t * k + rd.len()].copy_from_slice(&rd);
            }
        }
    };
    let threads = threads.max(1).min(nq.div_ceil(CHUNK));
    if threads == 1 {
        run();
    } else {
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(run);
            }
        });
    }
    (ids, dists, counts)
}

/// `FindNeighbors(object, dims, k.param, prune.SNN, n.trees, annoy.metric)` on the embeddings
/// (`n` cells by `ncol` columns, column-major): the forest, every cell's neighbours and both
/// graphs. The library entry point behind `scrust neighbors` (and, later, the clustering chain).
pub fn find_neighbors(emb: &[f64], n: usize, ncol: usize, params: &NeighborParams, threads: usize) -> Result<Neighbors, NeighborsError> {
    params.check(ncol).map_err(NeighborsError::Params)?;
    if n == 0 {
        return Err(NeighborsError::Params("no cells".into()));
    }
    if emb.len() != n * ncol {
        return Err(NeighborsError::Params(format!("embeddings: {} values for {n} x {ncol}", emb.len())));
    }
    // if (n.cells < k.param) k.param <- n.cells - 1, with a warning
    let k = if n < params.k_param { n - 1 } else { params.k_param };
    if k == 0 {
        // one cell: R's search returns no neighbours, then building SCT_nn from them fails
        return Err(NeighborsError::Params("a single cell: k.param becomes 0, and FindNeighbors cannot build its graphs".into()));
    }
    let f = params.dims.len();
    let items = items(emb, n, ncol, &params.dims);
    let t0 = Instant::now();
    let index = build_index(&items, n, f, params.metric, params.n_trees).map_err(NeighborsError::Params)?;
    let t1 = Instant::now();
    // the queries are the same rows, rounded to float again in getNNsByVectorList: the same
    // values as the items
    let (nn_idx, nn_dist, counts) = search(&index, &items, k, params.search_k, threads);
    if let Some((cell, &found)) = counts.iter().enumerate().find(|(_, &c)| c != k) {
        return Err(NeighborsError::TooFewNeighbours { cell, found, k });
    }
    let t2 = Instant::now();
    let nn = graph::nn_graph(&nn_idx, n, k);
    let snn = graph::snn_graph(&nn_idx, &nn, n, k, params.prune_snn, threads);
    let t3 = Instant::now();
    let timings = Timings { build: (t1 - t0).as_secs_f64(), search: (t2 - t1).as_secs_f64(), graphs: (t3 - t2).as_secs_f64() };
    Ok(Neighbors { n, k, index, nn_idx, nn_dist, nn, snn, timings })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn items_are_columns_rounded_to_float() {
        // 2 cells x 3 columns, column-major
        let emb = [0.1, 0.2, 1.0 / 3.0, 2.0, 5.0, 6.0];
        let it = items(&emb, 2, 3, &[3, 1]);
        assert_eq!(it, vec![5.0f32, 0.1f64 as f32, 6.0, 0.2f64 as f32]);
    }

    #[test]
    fn a_single_cell_is_refused() {
        let p = NeighborParams { dims: vec![1, 2], ..NeighborParams::default() };
        assert!(matches!(find_neighbors(&[0.5, 1.5], 1, 2, &p, 1), Err(NeighborsError::Params(_))));
    }

    #[test]
    fn fewer_cells_than_k() {
        // 5 cells: k = 4, every cell finds every other cell and itself... minus one
        let emb: Vec<f64> = (0..10).map(|v| v as f64 * 0.37).collect();
        let p = NeighborParams { dims: vec![1, 2], ..NeighborParams::default() };
        let r = find_neighbors(&emb, 5, 2, &p, 1).unwrap();
        assert_eq!(r.k, 4);
        assert_eq!(r.nn_idx.len(), 20);
        assert_eq!(r.nn.nnz(), 20);
        // every cell is its own nearest neighbour here (no duplicates)
        for i in 0..5 {
            assert_eq!(r.nn_idx[i * 4], i as i32);
        }
    }

    #[test]
    fn independent_of_thread_count() {
        let n = 700;
        let mut state = 5u64;
        let emb: Vec<f64> = (0..n * 4)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((state >> 11) as f64 / (1u64 << 53) as f64) * 10.0 - 5.0
            })
            .collect();
        let p = NeighborParams { dims: vec![1, 2, 3, 4], n_trees: 7, ..NeighborParams::default() };
        let one = find_neighbors(&emb, n, 4, &p, 1).unwrap();
        for t in [2, 5, 16] {
            let many = find_neighbors(&emb, n, 4, &p, t).unwrap();
            assert_eq!(many.nn_idx, one.nn_idx, "{t} threads");
            assert_eq!(many.nn_dist.iter().map(|d| d.to_bits()).collect::<Vec<_>>(), one.nn_dist.iter().map(|d| d.to_bits()).collect::<Vec<_>>());
            assert_eq!(many.snn, one.snn);
            assert_eq!(many.index.to_bytes(), one.index.to_bytes());
        }
    }
}
