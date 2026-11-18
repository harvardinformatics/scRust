//! The workflow's clustering chain as one library call (docs/PLAN.md C1, C3 and C3h; milestone
//! M11; purity: M12): SCTransform, the workflow's PCA guard, RunPCA, FindNeighbors, FindClusters,
//! the silhouette widths and the neighbourhood purities, from the counts, or from any stage's
//! output (`scrust cluster --sct`, `--pca`, `--snn`).
//!
//! ```text
//! obj <- SCTransform(obj, seed.use = 12345)            sctransform::start / finish
//! require_min_cells_for_pca(obj)                         the guard: stop at 50 cells or fewer
//! obj <- RunPCA(obj, seed.use = 12345)                   pca::run_pca (npcs 50)
//! obj <- FindNeighbors(obj, dims = 1:30)                 neighbors::find_neighbors
//! obj <- FindClusters(obj, random.seed = 12345)          crate::find_clusters
//! obj <- add_silhouette_to_metadata(obj)                 crate::silhouette, crate::purity
//! ```
//!
//! RunUMAP, which the workflow runs between RunPCA and FindNeighbors, is left out: it changes
//! nothing downstream (C1; M8 checked it on all 34 reference objects).
//!
//! **Memory.** Each stage is freed once the next has what it needs: the counts once
//! SCTransform's setup has kept the kept genes' rows (their per-cell totals are taken first),
//! SCTransform's working data once `scale.data` is made, `scale.data` once the PCA is done.
//! [`Observer`]s see each stage's result by reference as it is made (to write it, or record its
//! stage dumps) and keep nothing.
//!
//! **Random numbers.** R's stream (`&mut RRng`) is reseeded by SCTransform and by RunPCA, as in
//! R; Annoy and the modularity optimizer have their own generators; GroupSingletons reseeds R's
//! with `set.seed(1)`. The thread count changes no result.

use crate::chain_params::ChainParams;
use crate::find_clusters::{find_clusters, ClusterError, FindClusters};
use crate::purity::{neighbor_purity, Purity};
use crate::silhouette::{approx_silhouette, Silhouette};
use neighbors::{find_neighbors, Neighbors, NeighborsError};
use pca::{Pca, PcaError};
use rcompat::RRng;
use scio::CscMatrix;
use sctransform::{Sct, SctError};
use std::fmt;
use std::io;
use std::time::Instant;

/// Why the chain stops.
#[derive(Debug)]
pub enum ChainError {
    Params(String),
    Sct(SctError),
    /// The workflow's `require_min_cells_for_pca()`, with its message.
    Guard(String),
    Pca(PcaError),
    Neighbors(NeighborsError),
    Cluster(ClusterError),
    Silhouette(String),
    Purity(String),
    Io(io::Error),
}

impl fmt::Display for ChainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChainError::Params(m) | ChainError::Guard(m) => f.write_str(m),
            ChainError::Sct(e) => write!(f, "SCTransform: {e}"),
            ChainError::Pca(e) => write!(f, "RunPCA: {e}"),
            ChainError::Neighbors(e) => write!(f, "FindNeighbors: {e}"),
            ChainError::Cluster(e) => write!(f, "FindClusters: {e}"),
            ChainError::Silhouette(m) | ChainError::Purity(m) => write!(f, "add_silhouette_to_metadata: {m}"),
            ChainError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for ChainError {}

impl From<io::Error> for ChainError {
    fn from(e: io::Error) -> Self {
        ChainError::Io(e)
    }
}

/// The workflow's `require_min_cells_for_pca(seurat_obj, context = "", npcs = 50L)`: R stops
/// with this message when there are `npcs` cells or fewer.
pub fn pca_guard(n_cells: usize, npcs: usize) -> Result<(), ChainError> {
    if n_cells <= npcs {
        return Err(ChainError::Guard(format!(
            "[LOW_QUALITY_SAMPLE] only {n_cells} cell(s) remain - too few to compute {npcs} principal components (RunPCA and downstream UMAP/clustering require more cells than PCs). This usually means upstream QC/filtering removed nearly all cells for this sample; consider dataset-specific thresholds or excluding this sample."
        )));
    }
    Ok(())
}

/// Sees each stage's result as the chain makes it. Every method does nothing by default.
pub trait Observer {
    /// SCTransform's result (its RNG states at `vst`'s entry and SCTransform's exit included).
    fn sct(&mut self, _sct: &Sct) -> io::Result<()> {
        Ok(())
    }
    /// RunPCA's result: the input's rows (rows of the counts, in `VariableFeatures` order), the
    /// input (row-major, rows x cells), the generator's states at RunPCA's entry (`None` when the
    /// chain started from SCTransform's output, where R's state is not known) and exit.
    fn pca(&mut self, _input_rows: &[usize], _input: &[f64], _p: &Pca, _entry: Option<&[i32]>, _exit: &[i32]) -> io::Result<()> {
        Ok(())
    }
    /// FindNeighbors' result, on embeddings of `ncol` columns.
    fn neighbors(&mut self, _r: &Neighbors, _ncol: usize) -> io::Result<()> {
        Ok(())
    }
    fn clusters(&mut self, _r: &FindClusters) -> io::Result<()> {
        Ok(())
    }
    fn silhouette(&mut self, _s: &Silhouette) -> io::Result<()> {
        Ok(())
    }
    fn purity(&mut self, _p: &Purity) -> io::Result<()> {
        Ok(())
    }
    /// A stage ended (`sctransform`, `pca`, `neighbors`, `cluster`, `silhouette`, `purity`).
    fn phase(&mut self, _name: &str) {}
}

/// An observer that keeps nothing.
pub struct NoObserver;

impl Observer for NoObserver {}

/// What a run keeps beyond the result: every random start's partition (FindClusters), the full
/// residual matrix (for the `full` stage dumps; small inputs only), and purity's neighbour lists
/// (for the `standard` and `full` stage dumps).
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    pub keep_partitions: bool,
    pub keep_residuals: bool,
    pub keep_purity_lists: bool,
}

/// Wall-clock seconds of each stage the run made (0 for stages it started after).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Timings {
    pub sctransform: f64,
    pub pca: f64,
    pub neighbors: f64,
    pub cluster: f64,
    pub silhouette: f64,
    pub purity: f64,
}

/// The chain's result: the clustering and the scores. The stages before it were seen by the
/// observer and are gone.
#[derive(Debug, Clone)]
pub struct Chain {
    pub n_cells: usize,
    pub clusters: FindClusters,
    /// `None` when not asked for, or when the chain started from the SNN graph without
    /// embeddings.
    pub silhouette: Option<Silhouette>,
    /// Likewise.
    pub purity: Option<Purity>,
    /// The embeddings' columns (the PCA's components).
    pub n_components: Option<usize>,
    pub timings: Timings,
}

/// From the counts (genes x cells, every stored entry as read): the whole chain.
pub fn from_counts(counts: CscMatrix, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.check().map_err(ChainError::Params)?;
    let t0 = Instant::now();
    // SCTransform.Seurat calls set.seed(seed.use) first: R's state before it does not matter
    let mut rng = RRng::new(params.sct.seed_use.expect("checked"));
    let started = sctransform::start(&counts, &params.sct, &mut rng).map_err(ChainError::Sct)?;
    drop(counts);
    let mut sct = sctransform::finish(started, &params.sct, &mut rng, opts.keep_residuals, threads).map_err(ChainError::Sct)?;
    obs.sct(&sct)?;
    let sct_s = t0.elapsed().as_secs_f64();
    obs.phase("sctransform");
    let rows: Vec<usize> = sct.residuals.top_features.iter().map(|&k| sct.setup.genes[k]).collect();
    let n = sct.residuals.ncells;
    // only scale.data goes on: the kept genes' counts and the fits are freed here
    let scale = std::mem::take(&mut sct.residuals.scale_data);
    drop(sct);
    let entry = rng.random_seed();
    let mut chain = from_scale_data_with(scale, rows, n, Some((rng, entry)), params, threads, opts, obs)?;
    chain.timings.sctransform = sct_s;
    Ok(chain)
}

/// From SCTransform's `scale.data` (`scrust cluster --sct`): `scale` holds the variable features'
/// rows, row-major (`rows.len() x n`), `rows` their rows in the counts.
pub fn from_scale_data(scale: Vec<f64>, rows: Vec<usize>, n: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.check().map_err(ChainError::Params)?;
    from_scale_data_with(scale, rows, n, None, params, threads, opts, obs)
}

#[allow(clippy::too_many_arguments)]
fn from_scale_data_with(scale: Vec<f64>, rows: Vec<usize>, n: usize, rng: Option<(RRng, Vec<i32>)>, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    if let Some(npcs) = params.pca_guard {
        pca_guard(n, npcs)?;
    }
    let t0 = Instant::now();
    // RunPCA.default calls set.seed(seed.use) itself
    let (mut rng, entry) = match rng {
        Some((g, e)) => (g, Some(e)),
        None => (RRng::new(params.pca.seed_use.expect("checked")), None),
    };
    let nf = rows.len();
    let p = pca::run_pca(&scale, nf, n, &params.pca, &mut rng, threads).map_err(ChainError::Pca)?;
    obs.pca(&rows, &scale, &p, entry.as_deref(), &rng.random_seed())?;
    drop(scale);
    let pca_s = t0.elapsed().as_secs_f64();
    obs.phase("pca");
    let ncol = p.npcs;
    let emb = p.embeddings;
    let mut chain = from_embeddings_inner(emb, n, ncol, params, threads, opts, obs)?;
    chain.timings.pca = pca_s;
    Ok(chain)
}

/// From the PCA's embeddings (`scrust cluster --pca`): `n x ncol`, column-major.
pub fn from_embeddings(emb: Vec<f64>, n: usize, ncol: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.neighbors.check(ncol).map_err(|e| ChainError::Params(format!("FindNeighbors: {e}")))?;
    params.cluster.check().map_err(|e| ChainError::Params(format!("FindClusters: {e}")))?;
    from_embeddings_inner(emb, n, ncol, params, threads, opts, obs)
}

fn from_embeddings_inner(emb: Vec<f64>, n: usize, ncol: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    let t0 = Instant::now();
    let nb = find_neighbors(&emb, n, ncol, &params.neighbors, threads).map_err(ChainError::Neighbors)?;
    obs.neighbors(&nb, ncol)?;
    let neighbors_s = t0.elapsed().as_secs_f64();
    obs.phase("neighbors");
    let snn = nb.snn;
    let mut chain = from_snn(&snn, Some((&emb, ncol)), params, threads, opts, obs)?;
    chain.timings.neighbors = neighbors_s;
    Ok(chain)
}

/// From the SNN graph (`scrust cluster --snn`), with the embeddings (`n x ncol`, column-major)
/// when the scores are to be computed.
pub fn from_snn(snn: &CscMatrix, emb: Option<(&[f64], usize)>, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.cluster.check().map_err(|e| ChainError::Params(format!("FindClusters: {e}")))?;
    let t0 = Instant::now();
    let clusters = find_clusters(snn, &params.cluster, threads, opts.keep_partitions).map_err(ChainError::Cluster)?;
    obs.clusters(&clusters)?;
    let cluster_s = t0.elapsed().as_secs_f64();
    obs.phase("cluster");
    let n = clusters.n_cells;
    let t1 = Instant::now();
    let silhouette = match (params.stats.silhouette, emb) {
        (true, Some((e, ncol))) => {
            if e.len() != n * ncol {
                return Err(ChainError::Silhouette(format!("embeddings of {} values for {n} cells x {ncol} columns", e.len())));
            }
            let s = approx_silhouette(e, n, ncol, &params.neighbors.dims, &clusters.seurat_clusters).map_err(ChainError::Silhouette)?;
            obs.silhouette(&s)?;
            obs.phase("silhouette");
            Some(s)
        }
        _ => None,
    };
    let silhouette_s = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let purity = match (params.stats.purity, emb) {
        (true, Some((e, ncol))) => {
            if e.len() != n * ncol {
                return Err(ChainError::Purity(format!("embeddings of {} values for {n} cells x {ncol} columns", e.len())));
            }
            let st = &params.stats;
            let p = neighbor_purity(e, n, ncol, &params.neighbors.dims, &clusters.seurat_clusters, st.purity_k, st.purity_weighted, opts.keep_purity_lists, threads).map_err(ChainError::Purity)?;
            obs.purity(&p)?;
            obs.phase("purity");
            Some(p)
        }
        _ => None,
    };
    let purity_s = t2.elapsed().as_secs_f64();
    Ok(Chain {
        n_cells: n,
        clusters,
        silhouette,
        purity,
        n_components: emb.map(|(_, c)| c),
        timings: Timings { cluster: cluster_s, silhouette: silhouette_s, purity: purity_s, ..Timings::default() },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_guard_stops_at_50_cells_or_fewer() {
        assert!(pca_guard(51, 50).is_ok());
        let e = pca_guard(50, 50).unwrap_err().to_string();
        assert!(e.starts_with("[LOW_QUALITY_SAMPLE] only 50 cell(s) remain - too few to compute 50 principal components"), "{e}");
    }
}
