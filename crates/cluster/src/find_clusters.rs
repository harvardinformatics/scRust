//! FindClusters' result (docs/PLAN.md C1 step 5 and C3d): for each resolution, the optimizer
//! and GroupSingletons, then the metadata columns and their levels as `FindClusters.Seurat`
//! makes them, and `seurat_clusters`.
//!
//! *Ported from Seurat 5.5.1's `FindClusters.default` and `FindClusters.Seurat`
//! (`R/clustering.R`, Copyright 2021 Seurat authors, MIT), translated into Rust. NOTICE lists the
//! changes.*
//!
//! - **One optimizer run per resolution,** each with its own generator seeded from
//!   `random.seed`, and GroupSingletons after it (which seeds R's generator itself). The
//!   resolutions are therefore independent, and run in parallel here; Seurat runs them in turn.
//!   The network is the same for every resolution and is built once.
//! - **Columns:** `paste0("res.", r)` (R's 15 significant digits), prefixed with the graph name.
//! - **Each column's levels:** the labels that occur, the numeric ones (`^[0-9]+$`) sorted by
//!   value first, then the others sorted (only `"singleton"` occurs).
//! - **`seurat_clusters`:** the last column's labels. Its levels are the same levels sorted as
//!   numbers when they all are numbers, else as text (R's `as.numeric()` warns on
//!   `"singleton"`, and Seurat then keeps them as text: `"10"` sorts before `"2"`).

use crate::network::Network;
use crate::optimizer::{self, OptimizerError, OptimizerResult};
use crate::params::ClusterParams;
use crate::singletons::{self, GroupError, Grouped};
use scio::CscMatrix;
use std::fmt;
use std::sync::Mutex;
use std::time::Instant;

/// Why FindClusters stops.
#[derive(Debug, Clone, PartialEq)]
pub enum ClusterError {
    /// The arguments (refused here, or by Seurat's own checks).
    Params(String),
    /// The SNN graph cannot be used (shape, order, values).
    Graph(String),
    Optimizer(OptimizerError),
    Group(GroupError),
}

impl fmt::Display for ClusterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClusterError::Params(m) | ClusterError::Graph(m) => f.write_str(m),
            ClusterError::Optimizer(e) => write!(f, "RunModularityClusteringCpp: {e}"),
            ClusterError::Group(e) => write!(f, "GroupSingletons: {e}"),
        }
    }
}

impl std::error::Error for ClusterError {}

/// A factor: labels as codes into its levels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Factor {
    pub name: String,
    pub levels: Vec<String>,
    /// 0-based codes into `levels`, one per cell.
    pub codes: Vec<u32>,
}

impl Factor {
    /// The cells' labels as text.
    pub fn labels(&self) -> Vec<&str> {
        self.codes.iter().map(|&c| self.levels[c as usize].as_str()).collect()
    }

    /// A factor of `labels` with the given levels (every label must be one of them).
    fn new(name: &str, labels: &[String], levels: Vec<String>) -> Factor {
        let codes = labels.iter().map(|l| levels.iter().position(|v| v == l).expect("a label outside the levels") as u32).collect();
        Factor { name: name.to_string(), levels, codes }
    }
}

/// One resolution of a FindClusters call.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolutionResult {
    pub resolution: f64,
    /// `res.<r>`.
    pub res_name: String,
    pub optimizer: OptimizerResult,
    pub grouped: Grouped,
    /// The metadata column, with its levels.
    pub column: Factor,
    /// Wall-clock seconds of GroupSingletons.
    pub group_seconds: f64,
}

impl ResolutionResult {
    /// What FindClusters prints for this resolution (the optimizer's lines), and its message
    /// (GroupSingletons').
    pub fn printed(&self) -> Vec<String> {
        self.optimizer.printed()
    }

    pub fn message(&self) -> Option<String> {
        self.grouped.message()
    }
}

/// FindClusters' result.
#[derive(Debug, Clone, PartialEq)]
pub struct FindClusters {
    pub params: ClusterParams,
    pub n_cells: usize,
    /// The network, as `Number of nodes` and `Number of edges` print it.
    pub n_nodes: usize,
    pub n_edges: usize,
    /// The network's total edge weight (`getTotalEdgeWeight()`).
    pub total_edge_weight: f64,
    /// The resolutions, in the order given.
    pub resolutions: Vec<ResolutionResult>,
    /// `seurat_clusters` (and the identities).
    pub seurat_clusters: Factor,
    /// Wall-clock seconds: building the network, and all resolutions.
    pub network_seconds: f64,
    pub resolutions_seconds: f64,
}

impl FindClusters {
    /// R's `.Random.seed` after FindClusters, as far as FindClusters sets it: after the last
    /// resolution whose GroupSingletons drew a number (`None` if none did, when the state is
    /// whatever it was before).
    pub fn rng_after(&self) -> Option<&Vec<i32>> {
        self.resolutions.iter().rev().find_map(|r| r.grouped.rng_after.as_ref())
    }

    /// The cluster ids `write_cluster_metadata` writes: the identities' levels.
    pub fn cluster_ids(&self) -> &[String] {
        &self.seurat_clusters.levels
    }
}

/// Is the label made of digits only (`grepl("^[0-9]+$", x)`)?
fn is_int_label(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// The value of a label of digits (R's `as.integer`/`as.numeric`; labels are cluster ids).
fn label_value(s: &str) -> u64 {
    s.parse().unwrap_or(u64::MAX)
}

/// FindClusters.Seurat's levels of a `res.` column: the numeric labels by value, then the others
/// sorted (as text; only "singleton" occurs, so the locale does not matter).
pub fn res_levels(labels: &[String]) -> Vec<String> {
    let mut uniq: Vec<String> = labels.to_vec();
    uniq.sort_unstable();
    uniq.dedup();
    let (mut num, mut other): (Vec<String>, Vec<String>) = uniq.into_iter().partition(|s| is_int_label(s));
    num.sort_by_key(|s| label_value(s));
    other.sort();
    num.extend(other);
    num
}

/// The levels of `seurat_clusters`: `sort(as.numeric(levels))` when every level is a number,
/// else `sort(levels)` as text (R's `as.numeric()` warns, and Seurat keeps the text).
pub fn seurat_clusters_levels(levels: &[String]) -> Vec<String> {
    let mut lv = levels.to_vec();
    if lv.iter().all(|s| is_int_label(s)) {
        lv.sort_by_key(|s| label_value(s));
    } else {
        // character sort; the labels are digits and "singleton", which every locale's collation
        // orders as their bytes do
        lv.sort();
    }
    lv
}

/// Checks the SNN graph the optimizer and GroupSingletons read: square, in dgCMatrix order
/// (columns in turn, rows strictly increasing within each), every value finite, and the
/// strictly lower triangle (the network's edges) positive. Seurat's graphs always are; the
/// optimizer's bookkeeping (`edgeWeightPerCluster[l] == 0` marks a cluster not yet seen) is
/// only defined for positive weights.
pub fn check_graph(snn: &CscMatrix) -> Result<(), ClusterError> {
    if snn.nrow != snn.ncol {
        return Err(ClusterError::Graph(format!("the SNN graph is {} x {}, not square", snn.nrow, snn.ncol)));
    }
    snn.validate().map_err(|e| ClusterError::Graph(format!("the SNN graph is not a valid dgCMatrix: {e}")))?;
    for col in 0..snn.ncol {
        let (r, x) = snn.col(col);
        for (&ri, &xi) in r.iter().zip(x) {
            if !xi.is_finite() {
                return Err(ClusterError::Graph(format!("SNN entry ({}, {}) is {xi}", ri + 1, col + 1)));
            }
            if (ri as usize) > col && xi <= 0.0 {
                return Err(ClusterError::Graph(format!("SNN entry ({}, {}) is {xi}: the optimizer needs positive edge weights", ri + 1, col + 1)));
            }
        }
    }
    Ok(())
}

/// `FindClusters(object, graph.name, ...)` on the SNN graph `snn` (cells × cells, in cell
/// order): the library entry point behind `scrust cluster` and the clustering chain. With
/// `keep_partitions`, every start's own partition is kept. Resolutions run on up to `threads`
/// threads; the result does not depend on it.
pub fn find_clusters(snn: &CscMatrix, params: &ClusterParams, threads: usize, keep_partitions: bool) -> Result<FindClusters, ClusterError> {
    params.check().map_err(ClusterError::Params)?;
    check_graph(snn)?;
    let n = snn.ncol;
    let t0 = Instant::now();
    let network = Network::from_snn(snn, params.modularity_fxn).ok_or(ClusterError::Optimizer(OptimizerError::NoNetwork))?;
    let network_seconds = t0.elapsed().as_secs_f64();

    let t1 = Instant::now();
    let nres = params.resolutions.len();
    let one = |k: usize| -> Result<ResolutionResult, ClusterError> {
        let r = params.resolutions[k];
        let opt = optimizer::run(&network, &params.optimizer(r), keep_partitions).map_err(ClusterError::Optimizer)?;
        let tg = Instant::now();
        let grouped = singletons::group_singletons(&opt.ids, snn, params.group_singletons).map_err(ClusterError::Group)?;
        let group_seconds = tg.elapsed().as_secs_f64();
        let labels = grouped.label_strings();
        let column = Factor::new(&params.column_name(r), &labels, res_levels(&labels));
        Ok(ResolutionResult { resolution: r, res_name: ClusterParams::res_name(r), optimizer: opt, grouped, column, group_seconds })
    };
    let mut results: Vec<Option<Result<ResolutionResult, ClusterError>>> = (0..nres).map(|_| None).collect();
    let threads = threads.max(1).min(nres);
    if threads == 1 {
        for (k, slot) in results.iter_mut().enumerate() {
            *slot = Some(one(k));
        }
    } else {
        let next = Mutex::new(0usize);
        let slots: Vec<Mutex<&mut Option<Result<ResolutionResult, ClusterError>>>> = results.iter_mut().map(Mutex::new).collect();
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| loop {
                    let k = {
                        let mut g = next.lock().unwrap();
                        let k = *g;
                        *g += 1;
                        k
                    };
                    if k >= nres {
                        break;
                    }
                    let r = one(k);
                    **slots[k].lock().unwrap() = Some(r);
                });
            }
        });
    }
    // the first failing resolution in call order, as R would stop there
    let mut resolutions = Vec::with_capacity(nres);
    for r in results {
        resolutions.push(r.expect("every resolution ran")?);
    }
    let resolutions_seconds = t1.elapsed().as_secs_f64();

    let last = &resolutions.last().expect("at least one resolution").column;
    let labels: Vec<String> = last.labels().into_iter().map(String::from).collect();
    let seurat_clusters = Factor::new("seurat_clusters", &labels, seurat_clusters_levels(&last.levels));
    Ok(FindClusters {
        params: params.clone(),
        n_cells: n,
        n_nodes: network.n_nodes,
        n_edges: network.n_edges_undirected(),
        total_edge_weight: network.total_edge_weight(),
        resolutions,
        seurat_clusters,
        network_seconds,
        resolutions_seconds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn levels_as_seurat_sorts_them() {
        let labels = s(&["10", "2", "0", "singleton", "1", "2", "10"]);
        assert_eq!(res_levels(&labels), s(&["0", "1", "2", "10", "singleton"]));
        assert_eq!(seurat_clusters_levels(&res_levels(&labels)), s(&["0", "1", "10", "2", "singleton"]));
        let numeric = s(&["0", "1", "2", "10"]);
        assert_eq!(seurat_clusters_levels(&numeric), numeric);
    }

    #[test]
    fn graph_checks() {
        let ok = CscMatrix { nrow: 2, ncol: 2, p: vec![0, 2, 4], i: vec![0, 1, 0, 1], x: vec![1.0, 0.5, 0.5, 1.0] };
        assert!(check_graph(&ok).is_ok());
        let zero = CscMatrix { x: vec![1.0, 0.0, 0.0, 1.0], ..ok.clone() };
        assert!(check_graph(&zero).is_err());
        let nan = CscMatrix { x: vec![1.0, 0.5, f64::NAN, 1.0], ..ok.clone() };
        assert!(check_graph(&nan).is_err());
        let rect = CscMatrix { nrow: 3, ..ok.clone() };
        assert!(check_graph(&rect).is_err());
        // no edge below the diagonal: R stops
        let diag = CscMatrix { nrow: 2, ncol: 2, p: vec![0, 1, 2], i: vec![0, 1], x: vec![1.0, 1.0] };
        let e = find_clusters(&diag, &ClusterParams::workflow(), 1, false).unwrap_err();
        assert_eq!(e.to_string(), "RunModularityClusteringCpp: Matrix contained no network data.  Check format.");
    }
}
