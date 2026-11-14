//! `RunModularityClusteringCpp`: the random starts, their iterations, the best start and the
//! numbering of its clusters (docs/PLAN.md C1 step 5 and C3c), with what Seurat prints.
//!
//! *Ported from Seurat 5.5.1's `src/RModularityOptimizer.cpp` (Copyright 2021 Seurat authors,
//! MIT), translated into Rust statement by statement. NOTICE lists the changes.*
//!
//! - The resolution is rescaled to `γ / (2·W + self links)` for modularity function 1 (W the
//!   total edge weight), and used as given for function 2.
//! - One [`JavaRandom`] seeded once serves every start; each start begins from singletons and
//!   runs up to `n_iterations` iterations, stopping after one that changes nothing (the smart
//!   local moving algorithm never reports a change, so it always runs them all).
//! - The start with the highest modularity wins; the comparison is a strict `>`, so the first of
//!   equal starts wins.
//! - Its clusters are renumbered by decreasing size (`orderClustersByNNodes`).

use crate::clustering::Clustering;
use crate::java_random::JavaRandom;
use crate::network::Network;
use crate::vos;
use std::fmt;
use std::time::Instant;

/// The arguments of `RunModularityClusteringCpp` that FindClusters passes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizerParams {
    /// `modularityFunction`: 1 (standard) or 2 (alternative).
    pub modularity_function: i32,
    pub resolution: f64,
    /// 1 (Louvain), 2 (Louvain with multilevel refinement) or 3 (smart local moving).
    pub algorithm: i32,
    pub n_random_starts: i32,
    pub n_iterations: i32,
    /// `random.seed`, as the C++ `int` R passes.
    pub random_seed: i32,
}

/// Why `RunModularityClusteringCpp` stops, with Seurat's message.
#[derive(Debug, Clone, PartialEq)]
pub enum OptimizerError {
    /// One of the argument checks at the top of `RunModularityClusteringCpp`.
    Argument(String),
    /// "Matrix contained no network data.  Check format.": no entry below the diagonal.
    NoNetwork,
    /// "Clustering step failed.": no start's modularity exceeded -Inf (all NaN).
    Failed,
}

impl fmt::Display for OptimizerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptimizerError::Argument(m) => f.write_str(m),
            OptimizerError::NoNetwork => f.write_str("Matrix contained no network data.  Check format."),
            OptimizerError::Failed => f.write_str("Clustering step failed."),
        }
    }
}

impl std::error::Error for OptimizerError {}

impl OptimizerParams {
    /// The checks at the top of `RunModularityClusteringCpp`, with its messages. Algorithm 4
    /// (Leiden) is FindClusters' route to leidenbase, not to this optimizer; it is refused
    /// earlier ([`crate::ClusterParams::check`]).
    pub fn check(&self) -> Result<(), OptimizerError> {
        if self.modularity_function != 1 && self.modularity_function != 2 {
            return Err(OptimizerError::Argument("Modularity parameter must be equal to 1 or 2.".into()));
        }
        if !(1..=4).contains(&self.algorithm) {
            return Err(OptimizerError::Argument("Algorithm for modularity optimization must be 1, 2, 3, or 4".into()));
        }
        if self.n_random_starts < 1 {
            return Err(OptimizerError::Argument("Have to have at least one start".into()));
        }
        if self.n_iterations < 1 {
            return Err(OptimizerError::Argument("Need at least one interation".into()));
        }
        if self.modularity_function == 2 && self.resolution > 1.0 {
            return Err(OptimizerError::Argument("error: resolution<1 for alternative modularity".into()));
        }
        Ok(())
    }
}

/// One random start.
#[derive(Debug, Clone, PartialEq)]
pub struct Start {
    /// Iterations run (`j` when the loop ends).
    pub iterations: usize,
    /// The quality function after the last iteration, as the C++ computes it.
    pub modularity: f64,
    /// The start's own partition (cluster ids before `orderClustersByNNodes`), if asked for.
    pub partition: Option<Vec<u32>>,
}

/// `RunModularityClusteringCpp`'s result.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizerResult {
    pub params: OptimizerParams,
    /// `getNNodes()` and `getNEdges()`, as printed.
    pub n_nodes: usize,
    pub n_edges: usize,
    /// The resolution the algorithms use (`resolution2`).
    pub resolution2: f64,
    pub starts: Vec<Start>,
    /// The first start reaching the maximum, 0-based.
    pub best_start: usize,
    pub max_modularity: f64,
    /// `getNClusters()` of the best start, as printed ("Number of communities").
    pub n_communities: usize,
    /// The result: the best start's clusters renumbered by decreasing size, 0-based.
    pub ids: Vec<u32>,
    /// The number of clusters after renumbering.
    pub n_clusters: usize,
    /// Wall-clock seconds of the starts (the C++'s "Elapsed time").
    pub seconds: f64,
    /// `nextInt` calls over all starts, and those whose bound was a power of two (a network,
    /// reduced network or subnetwork of 2^k nodes).
    pub draws: u64,
    pub pow2_draws: u64,
}

impl OptimizerResult {
    /// What `RunModularityClusteringCpp` prints with `printOutput` (FindClusters'
    /// `verbose = TRUE`), line by line. The last line, the elapsed time, is this run's own.
    pub fn printed(&self) -> Vec<String> {
        let algo = match self.params.algorithm {
            1 => "Louvain algorithm",
            2 => "Louvain algorithm with multilevel refinement",
            _ => "smart local moving algorithm",
        };
        let mut out = vec![
            "Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck".to_string(),
            String::new(),
            format!("Number of nodes: {}", self.n_nodes),
            format!("Number of edges: {}", self.n_edges),
            String::new(),
            format!("Running {algo}..."),
        ];
        if self.params.n_random_starts == 1 {
            if self.params.n_iterations > 1 {
                out.push(String::new());
            }
            out.push(format!("Modularity: {}", fmt4(self.max_modularity)));
        } else {
            out.push(format!("Maximum modularity in {} random starts: {}", self.params.n_random_starts, fmt4(self.max_modularity)));
        }
        out.push(format!("Number of communities: {}", self.n_communities));
        out.push(format!("Elapsed time: {} seconds", self.seconds as i64));
        out
    }

    /// The maximum modularity as Seurat prints it (`%.4f`).
    pub fn max_modularity_printed(&self) -> String {
        fmt4(self.max_modularity)
    }
}

/// C's `%.4f` (Rust's fixed-precision formatting rounds the exact binary value as glibc's
/// printf does; rcompat's tests check this).
pub fn fmt4(x: f64) -> String {
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan".into() } else { "nan".into() };
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    format!("{x:.4}")
}

/// `RunModularityClusteringCpp(SNN, ...)` on the network built from the SNN graph
/// ([`Network::from_snn`]). With `keep_partitions`, each start's own partition is kept.
pub fn run(network: &Network, params: &OptimizerParams, keep_partitions: bool) -> Result<OptimizerResult, OptimizerError> {
    params.check()?;
    let resolution2 = if params.modularity_function == 1 {
        params.resolution / (2.0 * network.total_edge_weight() + network.total_edge_weight_self_links)
    } else {
        params.resolution
    };
    let begin = Instant::now();
    let mut best: Option<Clustering> = None;
    let mut best_start = 0usize;
    let mut max_modularity = f64::NEG_INFINITY;
    let mut random = JavaRandom::from_r_seed(params.random_seed);
    let mut starts = Vec::with_capacity(params.n_random_starts as usize);
    for i in 0..params.n_random_starts as usize {
        let mut clustering = Clustering::singletons(network.n_nodes);
        let mut j = 0usize;
        let mut update;
        let mut modularity;
        loop {
            update = match params.algorithm {
                1 => vos::run_louvain(network, &mut clustering, resolution2, &mut random),
                2 => vos::run_louvain_multilevel(network, &mut clustering, resolution2, &mut random),
                3 => {
                    vos::run_smart_local_moving(network, &mut clustering, resolution2, &mut random);
                    true
                }
                // algorithm 4 runs no algorithm in the C++ (FindClusters never sends it here)
                _ => true,
            };
            j += 1;
            modularity = vos::quality(network, &clustering, resolution2);
            if !(j < params.n_iterations as usize && update) {
                break;
            }
        }
        starts.push(Start { iterations: j, modularity, partition: keep_partitions.then(|| clustering.cluster.clone()) });
        if modularity > max_modularity {
            best = Some(clustering);
            best_start = i;
            max_modularity = modularity;
        }
    }
    let seconds = begin.elapsed().as_secs_f64();
    let Some(mut clustering) = best else {
        return Err(OptimizerError::Failed);
    };
    let n_communities = clustering.n_clusters;
    clustering.order_clusters_by_n_nodes();
    let (draws, pow2_draws) = random.draws();
    Ok(OptimizerResult {
        params: *params,
        n_nodes: network.n_nodes,
        n_edges: network.n_edges_undirected(),
        resolution2,
        starts,
        best_start,
        max_modularity,
        n_communities,
        n_clusters: clustering.n_clusters,
        ids: clustering.cluster,
        seconds,
        draws,
        pow2_draws,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> OptimizerParams {
        OptimizerParams { modularity_function: 1, resolution: 0.8, algorithm: 1, n_random_starts: 10, n_iterations: 10, random_seed: 12345 }
    }

    #[test]
    fn checks_with_seurats_messages() {
        let ok = params();
        assert!(ok.check().is_ok());
        let bad = |p: OptimizerParams| p.check().unwrap_err().to_string();
        assert_eq!(bad(OptimizerParams { modularity_function: 3, ..ok }), "Modularity parameter must be equal to 1 or 2.");
        assert_eq!(bad(OptimizerParams { algorithm: 5, ..ok }), "Algorithm for modularity optimization must be 1, 2, 3, or 4");
        assert_eq!(bad(OptimizerParams { n_random_starts: 0, ..ok }), "Have to have at least one start");
        assert_eq!(bad(OptimizerParams { n_iterations: 0, ..ok }), "Need at least one interation");
        assert_eq!(bad(OptimizerParams { modularity_function: 2, resolution: 1.2, ..ok }), "error: resolution<1 for alternative modularity");
        assert!(OptimizerParams { modularity_function: 2, resolution: 1.0, ..ok }.check().is_ok());
    }

    #[test]
    fn printed_like_seurat() {
        let net = Network::from_pairs(4, &[0, 2], &[1, 3], &[1.0, 1.0], 1);
        let r = run(&net, &params(), true).unwrap();
        let p = r.printed();
        assert_eq!(p[2], "Number of nodes: 4");
        assert_eq!(p[3], "Number of edges: 2");
        assert_eq!(p[5], "Running Louvain algorithm...");
        assert_eq!(p[6], "Maximum modularity in 10 random starts: 0.6000");
        assert_eq!(p[7], "Number of communities: 2");
        assert_eq!(r.ids, vec![0, 0, 1, 1]);
        assert_eq!(r.starts.len(), 10);
        let one = run(&net, &OptimizerParams { n_random_starts: 1, ..params() }, false).unwrap();
        assert_eq!(&one.printed()[6..8], &["".to_string(), "Modularity: 0.6000".to_string()]);
        let one1 = run(&net, &OptimizerParams { n_random_starts: 1, n_iterations: 1, ..params() }, false).unwrap();
        assert_eq!(one1.printed()[6], "Modularity: 0.6000");
        assert_eq!(fmt4(-0.00001), "-0.0000");
        assert_eq!(fmt4(0.33019971131680448), "0.3302");
    }

    #[test]
    fn smart_local_moving_runs_every_iteration() {
        let net = Network::from_pairs(4, &[0, 2], &[1, 3], &[1.0, 1.0], 1);
        let r = run(&net, &OptimizerParams { algorithm: 3, n_iterations: 7, ..params() }, false).unwrap();
        assert!(r.starts.iter().all(|s| s.iterations == 7));
        let r = run(&net, &OptimizerParams { algorithm: 1, n_iterations: 7, ..params() }, false).unwrap();
        assert!(r.starts.iter().all(|s| s.iterations == 2));
    }
}
