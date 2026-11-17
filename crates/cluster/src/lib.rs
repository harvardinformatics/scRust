//! The workflow's clustering without Seurat (docs/PLAN.md A7b, C1, C2, C3; milestones M10, M11
//! and M12): Seurat 5.5.1's `FindClusters` (the modularity optimizer that
//! `RunModularityClusteringCpp` runs, with its Java-style generator, then `GroupSingletons`, and the
//! metadata columns and levels `FindClusters.Seurat` makes), bluster's `approxSilhouette`
//! ([`silhouette`]) and `neighborPurity` ([`purity`], on `knn::purity`'s exhaustive search with its
//! certificate), and the whole chain from the counts as one library entry point ([`chain`]:
//! SCTransform, the workflow's PCA guard, RunPCA, FindNeighbors, FindClusters, the silhouette
//! widths and the purities, through the `sctransform`, `pca`, `neighbors` and `knn` crates), with
//! its parameters ([`chain_params`]), its outputs ([`chain_output`]) and its comparison with R
//! ([`validate`], kernel-gated).
//!
//! FindClusters, as `FindClusters(obj, random.seed = 12345)` runs it on the `SCT_snn` graph with
//! its defaults (modularity function 1, resolution 0.8, algorithm 1, 10 starts, 10 iterations):
//!
//! 1. **The network** ([`Network::from_snn`]): the SNN graph's strictly lower triangle, walked
//!    in stored order; each node weighs its total edge weight (function 1) or 1 (function 2).
//! 2. **The optimizer** ([`optimizer::run`]): one [`JavaRandom`] from the seed serves every
//!    start; each start runs Louvain iterations from singletons ([`vos`]); the first start with
//!    the highest modularity wins, and its clusters are renumbered by decreasing size.
//! 3. **GroupSingletons** ([`singletons::group_singletons`]): each cluster of one cell joins
//!    the cluster with the highest mean SNN weight to it, ties drawn by R's `sample()` after
//!    `set.seed(1)`.
//! 4. **The result** ([`find_clusters`]): one column per resolution (`SCT_snn_res.0.8`) with
//!    Seurat's level order, and `seurat_clusters`.
//!
//! FindClusters, the silhouette and the purity use no BLAS: their results are the same on any
//! machine and for any thread count, given the same SNN graph and embeddings. SCTransform and
//! RunPCA call OpenBLAS, so the chain from the counts reproduces R bit for bit under R's kernel.
//! The standard library is the only dependency besides the workspace's crates.

pub mod chain;
pub mod chain_output;
pub mod chain_params;
pub mod clustering;
pub mod find_clusters;
pub mod input;
pub mod java_random;
pub mod metadata;
pub mod network;
pub mod optimizer;
pub mod output;
pub mod params;
pub mod purity;
pub mod silhouette;
pub mod singletons;
pub mod validate;
pub mod vos;

pub use clustering::Clustering;
pub use find_clusters::{find_clusters, ClusterError, Factor, FindClusters, ResolutionResult};
pub use java_random::JavaRandom;
pub use network::Network;
pub use optimizer::{OptimizerParams, OptimizerResult};
pub use chain::{Chain, ChainError, Observer};
pub use chain_params::{ChainParams, Stats};
pub use params::{parse_resolutions, ClusterParams};
pub use silhouette::{approx_silhouette, Silhouette};
pub use purity::{neighbor_purity, Purity};

/// What this crate reproduces, for run records.
pub const REFERENCE: &str = "Seurat 5.5.1 FindClusters: RunModularityClusteringCpp (ModularityOptimizer 1.3.0 by Waltman and van Eck, Seurat's C++ translation, JavaRandom), GroupSingletons, FindClusters.Seurat's columns and levels";
