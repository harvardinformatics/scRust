//! Seurat 5.5.1's `FindClusters` without Seurat (docs/PLAN.md A7b, C1 step 5, C2, C3, C3c, C3d;
//! milestone M10): the modularity optimizer that `RunModularityClusteringCpp` runs, with its
//! Java-style generator, then `GroupSingletons`, and the metadata columns and levels
//! `FindClusters.Seurat` makes.
//!
//! The chain, as `FindClusters(obj, random.seed = 12345)` runs it on the `SCT_snn` graph with
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
//! None of it uses BLAS: the result is the same on any machine and for any thread count, given
//! the same SNN graph. Silhouette widths (M11) and neighbourhood purity (M12) are added to this
//! crate later. The standard library is the only dependency besides the workspace's crates.

pub mod clustering;
pub mod find_clusters;
pub mod input;
pub mod java_random;
pub mod network;
pub mod optimizer;
pub mod output;
pub mod params;
pub mod singletons;
pub mod validate;
pub mod vos;

pub use clustering::Clustering;
pub use find_clusters::{find_clusters, ClusterError, Factor, FindClusters, ResolutionResult};
pub use java_random::JavaRandom;
pub use network::Network;
pub use optimizer::{OptimizerParams, OptimizerResult};
pub use params::{parse_resolutions, ClusterParams};

/// What this crate reproduces, for run records.
pub const REFERENCE: &str = "Seurat 5.5.1 FindClusters: RunModularityClusteringCpp (ModularityOptimizer 1.3.0 by Waltman and van Eck, Seurat's C++ translation, JavaRandom), GroupSingletons, FindClusters.Seurat's columns and levels";
