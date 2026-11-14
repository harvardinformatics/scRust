//! `ModularityOptimizer::Clustering`: a partition of a network's nodes (docs/PLAN.md C3c).
//!
//! *Ported from Seurat 5.5.1's `src/ModularityOptimizer.cpp` (Copyright 2021 Seurat authors, MIT),
//! translated into Rust statement by statement. NOTICE lists the changes.*

/// The cluster of each node. As in the C++, `n_clusters` is a bound on the cluster ids that the
/// algorithms maintain themselves: after a local-moving pass it is the number of non-empty
/// clusters, numbered 0 to `n_clusters - 1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Clustering {
    pub n_nodes: usize,
    pub n_clusters: usize,
    pub cluster: Vec<u32>,
}

impl Clustering {
    /// `Clustering(int nNodes)`: every node in cluster 0, `nClusters = 1`.
    pub fn new(n_nodes: usize) -> Clustering {
        Clustering { n_nodes, n_clusters: 1, cluster: vec![0; n_nodes] }
    }

    /// `Clustering(nNodes)` followed by `initSingletonClusters()`, as the constructor of
    /// `VOSClusteringTechnique` makes it: node i in cluster i.
    pub fn singletons(n_nodes: usize) -> Clustering {
        let mut c = Clustering::new(n_nodes);
        c.init_singleton_clusters();
        c
    }

    /// `initSingletonClusters()`.
    pub fn init_singleton_clusters(&mut self) {
        for (i, c) in self.cluster.iter_mut().enumerate() {
            *c = i as u32;
        }
        self.n_clusters = self.n_nodes;
    }

    /// `getNNodesPerCluster()`.
    pub fn n_nodes_per_cluster(&self) -> Vec<usize> {
        let mut n = vec![0usize; self.n_clusters];
        for &c in &self.cluster {
            // .at(): the C++ throws std::out_of_range here; the algorithms never let it happen
            n[c as usize] += 1;
        }
        n
    }

    /// `getNodesPerCluster()`: each cluster's nodes, in increasing order.
    pub fn nodes_per_cluster(&self) -> Vec<Vec<u32>> {
        let counts = self.n_nodes_per_cluster();
        let mut out: Vec<Vec<u32>> = counts.iter().map(|&c| Vec::with_capacity(c)).collect();
        for (i, &c) in self.cluster.iter().enumerate() {
            out[c as usize].push(i as u32);
        }
        out
    }

    /// `orderClustersByNNodes()`: renumbers the clusters by decreasing size. `std::stable_sort`
    /// keeps clusters of equal size in their previous order (Rust's `sort_by` is stable too, and
    /// a stable sort's result is unique). Empty clusters are dropped: `nClusters` becomes the
    /// number of non-empty ones (but at least one).
    pub fn order_clusters_by_n_nodes(&mut self) {
        let counts = self.n_nodes_per_cluster();
        let mut pairs: Vec<(usize, usize)> = counts.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        // comparator: b.first < a.first (descending by size)
        pairs.sort_by(|a, b| b.0.cmp(&a.0));
        let mut new_cluster = vec![0u32; self.n_clusters];
        let mut i = 0usize;
        loop {
            new_cluster[pairs[i].1] = i as u32;
            i += 1;
            if !(i < self.n_clusters && pairs[i].0 > 0) {
                break;
            }
        }
        self.n_clusters = i;
        for c in self.cluster.iter_mut() {
            *c = new_cluster[*c as usize];
        }
    }

    /// `mergeClusters(clustering)`: node i moves to `clustering.cluster[cluster[i]]`, the
    /// partition of a reduced network whose nodes are this clustering's clusters.
    pub fn merge_clusters(&mut self, other: &Clustering) {
        for c in self.cluster.iter_mut() {
            *c = other.cluster[*c as usize];
        }
        self.n_clusters = other.n_clusters;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_by_size_is_stable_and_drops_empty_clusters() {
        // sizes: cluster 0 -> 1, 1 -> 3, 2 -> 0, 3 -> 3, 4 -> 1
        let mut c = Clustering { n_nodes: 8, n_clusters: 5, cluster: vec![1, 3, 0, 1, 3, 4, 1, 3] };
        c.order_clusters_by_n_nodes();
        // 1 (3 nodes) -> 0, 3 (3) -> 1, 0 (1) -> 2, 4 (1) -> 3; 2 is empty
        assert_eq!(c.cluster, vec![0, 1, 2, 0, 1, 3, 0, 1]);
        assert_eq!(c.n_clusters, 4);
    }

    #[test]
    fn merge_and_lists() {
        let mut c = Clustering::singletons(4);
        let reduced = Clustering { n_nodes: 4, n_clusters: 2, cluster: vec![1, 0, 1, 0] };
        c.merge_clusters(&reduced);
        assert_eq!(c.cluster, vec![1, 0, 1, 0]);
        assert_eq!(c.n_clusters, 2);
        assert_eq!(c.nodes_per_cluster(), vec![vec![1, 3], vec![0, 2]]);
        assert_eq!(c.n_nodes_per_cluster(), vec![2, 2]);
    }
}
