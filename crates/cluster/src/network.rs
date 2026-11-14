//! `ModularityOptimizer::Network`: the weighted graph the optimizer works on, its reduced
//! networks and subnetworks (docs/PLAN.md C3c), and `RunModularityClusteringCpp`'s conversion of
//! the SNN matrix (`matrixToNetwork`).
//!
//! *Ported from Seurat 5.5.1's `src/ModularityOptimizer.cpp` and `src/RModularityOptimizer.cpp`
//! (Copyright 2021 Seurat authors, MIT), translated into Rust statement by statement. NOTICE
//! lists the changes.*
//!
//! Every sum is accumulated in the C++'s order, left to right in double (`std::accumulate`
//! from `0.0`, and `+=` in loop order): a different order can change a node's or a cluster's
//! weight in the last bit, and with it a near-tie (C3i).

use crate::clustering::Clustering;
use scio::CscMatrix;

/// A network in adjacency-list form: node i's neighbours are
/// `neighbor[first_neighbor_index[i]..first_neighbor_index[i + 1]]`, with their edge weights in
/// `edge_weight`. Every edge appears once in each of its nodes' lists, so `n_edges` counts
/// each edge twice, as the C++'s `nEdges` does.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Network {
    pub n_nodes: usize,
    pub n_edges: usize,
    pub node_weight: Vec<f64>,
    pub first_neighbor_index: Vec<usize>,
    pub neighbor: Vec<u32>,
    pub edge_weight: Vec<f64>,
    pub total_edge_weight_self_links: f64,
}

impl Network {
    /// `Network(nNodes, nodeWeight, firstNeighborIndex, neighbor, edgeWeight)`: without node
    /// weights, each node weighs its total edge weight (`getTotalEdgeWeightPerNode`).
    pub fn new(n_nodes: usize, node_weight: Option<Vec<f64>>, first_neighbor_index: Vec<usize>, neighbor: Vec<u32>, edge_weight: Vec<f64>) -> Network {
        let n_edges = neighbor.len();
        let mut net = Network { n_nodes, n_edges, node_weight: Vec::new(), first_neighbor_index, neighbor, edge_weight, total_edge_weight_self_links: 0.0 };
        net.node_weight = match node_weight {
            Some(w) => w,
            None => net.total_edge_weight_per_node(),
        };
        net
    }

    /// `matrixToNetwork(node1, node2, edgeWeight1, modularityFunction, nNodes)`: the pairs with
    /// `node1 < node2` become edges, each appended to both nodes' lists in the order given.
    /// Modularity function 1 weighs each node by its total edge weight; function 2 by 1.
    pub fn from_pairs(n_nodes: usize, node1: &[u32], node2: &[u32], weight: &[f64], modularity_function: i32) -> Network {
        assert!(node1.len() == node2.len() && node1.len() == weight.len());
        let mut n_neighbors = vec![0usize; n_nodes];
        for (&a, &b) in node1.iter().zip(node2) {
            if a < b {
                n_neighbors[a as usize] += 1;
                n_neighbors[b as usize] += 1;
            }
        }
        let mut first_neighbor_index = vec![0usize; n_nodes + 1];
        let mut n_edges = 0usize;
        for i in 0..n_nodes {
            first_neighbor_index[i] = n_edges;
            n_edges += n_neighbors[i];
        }
        first_neighbor_index[n_nodes] = n_edges;
        let mut neighbor = vec![0u32; n_edges];
        let mut edge_weight2 = vec![0f64; n_edges];
        n_neighbors.fill(0);
        for ((&a, &b), &w) in node1.iter().zip(node2).zip(weight) {
            if a < b {
                let (a, b) = (a as usize, b as usize);
                let j = first_neighbor_index[a] + n_neighbors[a];
                neighbor[j] = b as u32;
                edge_weight2[j] = w;
                n_neighbors[a] += 1;
                let j = first_neighbor_index[b] + n_neighbors[b];
                neighbor[j] = a as u32;
                edge_weight2[j] = w;
                n_neighbors[b] += 1;
            }
        }
        if modularity_function == 1 {
            Network::new(n_nodes, None, first_neighbor_index, neighbor, edge_weight2)
        } else {
            Network::new(n_nodes, Some(vec![1.0; n_nodes]), first_neighbor_index, neighbor, edge_weight2)
        }
    }

    /// The network `RunModularityClusteringCpp` builds from the SNN matrix: the strictly lower
    /// triangle (`it.col() < it.row()`), walked as Eigen walks a column-major sparse matrix
    /// (column by column, entries in stored order), each entry giving `node1 = col`,
    /// `node2 = row`. The number of nodes is the larger dimension. Returns `None` when the
    /// triangle is empty (R stops: "Matrix contained no network data.  Check format.").
    pub fn from_snn(snn: &CscMatrix, modularity_function: i32) -> Option<Network> {
        let mut node1 = Vec::with_capacity(snn.nnz() / 2 + 3);
        let mut node2 = Vec::with_capacity(snn.nnz() / 2 + 3);
        let mut w = Vec::with_capacity(snn.nnz() / 2 + 3);
        for k in 0..snn.ncol {
            let (rows, vals) = snn.col(k);
            for (&r, &v) in rows.iter().zip(vals) {
                if k as u32 >= r {
                    continue;
                }
                node1.push(k as u32);
                node2.push(r);
                w.push(v);
            }
        }
        if node1.is_empty() {
            return None;
        }
        let n_nodes = snn.ncol.max(snn.nrow);
        Some(Network::from_pairs(n_nodes, &node1, &node2, &w, modularity_function))
    }

    /// `getNEdges()`: each edge once.
    pub fn n_edges_undirected(&self) -> usize {
        self.n_edges / 2
    }

    /// `getTotalEdgeWeight()`: the sum of every list entry's weight, in storage order, halved.
    pub fn total_edge_weight(&self) -> f64 {
        let mut s = 0.0;
        for &w in &self.edge_weight {
            s += w;
        }
        s / 2.0
    }

    /// `getTotalEdgeWeight(node)`.
    pub fn total_edge_weight_of(&self, node: usize) -> f64 {
        let mut s = 0.0;
        for &w in &self.edge_weight[self.first_neighbor_index[node]..self.first_neighbor_index[node + 1]] {
            s += w;
        }
        s
    }

    /// `getTotalEdgeWeightPerNode()`.
    pub fn total_edge_weight_per_node(&self) -> Vec<f64> {
        (0..self.n_nodes).map(|i| self.total_edge_weight_of(i)).collect()
    }

    /// `createReducedNetwork(clustering)`: one node per cluster, weighing the sum of its nodes'
    /// weights (added in node order); an edge between two clusters weighs the sum of the edges
    /// between their nodes (added in node order, then list order), and edges within a cluster
    /// are added to the self links. Each cluster's neighbours are listed in order of first
    /// appearance.
    pub fn create_reduced_network(&self, clustering: &Clustering) -> Network {
        let nc = clustering.n_clusters;
        let mut reduced = Network {
            n_nodes: nc,
            n_edges: 0,
            node_weight: vec![0.0; nc],
            first_neighbor_index: vec![0; nc + 1],
            neighbor: Vec::new(),
            edge_weight: Vec::new(),
            total_edge_weight_self_links: self.total_edge_weight_self_links,
        };
        let mut neighbor1 = vec![0u32; self.n_edges];
        let mut edge_weight1 = vec![0f64; self.n_edges];
        let mut neighbor2 = vec![0u32; nc.saturating_sub(1)];
        let mut edge_weight2 = vec![0f64; nc];
        let node_per_cluster = clustering.nodes_per_cluster();
        for i in 0..nc {
            let mut j = 0usize;
            for &l in &node_per_cluster[i] {
                let l = l as usize;
                reduced.node_weight[i] += self.node_weight[l];
                for m in self.first_neighbor_index[l]..self.first_neighbor_index[l + 1] {
                    let n = clustering.cluster[self.neighbor[m] as usize] as usize;
                    if n != i {
                        if edge_weight2[n] == 0.0 {
                            neighbor2[j] = n as u32;
                            j += 1;
                        }
                        edge_weight2[n] += self.edge_weight[m];
                    } else {
                        reduced.total_edge_weight_self_links += self.edge_weight[m];
                    }
                }
            }
            for k in 0..j {
                let n2 = neighbor2[k] as usize;
                neighbor1[reduced.n_edges + k] = n2 as u32;
                edge_weight1[reduced.n_edges + k] = edge_weight2[n2];
                edge_weight2[n2] = 0.0;
            }
            reduced.n_edges += j;
            reduced.first_neighbor_index[i + 1] = reduced.n_edges;
        }
        neighbor1.truncate(reduced.n_edges);
        edge_weight1.truncate(reduced.n_edges);
        reduced.neighbor = neighbor1;
        reduced.edge_weight = edge_weight1;
        reduced
    }

    /// `createSubnetworks(clustering)`: for each cluster, the network of its nodes and the edges
    /// between them (each node keeps its weight in this network; no self links).
    pub fn create_subnetworks(&self, clustering: &Clustering) -> Vec<Network> {
        let node_per_cluster = clustering.nodes_per_cluster();
        let mut subnetwork_node = vec![0u32; self.n_nodes];
        let mut subnetwork_neighbor = vec![0u32; self.n_edges];
        let mut subnetwork_edge_weight = vec![0f64; self.n_edges];
        (0..clustering.n_clusters)
            .map(|i| self.create_subnetwork(clustering, i, &node_per_cluster[i], &mut subnetwork_node, &mut subnetwork_neighbor, &mut subnetwork_edge_weight))
            .collect()
    }

    /// `createSubnetwork(clustering, cluster, node, ...)`, with the C++'s scratch arrays.
    fn create_subnetwork(&self, clustering: &Clustering, cluster: usize, node: &[u32], subnetwork_node: &mut [u32], subnetwork_neighbor: &mut [u32], subnetwork_edge_weight: &mut [f64]) -> Network {
        let n = node.len();
        let mut sub = Network { n_nodes: n, ..Network::default() };
        if n == 1 {
            sub.n_edges = 0;
            sub.node_weight = vec![self.node_weight[node[0] as usize]];
            sub.first_neighbor_index = vec![0; 2];
        } else {
            for (i, &v) in node.iter().enumerate() {
                subnetwork_node[v as usize] = i as u32;
            }
            sub.n_edges = 0;
            sub.node_weight = vec![0.0; n];
            sub.first_neighbor_index = vec![0; n + 1];
            for i in 0..n {
                let j = node[i] as usize;
                sub.node_weight[i] = self.node_weight[j];
                for k in self.first_neighbor_index[j]..self.first_neighbor_index[j + 1] {
                    if clustering.cluster[self.neighbor[k] as usize] as usize == cluster {
                        subnetwork_neighbor[sub.n_edges] = subnetwork_node[self.neighbor[k] as usize];
                        subnetwork_edge_weight[sub.n_edges] = self.edge_weight[k];
                        sub.n_edges += 1;
                    }
                }
                sub.first_neighbor_index[i + 1] = sub.n_edges;
            }
            sub.neighbor = subnetwork_neighbor[..sub.n_edges].to_vec();
            sub.edge_weight = subnetwork_edge_weight[..sub.n_edges].to_vec();
        }
        sub.total_edge_weight_self_links = 0.0;
        sub
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A path 0 - 1 - 2 plus an edge 0 - 3, given as (col, row) pairs.
    fn small() -> Network {
        Network::from_pairs(4, &[0, 0, 1], &[1, 3, 2], &[0.5, 0.25, 1.0], 1)
    }

    #[test]
    fn adjacency_lists_and_weights() {
        let n = small();
        assert_eq!(n.first_neighbor_index, vec![0, 2, 4, 5, 6]);
        assert_eq!(n.neighbor, vec![1, 3, 0, 2, 1, 0]);
        assert_eq!(n.edge_weight, vec![0.5, 0.25, 0.5, 1.0, 1.0, 0.25]);
        assert_eq!(n.node_weight, vec![0.75, 1.5, 1.0, 0.25]);
        assert_eq!(n.total_edge_weight(), 1.75);
        assert_eq!(n.n_edges_undirected(), 3);
        let m2 = Network::from_pairs(4, &[0, 0, 1], &[1, 3, 2], &[0.5, 0.25, 1.0], 2);
        assert_eq!(m2.node_weight, vec![1.0; 4]);
        // pairs that are not node1 < node2 are skipped
        let skip = Network::from_pairs(4, &[0, 0, 1, 2, 3], &[1, 3, 2, 2, 0], &[0.5, 0.25, 1.0, 9.0, 9.0], 1);
        assert_eq!(skip, n);
    }

    #[test]
    fn snn_lower_triangle_in_eigen_order() {
        // symmetric 3 x 3 with diagonal: entries (row, col)
        let m = CscMatrix { nrow: 3, ncol: 3, p: vec![0, 3, 5, 7], i: vec![0, 1, 2, 0, 1, 0, 2], x: vec![1.0, 0.5, 0.25, 0.5, 1.0, 0.25, 1.0] };
        let n = Network::from_snn(&m, 1).unwrap();
        assert_eq!(n.neighbor, vec![1, 2, 0, 0]);
        assert_eq!(n.edge_weight, vec![0.5, 0.25, 0.5, 0.25]);
        let diag = CscMatrix { nrow: 2, ncol: 2, p: vec![0, 1, 2], i: vec![0, 1], x: vec![1.0, 1.0] };
        assert!(Network::from_snn(&diag, 1).is_none());
    }

    #[test]
    fn reduced_network() {
        let n = small();
        // clusters {0, 1} and {2, 3}
        let c = Clustering { n_nodes: 4, n_clusters: 2, cluster: vec![0, 0, 1, 1] };
        let r = n.create_reduced_network(&c);
        assert_eq!(r.n_nodes, 2);
        assert_eq!(r.node_weight, vec![2.25, 1.25]);
        assert_eq!(r.neighbor, vec![1, 0]);
        // 0-3 (0.25) and 1-2 (1.0), seen from cluster 0: 0.25 first, then 1.0
        assert_eq!(r.edge_weight, vec![1.25, 1.25]);
        // the edge 0-1 twice
        assert_eq!(r.total_edge_weight_self_links, 1.0);
        let subs = n.create_subnetworks(&c);
        assert_eq!(subs.len(), 2);
        assert_eq!(subs[0].neighbor, vec![1, 0]);
        assert_eq!(subs[0].node_weight, vec![0.75, 1.5]);
        assert_eq!(subs[1].n_edges, 0);
        let one = Clustering { n_nodes: 4, n_clusters: 2, cluster: vec![0, 0, 0, 1] };
        let subs = n.create_subnetworks(&one);
        assert_eq!(subs[1].n_nodes, 1);
        assert_eq!(subs[1].node_weight, vec![0.25]);
        assert_eq!(subs[1].first_neighbor_index, vec![0, 0]);
    }
}
