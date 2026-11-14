//! `ModularityOptimizer::VOSClusteringTechnique`: the quality function and the three
//! algorithms FindClusters offers (docs/PLAN.md C3c): local moving, the Louvain algorithm
//! (algorithm 1), Louvain with multilevel refinement (2) and smart local moving (3).
//!
//! *Ported from Seurat 5.5.1's `src/ModularityOptimizer.cpp` (Copyright 2021 Seurat authors, MIT),
//! itself a translation of Ludo Waltman and Nees Jan van Eck's ModularityOptimizer 1.3.0 (Java);
//! translated into Rust statement by statement. NOTICE lists the changes.*
//!
//! The C++ keeps a technique's network and clustering behind `shared_ptr`s; here the functions
//! take the network by reference and the clustering by mutable reference, and the nested
//! techniques the C++ creates (on a reduced network or a subnetwork, starting from singletons)
//! are local values. Every comparison and every sum is the C++'s, in its order:
//!
//! - a node joins the neighbouring cluster with the largest gain `w(j→l) − k_j·K_l·γ`, the
//!   clusters taken in order of first appearance in its neighbour list; a strictly larger gain
//!   wins, and an equal gain wins only for a lower cluster index;
//! - when no gain is positive, the node takes the last unused cluster (its own, if it just
//!   emptied it);
//! - cluster weights change by `+=` and `−=` in visit order;
//! - the quality function adds the internal edge weights node by node, then subtracts each
//!   cluster's `K_c·K_c·γ` in cluster order. Two partitions that differ only in their numbering
//!   therefore give modularities that can differ in the last bits (M8 found such starts).

use crate::clustering::Clustering;
use crate::java_random::JavaRandom;
use crate::network::Network;

/// `calcQualityFunction()`: the modularity of `clustering` on `network` with resolution
/// `resolution` (already scaled by `RunModularityClusteringCpp` for modularity function 1).
pub fn quality(network: &Network, clustering: &Clustering, resolution: f64) -> f64 {
    let mut quality_function = 0.0f64;
    for i in 0..network.n_nodes {
        let j = clustering.cluster[i];
        let (a, b) = (network.first_neighbor_index[i], network.first_neighbor_index[i + 1]);
        for (&nb, &w) in network.neighbor[a..b].iter().zip(&network.edge_weight[a..b]) {
            if clustering.cluster[nb as usize] == j {
                quality_function += w;
            }
        }
    }
    quality_function += network.total_edge_weight_self_links;
    let mut cluster_weight = vec![0.0f64; clustering.n_clusters];
    for i in 0..network.n_nodes {
        cluster_weight[clustering.cluster[i] as usize] += network.node_weight[i];
    }
    for &w in &cluster_weight {
        quality_function -= w * w * resolution;
    }
    quality_function /= 2.0 * network.total_edge_weight() + network.total_edge_weight_self_links;
    quality_function
}

/// `runLocalMovingAlgorithm(random)`: visits the nodes in a random permutation (n draws of
/// `nextInt(n)`), cycling until every node has been visited without moving since the last
/// move; then drops the empty clusters and renumbers the others in increasing order of their
/// ids. Returns whether any node moved. A network of one node draws nothing.
pub fn run_local_moving(network: &Network, clustering: &mut Clustering, resolution: f64, random: &mut JavaRandom) -> bool {
    let n = network.n_nodes;
    if n == 1 {
        return false;
    }
    let mut update = false;
    let mut cluster_weight = vec![0.0f64; n];
    let mut n_nodes_per_cluster = vec![0usize; n];
    for i in 0..n {
        let c = clustering.cluster[i] as usize;
        cluster_weight[c] += network.node_weight[i];
        n_nodes_per_cluster[c] += 1;
    }
    let mut n_unused_clusters = 0usize;
    let mut unused_cluster = vec![0u32; n];
    for (i, &count) in n_nodes_per_cluster.iter().enumerate() {
        if count == 0 {
            unused_cluster[n_unused_clusters] = i as u32;
            n_unused_clusters += 1;
        }
    }
    let node_permutation = random.permutation(n);
    let mut edge_weight_per_cluster = vec![0.0f64; n];
    let mut neighboring_cluster = vec![0u32; n - 1];
    let mut n_stable_nodes = 0usize;
    let mut i = 0usize;
    loop {
        let j = node_permutation[i] as usize;
        let mut n_neighboring_clusters = 0usize;
        let (a, b) = (network.first_neighbor_index[j], network.first_neighbor_index[j + 1]);
        for (&nb, &w) in network.neighbor[a..b].iter().zip(&network.edge_weight[a..b]) {
            let l = clustering.cluster[nb as usize] as usize;
            if edge_weight_per_cluster[l] == 0.0 {
                neighboring_cluster[n_neighboring_clusters] = l as u32;
                n_neighboring_clusters += 1;
            }
            edge_weight_per_cluster[l] += w;
        }

        let cj = clustering.cluster[j] as usize;
        cluster_weight[cj] -= network.node_weight[j];
        n_nodes_per_cluster[cj] -= 1;
        if n_nodes_per_cluster[cj] == 0 {
            unused_cluster[n_unused_clusters] = cj as u32;
            n_unused_clusters += 1;
        }

        // bestCluster = -1, compared as an int
        let mut best_cluster: i64 = -1;
        let mut max_quality_function = 0.0f64;
        for &l in &neighboring_cluster[..n_neighboring_clusters] {
            let l = l as usize;
            let quality_function = edge_weight_per_cluster[l] - network.node_weight[j] * cluster_weight[l] * resolution;
            if quality_function > max_quality_function || (quality_function == max_quality_function && (l as i64) < best_cluster) {
                best_cluster = l as i64;
                max_quality_function = quality_function;
            }
            edge_weight_per_cluster[l] = 0.0;
        }
        if max_quality_function == 0.0 {
            best_cluster = unused_cluster[n_unused_clusters - 1] as i64;
            n_unused_clusters -= 1;
        }
        let best = best_cluster as usize;

        cluster_weight[best] += network.node_weight[j];
        n_nodes_per_cluster[best] += 1;
        if best == cj {
            n_stable_nodes += 1;
        } else {
            clustering.cluster[j] = best as u32;
            n_stable_nodes = 1;
            update = true;
        }

        i = if i < n - 1 { i + 1 } else { 0 };
        if n_stable_nodes >= n {
            break;
        }
    }

    let mut new_cluster = vec![0u32; n];
    clustering.n_clusters = 0;
    for (i, &count) in n_nodes_per_cluster.iter().enumerate() {
        if count > 0 {
            new_cluster[i] = clustering.n_clusters as u32;
            clustering.n_clusters += 1;
        }
    }
    for c in clustering.cluster.iter_mut() {
        *c = new_cluster[*c as usize];
    }
    update
}

/// `runLouvainAlgorithm(random)`: local moving, then, if some nodes merged, the Louvain
/// algorithm on the reduced network (its nodes starting as singletons), whose result is merged
/// back when it changed anything.
pub fn run_louvain(network: &Network, clustering: &mut Clustering, resolution: f64, random: &mut JavaRandom) -> bool {
    if network.n_nodes == 1 {
        return false;
    }
    let mut update = run_local_moving(network, clustering, resolution, random);
    if clustering.n_clusters < network.n_nodes {
        let reduced = network.create_reduced_network(clustering);
        let mut reduced_clustering = Clustering::singletons(reduced.n_nodes);
        let update2 = run_louvain(&reduced, &mut reduced_clustering, resolution, random);
        if update2 {
            update = true;
            clustering.merge_clusters(&reduced_clustering);
        }
    }
    update
}

/// `runLouvainAlgorithmWithMultilevelRefinement(random)`: as [`run_louvain`], with one more
/// local-moving pass at this level after merging the reduced network's result (its return
/// value is ignored, as in the C++).
pub fn run_louvain_multilevel(network: &Network, clustering: &mut Clustering, resolution: f64, random: &mut JavaRandom) -> bool {
    if network.n_nodes == 1 {
        return false;
    }
    let mut update = run_local_moving(network, clustering, resolution, random);
    if clustering.n_clusters < network.n_nodes {
        let reduced = network.create_reduced_network(clustering);
        let mut reduced_clustering = Clustering::singletons(reduced.n_nodes);
        let update2 = run_louvain_multilevel(&reduced, &mut reduced_clustering, resolution, random);
        if update2 {
            update = true;
            clustering.merge_clusters(&reduced_clustering);
            run_local_moving(network, clustering, resolution, random);
        }
    }
    update
}

/// `runSmartLocalMovingAlgorithm(random)`: local moving; then each cluster's subnetwork is
/// split by local moving from singletons, the network is reduced over those subclusters with
/// each subcluster starting in its parent cluster, and the algorithm recurses on the reduced
/// network, whose result is merged back.
pub fn run_smart_local_moving(network: &Network, clustering: &mut Clustering, resolution: f64, random: &mut JavaRandom) -> bool {
    if network.n_nodes == 1 {
        return false;
    }
    let mut update = run_local_moving(network, clustering, resolution, random);
    if clustering.n_clusters < network.n_nodes {
        let subnetwork = network.create_subnetworks(clustering);
        let node_per_cluster = clustering.nodes_per_cluster();
        clustering.n_clusters = 0;
        let mut n_nodes_per_cluster_reduced_network = vec![0usize; subnetwork.len()];
        for (i, sub) in subnetwork.iter().enumerate() {
            let mut sub_clustering = Clustering::singletons(sub.n_nodes);
            run_local_moving(sub, &mut sub_clustering, resolution, random);
            for j in 0..sub.n_nodes {
                clustering.cluster[node_per_cluster[i][j] as usize] = (clustering.n_clusters + sub_clustering.cluster[j] as usize) as u32;
            }
            clustering.n_clusters += sub_clustering.n_clusters;
            n_nodes_per_cluster_reduced_network[i] = sub_clustering.n_clusters;
        }

        let reduced = network.create_reduced_network(clustering);
        let mut clustering2 = Clustering::singletons(reduced.n_nodes);
        let mut i = 0usize;
        for (j, &count) in n_nodes_per_cluster_reduced_network.iter().enumerate() {
            for _ in 0..count {
                clustering2.cluster[i] = j as u32;
                i += 1;
            }
        }
        clustering2.n_clusters = n_nodes_per_cluster_reduced_network.len();

        update |= run_smart_local_moving(&reduced, &mut clustering2, resolution, random);

        clustering.merge_clusters(&clustering2);
    }
    update
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles joined by one weak edge.
    fn two_triangles() -> Network {
        let n1 = [0, 0, 1, 3, 3, 4, 2];
        let n2 = [1, 2, 2, 4, 5, 5, 3];
        let w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1];
        Network::from_pairs(6, &n1, &n2, &w, 1)
    }

    #[test]
    fn finds_the_two_triangles_with_every_algorithm() {
        let net = two_triangles();
        let res = 1.0 / (2.0 * net.total_edge_weight());
        type Algo = fn(&Network, &mut Clustering, f64, &mut JavaRandom) -> bool;
        let algos: [Algo; 3] = [run_louvain, run_louvain_multilevel, run_smart_local_moving];
        for (a, f) in algos.iter().enumerate() {
            let mut c = Clustering::singletons(6);
            let mut r = JavaRandom::new(12345);
            assert!(f(&net, &mut c, res, &mut r));
            assert_eq!(c.n_clusters, 2, "algorithm {}", a + 1);
            assert_eq!(c.cluster[0], c.cluster[1]);
            assert_eq!(c.cluster[1], c.cluster[2]);
            assert_eq!(c.cluster[3], c.cluster[4]);
            assert_ne!(c.cluster[0], c.cluster[3]);
            let q = quality(&net, &c, res);
            assert!(q > 0.3 && q < 0.5, "{q}");
        }
    }

    #[test]
    fn a_single_node_draws_nothing() {
        let net = Network::from_pairs(1, &[], &[], &[], 1);
        let mut c = Clustering::singletons(1);
        let mut r = JavaRandom::new(1);
        let before = r.clone();
        assert!(!run_louvain(&net, &mut c, 1.0, &mut r));
        assert!(!run_local_moving(&net, &mut c, 1.0, &mut r));
        assert_eq!(r, before);
    }

    #[test]
    fn isolated_nodes_stay_alone() {
        // a pair and two isolated nodes
        let net = Network::from_pairs(4, &[0], &[1], &[1.0], 1);
        let mut c = Clustering::singletons(4);
        let mut r = JavaRandom::new(3);
        run_louvain(&net, &mut c, 0.8 / (2.0 * net.total_edge_weight()), &mut r);
        assert_eq!(c.n_clusters, 3);
        assert_eq!(c.cluster[0], c.cluster[1]);
        assert!(c.cluster[2] != c.cluster[3] && c.cluster[2] != c.cluster[0] && c.cluster[3] != c.cluster[0]);
    }
}
