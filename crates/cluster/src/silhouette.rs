//! `bluster::approxSilhouette` (bluster 1.20.0), as the workflow's `add_silhouette_to_metadata`
//! calls it (docs/PLAN.md C1 step 6 and C3d): on the PCA embeddings' columns `dims` (1:30) and
//! `seurat_clusters`, giving the `silhouette_width` column, and the closest other cluster, which
//! the workflow drops and the port writes as `silhouette_other`.
//!
//! *Ported from bluster 1.20.0's `R/approxSilhouette.R` (Copyright Aaron Lun, GPL-3), translated
//! into Rust. NOTICE lists the changes.*
//!
//! Each cluster's centroid and within-cluster variance, then each cell's root-mean-square
//! distance to every cluster, which approximates the average distance:
//!
//! ```r
//! uclust <- sort(unique(clusters))                         # the levels present, in level order
//! centroid <- colMeans(xcurrent)                           # long double sum, / n, to double
//! clust.var[i] <- sum(colMeans(sweep(xcurrent, 2, centroid)^2))   # sum(): long double
//! D <- sqrt(colSums((tx - averaged[[i]])^2) + clust.var[i])        # colSums: long double
//! better <- other.D < other.dist[is.other]                 # strict <, clusters in level order
//! width <- (other.dist - self.dist) / pmax(other.dist, self.dist)
//! ```
//!
//! R's `colMeans` and `colSums` (`do_colsum`) and `sum` (`rsum`) accumulate in 80-bit long
//! double, emulated exactly by [`rcompat::F80`]; `x^2` is `x * x`. No BLAS is involved, so the
//! result is the same on any machine, given the same embeddings.
//!
//! The workflow's own guard: with fewer than 2 cells or fewer than 2 clusters, every width is NA
//! (`add_silhouette_to_metadata` returns before calling bluster).

use crate::find_clusters::Factor;
use rcompat::F80;

/// approxSilhouette's result, and its intermediates (what M8's instrumented runs dump as `K8_*`).
#[derive(Debug, Clone, PartialEq)]
pub struct Silhouette {
    /// The embedding columns used (1-based), after `dims[dims <= ncol(emb)]`.
    pub dims: Vec<usize>,
    /// `uclust`: the factor's codes present, in level order.
    pub uclust: Vec<u32>,
    /// Per cell: the width (NaN where not computed), and the closest other cluster as a position
    /// in `uclust` (`None` where not computed).
    pub width: Vec<f64>,
    pub other: Vec<Option<u32>>,
    /// The centroids, `uclust.len() x dims.len()`, row-major (R's `K8_centroids` is that matrix).
    pub centroids: Vec<f64>,
    /// `clust.var`, per cluster of `uclust`.
    pub cluster_var: Vec<f64>,
    /// False when the workflow's guard returned NA (fewer than 2 cells or clusters).
    pub computed: bool,
}

impl Silhouette {
    /// The closest other cluster's label for each cell (`NA` where not computed).
    pub fn other_labels(&self, f: &Factor) -> Vec<String> {
        self.other.iter().map(|o| o.map(|k| f.levels[self.uclust[k as usize] as usize].clone()).unwrap_or_else(|| "NA".into())).collect()
    }
}

/// `dims[dims <= ncol(emb)]`, refusing an empty result as the workflow does.
pub fn usable_dims(dims: &[usize], ncol: usize) -> Result<Vec<usize>, String> {
    let d: Vec<usize> = dims.iter().copied().filter(|&d| d >= 1 && d <= ncol).collect();
    if d.is_empty() {
        return Err("No pca dimensions available for silhouette calculation".into());
    }
    Ok(d)
}

/// `approxSilhouette(emb[, dims], clusters)` as `add_silhouette_to_metadata` calls it: `emb` is
/// `n x ncol`, column-major (R's `cell.embeddings`), `dims` 1-based, `clusters` the cells'
/// labels with their levels.
pub fn approx_silhouette(emb: &[f64], n: usize, ncol: usize, dims: &[usize], clusters: &Factor) -> Result<Silhouette, String> {
    assert_eq!(emb.len(), n * ncol, "embeddings: {} values for {n} x {ncol}", emb.len());
    assert_eq!(clusters.codes.len(), n, "one label per cell");
    let dims = usable_dims(dims, ncol)?;
    let nd = dims.len();
    let col = |d: usize| &emb[(dims[d] - 1) * n..dims[d] * n];
    let mut present = vec![false; clusters.levels.len()];
    for &c in &clusters.codes {
        present[c as usize] = true;
    }
    let uclust: Vec<u32> = (0..clusters.levels.len() as u32).filter(|&c| present[c as usize]).collect();
    if n < 2 || uclust.len() < 2 {
        return Ok(Silhouette { dims, uclust, width: vec![f64::NAN; n], other: vec![None; n], centroids: Vec::new(), cluster_var: Vec::new(), computed: false });
    }
    let nk = uclust.len();
    let mut pos = vec![u32::MAX; clusters.levels.len()];
    for (k, &c) in uclust.iter().enumerate() {
        pos[c as usize] = k as u32;
    }
    let cell_k: Vec<usize> = clusters.codes.iter().map(|&c| pos[c as usize] as usize).collect();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); nk];
    for (i, &k) in cell_k.iter().enumerate() {
        members[k].push(i);
    }

    // centroids and within-cluster variances
    let mut centroids = vec![0.0; nk * nd];
    let mut cluster_var = vec![0.0; nk];
    for k in 0..nk {
        let m = &members[k];
        let cnt = F80::from_u64(m.len() as u64);
        let mut colmeans_sq = Vec::with_capacity(nd);
        for d in 0..nd {
            let x = col(d);
            // colMeans(xcurrent)
            let mut s = F80::ZERO;
            for &i in m {
                s = s.add_f64(x[i]);
            }
            let c = (s / cnt).to_f64();
            centroids[k * nd + d] = c;
            // colMeans(sweep(xcurrent, 2, centroid)^2)
            let mut s2 = F80::ZERO;
            for &i in m {
                let v = x[i] - c;
                s2 = s2.add_f64(v * v);
            }
            colmeans_sq.push((s2 / cnt).to_f64());
        }
        // sum(): long double, then double (beyond the double range, +-Inf)
        cluster_var[k] = rcompat::sum(&colmeans_sq, false);
    }

    // each cell's distance to every cluster, clusters in level order
    let mut self_dist = vec![f64::INFINITY; n];
    let mut other_dist = vec![f64::INFINITY; n];
    let mut other: Vec<u32> = vec![0; n];
    for k in 0..nk {
        let cent = &centroids[k * nd..(k + 1) * nd];
        for i in 0..n {
            // colSums((tx - averaged[[i]])^2), then + clust.var[i], then sqrt
            let mut s = F80::ZERO;
            for (d, &c) in cent.iter().enumerate() {
                let v = col(d)[i] - c;
                s = s.add_f64(v * v);
            }
            let dist = (s.to_f64() + cluster_var[k]).sqrt();
            if cell_k[i] == k {
                self_dist[i] = dist;
            } else if dist < other_dist[i] {
                other_dist[i] = dist;
                other[i] = k as u32;
            }
        }
    }
    let width = (0..n)
        .map(|i| {
            let (o, s) = (other_dist[i], self_dist[i]);
            // pmax(other.dist, self.dist): the first, replaced by the second when it is larger
            let m = if s > o { s } else { o };
            (o - s) / m
        })
        .collect();
    Ok(Silhouette { dims, uclust, width, other: other.into_iter().map(Some).collect(), centroids, cluster_var, computed: true })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn factor(codes: Vec<u32>, nlev: usize) -> Factor {
        Factor { name: "seurat_clusters".into(), levels: (0..nlev).map(|v| v.to_string()).collect(), codes }
    }

    #[test]
    fn small_example_by_hand() {
        // 4 cells in 1 dimension: cluster 0 = {0, 2}, cluster 1 = {10, 12}
        let emb = [0.0, 2.0, 10.0, 12.0];
        let f = factor(vec![0, 0, 1, 1], 2);
        let s = approx_silhouette(&emb, 4, 1, &[1, 5], &f).unwrap();
        assert_eq!(s.dims, vec![1]);
        assert_eq!(s.centroids, vec![1.0, 11.0]);
        assert_eq!(s.cluster_var, vec![1.0, 1.0]);
        // cell 0: self sqrt(1 + 1), other sqrt(121 + 1)
        let (a, b) = (2f64.sqrt(), 122f64.sqrt());
        assert_eq!(s.width[0], (b - a) / b);
        assert_eq!(s.other, vec![Some(1), Some(1), Some(0), Some(0)]);
        assert_eq!(s.other_labels(&f), vec!["1", "1", "0", "0"]);
    }

    #[test]
    fn one_cluster_gives_na() {
        let f = factor(vec![1, 1, 1], 3);
        let s = approx_silhouette(&[1.0, 2.0, 3.0], 3, 1, &[1], &f).unwrap();
        assert!(!s.computed && s.width.iter().all(|w| w.is_nan()));
        assert_eq!(s.uclust, vec![1]);
        assert_eq!(s.other_labels(&f), vec!["NA"; 3]);
        assert!(approx_silhouette(&[1.0, 2.0], 2, 1, &[2], &factor(vec![0, 1], 2)).is_err());
    }

    #[test]
    fn levels_absent_from_the_cells_are_skipped_and_ties_go_to_the_first() {
        // levels 0..3, level 1 absent; cell 1 (cluster 0) is equally far from clusters 2 and 3
        let emb = [0.0, 5.0, 4.0, 6.0];
        let f = factor(vec![0, 0, 2, 3], 4);
        let s = approx_silhouette(&emb, 4, 1, &[1], &f).unwrap();
        assert_eq!(s.uclust, vec![0, 2, 3]);
        // cluster 0's centroid 2.5, var 6.25; clusters 2 and 3 single cells at 4 and 6, var 0
        assert_eq!(s.other[1], Some(1), "a tie goes to the first cluster in level order");
        assert_eq!(s.other_labels(&f)[1], "2");
    }

    #[test]
    fn long_double_sums() {
        // values whose double sum loses a unit that the long double sum keeps
        let big = 2f64.powi(53);
        let emb = [big, 1.0, 1.0, -big, 0.0, 0.0];
        let f = factor(vec![0, 0, 0, 0, 1, 1], 2);
        let s = approx_silhouette(&emb, 6, 1, &[1], &f).unwrap();
        // long double: (2^53 + 1 + 1 - 2^53) / 4 = 0.5; in double the sum would be 0
        assert_eq!(s.centroids[0], 0.5);
    }
}
