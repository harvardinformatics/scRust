//! `bluster::neighborPurity` (bluster 1.20.0), as the workflow's `add_silhouette_to_metadata`
//! calls it (docs/PLAN.md C1 step 6 and C3e; milestone M12): on the PCA embeddings' columns
//! `dims` (1:30) and `seurat_clusters`, with `k = 50`, `weighted = TRUE` and BiocNeighbors 2.4.0's
//! default index (`BNPARAM = NULL`, KMKNN), giving the `neighborhood_purity` column, and the
//! cluster with the largest share, which the workflow drops and the port writes as
//! `purity_maximum`.
//!
//! *Ported from bluster 1.20.0's `R/neighborPurity.R` and `src/sum_neighbor_weights.cpp`
//! (Copyright Aaron Lun, GPL-3), translated into Rust. NOTICE lists the changes.*
//!
//! ```r
//! dist <- median(findDistance(idx, k = k))                         # knn::purity
//! nout <- findNeighbors(idx, threshold = dist, get.distance = FALSE)$index
//! w <- 1/table(clusters); w <- as.numeric(w[as.character(clusters)])  # 1/(cluster size)
//! uclust <- sort(unique(clusters))                                 # the levels present
//! aggregated <- sum_neighbor_weights(length(uclust), nout, m - 1L, w)
//! targets <- t(aggregated)
//! totals <- rowSums(targets)                                       # long double
//! purity <- targets[cbind(seq_along(m), m)]/totals
//! maximum <- uclust[max.col(targets, ties.method = "first")]
//! ```
//!
//! - **The search** is [`knn::purity`]'s: knncolle's arithmetic, searched exhaustively, with its
//!   certificate that KMKNN finds the same.
//! - **The sums** (`sum_neighbor_weights`): for each cell, its own cluster's entry is set to its
//!   weight, then each neighbour's weight is added to its cluster's entry, in double, in the
//!   order `findNeighbors` lists them. Every cell of a cluster has the same weight, so each
//!   entry is that weight added to 0 once per cell (the cell itself included), whatever the
//!   order: the entry for a cluster with `m` such cells is the `m`-th repeated sum, taken from a
//!   table per cluster.
//! - **The totals** are R's `rowSums`: an 80-bit long-double accumulator per row, the clusters
//!   added in level order, then rounded to double ([`rcompat::F80`]).
//! - **The maximum** is `max.col(ties.method = "first")`: the first cluster whose entry is
//!   strictly larger than every earlier one.
//!
//! The workflow's own guard: with fewer than 2 cells or fewer than 2 clusters, every purity is NA
//! (`add_silhouette_to_metadata` returns before calling bluster). No BLAS is involved, so given
//! the same embeddings the result is the same on any machine and for any thread count.

use crate::find_clusters::Factor;
use crate::silhouette::usable_dims;
use knn::purity::{Certificate, KthDistances};
use knn::Points;
use rcompat::F80;

/// neighborPurity's result, and its intermediates (what M8's instrumented runs dump as `K9_*`).
#[derive(Debug, Clone, PartialEq)]
pub struct Purity {
    /// The embedding columns used (1-based), after `dims[dims <= ncol(emb)]`.
    pub dims: Vec<usize>,
    /// `uclust`: the factor's codes present, in level order.
    pub uclust: Vec<u32>,
    /// `k` as asked for, and after BiocNeighbors' cap at n - 1.
    pub k_requested: usize,
    pub k: usize,
    pub weighted: bool,
    /// `findDistance`: each cell's distance to its k-th nearest other cell.
    pub dist_k: Vec<f64>,
    /// Their median, the radius.
    pub radius: f64,
    /// Each cell's neighbours within the radius (itself excluded).
    pub neighbor_count: Vec<u32>,
    /// The neighbours (0-based) in BiocNeighbors' order, by raw distance, ties by index, as
    /// pointers and indices; kept when asked for (the stage dumps).
    pub neighbor_lists: Option<(Vec<u64>, Vec<u32>)>,
    /// `w`: each cell's weight.
    pub weights: Vec<f64>,
    /// `targets`: cells x clusters of `uclust`, row-major.
    pub sums: Vec<f64>,
    /// `rowSums(targets)`.
    pub totals: Vec<f64>,
    /// Per cell: the purity (NaN where not computed), and the maximum as a position in `uclust`
    /// (`None` where not computed).
    pub purity: Vec<f64>,
    pub maximum: Vec<Option<u32>>,
    /// The search's certificate.
    pub certificate: Certificate,
    /// The search's margin (the certificate's constants).
    pub margin: Option<knn::purity::Margin>,
    /// False when the workflow's guard returned NA (fewer than 2 cells or clusters).
    pub computed: bool,
}

impl Purity {
    /// The maximum's label for each cell (`NA` where not computed).
    pub fn maximum_labels(&self, f: &Factor) -> Vec<String> {
        self.maximum.iter().map(|o| o.map(|k| f.levels[self.uclust[k as usize] as usize].clone()).unwrap_or_else(|| "NA".into())).collect()
    }

    /// R's `mean(lengths(nout))`: the long-double mean of the neighbour counts.
    pub fn mean_neighbors(&self) -> f64 {
        let n = self.neighbor_count.len();
        if n == 0 {
            return f64::NAN;
        }
        let s: u64 = self.neighbor_count.iter().map(|&c| c as u64).sum();
        (F80::from_u64(s) / F80::from_u64(n as u64)).to_f64()
    }
}

/// `w`'s value for a cell of a cluster of `size` cells: `1/table(clusters)` (a double division
/// of 1 by the integer count), or 1 unweighted.
fn weight(size: usize, weighted: bool) -> f64 {
    if weighted {
        1.0 / size as f64
    } else {
        1.0
    }
}

/// One row of R's `rowSums`: a long-double accumulator from 0, the columns added in order, then
/// rounded to double (`do_colsum` with `keepNA`).
pub fn row_sum(row: &[f64]) -> f64 {
    let mut acc = F80::ZERO;
    for &v in row {
        acc = acc.add_f64(v);
    }
    acc.to_f64()
}

/// `neighborPurity(emb[, dims], clusters, k, weighted)` as `add_silhouette_to_metadata` calls it:
/// `emb` is `n x ncol`, column-major (R's `cell.embeddings`), `dims` 1-based, `clusters` the
/// cells' labels with their levels. With `keep_lists`, the neighbour lists are kept (sorted as
/// BiocNeighbors sorts them).
pub fn neighbor_purity(emb: &[f64], n: usize, ncol: usize, dims: &[usize], clusters: &Factor, k: usize, weighted: bool, keep_lists: bool, threads: usize) -> Result<Purity, String> {
    assert_eq!(emb.len(), n * ncol, "embeddings: {} values for {n} x {ncol}", emb.len());
    assert_eq!(clusters.codes.len(), n, "one label per cell");
    let dims = usable_dims(dims, ncol)?;
    let mut present = vec![0usize; clusters.levels.len()];
    for &c in &clusters.codes {
        present[c as usize] += 1;
    }
    let uclust: Vec<u32> = (0..clusters.levels.len() as u32).filter(|&c| present[c as usize] > 0).collect();
    if n < 2 || uclust.len() < 2 {
        return Ok(Purity {
            dims,
            uclust,
            k_requested: k,
            k: 0,
            weighted,
            dist_k: Vec::new(),
            radius: f64::NAN,
            neighbor_count: Vec::new(),
            neighbor_lists: None,
            weights: Vec::new(),
            sums: Vec::new(),
            totals: Vec::new(),
            purity: vec![f64::NAN; n],
            maximum: vec![None; n],
            certificate: Certificate::default(),
            margin: None,
            computed: false,
        });
    }
    if emb.iter().any(|v| !v.is_finite()) {
        return Err("neighborPurity: the embeddings hold non-finite values".into());
    }
    let nk = uclust.len();
    let mut pos = vec![u32::MAX; clusters.levels.len()];
    for (p, &c) in uclust.iter().enumerate() {
        pos[c as usize] = p as u32;
    }
    let m: Vec<u32> = clusters.codes.iter().map(|&c| pos[c as usize]).collect();

    // emb[, dims], column-major
    let mut x = Vec::with_capacity(n * dims.len());
    for &d in &dims {
        x.extend_from_slice(&emb[(d - 1) * n..d * n]);
    }
    let points = Points::new(&x, n, dims.len());
    let search = knn::purity::kth_distances(&points, k, threads);

    // the weights, and each cluster's repeated sums: sum_c[m] = w_c added m times to 0
    let size_of = |p: usize| present[uclust[p] as usize];
    let wk: Vec<f64> = (0..nk).map(|p| weight(size_of(p), weighted)).collect();
    let tables: Vec<Vec<f64>> = (0..nk)
        .map(|p| {
            let mut t = Vec::with_capacity(size_of(p) + 1);
            let mut s = 0.0f64;
            t.push(s);
            for _ in 0..size_of(p) {
                s += wk[p];
                t.push(s);
            }
            t
        })
        .collect();
    let weights: Vec<f64> = m.iter().map(|&p| wk[p as usize]).collect();

    // the neighbours, counted by cluster (and kept, if asked)
    let (per_cell, certificate) = search.within(&points, keep_lists, threads, |i, nb| {
        let mut counts = vec![0u32; nk];
        counts[m[i] as usize] = 1; // the cell's own weight
        for &j in nb {
            counts[m[j as usize] as usize] += 1;
        }
        (nb.len() as u32, counts, if keep_lists { nb.to_vec() } else { Vec::new() })
    });
    let mut neighbor_count = Vec::with_capacity(n);
    let mut sums = Vec::with_capacity(n * nk);
    let mut lists: Option<(Vec<u64>, Vec<u32>)> = if keep_lists { Some((vec![0], Vec::new())) } else { None };
    for (cnt, counts, list) in per_cell {
        neighbor_count.push(cnt);
        for (p, &c) in counts.iter().enumerate() {
            sums.push(tables[p][c as usize]);
        }
        if let Some((ptr, idx)) = lists.as_mut() {
            idx.extend_from_slice(&list);
            ptr.push(idx.len() as u64);
        }
    }

    // rowSums in long double, the purity, max.col(ties.method = "first")
    let mut totals = Vec::with_capacity(n);
    let mut purity = Vec::with_capacity(n);
    let mut maximum = Vec::with_capacity(n);
    for i in 0..n {
        let row = &sums[i * nk..(i + 1) * nk];
        let total = row_sum(row);
        totals.push(total);
        purity.push(row[m[i] as usize] / total);
        let (mut a, mut best) = (row[0], 0usize);
        for (c, &b) in row.iter().enumerate().skip(1) {
            if a < b {
                a = b;
                best = c;
            }
        }
        maximum.push(Some(best as u32));
    }
    Ok(Purity {
        dims,
        uclust,
        k_requested: search.k_requested,
        k: search.k,
        weighted,
        dist_k: search.dist_k.clone(),
        radius: search.radius,
        neighbor_count,
        neighbor_lists: lists,
        weights,
        sums,
        totals,
        purity,
        maximum,
        certificate,
        margin: Some(search.margin.clone()),
        computed: true,
    })
}

/// The search's first two steps alone (the k-th distances and the radius), for callers that
/// need nothing more.
pub fn kth_distances(emb: &[f64], n: usize, ncol: usize, dims: &[usize], k: usize, threads: usize) -> Result<KthDistances, String> {
    let dims = usable_dims(dims, ncol)?;
    let mut x = Vec::with_capacity(n * dims.len());
    for &d in &dims {
        x.extend_from_slice(&emb[(d - 1) * n..d * n]);
    }
    Ok(knn::purity::kth_distances(&Points::new(&x, n, dims.len()), k, threads))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn factor(codes: Vec<u32>, nlev: usize) -> Factor {
        Factor { name: "seurat_clusters".into(), levels: (0..nlev).map(|v| v.to_string()).collect(), codes }
    }

    #[test]
    fn small_example_by_hand() {
        // 6 cells on a line: 0, 1, 2 (cluster 0) and 10, 11, 13 (cluster 1); k = 1
        let emb = [0.0, 1.0, 2.0, 10.0, 11.0, 13.0];
        let f = factor(vec![0, 0, 0, 1, 1, 1], 2);
        let p = neighbor_purity(&emb, 6, 1, &[1], &f, 1, true, true, 1).unwrap();
        // nearest other cell: 1, 1, 1, 1, 1, 2; median of (1,1,1,1,1,2) = 1
        assert_eq!(p.dist_k, vec![1.0, 1.0, 1.0, 1.0, 1.0, 2.0]);
        assert_eq!(p.radius, 1.0);
        // within 1: cell 0 {1}, 1 {0, 2}, 2 {1}, 3 {4}, 4 {3}, 5 {}
        assert_eq!(p.neighbor_count, vec![1, 2, 1, 1, 1, 0]);
        assert_eq!(p.neighbor_lists.as_ref().unwrap().1, vec![1, 0, 2, 1, 4, 3]);
        let w = 1.0 / 3.0;
        assert_eq!(p.weights, vec![w; 6]);
        assert_eq!(&p.sums[0..2], &[w + w, 0.0]);
        assert_eq!(&p.sums[2..4], &[w + w + w, 0.0]);
        assert_eq!(&p.sums[10..12], &[0.0, w]);
        assert_eq!(p.purity, vec![1.0; 6]);
        assert_eq!(p.maximum_labels(&f), vec!["0", "0", "0", "1", "1", "1"]);
        assert_eq!(p.mean_neighbors(), 1.0);
    }

    #[test]
    fn repeated_sums_are_r_s_and_ties_go_to_the_first_cluster() {
        // cluster 0 has 3 cells (w = 1/3), cluster 1 has 3 (w = 1/3): a cell with as many
        // neighbours in each ties, and max.col takes the first
        let emb = [0.0, 0.1, 0.2, 0.05, 0.15, 0.25];
        let f = factor(vec![0, 0, 0, 1, 1, 1], 2);
        let p = neighbor_purity(&emb, 6, 1, &[1], &f, 5, true, false, 2).unwrap();
        // k = 5 = n - 1: every cell's farthest; the radius covers everything for some cells
        assert_eq!(p.k, 5);
        let w = 1.0f64 / 3.0;
        let s3 = w + w + w;
        for i in 0..6 {
            let row = &p.sums[i * 2..i * 2 + 2];
            if p.neighbor_count[i] == 5 {
                assert_eq!(row, &[s3, s3]);
                assert_eq!(p.maximum[i], Some(0), "a tie goes to the first cluster");
                assert_eq!(p.purity[i], s3 / (F80::ZERO.add_f64(s3).add_f64(s3)).to_f64());
            }
        }
        // unweighted: the sums are counts
        let q = neighbor_purity(&emb, 6, 1, &[1], &f, 5, false, false, 1).unwrap();
        assert!(q.sums.iter().all(|v| v.fract() == 0.0));
        assert_eq!(q.weights, vec![1.0; 6]);
    }

    #[test]
    fn one_cluster_gives_na_and_absent_levels_are_skipped() {
        let f = factor(vec![1, 1, 1], 3);
        let p = neighbor_purity(&[1.0, 2.0, 3.0], 3, 1, &[1], &f, 50, true, false, 1).unwrap();
        assert!(!p.computed && p.purity.iter().all(|v| v.is_nan()));
        assert_eq!(p.maximum_labels(&f), vec!["NA"; 3]);
        // levels 0..3 with level 1 absent: three columns
        let f = factor(vec![0, 2, 3, 0, 2, 3], 4);
        let p = neighbor_purity(&[0.0, 1.0, 2.0, 0.5, 1.5, 2.5], 6, 1, &[1], &f, 2, true, false, 1).unwrap();
        assert_eq!(p.uclust, vec![0, 2, 3]);
        assert_eq!(p.sums.len(), 18);
        assert!(p.maximum_labels(&f).iter().all(|l| l != "1"));
    }

    #[test]
    fn long_double_totals() {
        // R's rowSums accumulates in long double: 1 + 2^-53 + 2^-53 is 1 + 2^-52 there, 1 in double
        let e = 2f64.powi(-53);
        assert_eq!(row_sum(&[1.0, e, e]), 1.0 + 2.0 * e);
        assert_eq!((1.0 + e) + e, 1.0);
        let f = factor((0..12).map(|i| (i % 3) as u32).collect(), 3);
        let emb: Vec<f64> = (0..12).map(|i| i as f64 * 1e-3).collect();
        let p = neighbor_purity(&emb, 12, 1, &[1], &f, 11, true, false, 1).unwrap();
        for i in 0..12 {
            assert_eq!(p.totals[i].to_bits(), row_sum(&p.sums[i * 3..i * 3 + 3]).to_bits());
        }
    }
}
