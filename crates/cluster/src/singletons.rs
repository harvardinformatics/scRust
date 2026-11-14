//! Seurat's `GroupSingletons` (R code, `R/clustering.R`; docs/PLAN.md C2 and C3d): each
//! cluster of one cell joins the cluster with the highest mean SNN weight to that cell, ties
//! broken by R's random draw.
//!
//! *Ported from Seurat 5.5.1's `GroupSingletons` (Copyright 2021 Seurat authors, MIT), translated
//! into Rust statement by statement. NOTICE lists the changes.*
//!
//! As R runs it on the optimizer's integer ids:
//!
//! - **the singletons** are the clusters of one cell, in order of their first appearance among
//!   the cells (`intersect(unique(ids), names(which(table(ids) == 1)))`), and the candidate
//!   clusters the others, in the same order (`setdiff`);
//! - **a candidate's connectivity** is `mean(SNN[i.cells, j.cells])`: R's `mean`, in long double
//!   with its second pass ([`rcompat::mean`]), over the singleton's SNN row at the candidate's
//!   *current* cells in cell order, zeros included. A singleton merged earlier is one of those
//!   cells for the singletons after it;
//! - **the highest connectivity wins**; the tied candidates go to `sample(names, 1)` after
//!   `set.seed(1)`, which R calls before each candidate's connectivity, so every singleton's
//!   draw starts from the same state. A cell without SNN edges ties at 0 with every candidate;
//! - **with `group.singletons = FALSE`** the singletons are labelled `"singleton"` and nothing
//!   is drawn.
//!
//! The labels become character once a singleton is merged (`ids[i.cells] <- closest_cluster`);
//! the levels Seurat gives them do not depend on that ([`crate::find_clusters`]).

use rcompat::RRng;
use scio::CscMatrix;
use std::collections::HashMap;
use std::fmt;

/// GroupSingletons' result for one resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct Grouped {
    /// The singleton clusters (optimizer ids), in order of first appearance.
    pub singletons: Vec<u32>,
    /// Their cells, in the same order.
    pub singleton_cells: Vec<usize>,
    /// For each singleton: how many candidates tied at the highest connectivity, and that
    /// connectivity (empty with `group.singletons = FALSE`).
    pub n_tied: Vec<usize>,
    pub max_connectivity: Vec<f64>,
    /// For each singleton: the cluster it joined.
    pub joined: Vec<u32>,
    /// Each cell's final label: an optimizer id, or `None` for `"singleton"`.
    pub labels: Vec<Option<u32>>,
    /// R's `.Random.seed` after the last draw (`None` when nothing was drawn).
    pub rng_after: Option<Vec<i32>>,
    /// `group.singletons`.
    pub group: bool,
}

impl Grouped {
    /// The labels as R holds them: the id as text, or `"singleton"`.
    pub fn label_strings(&self) -> Vec<String> {
        self.labels.iter().map(|l| l.map_or_else(|| "singleton".to_string(), |v| v.to_string())).collect()
    }

    /// The number of distinct final labels (`length(unique(ids))`).
    pub fn n_final(&self) -> usize {
        let mut seen: Vec<Option<u32>> = self.labels.clone();
        seen.sort_unstable();
        seen.dedup();
        seen.len()
    }

    /// The message GroupSingletons gives (`verbose = TRUE`) when it grouped singletons. With
    /// `group.singletons = FALSE` it returns before its message, and says nothing.
    pub fn message(&self) -> Option<String> {
        (self.group && !self.singletons.is_empty()).then(|| format!("{} singletons identified. {} final clusters.", self.singletons.len(), self.n_final()))
    }
}

/// Where R's GroupSingletons stops.
#[derive(Debug, Clone, PartialEq)]
pub enum GroupError {
    /// Every cluster is a singleton: `sample()` of an empty candidate list fails in R
    /// ("invalid first argument").
    NoCandidates,
    /// The SNN matrix does not fit the ids.
    Shape(String),
}

impl fmt::Display for GroupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GroupError::NoCandidates => f.write_str("every cluster is a singleton, so GroupSingletons has no cluster to join them to (R's sample() fails: invalid first argument)"),
            GroupError::Shape(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for GroupError {}

/// `GroupSingletons(ids, SNN, group.singletons)` on the optimizer's ids (one per cell, the
/// SNN's row and column order).
pub fn group_singletons(ids: &[u32], snn: &CscMatrix, group: bool) -> Result<Grouped, GroupError> {
    let n = ids.len();
    if snn.nrow != n || snn.ncol != n {
        return Err(GroupError::Shape(format!("SNN graph is {} x {}, but there are {n} cells", snn.nrow, snn.ncol)));
    }
    let n_ids = ids.iter().map(|&v| v as usize + 1).max().unwrap_or(0);
    let mut size = vec![0usize; n_ids];
    for &v in ids {
        size[v as usize] += 1;
    }
    // unique(ids): order of first appearance
    let mut seen = vec![false; n_ids];
    let mut unique = Vec::new();
    let mut first_cell = vec![0usize; n_ids];
    for (c, &v) in ids.iter().enumerate() {
        if !seen[v as usize] {
            seen[v as usize] = true;
            unique.push(v);
            first_cell[v as usize] = c;
        }
    }
    let singletons: Vec<u32> = unique.iter().copied().filter(|&v| size[v as usize] == 1).collect();
    let singleton_cells: Vec<usize> = singletons.iter().map(|&v| first_cell[v as usize]).collect();
    let mut out = Grouped {
        singletons: singletons.clone(),
        singleton_cells: singleton_cells.clone(),
        n_tied: Vec::new(),
        max_connectivity: Vec::new(),
        joined: Vec::new(),
        labels: ids.iter().map(|&v| Some(v)).collect(),
        rng_after: None,
        group,
    };
    if !group {
        for &c in &singleton_cells {
            out.labels[c] = None;
        }
        return Ok(out);
    }
    if singletons.is_empty() {
        return Ok(out);
    }
    let cluster_names: Vec<u32> = unique.iter().copied().filter(|&v| size[v as usize] != 1).collect();
    if cluster_names.is_empty() {
        return Err(GroupError::NoCandidates);
    }
    // each candidate's current cells, in cell order
    let mut position = vec![usize::MAX; n_ids];
    for (k, &v) in cluster_names.iter().enumerate() {
        position[v as usize] = k;
    }
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); cluster_names.len()];
    for (c, &v) in ids.iter().enumerate() {
        if position[v as usize] != usize::MAX {
            members[position[v as usize]].push(c);
        }
    }
    // the singletons' SNN rows, from one pass over the matrix
    let mut rows: HashMap<usize, Vec<(usize, f64)>> = singleton_cells.iter().map(|&c| (c, Vec::new())).collect();
    for col in 0..snn.ncol {
        let (r, x) = snn.col(col);
        for (&ri, &xi) in r.iter().zip(x) {
            if let Some(row) = rows.get_mut(&(ri as usize)) {
                row.push((col, xi));
            }
        }
    }
    let mut dense = vec![0.0f64; n];
    let mut connectivity = vec![0.0f64; cluster_names.len()];
    let mut values = Vec::new();
    let mut rng = RRng::new(1);
    for &c in &singleton_cells {
        for &(col, v) in &rows[&c] {
            dense[col] = v;
        }
        for (k, cells) in members.iter().enumerate() {
            values.clear();
            values.extend(cells.iter().map(|&j| dense[j]));
            connectivity[k] = rcompat::mean(&values, false);
        }
        for &(col, _) in &rows[&c] {
            dense[col] = 0.0;
        }
        // m <- max(connectivity, na.rm = TRUE); which(connectivity == m)
        let m = connectivity.iter().copied().filter(|v| !v.is_nan()).fold(f64::NEG_INFINITY, f64::max);
        let tied: Vec<usize> = (0..connectivity.len()).filter(|&k| connectivity[k] == m).collect();
        // set.seed(1); sample(names(connectivity[mi]), 1)
        rng.set_seed(1);
        let pick = rcompat::sample::sample(&mut rng, &tied, 1, false).map_err(|_| GroupError::NoCandidates)?[0];
        let target = cluster_names[pick];
        out.labels[c] = Some(target);
        let list = &mut members[pick];
        let at = list.partition_point(|&j| j < c);
        list.insert(at, c);
        out.n_tied.push(tied.len());
        out.max_connectivity.push(m);
        out.joined.push(target);
    }
    out.rng_after = Some(rng.random_seed());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A symmetric matrix from (row, col, value) entries, both triangles given.
    fn csc(n: usize, entries: &[(usize, usize, f64)]) -> CscMatrix {
        let mut e: Vec<(usize, usize, f64)> = entries.to_vec();
        e.sort_by_key(|&(r, c, _)| (c, r));
        let mut m = CscMatrix::zeros(n, n);
        m.p.clear();
        m.p.push(0);
        let mut k = 0;
        for c in 0..n {
            while k < e.len() && e[k].1 == c {
                m.i.push(e[k].0 as u32);
                m.x.push(e[k].2);
                k += 1;
            }
            m.p.push(m.i.len());
        }
        m
    }

    #[test]
    fn no_singletons_changes_nothing() {
        let m = csc(4, &[(0, 1, 1.0), (1, 0, 1.0), (2, 3, 1.0), (3, 2, 1.0)]);
        let g = group_singletons(&[0, 0, 1, 1], &m, true).unwrap();
        assert!(g.singletons.is_empty() && g.rng_after.is_none());
        assert_eq!(g.labels, vec![Some(0), Some(0), Some(1), Some(1)]);
        assert_eq!(g.message(), None);
    }

    #[test]
    fn the_most_connected_cluster_wins_and_earlier_merges_count() {
        // clusters 0 = {0, 1}, 1 = {2, 3}; singletons 2 (cell 4) and 3 (cell 5)
        // cell 4 is tied to both pairs by its edges to cells 1 and 2, and cell 5 is joined to
        // cell 4 and cell 3; the order of appearance puts cell 5's cluster first
        let ids = [0, 0, 1, 1, 3, 2];
        let m = csc(6, &[(4, 1, 0.5), (1, 4, 0.5), (4, 3, 0.25), (3, 4, 0.25), (5, 4, 0.75), (4, 5, 0.75), (5, 0, 0.1), (0, 5, 0.1)]);
        let g = group_singletons(&ids, &m, true).unwrap();
        // singletons in order of first appearance: cluster 3 (cell 4), then 2 (cell 5)
        assert_eq!(g.singletons, vec![3, 2]);
        assert_eq!(g.singleton_cells, vec![4, 5]);
        // cell 4: mean 0.25 with {0, 1}, 0.125 with {2, 3}: joins 0
        assert_eq!(g.joined[0], 0);
        // cell 5 then sees {0, 1, 4}: (0.1 + 0 + 0.75) / 3, against 0 for {2, 3}
        assert_eq!(g.joined[1], 0);
        assert_eq!(g.max_connectivity[1], rcompat::mean(&[0.1, 0.0, 0.75], false));
        assert_eq!(g.labels, vec![Some(0), Some(0), Some(1), Some(1), Some(0), Some(0)]);
        assert_eq!(g.message().unwrap(), "2 singletons identified. 2 final clusters.");
        // the state after set.seed(1) and one draw
        let mut r = RRng::new(1);
        r.unif_rand();
        assert_eq!(g.rng_after.unwrap(), r.random_seed());
    }

    #[test]
    fn isolated_cells_tie_with_every_cluster() {
        let ids = [0, 0, 1, 1, 2];
        let m = csc(5, &[(0, 1, 1.0), (1, 0, 1.0), (2, 3, 1.0), (3, 2, 1.0), (4, 4, 1.0)]);
        let g = group_singletons(&ids, &m, true).unwrap();
        assert_eq!(g.n_tied, vec![2]);
        assert_eq!(g.max_connectivity, vec![0.0]);
        // set.seed(1); sample(c("0", "1"), 1) gives "0" (M8's K0 vectors: position 1 for n = 2)
        assert_eq!(g.joined, vec![0]);
    }

    #[test]
    fn without_grouping_singletons_are_labelled() {
        let ids = [0, 0, 1];
        let m = csc(3, &[(0, 1, 1.0), (1, 0, 1.0)]);
        let g = group_singletons(&ids, &m, false).unwrap();
        assert_eq!(g.label_strings(), vec!["0", "0", "singleton"]);
        assert!(g.rng_after.is_none());
        // GroupSingletons returns before its message
        assert_eq!(g.message(), None);
        assert_eq!(g.singletons, vec![1]);
        assert_eq!(group_singletons(&[0, 1], &CscMatrix::zeros(2, 2), true).unwrap_err(), GroupError::NoCandidates);
    }
}
