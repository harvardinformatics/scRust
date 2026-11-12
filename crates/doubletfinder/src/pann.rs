//! Stages S4 and S5: neighbours and pANN, in the sweep (`parallel_paramSweep`) and the final run
//! (`doubletFinder`), and stage S7's classification.
//!
//! R ranks every cell by distance with `order(fields::rdist(pca)[, i])` and counts the
//! artificial doublets (indices above the number of real cells) among order positions
//! `2:(k + 1)`, whatever cell position 1 holds. The [`knn`] crate gives the same positions
//! without the distance matrix.

use crate::DfError;
use knn::{neighbours, Digest, Order, Points};
use rcompat::arith::round;

/// Cell embeddings: `n` cells (real cells first, then the artificial doublets) by `dim`
/// components, column-major as R stores `cell.embeddings`.
#[derive(Debug, Clone)]
pub struct Embeddings {
    pub data: Vec<f64>,
    pub n: usize,
    pub dim: usize,
}

impl Embeddings {
    pub fn new(data: Vec<f64>, n: usize, dim: usize) -> Result<Embeddings, DfError> {
        if data.len() != n * dim || n == 0 || dim == 0 {
            return Err(DfError::Input(format!("embeddings: {} values for {n} cells x {dim} components", data.len())));
        }
        Ok(Embeddings { data, n, dim })
    }

    /// The first `pcs` components (R's `cell.embeddings[, PCs]` with `PCs = 1:pcs`).
    pub fn points(&self, pcs: usize) -> Result<Points<'_>, DfError> {
        if pcs == 0 || pcs > self.dim {
            return Err(DfError::Input(format!("{pcs} PCs requested, embeddings have {}", self.dim)));
        }
        Ok(Points::new(&self.data[..pcs * self.n], self.n, pcs))
    }
}

/// The 0-based order positions R reads with `[2:(k + 1)]`, as an inclusive range: positions 2
/// to k + 1 for k >= 1, and `2:1`, i.e. positions 2 and 1, for k = 0.
pub fn r_neighbour_range(k: usize) -> (usize, usize) {
    if k >= 1 {
        (1, k)
    } else {
        (0, 1)
    }
}

/// `k` from `round(nCells * pK)`, checked to be a count.
fn count_of(x: f64, what: &str) -> Result<usize, DfError> {
    if x.is_finite() && x >= 0.0 && x <= u32::MAX as f64 {
        Ok(x as usize)
    } else {
        Err(DfError::Input(format!("{what} = {x} is not a valid neighbour count")))
    }
}

/// Optional per-cell outputs for validation: R's S4 digest and the order positions themselves.
#[derive(Debug, Clone, Copy, Default)]
pub struct NeighbourDumps {
    pub digests: bool,
    pub positions: bool,
}

/// One iteration of `parallel_paramSweep`: pANN of every real cell at every pK.
#[derive(Debug, Clone)]
pub struct SweepIteration {
    /// Real cells (the first rows of the embeddings).
    pub n_real: usize,
    /// `nCells`: real cells plus artificial doublets.
    pub n_cells: usize,
    /// `ind = round(nCells * max(pK)) + 5`: order positions kept per cell.
    pub ind: usize,
    /// `pk.temp = round(nCells * pK)` for each pK.
    pub pk_temp: Vec<usize>,
    /// pANN for each pK (outer) and real cell (inner), as `sweep.res.list`.
    pub pann: Vec<Vec<f64>>,
    /// S4 digest `(h1, h2)` per real cell, if requested.
    pub digests: Option<Vec<(u64, u64)>>,
    /// The first `ind` order positions per real cell (0-based), if requested.
    pub positions: Option<Vec<Vec<u32>>>,
    /// Real cells whose order position 1 is the cell itself (R's dumps record it; a
    /// lower-index twin at distance 0 comes first otherwise).
    pub first_is_self: usize,
}

/// `parallel_paramSweep`'s neighbour ranking and pANN for one pN, from the embeddings of the
/// real cells and doublets: `ind = round(nCells * max(pK)) + 5` positions per real cell, then
/// for each pK, `pk.temp = round(nCells * pK)` and pANN = (doublets among positions
/// `2:(pk.temp + 1)`) / pk.temp.
pub fn sweep_iteration(emb: &Embeddings, pcs: usize, n_real: usize, pk: &[f64], threads: usize, dumps: NeighbourDumps) -> Result<SweepIteration, DfError> {
    let points = emb.points(pcs)?;
    let n_cells = emb.n;
    if n_real == 0 || n_real > n_cells || pk.is_empty() {
        return Err(DfError::Input(format!("sweep: {n_real} real cells of {n_cells}, {} pK values", pk.len())));
    }
    let max_pk = pk.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let ind = count_of(round(n_cells as f64 * max_pk) + 5.0, "ind")?;
    if ind > n_cells {
        // dist.mat[1:ind, ]: "subscript out of bounds"
        return Err(DfError::Input(format!("sweep: ind = {ind} exceeds the {n_cells} cells (subscript out of bounds in R)")));
    }
    let pk_temp: Vec<usize> = pk.iter().map(|&p| count_of(round(n_cells as f64 * p), "pk.temp")).collect::<Result<_, _>>()?;
    let ranges: Vec<(usize, usize)> = pk_temp.iter().map(|&k| r_neighbour_range(k)).collect();
    let digest = dumps.digests.then(|| Digest::new(ind));
    let per_cell = neighbours(&points, n_real, ind, Order::Sorted, threads, |q, pos| {
        // cum[m]: doublets among the first m positions
        let mut cum = Vec::with_capacity(pos.len() + 1);
        cum.push(0u32);
        let mut c = 0u32;
        for &p in pos {
            c += (p as usize >= n_real) as u32;
            cum.push(c);
        }
        let counts: Vec<u32> = ranges.iter().map(|&(lo, hi)| cum[hi + 1] - cum[lo]).collect();
        (counts, digest.as_ref().map(|d| d.of(pos)), dumps.positions.then(|| pos.to_vec()), pos.first() == Some(&(q as u32)))
    });
    let mut pann = vec![Vec::with_capacity(n_real); pk.len()];
    let mut digests = dumps.digests.then(|| Vec::with_capacity(n_real));
    let mut positions = dumps.positions.then(|| Vec::with_capacity(n_real));
    let mut first_is_self = 0;
    for (counts, dg, ps, own) in per_cell {
        first_is_self += own as usize;
        for (j, &cnt) in counts.iter().enumerate() {
            pann[j].push(cnt as f64 / pk_temp[j] as f64);
        }
        if let (Some(d), Some(v)) = (dg, digests.as_mut()) {
            v.push(d);
        }
        if let (Some(p), Some(v)) = (ps, positions.as_mut()) {
            v.push(p);
        }
    }
    Ok(SweepIteration { n_real, n_cells, ind, pk_temp, pann, digests, positions, first_is_self })
}

/// The final run of `doubletFinder`: pANN of every real cell.
#[derive(Debug, Clone)]
pub struct FinalRun {
    pub n_real: usize,
    pub n_cells: usize,
    /// `k = round(nCells * pK)`.
    pub k: f64,
    pub pann: Vec<f64>,
    /// S4 digest of the first k + 1 positions per cell, if requested.
    pub digests: Option<Vec<(u64, u64)>>,
    /// The first k + 1 order positions per cell (0-based), if requested.
    pub positions: Option<Vec<Vec<u32>>>,
    /// Real cells whose order position 1 is the cell itself.
    pub first_is_self: usize,
}

/// `doubletFinder`'s pANN: `k = round(nCells * pK)`, then for each real cell the doublets among
/// order positions `2:(k + 1)`, over k. Without dumps, only the set of the first k + 1 positions
/// and position 1 are needed, so no sort is done.
pub fn final_pann(emb: &Embeddings, pcs: usize, n_real: usize, pk: f64, threads: usize, dumps: NeighbourDumps) -> Result<FinalRun, DfError> {
    let points = emb.points(pcs)?;
    let n_cells = emb.n;
    if n_real == 0 || n_real > n_cells {
        return Err(DfError::Input(format!("final run: {n_real} real cells of {n_cells}")));
    }
    let k = round(n_cells as f64 * pk);
    let ku = count_of(k, "k")?;
    let (lo, hi) = r_neighbour_range(ku);
    // positions beyond nCells are NA in R and count nothing
    let l = (hi + 1).max(ku + 1).min(n_cells);
    let sorted = dumps.digests || dumps.positions;
    let digest = dumps.digests.then(|| Digest::new((ku + 1).min(n_cells)));
    let per_cell = neighbours(&points, n_real, l, if sorted { Order::Sorted } else { Order::FirstThenSet }, threads, |q, pos| {
        let is_doublet = |p: u32| (p as usize >= n_real) as u32;
        let count = if sorted {
            pos[lo.min(pos.len())..(hi + 1).min(pos.len())].iter().map(|&p| is_doublet(p)).sum::<u32>()
        } else {
            // position 1 first, then the rest of the first l positions
            let all: u32 = pos[..(hi + 1).min(pos.len())].iter().map(|&p| is_doublet(p)).sum();
            if lo == 1 {
                all - is_doublet(pos[0])
            } else {
                all
            }
        };
        let take = (ku + 1).min(pos.len());
        // either order puts position 1 first
        (count, digest.as_ref().map(|d| d.of(&pos[..take])), dumps.positions.then(|| pos[..take].to_vec()), pos.first() == Some(&(q as u32)))
    });
    let mut pann = Vec::with_capacity(n_real);
    let mut digests = dumps.digests.then(|| Vec::with_capacity(n_real));
    let mut positions = dumps.positions.then(|| Vec::with_capacity(n_real));
    let mut first_is_self = 0;
    for (count, dg, ps, own) in per_cell {
        first_is_self += own as usize;
        pann.push(count as f64 / k);
        if let (Some(d), Some(v)) = (dg, digests.as_mut()) {
            v.push(d);
        }
        if let (Some(p), Some(v)) = (ps, positions.as_mut()) {
            v.push(p);
        }
    }
    Ok(FinalRun { n_real, n_cells, k, pann, digests, positions, first_is_self })
}

/// Stage S7: `classifications[order(pANN, decreasing = TRUE)[1:nExp]] <- "Doublet"`, a stable
/// descending order (ties to the lower index). R's `1:nExp` gives one cell for 0 <= nExp < 1
/// (`1:0` is `c(1, 0)`) and stops for nExp <= -1; indices beyond the cells are NA and ignored.
pub fn classify(pann: &[f64], nexp: f64) -> Result<Vec<bool>, DfError> {
    let top = if nexp >= 1.0 {
        nexp.floor()
    } else if nexp > -1.0 {
        1.0
    } else {
        return Err(DfError::Input(format!("nExp = {nexp}: only 0's may be mixed with negative subscripts")));
    };
    let top = (top.min(pann.len() as f64)) as usize;
    let order = rcompat::order(pann, true);
    let mut doublet = vec![false; pann.len()];
    for &i in &order[..top] {
        doublet[i] = true;
    }
    Ok(doublet)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_positions() {
        assert_eq!(r_neighbour_range(1), (1, 1));
        assert_eq!(r_neighbour_range(5), (1, 5));
        assert_eq!(r_neighbour_range(0), (0, 1));
    }

    #[test]
    fn classification_ties_and_edges() {
        let p = [0.2, 0.5, 0.5, 0.1, 0.5];
        assert_eq!(classify(&p, 2.0).unwrap(), [false, true, true, false, false]);
        assert_eq!(classify(&p, 0.0).unwrap(), [false, true, false, false, false]);
        assert_eq!(classify(&p, 9.0).unwrap(), [true; 5]);
        assert!(classify(&p, -1.0).is_err());
    }

    #[test]
    fn sweep_and_final_agree_with_brute_force() {
        // 40 real cells and 20 doublets in 3 dimensions, with duplicated coordinates
        let (n_real, n) = (40usize, 60usize);
        let mut s = 99u64;
        let mut next = move || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 33) % 7) as f64
        };
        let mut data = vec![0.0; n * 3];
        for v in data.iter_mut() {
            *v = next();
        }
        let emb = Embeddings::new(data, n, 3).unwrap();
        let pts = emb.points(3).unwrap();
        let brute = |q: usize| -> Vec<usize> {
            let mut d = vec![0.0; n];
            knn::distances(&pts, q, &mut d);
            rcompat::order(&d, false)
        };
        let pk = [0.05, 0.1, 0.2];
        let it = sweep_iteration(&emb, 3, n_real, &pk, 2, NeighbourDumps::default()).unwrap();
        assert_eq!(it.pk_temp, [3, 6, 12]);
        assert_eq!(it.ind, 17);
        for q in 0..n_real {
            let o = brute(q);
            for (j, &k) in it.pk_temp.iter().enumerate() {
                let c = o[1..=k].iter().filter(|&&i| i >= n_real).count();
                assert_eq!(it.pann[j][q], c as f64 / k as f64);
            }
        }
        for dumps in [NeighbourDumps::default(), NeighbourDumps { digests: true, positions: true }] {
            let f = final_pann(&emb, 3, n_real, 0.15, 3, dumps).unwrap();
            assert_eq!(f.k, 9.0);
            for q in 0..n_real {
                let c = brute(q)[1..=9].iter().filter(|&&i| i >= n_real).count();
                assert_eq!(f.pann[q], c as f64 / 9.0);
            }
        }
    }
}
