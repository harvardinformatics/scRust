//! FindNeighbors' two graphs from the neighbour lists (docs/PLAN.md C1 step 4 and C3b):
//!
//! - **`SCT_nn`**: `sparseMatrix(i = rep(cells, each = k), j = neighbours, x = 1)`, the 0/1
//!   matrix with a 1 at (cell, neighbour); Matrix stores it as compressed columns with the rows
//!   of each column increasing.
//! - **`SCT_snn`**: Seurat's `ComputeSNN` (src/snn.cpp): with A that matrix, S = A Aᵀ counts the
//!   neighbours s that two cells share, an exact integer; each stored entry becomes
//!   `s / (k + (k - s))` in double, in that order (`int k` against the double s), and is removed
//!   when it is `< prune` (`SNN.prune(0.0)` drops exactly the zeroed entries). The diagonal stays
//!   (s = k, value 1). Eigen's product leaves the rows of each column increasing, which the
//!   reference graphs confirm (M8, every object).
//!
//! Both are built from the neighbour lists without forming A Aᵀ: column j of S holds the cells
//! i that share a neighbour c with j, found through the lists of cells having c as a neighbour
//! (the columns of A). Columns are independent, so they are computed in parallel; the result
//! does not depend on the thread count.

use scio::CscMatrix;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// The `SCT_nn` graph of `n` cells from their `k` neighbours each (`nn[i * k + r]`, 0-based):
/// value 1 at (i, nn[i][r]); each column's rows increasing. Lists must not repeat an id.
pub fn nn_graph(nn: &[i32], n: usize, k: usize) -> CscMatrix {
    assert_eq!(nn.len(), n * k);
    let mut p = vec![0usize; n + 1];
    for &c in nn {
        p[c as usize + 1] += 1;
    }
    for j in 0..n {
        p[j + 1] += p[j];
    }
    let mut next = p.clone();
    let mut rows = vec![0u32; nn.len()];
    // rows in increasing order: cells are visited in order
    for i in 0..n {
        for &c in &nn[i * k..(i + 1) * k] {
            let c = c as usize;
            rows[next[c]] = i as u32;
            next[c] += 1;
        }
    }
    CscMatrix { nrow: n, ncol: n, p, i: rows, x: vec![1.0; nn.len()] }
}

/// `ComputeSNN(nn_ranked, prune)`: the `SCT_snn` graph from the neighbour lists and the
/// `SCT_nn` graph built from them (whose column c lists the cells having c as a neighbour).
pub fn snn_graph(nn: &[i32], nn_graph: &CscMatrix, n: usize, k: usize, prune: f64, threads: usize) -> CscMatrix {
    assert_eq!(nn.len(), n * k);
    let kf = k as f64;
    let threads = threads.max(1).min(n.max(1));
    const CHUNK: usize = 256;
    let next = AtomicUsize::new(0);
    let parts: Mutex<Vec<(usize, Vec<usize>, Vec<u32>, Vec<f64>)>> = Mutex::new(Vec::new());
    let work = || {
        let mut count = vec![0u32; n];
        let mut touched: Vec<u32> = Vec::new();
        loop {
            let start = next.fetch_add(CHUNK, Ordering::Relaxed);
            if start >= n {
                break;
            }
            let end = (start + CHUNK).min(n);
            let mut lens = Vec::with_capacity(end - start);
            let (mut rows, mut vals) = (Vec::new(), Vec::new());
            for j in start..end {
                for &c in &nn[j * k..(j + 1) * k] {
                    let (r, _) = nn_graph.col(c as usize);
                    for &i in r {
                        let slot = &mut count[i as usize];
                        if *slot == 0 {
                            touched.push(i);
                        }
                        *slot += 1;
                    }
                }
                touched.sort_unstable();
                let before = rows.len();
                for &i in &touched {
                    let s = count[i as usize] as f64;
                    count[i as usize] = 0;
                    // it.value() / (k + (k - it.value())); pruned when < prune
                    let v = s / (kf + (kf - s));
                    if !(v < prune) {
                        rows.push(i);
                        vals.push(v);
                    }
                }
                touched.clear();
                lens.push(rows.len() - before);
            }
            parts.lock().unwrap().push((start, lens, rows, vals));
        }
    };
    if threads == 1 {
        work();
    } else {
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(work);
            }
        });
    }
    let mut parts = parts.into_inner().unwrap();
    parts.sort_unstable_by_key(|p| p.0);
    let nnz: usize = parts.iter().map(|p| p.2.len()).sum();
    let mut out = CscMatrix { nrow: n, ncol: n, p: Vec::with_capacity(n + 1), i: Vec::with_capacity(nnz), x: Vec::with_capacity(nnz) };
    out.p.push(0);
    for (_, lens, rows, vals) in parts {
        for l in lens {
            let last = *out.p.last().unwrap();
            out.p.push(last + l);
        }
        out.i.extend_from_slice(&rows);
        out.x.extend_from_slice(&vals);
    }
    debug_assert_eq!(out.p.len(), n + 1);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// S = A Aᵀ by brute force, then Seurat's transformation and pruning.
    fn brute(nn: &[i32], n: usize, k: usize, prune: f64) -> CscMatrix {
        let mut a = vec![vec![0u32; n]; n];
        for i in 0..n {
            for &c in &nn[i * k..(i + 1) * k] {
                a[i][c as usize] += 1;
            }
        }
        let mut m = CscMatrix::zeros(n, n);
        m.p.clear();
        m.p.push(0);
        for j in 0..n {
            for i in 0..n {
                let s: u32 = (0..n).map(|c| a[i][c] * a[j][c]).sum();
                if s > 0 {
                    let s = s as f64;
                    let v = s / (k as f64 + (k as f64 - s));
                    if !(v < prune) {
                        m.i.push(i as u32);
                        m.x.push(v);
                    }
                }
            }
            m.p.push(m.i.len());
        }
        m
    }

    #[test]
    fn snn_equals_the_matrix_product() {
        let (n, k) = (60usize, 7usize);
        let mut state = 17u64;
        let mut nn = Vec::with_capacity(n * k);
        for i in 0..n {
            let mut row = vec![i as i32];
            while row.len() < k {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let c = ((state >> 33) % (n as u64)) as i32;
                if !row.contains(&c) {
                    row.push(c);
                }
            }
            nn.extend(row);
        }
        let g = nn_graph(&nn, n, k);
        g.validate().unwrap();
        for prune in [0.0, 1.0 / 15.0, 0.2, 1.0] {
            let want = brute(&nn, n, k, prune);
            for t in [1, 3] {
                let s = snn_graph(&nn, &g, n, k, prune, t);
                s.validate().unwrap();
                assert_eq!(s, want, "prune {prune}, {t} threads");
            }
        }
    }

    #[test]
    fn values_are_jaccard_in_seurats_order() {
        // k = 20: s shared neighbours give s / (40 - s); s <= 2 is pruned at 1/15
        for s in 1..=20u32 {
            let v = s as f64 / (20.0 + (20.0 - s as f64));
            assert_eq!(v < 1.0 / 15.0, s <= 2, "s = {s}");
        }
    }
}
