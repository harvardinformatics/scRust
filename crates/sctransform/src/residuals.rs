//! Stage S2d: the Pearson residuals `vst` computes with `residual_type = "pearson"` (Seurat 5.5.1's
//! default `sct.method`), and what `SCTransform.default` makes of them:
//!
//! 1. `min_var = (median(umi@x) / 5)^2` (`min_variance = "umi_median"`), over the stored entries
//!    of the kept genes' counts;
//! 2. for each bin of `bin_size` kept genes, `mu = exp(tcrossprod(model_pars_fit[, -1],
//!    model.matrix(~ log_umi)))` (OpenBLAS `dgemm`, as R calls it), the residuals
//!    `(y - mu) / sqrt(max(mu + mu^2 / theta, min_var))`, clipped to `+-sqrt(ncells)`;
//! 3. per gene, `residual_mean = rowMeans(y)` (base R: long double) and `residual_variance =
//!    matrixStats::rowVars(y)` (double, refined);
//! 4. the variable features: `sort(residual_variance, decreasing = TRUE)` (stable; NA dropped),
//!    the first `variable_features_n`;
//! 5. `scale.data`: their residuals clipped to `+-sqrt(ncells / 30)`, then centred by
//!    `ScaleData` (`FastRowScale`: `mat - rowMeans2(mat, na.rm = TRUE)`).
//!
//! R holds the full residual matrix (genes x cells); here bins are processed in parallel and only
//! the rows of the current best `variable_features_n` genes are kept (a worker drops a row whose
//! variance is already below the k-th best so far), so memory stays near `(variable_features_n +
//! threads * bin_size) * ncells` doubles. Results do not depend on the number of threads.

use crate::matrix::GeneRows;
use crate::params::{MinVariance, SctParams};
use crate::setup::bins;
use crate::SctError;
use rcompat::blas::tcrossprod;
use rcompat::libm::exp;
use rcompat::matrixstats::{mean2_na_rm, var_refined};
use rcompat::F80;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc;

/// What S2d reads.
#[derive(Debug, Clone)]
pub struct ResidualInput<'a> {
    /// per kept gene: theta, `(Intercept)`, `log_umi` (`model_pars_fit`)
    pub model_pars_fit: &'a [[f64; 3]],
    /// `log10(colSums(umi))` per cell (`cell_attr$log_umi`)
    pub log_umi: &'a [f64],
}

/// The result of S2d.
#[derive(Debug, Clone)]
pub struct Residuals {
    pub ncells: usize,
    /// `min_var`, and the clips `sqrt(ncells)` and `sqrt(ncells / 30)` (upper bounds)
    pub min_var: f64,
    pub res_clip: (f64, f64),
    pub clip: (f64, f64),
    /// per kept gene
    pub residual_mean: Vec<f64>,
    pub residual_variance: Vec<f64>,
    /// the variable features, positions in the kept genes, highest residual variance first
    pub top_features: Vec<usize>,
    /// `scale.data`: row `k` (cells in order) belongs to `top_features[k]`; row-major
    pub scale_data: Vec<f64>,
    /// the clipped residuals of every kept gene (row-major), when asked for
    pub y: Option<Vec<f64>>,
}

impl Residuals {
    /// Row `k` of `scale.data`.
    pub fn scale_row(&self, k: usize) -> &[f64] {
        &self.scale_data[k * self.ncells..(k + 1) * self.ncells]
    }
}

/// `min_var` for the kept genes' counts.
pub fn min_var(rows: &GeneRows, params: &SctParams) -> f64 {
    match params.min_variance {
        MinVariance::UmiMedian => {
            let m = rows.stored_median() / 5.0;
            m * m // R's x^2 is x * x
        }
        MinVariance::Value(v) => v,
    }
}

/// A candidate variable feature: ordered by residual variance, ties by gene position (the earlier
/// gene first, as R's stable sort keeps it).
#[derive(Debug, Clone, Copy)]
struct Cand {
    var: f64,
    gene: usize,
}

impl Cand {
    /// "better" = earlier in R's decreasing order
    fn better(&self, other: &Cand) -> Ordering {
        match self.var.partial_cmp(&other.var).expect("NaN residual variance in the ranking") {
            Ordering::Equal => other.gene.cmp(&self.gene),
            o => o,
        }
    }
}

impl PartialEq for Cand {
    fn eq(&self, o: &Self) -> bool {
        self.better(o) == Ordering::Equal
    }
}
impl Eq for Cand {}
impl PartialOrd for Cand {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Cand {
    /// reversed, so that a max-heap keeps the worst candidate on top
    fn cmp(&self, o: &Self) -> Ordering {
        o.better(self)
    }
}

/// One bin's results: every gene's statistics, and the clipped residual rows that may still
/// enter the variable features (all of them when the full matrix is kept).
struct BinOut {
    first: usize,
    mean: Vec<f64>,
    var: Vec<f64>,
    rows: Vec<(usize, Vec<f64>)>,
}

/// The current entry bar of the variable features: a row whose residual variance is below it
/// can no longer be selected (the bar only rises), so the workers drop it.
struct Bar(AtomicU64);

impl Bar {
    fn get(&self) -> f64 {
        f64::from_bits(self.0.load(AtomicOrdering::Relaxed))
    }
    fn set(&self, v: f64) {
        self.0.store(v.to_bits(), AtomicOrdering::Relaxed);
    }
}

fn residual_bin(rows: &GeneRows, input: &ResidualInput, genes: std::ops::Range<usize>, xreg: &[f64], min_var: f64, clip: (f64, f64), keep_all: bool, bar: &Bar) -> BinOut {
    let nc = rows.ncol;
    let nb = genes.len();
    // model_pars_final[genes_bin, -1]: (Intercept) and log_umi, column-major nb x 2
    let mut p = vec![0.0f64; nb * 2];
    for (r, g) in genes.clone().enumerate() {
        p[r] = input.model_pars_fit[g][1];
        p[r + nb] = input.model_pars_fit[g][2];
    }
    // tcrossprod(nb x 2, ncells x 2): nb x ncells, column-major
    let eta = tcrossprod(&p, nb, 2, xreg, nc, 2);
    let mut row = vec![0.0f64; nc];
    let mut counts = vec![0.0f64; nc];
    let mut mean = Vec::with_capacity(nb);
    let mut var = Vec::with_capacity(nb);
    let mut kept = Vec::new();
    let nf = F80::from_u64(nc as u64);
    for (r, g) in genes.enumerate() {
        rows.dense_row(g, &mut counts);
        let theta = input.model_pars_fit[g][0];
        for c in 0..nc {
            let mu = exp(eta[r + c * nb]);
            // pearson_residual(): model_var <- mu + mu^2 / theta; floored at min_var
            let mut model_var = mu + mu * mu / theta;
            if model_var < min_var {
                model_var = min_var;
            }
            let mut v = (counts[c] - mu) / model_var.sqrt();
            // clip_matrix_values(y, res_clip_range)
            if v < clip.0 {
                v = clip.0;
            }
            if v > clip.1 {
                v = clip.1;
            }
            row[c] = v;
        }
        // rowMeans(y): long-double sum in column order, divided by ncol in long double
        let mut s = F80::ZERO;
        for &v in row.iter() {
            s = s.add_f64(v);
        }
        mean.push((s / nf).to_f64());
        let v = var_refined(&row);
        var.push(v);
        if keep_all || v >= bar.get() {
            kept.push((g, row.clone()));
        }
    }
    BinOut { first: 0, mean, var, rows: kept }
}

/// Stage S2d. `rows` are the kept genes' counts, `threads` the number of bins computed at once.
pub fn residuals(rows: &GeneRows, input: &ResidualInput, params: &SctParams, keep_y: bool, threads: usize) -> Result<Residuals, SctError> {
    let ng = rows.nrow();
    let nc = rows.ncol;
    if input.model_pars_fit.len() != ng || input.log_umi.len() != nc {
        return Err(SctError::Input("residuals: model_pars_fit or log_umi does not match the counts".into()));
    }
    let min_var = min_var(rows, params);
    let res_clip = params.res_clip(nc);
    let clip = params.sct_clip(nc);
    // model.matrix(~ log_umi, cell_attr): (Intercept), log_umi; column-major ncells x 2
    let mut xreg = vec![1.0f64; nc * 2];
    xreg[nc..].copy_from_slice(input.log_umi);
    let bin_list = bins(ng, params.bin_size);
    // variable.features.n, or every gene at or above variable.features.rv.th
    let (k_top, threshold) = match params.variable_features_n {
        Some(n) => (n.min(ng), f64::NEG_INFINITY),
        None => (usize::MAX, params.variable_features_rv_th),
    };
    let mut residual_mean = vec![0.0f64; ng];
    let mut residual_variance = vec![0.0f64; ng];
    let mut heap: BinaryHeap<Cand> = BinaryHeap::with_capacity(k_top.min(ng) + 1);
    let mut kept: HashMap<usize, Vec<f64>> = HashMap::with_capacity(k_top.min(ng) + 1);
    let mut all_y = if keep_y { Some(vec![0.0f64; ng * nc]) } else { None };
    let threads = threads.max(1).min(bin_list.len().max(1));
    let next = AtomicUsize::new(0);
    // rows below the bar cannot be selected: the threshold, or the k-th best variance so far
    let bar = Bar(AtomicU64::new(threshold.to_bits()));
    let (tx, rx) = mpsc::sync_channel::<BinOut>(threads);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next, bin_list, xreg, bar) = (&next, &bin_list, &xreg, &bar);
            scope.spawn(move || loop {
                let b = next.fetch_add(1, AtomicOrdering::Relaxed);
                if b >= bin_list.len() {
                    break;
                }
                let mut out = residual_bin(rows, input, bin_list[b].clone(), xreg, min_var, res_clip, keep_y, bar);
                out.first = bin_list[b].start;
                if tx.send(out).is_err() {
                    break;
                }
            });
        }
        drop(tx);
        for out in rx {
            for (r, (&m, &v)) in out.mean.iter().zip(&out.var).enumerate() {
                residual_mean[out.first + r] = m;
                residual_variance[out.first + r] = v;
            }
            for (g, row) in out.rows {
                let v = residual_variance[g];
                if let Some(all) = all_y.as_mut() {
                    all[g * nc..(g + 1) * nc].copy_from_slice(&row);
                }
                if v.is_nan() || k_top == 0 || !(v >= threshold) {
                    continue; // sort() drops NA
                }
                let cand = Cand { var: v, gene: g };
                if heap.len() < k_top {
                    heap.push(cand);
                    kept.insert(g, row);
                } else if cand.better(heap.peek().unwrap()) == Ordering::Greater {
                    let worst = heap.pop().unwrap();
                    kept.remove(&worst.gene);
                    heap.push(cand);
                    kept.insert(g, row);
                }
                if heap.len() == k_top {
                    bar.set(heap.peek().unwrap().var);
                }
            }
        }
    });
    // the variable features in R's order
    let mut top: Vec<Cand> = heap.into_vec();
    top.sort_by(|a, b| b.better(a));
    let top_features: Vec<usize> = top.iter().map(|c| c.gene).collect();
    // scale.data: clip to +-sqrt(ncells / 30), then (do.center) subtract rowMeans2(na.rm = TRUE)
    let mut scale_data = vec![0.0f64; top_features.len() * nc];
    for (k, &g) in top_features.iter().enumerate() {
        let src = &kept[&g];
        let row = &mut scale_data[k * nc..(k + 1) * nc];
        for (d, &v) in row.iter_mut().zip(src) {
            let mut v = v;
            if v < clip.0 {
                v = clip.0;
            }
            if v > clip.1 {
                v = clip.1;
            }
            *d = v;
        }
        if params.do_center {
            let rm = mean2_na_rm(row);
            for d in row.iter_mut() {
                *d -= rm;
            }
        }
    }
    Ok(Residuals { ncells: nc, min_var, res_clip, clip, residual_mean, residual_variance, top_features, scale_data, y: all_y })
}
