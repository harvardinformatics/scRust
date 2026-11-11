//! Stage S2a, `vst`'s setup (sctransform 0.4.3, `R/vst.R` up to `get_model_pars`): cell
//! attributes, the gene filter, geometric means, the fitting cells and genes, and the bins.
//!
//! Arithmetic, in R's order:
//! - `umi = colSums(umi)` (Matrix: a double sum of each column's stored entries, exact for the
//!   integer and half-integer counts here) and `log_umi = log10(umi)`;
//! - genes kept if `rowSums(umi >= 0.01) >= min_cells`;
//! - `row_gmean_dgcmatrix` (sctransform's C++): per row, `sum += log(x + eps)` over the stored
//!   entries in column order, then `exp((sum + log(eps) * n_zero) / ncol) - eps`, `log10` of it;
//! - fitting cells: `sample(colnames(umi), n_cells)` when `n_cells < ncol` (draw order kept);
//! - fitting genes: detected (`> 0`) in at least `min_cells` fitting cells, then only the
//!   overdispersed ones (`row_var - rowMeans > 0`, over all cells; `row_var_dgcmatrix` adds the
//!   stored entries' squared deviations in column order and the implicit zeros' at the end, so
//!   explicit zeros change the bits), then, if more than `n_genes` remain, `n_genes` of them drawn
//!   without replacement with probability `1 / (density(log10 gmean) + eps)`;
//! - bins of `bin_size` genes in draw order.

use crate::matrix::GeneRows;
use crate::params::SctParams;
use crate::SctError;
use rcompat::approx::{ApproxFun, Method};
use rcompat::density::density_nrd;
use rcompat::libm::{exp, log, log10};
use rcompat::{sample_int, sample_int_weighted, RRng};

/// `vst`'s cell attributes (`make_cell_attr` for `latent_var = "log_umi"`).
#[derive(Debug, Clone, PartialEq)]
pub struct CellAttr {
    /// `colSums(umi)`
    pub umi: Vec<f64>,
    /// `log10(umi)`
    pub log_umi: Vec<f64>,
}

/// The density-weighted gene sampling, when it happens.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneSampling {
    /// `density(log10 gmean of the pool, bw = "nrd")$x` and `$y`
    pub density_x: Vec<f64>,
    pub density_y: Vec<f64>,
    pub density_bw: f64,
    /// `1 / (approx(density, xout = pool)$y + .Machine$double.eps)`, one per pool gene
    pub sampling_prob: Vec<f64>,
}

/// The result of S2a. Gene indices are positions in `genes` (the filtered genes) unless noted;
/// cell indices are columns of the vst input.
#[derive(Debug, Clone, PartialEq)]
pub struct Setup {
    pub n_cells: usize,
    pub cell_attr: CellAttr,
    /// rows of the input kept by the `min_cells` filter, in input order
    pub genes: Vec<usize>,
    /// `rowSums(umi >= 0.01)` of the kept genes
    pub genes_cell_count: Vec<u32>,
    /// `log10(row_gmean(umi))` of the kept genes
    pub log_gmean: Vec<f64>,
    /// `rowMeans(umi)` and `row_var(umi)` of the kept genes
    pub amean: Vec<f64>,
    pub var: Vec<f64>,
    /// fitting cells, in draw order (all cells in order without cell sampling)
    pub cells_step1: Vec<usize>,
    /// fitting-gene candidates before the density sampling (after the Poisson exclusion)
    pub genes_step1_pool: Vec<usize>,
    /// fitting genes, in draw order
    pub genes_step1: Vec<usize>,
    /// `log10(row_gmean)` of the fitting genes
    pub log_gmean_step1: Vec<f64>,
    pub sampling: Option<GeneSampling>,
    /// `.Random.seed` after the cell sampling and after the gene sampling (R's checkpoints
    /// "after_cell_sampling" and "fits_entry")
    pub rng_after_cells: Vec<i32>,
    pub rng_after_genes: Vec<i32>,
}

impl Setup {
    /// `detection_rate = genes_cell_count / ncol(umi)` of the kept genes.
    pub fn detection_rate(&self) -> Vec<f64> {
        self.genes_cell_count.iter().map(|&c| c as f64 / self.n_cells as f64).collect()
    }

    /// `gmean = 10 ^ log_gmean` (R's `^`: glibc `pow`).
    pub fn gmean(&self) -> Vec<f64> {
        self.log_gmean.iter().map(|&l| rcompat::r_pow(10.0, l)).collect()
    }

    /// Bin (0-based) of each fitting gene: `ceiling(1:n / bin_size) - 1`.
    pub fn fit_bins(&self, bin_size: usize) -> Vec<std::ops::Range<usize>> {
        bins(self.genes_step1.len(), bin_size)
    }
}

/// Consecutive ranges of `bin_size` (the last one shorter): `bin_ind <- ceiling(1:n / bin_size)`.
pub fn bins(n: usize, bin_size: usize) -> Vec<std::ops::Range<usize>> {
    (0..n.div_ceil(bin_size)).map(|b| b * bin_size..((b + 1) * bin_size).min(n)).collect()
}

/// `colSums` of a counts matrix, as Matrix computes it for a `dgCMatrix`.
pub fn col_sums(m: &scio::CscMatrix) -> Vec<f64> {
    (0..m.ncol)
        .map(|j| {
            let mut s = 0.0;
            for &v in m.col(j).1 {
                s += v;
            }
            s
        })
        .collect()
}

/// `row_gmean_dgcmatrix(umi, eps)` for one row's stored values (in column order).
pub fn row_gmean(values: &[f64], ncol: usize, eps: f64) -> f64 {
    let mut ret = 0.0;
    for &v in values {
        ret += log(v + eps);
    }
    let nzero = (ncol - values.len()) as i32;
    let log_eps = log(eps);
    exp((ret + log_eps * nzero as f64) / ncol as i32 as f64) - eps
}

/// `row_var_dgcmatrix(x, i, rows, cols)` for one row's stored values (in column order).
pub fn row_var(values: &[f64], ncol: usize) -> f64 {
    let mut mean = 0.0;
    for &v in values {
        mean += v;
    }
    mean /= ncol as i32 as f64;
    let mut var = 0.0;
    for &v in values {
        let d = v - mean;
        var += d * d;
    }
    let nzero = (ncol - values.len()) as i32;
    (var + (mean * mean) * nzero as f64) / (ncol as i32 - 1) as f64
}

/// Matrix's `rowMeans` of a `dgCMatrix` row: the double sum of the stored values divided by `ncol`.
pub fn row_mean(values: &[f64], ncol: usize) -> f64 {
    let mut s = 0.0;
    for &v in values {
        s += v;
    }
    s / ncol as i32 as f64
}

/// Stage S2a on the vst input `counts` (genes x cells). `rng` must be in R's state at `vst`'s
/// entry (after `SCTransform`'s `set.seed` calls); it is advanced by the cell and gene sampling.
/// Also returns the kept genes' rows, which the later stages read.
pub fn setup(counts: &scio::CscMatrix, params: &SctParams, rng: &mut RRng) -> Result<(Setup, GeneRows), SctError> {
    params.check().map_err(SctError::Input)?;
    let n_cells = counts.ncol;
    if n_cells < 2 {
        return Err(SctError::Input("vst needs at least 2 cells".into()));
    }
    let umi = col_sums(counts);
    let log_umi: Vec<f64> = umi.iter().map(|&u| log10(u)).collect();
    if log_umi.iter().any(|v| !v.is_finite()) {
        return Err(SctError::Input("cell attribute \"log_umi\" contains NA, NaN, or infinite value".into()));
    }
    let csr = counts.to_csr();
    // genes_cell_count <- rowSums(umi >= 0.01); genes <- rownames(umi)[genes_cell_count >= min_cells]
    let mut genes = Vec::new();
    let mut genes_cell_count = Vec::new();
    for r in 0..counts.nrow {
        let c = csr.row(r).1.iter().filter(|&&v| v >= 0.01).count() as u32;
        if c as usize >= params.min_cells {
            genes.push(r);
            genes_cell_count.push(c);
        }
    }
    let ng = genes.len();
    let eps = params.gmean_eps;
    let mut log_gmean = Vec::with_capacity(ng);
    let mut amean = Vec::with_capacity(ng);
    let mut var = Vec::with_capacity(ng);
    for &g in &genes {
        let vals = csr.row(g).1;
        log_gmean.push(log10(row_gmean(vals, n_cells, eps)));
        amean.push(row_mean(vals, n_cells));
        var.push(row_var(vals, n_cells));
    }

    // fitting cells
    let n_cells_fit = params.n_cells_fit(n_cells);
    let (cells_step1, mut genes_step1): (Vec<usize>, Vec<usize>) = if n_cells_fit < n_cells {
        let cells = sample_int(rng, n_cells, n_cells_fit, false).map_err(|e| SctError::Input(e.to_string()))?;
        let mut in_step1 = vec![false; n_cells];
        for &c in &cells {
            in_step1[c] = true;
        }
        let g1: Vec<usize> = (0..ng)
            .filter(|&k| {
                let (cols, vals) = csr.row(genes[k]);
                let det = cols.iter().zip(vals).filter(|(&j, &v)| in_step1[j as usize] && v > 0.0).count();
                det >= params.min_cells
            })
            .collect();
        (cells, g1)
    } else {
        ((0..n_cells).collect(), (0..ng).collect())
    };
    let rng_after_cells = rng.random_seed();

    // exclude_poisson: keep the overdispersed fitting genes
    genes_step1.retain(|&k| var[k] - amean[k] > 0.0);
    let genes_step1_pool = genes_step1.clone();
    let mut log_gmean_step1: Vec<f64> = genes_step1.iter().map(|&k| log_gmean[k]).collect();

    let mut sampling = None;
    if let Some(n_genes) = params.n_genes {
        if n_genes < genes_step1.len() {
            let d = density_nrd(&log_gmean_step1, 1.0, 512).map_err(|e| SctError::Input(format!("density: {e}")))?;
            let f = ApproxFun::new(&d.x, &d.y, Method::Linear).map_err(|e| SctError::Input(format!("approx: {e}")))?;
            let prob: Vec<f64> = log_gmean_step1.iter().map(|&x| 1.0 / (f.eval(x) + f64::EPSILON)).collect();
            let draw = sample_int_weighted(rng, &prob, n_genes, false).map_err(|e| SctError::Input(e.to_string()))?;
            genes_step1 = draw.iter().map(|&k| genes_step1[k]).collect();
            log_gmean_step1 = genes_step1.iter().map(|&k| log_gmean[k]).collect();
            sampling = Some(GeneSampling { density_x: d.x, density_y: d.y, density_bw: d.bw, sampling_prob: prob });
        }
    }
    let rng_after_genes = rng.random_seed();
    let rows = GeneRows::from_csr(&csr, &genes);
    let setup = Setup {
        n_cells,
        cell_attr: CellAttr { umi, log_umi },
        genes,
        genes_cell_count,
        log_gmean,
        amean,
        var,
        cells_step1,
        genes_step1_pool,
        genes_step1,
        log_gmean_step1,
        sampling,
        rng_after_cells,
        rng_after_genes,
    };
    Ok((setup, rows))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_like_r() {
        assert_eq!(bins(707, 500), vec![0..500, 500..707]);
        assert_eq!(bins(2000, 500).len(), 4);
        assert_eq!(bins(500, 500), vec![0..500]);
    }

    #[test]
    fn row_statistics() {
        // row (1, 0, 0, 3) stored as 1, 0(explicit), 3 over 4 columns
        let v = [1.0, 0.0, 3.0];
        assert_eq!(row_mean(&v, 4), 1.0);
        // mean 1: deviations 0, 1 (explicit zero), 2 (in order), plus one implicit zero
        assert_eq!(row_var(&v, 4), (0.0 + 1.0 + 4.0 + 1.0) / 3.0);
    }
}
