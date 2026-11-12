//! Truncated PCA, a port of irlba 2.4.1 as Seurat 5.5.1's `RunPCA` calls it (docs/PLAN.md B3b;
//! milestone M5).
//!
//! The call path is `RunPCA.Seurat` → `RunPCA.Assay` → `PrepDR` → `RunPCA.default` (the SCT
//! assay is a v3 `Assay`). [`run_pca`] reproduces it on the rows of `scale.data` for the
//! variable features:
//!
//! - `PrepDR`: drops features whose `RowVar` is not positive ([`rowvar`], Seurat's Eigen code);
//! - `RunPCA.default`: `set.seed(seed.use)`, `total.variance <- sum(RowVar(object))`,
//!   `npcs <- min(npcs, nrow(object) - 1)`, `irlba(A = t(object), nv = npcs)` ([`irlba`]), then
//!   `sdev = d / sqrt(max(1, ncol(object) - 1))`, loadings `v` and embeddings `u %*% diag(d)`
//!   (`weight.by.var = TRUE`).
//!
//! Every product is the OpenBLAS call R makes, through FFI ([`rcompat::blas`]); the two large
//! ones are split across threads only where the kernels give the single call's bits
//! ([`products`]), so the results do not depend on the number of threads. irlba draws its start
//! vector, and in rare cases more vectors, from R's stream, so the generator continues as in R.
//!
//! Shared by DoubletFinder (7 calls per object, `npcs = 10`) and, later, the clustering crate.

pub mod irlba;
pub mod products;
pub mod rowvar;

pub use irlba::{irlba, Irlba, IrlbaParams};
pub use products::Products;

use rcompat::RRng;
use std::fmt;

/// The original implementation this crate reproduces.
pub const REFERENCE: &str = "irlba 2.4.1, as Seurat 5.5.1's RunPCA calls it";

/// Why a run stops.
#[derive(Debug)]
pub enum PcaError {
    /// Invalid input or parameters, or a condition under which R stops.
    Input(String),
}

impl fmt::Display for PcaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PcaError::Input(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for PcaError {}

/// `RunPCA()`'s arguments (Seurat's defaults; DoubletFinder passes `npcs = 10`).
#[derive(Debug, Clone, PartialEq)]
pub struct PcaParams {
    /// `npcs` (default 50)
    pub npcs: usize,
    /// `seed.use` (default 42; `None` sets no seed)
    pub seed_use: Option<i32>,
    /// `weight.by.var` (default `TRUE`): embeddings `u %*% diag(d)` rather than `u`
    pub weight_by_var: bool,
    /// irlba's own arguments (its defaults)
    pub irlba: IrlbaParams,
}

impl Default for PcaParams {
    fn default() -> Self {
        PcaParams { npcs: 50, seed_use: Some(42), weight_by_var: true, irlba: IrlbaParams::default() }
    }
}

/// `RunPCA`'s results, and what the instrumented R runs record of it.
#[derive(Debug, Clone)]
pub struct Pca {
    /// the input rows `PrepDR` keeps (positive `RowVar`), in order: the loadings' rows
    pub features: Vec<usize>,
    pub ncells: usize,
    /// `npcs` after `min(npcs, nrow(object) - 1)`
    pub npcs: usize,
    /// cell embeddings, `ncells x npcs`, column-major
    pub embeddings: Vec<f64>,
    /// feature loadings, `features.len() x npcs`, column-major
    pub loadings: Vec<f64>,
    /// standard deviations of the components
    pub sdev: Vec<f64>,
    /// `misc$total.variance`
    pub total_variance: f64,
    /// irlba's result (`d`, `iter`, `mprod`, ...); its `u` and `v` are moved into the embeddings
    /// and loadings
    pub irlba: Irlba,
    /// R's `.Random.seed` when irlba starts (after `set.seed(seed.use)`)
    pub rng_irlba_entry: Vec<i32>,
}

impl Pca {
    /// Component `c` (0-based) of the embeddings.
    pub fn embedding(&self, c: usize) -> &[f64] {
        &self.embeddings[c * self.ncells..(c + 1) * self.ncells]
    }
}

/// `RunPCA(object, npcs)` on an SCT assay: `rows` holds the `scale.data` rows of the variable
/// features in `VariableFeatures` order, each over the `ncells` cells (row-major
/// `nfeatures x ncells`, which is `t(scale.data)` column-major, irlba's `A`). Draws from `rng` as
/// R does; `threads` only changes the speed.
pub fn run_pca(rows: &[f64], nfeatures: usize, ncells: usize, params: &PcaParams, rng: &mut RRng, threads: usize) -> Result<Pca, PcaError> {
    if rows.len() != nfeatures * ncells {
        return Err(PcaError::Input(format!("run_pca: {} values for {nfeatures} features x {ncells} cells", rows.len())));
    }
    if nfeatures < 2 || ncells < 2 {
        return Err(PcaError::Input(format!("run_pca: needs at least 2 features and 2 cells, not {nfeatures} x {ncells}")));
    }
    // PrepDR: features.keep <- features[RowVar(data.use[features, ]) > 0] (NA dropped too)
    let var: Vec<f64> = (0..nfeatures).map(|f| rowvar::row_var(&rows[f * ncells..(f + 1) * ncells])).collect();
    let features: Vec<usize> = (0..nfeatures).filter(|&f| var[f] > 0.0).collect();
    if features.is_empty() {
        return Err(PcaError::Input("None of the requested features have any variance".into()));
    }
    let compact;
    let a: &[f64] = if features.len() == nfeatures {
        rows
    } else {
        compact = features.iter().flat_map(|&f| rows[f * ncells..(f + 1) * ncells].iter().copied()).collect::<Vec<f64>>();
        &compact
    };
    let nf = features.len();

    // RunPCA.default
    if let Some(seed) = params.seed_use {
        rng.set_seed(seed);
    }
    // total.variance <- sum(RowVar(object)): the kept rows' variances again, R's long-double sum
    let kept_var: Vec<f64> = features.iter().map(|&f| var[f]).collect();
    let total_variance = rcompat::sum(&kept_var, false);
    if params.npcs == 0 {
        return Err(PcaError::Input("npcs must be positive".into()));
    }
    let npcs = params.npcs.min(nf - 1);
    if params.weight_by_var && npcs < 2 {
        // diag(d) of a single d is an identity matrix of size d, so u %*% diag(d) stops in R
        // unless 1 <= d < 2; not reproduced
        return Err(PcaError::Input(format!("RunPCA with {npcs} component(s) and weight.by.var: u %*% diag(d) of a single d is not supported")));
    }
    let rng_irlba_entry = rng.random_seed();
    let prod = Products::new(a, ncells, nf, threads);
    let mut ir = irlba(&prod, npcs, &params.irlba, rng)?;
    // sdev <- pca.results$d / sqrt(max(1, ncol(object) - 1))
    let denom = (ncells as f64 - 1.0).max(1.0).sqrt();
    let sdev: Vec<f64> = ir.d.iter().map(|d| d / denom).collect();
    let u = std::mem::take(&mut ir.u);
    let embeddings = if params.weight_by_var {
        // cell.embeddings <- pca.results$u %*% diag(pca.results$d)
        let mut dg = vec![0.0; npcs * npcs];
        for i in 0..npcs {
            dg[i + i * npcs] = ir.d[i];
        }
        rcompat::blas::matprod(&u, ncells, npcs, &dg, npcs, npcs)
    } else {
        u
    };
    let loadings = std::mem::take(&mut ir.v);
    Ok(Pca { features, ncells, npcs, embeddings, loadings, sdev, total_variance, irlba: ir, rng_irlba_entry })
}
