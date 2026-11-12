//! `SctParams`: the arguments of Seurat 5.5.1's `SCTransform()` and of sctransform 0.4.3's
//! `vst()` that the port implements, with their defaults as DoubletFinder 2.0.6 calls
//! `SCTransform(seu)`.
//!
//! Only the defaults are validated against R (milestone M4). Other values are accepted where the
//! code path is the same one with a different number (cells and genes for the fits, bin size,
//! clipping, bandwidth adjustment, number of variable features, seed); options that would take R
//! down another path (other methods, batch variables, `do.scale`, extra regressors) are rejected
//! by [`SctParams::check`].

/// How `min_variance` is set (`vst(min_variance = )`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MinVariance {
    /// `"umi_median"`: `(median(umi@x) / 5)^2` over the stored entries of the filtered counts
    /// (the default of `vst.flavor = "v2"`).
    UmiMedian,
    /// a fixed value (`-Inf` disables the floor)
    Value(f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SctParams {
    // --- SCTransform (Seurat) ---
    /// `seed.use`: `set.seed()` at the start of `SCTransform.Seurat`, `.Assay` and `.default`
    /// (`None`: no reseeding).
    pub seed_use: Option<i32>,
    /// `ncells`: cells used to fit the per-gene models (`vst(n_cells = min(ncells, ncol))`).
    pub ncells: usize,
    /// `variable.features.n`: number of variable features (top residual variance).
    pub variable_features_n: Option<usize>,
    /// `variable.features.rv.th`: residual-variance threshold when `variable_features_n` is None.
    pub variable_features_rv_th: f64,
    /// `do.center`: centre `scale.data` rows (through `ScaleData`).
    pub do_center: bool,
    /// `clip.range`: `scale.data` clipping; `None` is `c(-sqrt(N/30), sqrt(N/30))`.
    pub clip_range: Option<(f64, f64)>,
    /// `return.only.var.genes`: `scale.data` holds the variable features only.
    pub return_only_var_genes: bool,

    // --- vst (sctransform, vst.flavor = "v2") ---
    /// `n_genes`: genes used to fit the models (density-weighted sample when more qualify).
    pub n_genes: Option<usize>,
    /// `bin_size`: genes per bin, for the fits and for the residuals.
    pub bin_size: usize,
    /// `min_cells`: genes detected in fewer cells are dropped.
    pub min_cells: usize,
    /// `res_clip_range`: residual clipping inside `vst`; `None` is `c(-sqrt(N), sqrt(N))`.
    pub res_clip_range: Option<(f64, f64)>,
    /// `min_variance`.
    pub min_variance: MinVariance,
    /// `bw_adjust`: multiplier of `bw.SJ` in the regularisation.
    pub bw_adjust: f64,
    /// `gmean_eps`: pseudocount of the geometric mean.
    pub gmean_eps: f64,
    /// `exclude_poisson` (TRUE in v2).
    pub exclude_poisson: bool,
    /// The `future` plan's number of workers (`future::nbrOfWorkers()`); sctransform splits each
    /// fitting bin into this many chunks. DoubletFinder runs with R's default sequential plan: 1.
    pub n_workers: usize,
}

impl Default for SctParams {
    fn default() -> SctParams {
        SctParams {
            seed_use: Some(1448145),
            ncells: 5000,
            variable_features_n: Some(3000),
            variable_features_rv_th: 1.3,
            do_center: true,
            clip_range: None,
            return_only_var_genes: true,
            n_genes: Some(2000),
            bin_size: 500,
            min_cells: 5,
            res_clip_range: None,
            min_variance: MinVariance::UmiMedian,
            bw_adjust: 3.0,
            gmean_eps: 1.0,
            exclude_poisson: true,
            n_workers: 1,
        }
    }
}

impl SctParams {
    /// Rejects values this port does not implement.
    pub fn check(&self) -> Result<(), String> {
        if self.bin_size == 0 {
            return Err("bin_size must be positive".into());
        }
        if !self.exclude_poisson {
            return Err("exclude_poisson = FALSE is not implemented (vst.flavor = \"v2\" sets it)".into());
        }
        if self.n_workers != 1 {
            return Err("only a sequential future plan (one worker) is implemented".into());
        }
        if self.gmean_eps != 1.0 {
            return Err("gmean_eps other than 1 is not implemented".into());
        }
        if !(self.bw_adjust > 0.0) {
            return Err("bw_adjust must be positive".into());
        }
        Ok(())
    }

    /// `vst(n_cells = min(ncells, ncol(umi)))`.
    pub fn n_cells_fit(&self, ncol: usize) -> usize {
        self.ncells.min(ncol)
    }

    /// vst's residual clipping range for `n` cells.
    pub fn res_clip(&self, n: usize) -> (f64, f64) {
        self.res_clip_range.unwrap_or_else(|| {
            let s = (n as f64).sqrt();
            (-s, s)
        })
    }

    /// `SCTransform`'s `scale.data` clipping range for `n` cells.
    pub fn sct_clip(&self, n: usize) -> (f64, f64) {
        self.clip_range.unwrap_or_else(|| {
            let s = (n as f64 / 30.0).sqrt();
            (-s, s)
        })
    }
}
