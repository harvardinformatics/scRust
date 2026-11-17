//! Stage dumps of a full Rust run (`scrust doubletfinder --dump-stages`, milestone M6): stages
//! S1 to S7, the `.Random.seed` checkpoints and the per-context scalars, in the layout of the
//! instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`), so that
//! `scrust validate doubletfinder` compares the two runs block by block, checkpoint by
//! checkpoint and scalar by scalar.
//!
//! Each R hook has its counterpart here:
//!
//! | R hook (trace) | what it records | here |
//! |---|---|---|
//! | `paramSweep` entry | RNG `paramSweep_entry` (context `sweep`) | [`RustSource`](crate::chain::RustSource) |
//! | `parallel_paramSweep` entry | RNG `before_doublets`; `S1_subsample`; `pN`, `n_real`, `subsampled` | idem, [`Recorder::finish`] |
//! | `SCTransform.Seurat` entry | RNG `after_doublets` | idem |
//! | `vst` (statement 2, before the gene sampling, exit) | RNG `vst_entry`, `after_cell_sampling`, `vst_exit`; `S2a_*`, `S2b_model_pars`, `S2c_model_pars_fit`, `S2c_outliers`, `S2d_residual_stats`, `S2d_y` (full) | [`Recorder::sct`] |
//! | `get_model_pars`, `future_lapply` | RNG `fits_entry`, `bin<k>_entry`, `fits_exit`; `S2b_model_pars_raw`, `S2b_glm` | idem |
//! | `reg_model_pars` exit | `S2c_poisson_genes`, `S2c_smoothing_genes`, bandwidth, counts | idem |
//! | `SCTransform.default`, `.Seurat` exits | `S2d_top_features`, `S2d_variable_features`; RNG `sctransform_exit` | idem |
//! | `RunPCA.default`, `irlba` | RNG `pca_entry`, `irlba_entry`, `irlba_exit`, `pca_exit`; `S2d_pca_features`, `S2d_scale_data` (standard, full), `S3_*` | [`Recorder::pca`] |
//! | `parallel_paramSweep` exit | `S1_doublets`, `S4_digest`, `S4_positions` (full), `S5_pann_sweep`; RNG `sweep_iteration_exit` | `RustSource`, [`Recorder::finish`] |
//! | `doubletFinder` entry and exit | RNG `before_doublets`, `doubletFinder_exit`; `S1_doublets`, `S4_*`, `S5_pann_final`, `S7_doublet` | idem |
//! | `bkde`, `bimodality_coefficient` | `S6_bkde`, `S6_gkde` | [`Recorder::finish`] |
//!
//! R's structural checks of its own containers (`vst_rows_are_input_features`,
//! `scale_data_rows_are_top_features`, `variable_features_equal_top_features`,
//! `scale_data_equals_sctransform_default_y`, `pca_cells_are_vst_cells`) have no Rust
//! counterpart: they hold by construction here, and are not written.
//!
//! Levels, as R's: `light` (no residual-sized matrix, no S4 positions), `standard` (also
//! `scale.data`, the PCA input, per context), `full` (also the residual matrix `y` of every kept
//! gene and the S4 order positions).

use crate::run::Output;
use crate::stages::add_outputs;
use scio::stages::StageRecorder;
use scio::Json;
use sctransform::{Fits, Reg, Residuals, Setup};
use std::io;
use std::path::Path;

pub use scio::stages::DumpLevel;

fn int(v: usize) -> Json {
    Json::Int(v as i64)
}

fn ints(v: &[usize]) -> Vec<i32> {
    v.iter().map(|&x| x as i32).collect()
}

/// Collects a Rust run's stage dumps and writes them into `<dir>` in R's layout. Per-context
/// matrices (`S2d_scale_data`, `S2d_y`) are written as they come; everything else is small and
/// kept until [`finish`](Self::finish). The layout's machinery is [`scio::stages`]'s; the S2 and
/// S3 dumps are [`sctransform::record`]'s and [`pca::record`]'s, shared with the clustering
/// chain.
pub struct Recorder {
    inner: StageRecorder,
}

impl Recorder {
    pub fn new(dir: &Path, level: DumpLevel) -> Recorder {
        Recorder { inner: StageRecorder::new(dir, level) }
    }

    pub fn level(&self) -> DumpLevel {
        self.inner.level
    }

    /// A `.Random.seed` checkpoint.
    pub fn rng(&mut self, ctx: &str, stage: &str, label: &str, state: Vec<i32>) {
        self.inner.rng(ctx, stage, label, state);
    }

    /// A context that called SCTransform (`sweep_pN<pN>` or `final`), in run order.
    pub fn context(&mut self, ctx: &str) {
        self.inner.context(ctx);
    }

    /// A scalar of context `ctx` (R's `dfi_scalar`).
    pub fn scalar(&mut self, ctx: &str, name: &str, value: Json) {
        self.inner.scalar(ctx, name, value);
    }

    /// The sweep subsample (0-based input cells, in draw order).
    pub fn s1_subsample(&mut self, cells: &[usize]) {
        self.inner.writer.add_i32("S1_subsample", "sweep", ints(cells), cells.len(), 1, "S1");
    }

    /// One context's doublet pairs (0-based positions in the real cells).
    pub fn s1_doublets(&mut self, ctx: &str, pairs: &[(usize, usize)]) {
        let mut v: Vec<i32> = pairs.iter().map(|p| p.0 as i32).collect();
        v.extend(pairs.iter().map(|p| p.1 as i32));
        self.inner.writer.add_i32("S1_doublets", ctx, v, pairs.len(), 2, "S1");
    }

    /// Stage S2 of one context ([`sctransform::record::record`]): what R's hooks in `vst`,
    /// `get_model_pars`, `future_lapply`, glmGamPoi, `reg_model_pars` and the `SCTransform`
    /// methods record. `rng_entry` and `rng_exit` are the states at `vst`'s entry and at
    /// `SCTransform`'s exit.
    #[allow(clippy::too_many_arguments)]
    pub fn sct(&mut self, ctx: &str, setup: &Setup, fits: &Fits, mean_cell_sum: f64, reg: &Reg, res: &Residuals, rng_entry: &[i32], rng_exit: &[i32]) -> io::Result<()> {
        sctransform::record::record(&mut self.inner, ctx, setup, fits, mean_cell_sum, reg, res, rng_entry, rng_exit)
    }

    /// Stage S3 of one context ([`pca::record::record`]): RunPCA on `scale.data` (rows
    /// `res.top_features`, whose input rows are `rows`), what R's hooks in `RunPCA.default` and
    /// `irlba` record. `entry` and `exit` are the generator's states at RunPCA's entry and exit.
    pub fn pca(&mut self, ctx: &str, rows: &[usize], res: &Residuals, p: &pca::Pca, entry: &[i32], exit: &[i32]) -> io::Result<()> {
        pca::record::record(&mut self.inner, ctx, rows, &res.scale_data, p, Some(entry), exit)
    }

    /// Adds stages S4-S7 and the sweep's and final run's scalars from the run's output, then
    /// writes everything: the arrays, `index.tsv`, `rng.npy` and `rng.tsv`, `stages.json`.
    pub fn finish(mut self, out: &Output) -> io::Result<()> {
        add_outputs(&mut self.inner.writer, out);
        if let Some(sweep) = &out.sweep {
            let subsampled = sweep.real_cells.len() < out.n_cells;
            for (it, pn) in sweep.iterations.iter().zip(&sweep.pn) {
                let ctx = format!("sweep_pN{}", pn.label);
                let pk: Vec<f64> = sweep.pk.iter().map(|p| p.value).collect();
                for (k, v) in [
                    ("pN", Json::Num(pn.value)),
                    ("n_real", int(it.n_real)),
                    ("subsampled", Json::Bool(subsampled)),
                    ("n_doublets", int(it.n_cells - it.n_real)),
                    ("nCells", int(it.n_cells)),
                    ("ind", int(it.ind)),
                    ("pK", Json::nums(&pk)),
                    ("pk_temp", Json::nums(&it.pk_temp.iter().map(|&k| k as f64).collect::<Vec<f64>>())),
                    ("first_position_is_self", int(it.first_is_self)),
                ] {
                    self.scalar(&ctx, k, v);
                }
            }
        }
        let f = &out.final_run;
        for (k, v) in [
            ("n_doublets", int(out.n_doublets)),
            ("nCells", int(f.n_cells)),
            ("k", Json::Num(f.k)),
            ("pK", Json::Num(out.pk)),
            ("nExp", Json::Num(out.nexp)),
            ("first_position_is_self", int(f.first_is_self)),
        ] {
            self.scalar("final", k, v);
        }
        self.inner.finish(&format!("scrust doubletfinder {}", env!("CARGO_PKG_VERSION")))
    }
}
