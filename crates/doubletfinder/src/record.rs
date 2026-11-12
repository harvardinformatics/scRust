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
use crate::stages::{add_outputs, StageWriter};
use scio::npy;
use scio::{Json, Table};
use sctransform::{Fits, Reg, Residuals, Setup};
use std::io;
use std::path::{Path, PathBuf};

/// What a run dumps, as the instrumented R runs' levels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DumpLevel {
    Light,
    Standard,
    Full,
}

impl DumpLevel {
    /// `none` (no dumps), `light`, `standard` or `full`.
    pub fn parse(s: &str) -> Result<Option<DumpLevel>, String> {
        match s {
            "none" => Ok(None),
            "light" => Ok(Some(DumpLevel::Light)),
            "standard" => Ok(Some(DumpLevel::Standard)),
            "full" => Ok(Some(DumpLevel::Full)),
            _ => Err(format!("unknown dump level '{s}' (none, light, standard or full)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            DumpLevel::Light => "light",
            DumpLevel::Standard => "standard",
            DumpLevel::Full => "full",
        }
    }
}

/// One `.Random.seed` checkpoint: context, stage, label, state.
#[derive(Debug, Clone)]
struct Checkpoint {
    ctx: String,
    stage: String,
    label: String,
    state: Vec<i32>,
}

fn int(v: usize) -> Json {
    Json::Int(v as i64)
}

fn ints(v: &[usize]) -> Vec<i32> {
    v.iter().map(|&x| x as i32).collect()
}

/// Collects a Rust run's stage dumps and writes them into `<dir>` in R's layout. Per-context
/// matrices (`S2d_scale_data`, `S2d_y`) are written as they come; everything else is small and
/// kept until [`finish`](Self::finish).
pub struct Recorder {
    level: DumpLevel,
    dir: PathBuf,
    writer: StageWriter,
    rng: Vec<Checkpoint>,
    /// per context (in order of first use): scalar name -> value
    scalars: Vec<(String, Vec<(String, Json)>)>,
    contexts: Vec<String>,
}

impl Recorder {
    pub fn new(dir: &Path, level: DumpLevel) -> Recorder {
        Recorder { level, dir: dir.to_path_buf(), writer: StageWriter::default(), rng: Vec::new(), scalars: Vec::new(), contexts: Vec::new() }
    }

    pub fn level(&self) -> DumpLevel {
        self.level
    }

    /// A `.Random.seed` checkpoint.
    pub fn rng(&mut self, ctx: &str, stage: &str, label: &str, state: Vec<i32>) {
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: stage.into(), label: label.into(), state });
    }

    /// A context that called SCTransform (`sweep_pN<pN>` or `final`), in run order.
    pub fn context(&mut self, ctx: &str) {
        if !self.contexts.iter().any(|c| c == ctx) {
            self.contexts.push(ctx.into());
        }
    }

    /// A scalar of context `ctx` (R's `dfi_scalar`).
    pub fn scalar(&mut self, ctx: &str, name: &str, value: Json) {
        let i = match self.scalars.iter().position(|(c, _)| c == ctx) {
            Some(i) => i,
            None => {
                self.scalars.push((ctx.into(), Vec::new()));
                self.scalars.len() - 1
            }
        };
        let m = &mut self.scalars[i].1;
        match m.iter_mut().find(|(k, _)| k == name) {
            Some(slot) => slot.1 = value,
            None => m.push((name.into(), value)),
        }
    }

    /// The sweep subsample (0-based input cells, in draw order).
    pub fn s1_subsample(&mut self, cells: &[usize]) {
        self.writer.add_i32("S1_subsample", "sweep", ints(cells), cells.len(), 1, "S1");
    }

    /// One context's doublet pairs (0-based positions in the real cells).
    pub fn s1_doublets(&mut self, ctx: &str, pairs: &[(usize, usize)]) {
        let mut v: Vec<i32> = pairs.iter().map(|p| p.0 as i32).collect();
        v.extend(pairs.iter().map(|p| p.1 as i32));
        self.writer.add_i32("S1_doublets", ctx, v, pairs.len(), 2, "S1");
    }

    /// Stage S2 of one context: what R's hooks in `vst`, `get_model_pars`, `future_lapply`,
    /// glmGamPoi, `reg_model_pars` and the `SCTransform` methods record. `rng_entry` and
    /// `rng_exit` are the states at `vst`'s entry and at `SCTransform`'s exit.
    #[allow(clippy::too_many_arguments)]
    pub fn sct(&mut self, ctx: &str, setup: &Setup, fits: &Fits, mean_cell_sum: f64, reg: &Reg, res: &Residuals, rng_entry: &[i32], rng_exit: &[i32]) -> io::Result<()> {
        let w = &mut self.writer;
        let (ncells, ng, n1) = (setup.n_cells, setup.genes.len(), setup.genes_step1.len());
        // RNG checkpoints in R's order
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2a".into(), label: "vst_entry".into(), state: rng_entry.to_vec() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2a".into(), label: "after_cell_sampling".into(), state: setup.rng_after_cells.clone() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2b".into(), label: "fits_entry".into(), state: setup.rng_after_genes.clone() });
        for (k, s) in fits.rng_bins.iter().enumerate() {
            self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2b".into(), label: format!("bin{}_entry", k + 1), state: s.clone() });
        }
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2b".into(), label: "fits_exit".into(), state: fits.rng_exit.clone() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2".into(), label: "vst_exit".into(), state: rng_exit.to_vec() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S2".into(), label: "sctransform_exit".into(), state: rng_exit.to_vec() });

        // S2a
        w.add_i32("S2a_genes_step1_pool", ctx, ints(&setup.genes_step1_pool), setup.genes_step1_pool.len(), 1, "S2a");
        w.add_i32("S2a_genes", ctx, ints(&setup.genes), ng, 1, "S2a");
        let mut cells = setup.cell_attr.umi.clone();
        cells.extend_from_slice(&setup.cell_attr.log_umi);
        w.add_f64("S2a_cells", ctx, cells, ncells, 2, "S2a");
        let mut gs = setup.log_gmean.clone();
        gs.extend_from_slice(&setup.amean);
        gs.extend_from_slice(&setup.var);
        gs.extend(setup.detection_rate());
        gs.extend(setup.gmean());
        w.add_f64("S2a_gene_stats", ctx, gs, ng, 5, "S2a");
        w.add_i32("S2a_cells_step1", ctx, ints(&setup.cells_step1), setup.cells_step1.len(), 1, "S2a");
        w.add_i32("S2a_genes_step1", ctx, ints(&setup.genes_step1), n1, 1, "S2a");
        w.add_f64("S2a_log_gmean_step1", ctx, setup.log_gmean_step1.clone(), n1, 1, "S2a");
        if let Some(s) = &setup.sampling {
            w.add_f64("S2a_sampling_prob", ctx, s.sampling_prob.clone(), s.sampling_prob.len(), 1, "S2a");
            let mut d = s.density_x.clone();
            d.extend_from_slice(&s.density_y);
            w.add_f64("S2a_density", ctx, d, s.density_x.len(), 2, "S2a");
        }
        // S2b
        let pars = |p: &[[f64; 3]]| -> Vec<f64> { (0..3).flat_map(|k| p.iter().map(move |x| x[k])).collect() };
        w.add_f64("S2b_model_pars", ctx, pars(&fits.model_pars), fits.model_pars.len(), 3, "S2b");
        w.add_f64("S2b_model_pars_raw", ctx, pars(&fits.model_pars_raw), fits.model_pars_raw.len(), 3, "S2b");
        let g = &fits.glm;
        let mut glm: Vec<f64> = Vec::with_capacity(9 * g.len());
        glm.extend(g.iter().map(|x| x.disp_init));
        glm.extend(g.iter().map(|x| x.beta_init));
        glm.extend(g.iter().map(|x| x.beta1));
        glm.extend(g.iter().map(|x| x.iter1 as f64));
        glm.extend(g.iter().map(|x| x.disp_mle));
        glm.extend(g.iter().map(|x| x.disp_trend));
        glm.extend(g.iter().map(|x| x.beta2));
        glm.extend(g.iter().map(|x| x.iter2 as f64));
        glm.extend(g.iter().map(|x| x.deviance2));
        w.add_f64("S2b_glm", ctx, glm, g.len(), 9, "S2b");
        // S2c
        w.add_f64("S2c_model_pars_fit", ctx, pars(&reg.model_pars_fit), reg.model_pars_fit.len(), 3, "S2c");
        w.add_i32("S2c_outliers", ctx, reg.outliers.iter().map(|&o| o as i32).collect(), reg.outliers.len(), 1, "S2c");
        w.add_i32("S2c_poisson_genes", ctx, ints(&reg.poisson_genes), reg.poisson_genes.len(), 1, "S2c");
        w.add_i32("S2c_smoothing_genes", ctx, ints(&reg.smoothing_genes), reg.smoothing_genes.len(), 1, "S2c");
        // S2d
        let mut rs = res.residual_mean.clone();
        rs.extend_from_slice(&res.residual_variance);
        w.add_f64("S2d_residual_stats", ctx, rs, ng, 2, "S2d");
        let top: Vec<i32> = res.top_features.iter().map(|&k| setup.genes[k] as i32).collect();
        w.add_i32("S2d_top_features", ctx, top.clone(), top.len(), 1, "S2d");
        w.add_i32("S2d_variable_features", ctx, top.clone(), top.len(), 1, "S2d");
        if self.level == DumpLevel::Full {
            if let Some(y) = &res.y {
                w.write_big_f64_rows(&self.dir, "S2d_y", ctx, y, ng, ncells, "S2d")?;
            }
        }

        // scalars, as the hooks set them
        let raw_theta = fits.model_pars_raw.iter().map(|p| p[0]);
        let sampled = setup.sampling.as_ref();
        let mut sc: Vec<(&str, Json)> = vec![("vst_n_cells", int(ncells))];
        if let Some(s) = sampled {
            sc.push(("density_bw", Json::Num(s.density_bw)));
            sc.push(("density_n", int(s.sampling_prob.len())));
        }
        sc.extend([
            ("min_variance", Json::Num(res.min_var)),
            ("n_genes", int(ng)),
            ("n_genes_step1", int(n1)),
            ("n_cells_step1", int(setup.cells_step1.len())),
            ("n_genes_step1_pool", int(setup.genes_step1_pool.len())),
            ("gene_sampling", Json::Bool(sampled.is_some())),
            ("n_cells", int(ncells)),
            ("res_clip_range", Json::nums(&[res.res_clip.0, res.res_clip.1])),
            ("n_bins_fit", int(fits.rng_bins.len())),
            ("n_theta_below_1e7", int(raw_theta.clone().filter(|&t| t < 1e-7).count())),
            ("n_theta_inf", int(raw_theta.filter(|t| !t.is_finite()).count())),
            ("reg_bw", Json::Num(reg.bw)),
            ("mean_cell_sum", Json::Num(mean_cell_sum)),
            ("n_outliers", int(reg.outliers.iter().filter(|&&o| o).count())),
            ("n_all_poisson_genes", int(reg.poisson_genes.len())),
            ("n_smoothing_genes", int(reg.smoothing_genes.len())),
            ("clip_range", Json::Num(res.clip.1)),
        ]);
        for (k, v) in sc {
            self.scalar(ctx, k, v);
        }
        Ok(())
    }

    /// Stage S3 of one context: RunPCA on `scale.data` (rows `res.top_features`, whose input
    /// rows are `rows`), what R's hooks in `RunPCA.default` and `irlba` record. `entry` and
    /// `exit` are the generator's states at RunPCA's entry and exit.
    pub fn pca(&mut self, ctx: &str, rows: &[usize], res: &Residuals, p: &pca::Pca, entry: &[i32], exit: &[i32]) -> io::Result<()> {
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S3".into(), label: "pca_entry".into(), state: entry.to_vec() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S3".into(), label: "irlba_entry".into(), state: p.rng_irlba_entry.clone() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S3".into(), label: "irlba_exit".into(), state: exit.to_vec() });
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: "S3".into(), label: "pca_exit".into(), state: exit.to_vec() });
        let w = &mut self.writer;
        let feats: Vec<i32> = p.features.iter().map(|&k| rows[k] as i32).collect();
        w.add_i32("S2d_pca_features", ctx, feats.clone(), feats.len(), 1, "S2d");
        if self.level >= DumpLevel::Standard {
            // the matrix RunPCA.default receives: scale.data's rows that PrepDR keeps
            let n = res.ncells;
            if p.features.len() == res.top_features.len() {
                w.write_big_f64_rows(&self.dir, "S2d_scale_data", ctx, &res.scale_data, feats.len(), n, "S2d")?;
            } else {
                let kept: Vec<f64> = p.features.iter().flat_map(|&k| res.scale_row(k).iter().copied()).collect();
                w.write_big_f64_rows(&self.dir, "S2d_scale_data", ctx, &kept, feats.len(), n, "S2d")?;
            }
        }
        w.add_f64("S3_embeddings", ctx, p.embeddings.clone(), p.ncells, p.npcs, "S3");
        w.add_f64("S3_loadings", ctx, p.loadings.clone(), p.features.len(), p.npcs, "S3");
        w.add_f64("S3_sdev", ctx, p.sdev.clone(), 1, p.npcs, "S3");
        w.add_f64("S3_irlba_d", ctx, p.irlba.d.clone(), 1, p.npcs, "S3");
        for (k, v) in [
            ("pca_n_features", int(feats.len())),
            ("total_variance", Json::Num(p.total_variance)),
            ("irlba_iter", int(p.irlba.iter)),
            ("irlba_mprod", int(p.irlba.mprod)),
        ] {
            self.scalar(ctx, k, v);
        }
        Ok(())
    }

    /// Adds stages S4-S7 and the sweep's and final run's scalars from the run's output, then
    /// writes everything: the arrays, `index.tsv`, `rng.npy` and `rng.tsv`, `stages.json`.
    pub fn finish(mut self, out: &Output) -> io::Result<()> {
        add_outputs(&mut self.writer, out);
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

        let dir = self.dir.clone();
        self.writer.write(&dir)?;
        // rng.npy: one row of 626 per checkpoint (column-major, as R writes the matrix)
        let len = self.rng.first().map_or(0, |c| c.state.len());
        if self.rng.iter().any(|c| c.state.len() != len) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "RNG checkpoints of different lengths"));
        }
        let n = self.rng.len();
        let mut m = Vec::with_capacity(n * len);
        for j in 0..len {
            m.extend(self.rng.iter().map(|c| c.state[j]));
        }
        npy::write_i32(&dir.join("rng.npy"), &m, &[n, len])?;
        Table {
            header: ["row", "ctx", "stage", "label"].iter().map(|s| s.to_string()).collect(),
            rows: self.rng.iter().enumerate().map(|(i, c)| vec![i.to_string(), c.ctx.clone(), c.stage.clone(), c.label.clone()]).collect(),
        }
        .write(&dir.join("rng.tsv"))?;
        let scalars = Json::Object(self.scalars.iter().map(|(c, m)| (c.clone(), Json::Object(m.clone()))).collect());
        let json = Json::obj(vec![
            ("level", Json::Str(self.level.name().into())),
            ("contexts", Json::strs(&self.contexts)),
            ("scalars", scalars),
            ("producer", Json::Str(format!("scrust doubletfinder {}", env!("CARGO_PKG_VERSION")))),
        ]);
        std::fs::write(dir.join("stages.json"), scio::json::to_file_text(&json))
    }
}
