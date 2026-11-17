//! Stage dumps of one `SCTransform()` call (stages S2a to S2d), in the layout of the
//! instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`, whose SCTransform hooks the
//! clustering runs reuse), so that a Rust run compares with an R run block by block. Shared by
//! DoubletFinder's runs (M6, one context per SCTransform call) and the clustering chain (M11,
//! context `chain`).
//!
//! What R's hooks in `vst`, `get_model_pars`, `future_lapply`, glmGamPoi, `reg_model_pars` and
//! the `SCTransform` methods record, here from the port's own results. Indices are R's 0-based
//! conventions: genes are rows of the input matrix, cells its columns.

use crate::{Fits, Reg, Residuals, Setup};
use scio::stages::{DumpLevel, StageRecorder};
use scio::Json;
use std::io;

fn int(v: usize) -> Json {
    Json::Int(v as i64)
}

fn ints(v: &[usize]) -> Vec<i32> {
    v.iter().map(|&x| x as i32).collect()
}

/// Records one SCTransform call in context `ctx`: the RNG checkpoints in R's order (`vst_entry`,
/// `after_cell_sampling`, `fits_entry`, `bin<k>_entry`, `fits_exit`, `vst_exit`,
/// `sctransform_exit`), the arrays `S2a_*`, `S2b_*`, `S2c_*`, `S2d_residual_stats`,
/// `S2d_top_features`, `S2d_variable_features` (and `S2d_y` at level full, when the residual
/// matrix was kept), and the scalars R's hooks set. `rng_entry` and `rng_exit` are the states at
/// `vst`'s entry and at `SCTransform`'s exit.
#[allow(clippy::too_many_arguments)]
pub fn record(rec: &mut StageRecorder, ctx: &str, setup: &Setup, fits: &Fits, mean_cell_sum: f64, reg: &Reg, res: &Residuals, rng_entry: &[i32], rng_exit: &[i32]) -> io::Result<()> {
    let (ncells, ng, n1) = (setup.n_cells, setup.genes.len(), setup.genes_step1.len());
    // RNG checkpoints in R's order
    rec.rng(ctx, "S2a", "vst_entry", rng_entry.to_vec());
    rec.rng(ctx, "S2a", "after_cell_sampling", setup.rng_after_cells.clone());
    rec.rng(ctx, "S2b", "fits_entry", setup.rng_after_genes.clone());
    for (k, s) in fits.rng_bins.iter().enumerate() {
        rec.rng(ctx, "S2b", &format!("bin{}_entry", k + 1), s.clone());
    }
    rec.rng(ctx, "S2b", "fits_exit", fits.rng_exit.clone());
    rec.rng(ctx, "S2", "vst_exit", rng_exit.to_vec());
    rec.rng(ctx, "S2", "sctransform_exit", rng_exit.to_vec());

    let full = rec.level == DumpLevel::Full;
    let dir = rec.dir.clone();
    let w = &mut rec.writer;
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
    if full {
        if let Some(y) = &res.y {
            w.write_big_f64_rows(&dir, "S2d_y", ctx, y, ng, ncells, "S2d")?;
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
        rec.scalar(ctx, k, v);
    }
    Ok(())
}
