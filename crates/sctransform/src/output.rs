//! `scrust sctransform`'s output layout (docs/PLAN.md C3), written by the command and by the
//! clustering chain (`scrust cluster --keep-stages`, M11), so that both give the same bytes:
//!
//! - `cell_attr.tsv`: barcode, umi, log_umi;
//! - `gene_attr.tsv`: the genes vst keeps: id, name, detection_rate, gmean, amean, variance,
//!   residual_mean, residual_variance;
//! - `model_pars.tsv`: the fitting genes' parameters (theta after the 1e-7 floor, (Intercept),
//!   log_umi), in sampling order;
//! - `model_pars_fit.tsv`: the regularised parameters of every kept gene;
//! - `variable_features.tsv`: id and name of the variable features, highest residual variance
//!   first;
//! - `scale_data.npy`: scale.data, variable features x cells, float64, column-major (R's layout).
//!
//! Doubles are written with 17 significant digits.

use crate::{SctParams, Sct};
use scio::{fmt17, Json};
use std::fmt::Write as _;
use std::io;
use std::path::Path;

fn write(path: &Path, text: &str) -> io::Result<()> {
    std::fs::write(path, text).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
}

/// A gene's id and name (the input's feature lines).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

/// Writes the layout into `dir` (which must exist): `barcodes` are the input's cells, `genes`
/// the input's rows. Returns the names written.
pub fn write_outputs(dir: &Path, barcodes: &[String], genes: &[Gene], sct: &Sct) -> io::Result<Vec<String>> {
    let s = &sct.setup;
    let res = &sct.residuals;
    let mut t = String::from("barcode\tumi\tlog_umi\n");
    for (c, b) in barcodes.iter().enumerate() {
        let _ = writeln!(t, "{b}\t{}\t{}", fmt17(s.cell_attr.umi[c]), fmt17(s.cell_attr.log_umi[c]));
    }
    write(&dir.join("cell_attr.tsv"), &t)?;
    let gene = |k: usize| &genes[s.genes[k]];
    let (detection, gmean) = (s.detection_rate(), s.gmean());
    let mut t = String::from("id\tname\tdetection_rate\tgmean\tamean\tvariance\tresidual_mean\tresidual_variance\n");
    for k in 0..s.genes.len() {
        let f = gene(k);
        let _ = writeln!(
            t,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.id,
            f.name,
            fmt17(detection[k]),
            fmt17(gmean[k]),
            fmt17(s.amean[k]),
            fmt17(s.var[k]),
            fmt17(res.residual_mean[k]),
            fmt17(res.residual_variance[k])
        );
    }
    write(&dir.join("gene_attr.tsv"), &t)?;
    let pars = |rows: &mut dyn Iterator<Item = (usize, [f64; 3])>| {
        let mut t = String::from("id\tname\ttheta\t(Intercept)\tlog_umi\n");
        for (k, p) in rows {
            let f = gene(k);
            let _ = writeln!(t, "{}\t{}\t{}\t{}\t{}", f.id, f.name, fmt17(p[0]), fmt17(p[1]), fmt17(p[2]));
        }
        t
    };
    write(&dir.join("model_pars.tsv"), &pars(&mut s.genes_step1.iter().copied().zip(sct.fits.model_pars.iter().copied())))?;
    write(&dir.join("model_pars_fit.tsv"), &pars(&mut sct.reg.model_pars_fit.iter().copied().enumerate()))?;
    let mut t = String::from("id\tname\n");
    for &k in &res.top_features {
        let _ = writeln!(t, "{}\t{}", gene(k).id, gene(k).name);
    }
    write(&dir.join("variable_features.tsv"), &t)?;
    // scale.data, column-major as R stores it
    let (nr, nc) = (res.top_features.len(), res.ncells);
    scio::npy::write_f64_rows(&dir.join("scale_data.npy"), &res.scale_data, nr, nc)?;
    Ok(["cell_attr.tsv", "gene_attr.tsv", "model_pars.tsv", "model_pars_fit.tsv", "variable_features.tsv", "scale_data.npy"].iter().map(|s| s.to_string()).collect())
}

/// The parameters as `summary.json` records them.
pub fn params_json(params: &SctParams) -> Json {
    Json::obj(vec![
        ("seed_use", Json::Int(params.seed_use.unwrap_or(0) as i64)),
        ("ncells", Json::Int(params.ncells as i64)),
        ("variable_features_n", Json::Int(params.variable_features_n.unwrap_or(0) as i64)),
        ("n_genes", Json::Int(params.n_genes.map_or(0, |v| v as i64))),
        ("bin_size", Json::Int(params.bin_size as i64)),
        ("min_cells", Json::Int(params.min_cells as i64)),
        ("bw_adjust", Json::Num(params.bw_adjust)),
    ])
}

/// The `result` section of `summary.json` (deterministic: equal between runs that agree).
pub fn result_json(sct: &Sct, n_features: usize) -> Json {
    let s = &sct.setup;
    let res = &sct.residuals;
    let n_out = sct.reg.outliers.iter().filter(|&&o| o).count();
    Json::obj(vec![
        ("n_cells", Json::Int(s.n_cells as i64)),
        ("n_features", Json::Int(n_features as i64)),
        ("n_genes_kept", Json::Int(s.genes.len() as i64)),
        ("n_cells_step1", Json::Int(s.cells_step1.len() as i64)),
        ("n_genes_step1", Json::Int(s.genes_step1.len() as i64)),
        ("gene_sampling", Json::Bool(s.sampling.is_some())),
        ("n_theta_inf", Json::Int(sct.fits.model_pars_raw.iter().filter(|p| !p[0].is_finite()).count() as i64)),
        ("mean_cell_sum", Json::Num(sct.mean_cell_sum)),
        ("reg_bw", Json::Num(sct.reg.bw)),
        ("n_outliers", Json::Int(n_out as i64)),
        ("n_poisson_genes", Json::Int(sct.reg.poisson_genes.len() as i64)),
        ("n_smoothing_genes", Json::Int(sct.reg.smoothing_genes.len() as i64)),
        ("min_variance", Json::Num(res.min_var)),
        ("res_clip_range", Json::nums(&[res.res_clip.0, res.res_clip.1])),
        ("clip_range", Json::nums(&[res.clip.0, res.clip.1])),
        ("n_variable_features", Json::Int(res.top_features.len() as i64)),
    ])
}

/// SHA-256 of `scale.data` as R's reference runs digest it (`raw_sha256(sd)` in
/// `validation/R/lib/clustering.R`): the doubles, little-endian, column-major, with the rows in
/// the SCT assay's gene order (the kept genes' order, which Seurat keeps in `scale.data`; not
/// `VariableFeatures()` order). `res.top_features` are positions in the kept genes.
pub fn scale_data_sha256_r_order(res: &crate::Residuals) -> String {
    let nr = res.top_features.len();
    let nc = res.ncells;
    // rows of scale_data (top-feature order) sorted by their kept-gene position
    let mut order: Vec<usize> = (0..nr).collect();
    order.sort_by_key(|&r| res.top_features[r]);
    let mut h = scio::sha256::Sha256::new();
    let mut buf: Vec<u8> = Vec::with_capacity(nr * 8);
    for c in 0..nc {
        buf.clear();
        for &r in &order {
            buf.extend_from_slice(&res.scale_data[r * nc + c].to_le_bytes());
        }
        h.update(&buf);
    }
    h.hex()
}
