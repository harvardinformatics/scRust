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
    let mut v = vec![
        ("seed_use", Json::Int(params.seed_use.unwrap_or(0) as i64)),
        ("ncells", Json::Int(params.ncells as i64)),
        ("variable_features_n", Json::Int(params.variable_features_n.unwrap_or(0) as i64)),
        ("n_genes", Json::Int(params.n_genes.map_or(0, |v| v as i64))),
        ("bin_size", Json::Int(params.bin_size as i64)),
        ("min_cells", Json::Int(params.min_cells as i64)),
        ("bw_adjust", Json::Num(params.bw_adjust)),
    ];
    // M15b: only when set, so that a run without it records what earlier runs recorded
    if !params.vars_to_regress.is_empty() {
        v.push(("vars_to_regress", Json::strs(&params.vars_to_regress)));
    }
    Json::obj(v)
}

/// `result`'s record of `vars.to.regress`' regression (M15b): the covariates in the formula's
/// order, and the QR `lm()` made of the model matrix (the intercept and the covariates): its rank,
/// the columns' order after `dqrdc2`'s pivoting (1-based) and its tolerance.
pub fn regression_json(r: &crate::Regression) -> Json {
    Json::obj(vec![
        ("vars_to_regress", Json::strs(&r.names)),
        ("model", Json::Str("linear: lm(GENE ~ vars)$qr once, qr.resid() per row (Seurat's RegressOutMatrix)".into())),
        ("n_columns", Json::Int(r.qr.p as i64)),
        ("rank", Json::Int(r.qr.rank as i64)),
        ("pivot", Json::Vector(r.qr.pivot.iter().map(|&v| Json::Int(v as i64)).collect())),
        ("tol", Json::Num(r.qr.tol)),
    ])
}

/// The `result` section of `summary.json` (deterministic: equal between runs that agree).
pub fn result_json(sct: &Sct, n_features: usize) -> Json {
    let s = &sct.setup;
    let res = &sct.residuals;
    let n_out = sct.reg.outliers.iter().filter(|&&o| o).count();
    let mut out = Json::obj(vec![
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
    ]);
    // the corrected counts (M13), when computed: the members R's reference runs record
    if let (Json::Object(m), Some(c)) = (&mut out, &sct.corrected) {
        m.push(("corrected".into(), corrected_json(c)));
    }
    // vars.to.regress (M15b), when regressed
    if let (Json::Object(m), Some(r)) = (&mut out, &sct.regression) {
        m.push(("regression".into(), regression_json(r)));
    }
    out
}

/// `result`'s summary of the corrected counts, with the digests R's reference runs record
/// (`sct_counts_x_sha256`, `sct_data_x_sha256`: the stored values as doubles, in stored order).
pub fn corrected_json(c: &crate::Corrected) -> Json {
    let nc = c.n_count();
    let total: f64 = nc.iter().sum();
    Json::obj(vec![
        ("n_features", Json::Int(c.ngenes as i64)),
        ("n_cells", Json::Int(c.ncells as i64)),
        ("nnz", Json::Int(c.nnz() as i64)),
        ("min", c.min().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
        ("explicit_zeros", Json::Int(0)),
        ("median_log_umi", Json::Num(c.median_log_umi)),
        ("residuals_beyond_clip", Json::Int(c.beyond_clip as i64)),
        ("total", Json::Num(total)),
        ("counts_x_sha256", Json::Str(c.x_sha256(false))),
        ("data_x_sha256", Json::Str(c.x_sha256(true))),
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

/// The Matrix Market comment lines of the SCT assay's two layers (docs/PLAN.md C3f): what the
/// matrix is, that it is not raw counts, how it was made (scrust's version and command) and the
/// input's checksums (`inputs`: SHA-256 and name of each input file, as `sha256sum` prints
/// them). Each line is written after `%`. Returns (counts layer, data layer).
pub fn corrected_comments(command: &str, inputs: &[(String, String)], c: &crate::Corrected) -> (Vec<String>, Vec<String>) {
    let version = env!("CARGO_PKG_VERSION");
    let tail = |v: &mut Vec<String>| {
        v.push(format!(" They are sctransform 0.4.3's correct() as Seurat 5.5.1's SCTransform() makes them: every cell's counts corrected to the median sequencing depth (log10 UMIs {}) and rounded.", scio::fmt17(c.median_log_umi)));
        v.push(" Rows: the model genes, sct_features.tsv.gz; columns: the cells, sct_cells.tsv.gz. The raw counts are elsewhere (an export's counts/).".into());
        v.push(format!(" Made by scrust {version}: {command}"));
        for (sha, name) in inputs {
            v.push(format!(" Input sha256: {sha}  {name}"));
        }
        v.push(" Described in summary.json (files) and README.txt.".into());
    };
    let mut counts = vec![" scrust: SCTransform's CORRECTED counts (Seurat's SCT assay, counts layer). NOT raw counts: do not use them as raw UMI counts.".to_string()];
    tail(&mut counts);
    let mut data = vec![" scrust: log1p of SCTransform's corrected counts (Seurat's SCT assay, data layer). NOT raw counts.".to_string()];
    tail(&mut data);
    (counts, data)
}

/// The manifest's entries for `scrust sctransform`'s layout (as [`write_outputs`] writes it), in
/// the directory `dir` relative to the output directory (`""` or `"sct"`): `ncells` cells,
/// `ngenes` kept genes, `n_step1` fitting genes, `nvf` variable features. With `corrected`, also
/// the SCT assay's `counts` and `data` layers and their rows and columns ([`crate::Corrected::write`]).
pub fn manifest_entries(dir: &str, ncells: usize, ngenes: usize, n_step1: usize, nvf: usize, corrected: Option<&crate::Corrected>) -> Vec<scio::manifest::Entry> {
    use scio::manifest::{Entry, Ids};
    let cells = Ids::named("cells", "cell_attr.tsv", "barcode");
    let mut v = vec![
        Entry::new("cell_attr.tsv", "sct_cell_attributes", "none", "tsv", "vst's cell attributes: barcode, umi (the cell's total count), log_umi (its log10); Seurat's SCTModel cell.attributes")
            .header(true)
            .shape(&[ncells, 3])
            .rows(cells.clone()),
        Entry::new("gene_attr.tsv", "sct_feature_attributes", "none", "tsv", "vst's attributes of the kept genes (detected in at least 5 cells), in the counts' order: id, name, detection_rate, gmean, amean, variance, residual_mean, residual_variance; Seurat's SCTModel feature.attributes")
            .header(true)
            .shape(&[ngenes, 8])
            .rows(Ids::named("features", "gene_attr.tsv", "name")),
        Entry::new("model_pars.tsv", "sct_model_pars", "none", "tsv", "the fitting genes' parameters (theta after the 1e-7 floor, (Intercept), log_umi), in sampling order")
            .header(true)
            .shape(&[n_step1, 5])
            .rows(Ids::named("features", "model_pars.tsv", "name")),
        Entry::new("model_pars_fit.tsv", "sct_model_pars_fit", "none", "tsv", "the regularised parameters (theta, (Intercept), log_umi) of every kept gene, in the counts' order; Seurat's SCTModel model.pars")
            .header(true)
            .shape(&[ngenes, 5])
            .rows(Ids::named("features", "model_pars_fit.tsv", "name")),
        Entry::new("variable_features.tsv", "sct_variable_features", "none", "tsv", "the variable features (id, name), highest residual variance first: Seurat's VariableFeatures(), the rows of scale_data.npy")
            .header(true)
            .shape(&[nvf, 2])
            .rows(Ids::named("variable_features", "variable_features.tsv", "name")),
        Entry::new("scale_data.npy", "sct_scale_data", "residuals", "npy", "Seurat's SCT scale.data: the variable features' Pearson residuals, clipped to +-sqrt(cells/30) and centred; float64, variable features x cells, column-major. NOT counts")
            .shape(&[nvf, ncells])
            .rows(Ids::named("variable_features", "variable_features.tsv", "name"))
            .columns(cells),
    ];
    if let Some(c) = corrected {
        let genes = Ids::numbered("features", "sct_features.tsv.gz", 1);
        let cols = Ids::numbered("cells", "sct_cells.tsv.gz", 1);
        v.push(
            Entry::new("corrected_counts.mtx.gz", "sct_counts", "corrected_counts", "matrix_market", "Seurat's SCT counts layer: sctransform's correct() counts, corrected to the median depth and rounded (integers); model genes x cells. NOT raw counts")
                .shape(&[c.ngenes, c.ncells])
                .nnz(c.nnz())
                .rows(genes.clone())
                .columns(cols.clone()),
        );
        v.push(
            Entry::new("log1p_corrected_counts.mtx.gz", "sct_data", "log1p_corrected_counts", "matrix_market", "Seurat's SCT data layer: log1p of the corrected counts; model genes x cells, the same entries. NOT raw counts")
                .shape(&[c.ngenes, c.ncells])
                .nnz(c.nnz())
                .rows(genes)
                .columns(cols),
        );
        v.push(Entry::new("sct_features.tsv.gz", "sct_features", "none", "tsv", "the model genes, one name per line, in the counts' order: the rows of the SCT counts and data layers").header(false).shape(&[c.ngenes, 1]));
        v.push(Entry::new("sct_cells.tsv.gz", "sct_cells", "none", "tsv", "the cells, one barcode per line, in the counts' order: the columns of the SCT counts and data layers").header(false).shape(&[c.ncells, 1]));
    }
    v.into_iter().map(|e| e.under(dir)).collect()
}
