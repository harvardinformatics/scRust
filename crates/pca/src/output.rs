//! `scrust pca`'s output layout (docs/PLAN.md C3), written by the command and by the clustering
//! chain (`scrust cluster`'s `pca/`, M11), so that both give the same bytes:
//!
//! - `embeddings.tsv`: barcode, PC_1 .. PC_<npcs>: the cell embeddings;
//! - `loadings.tsv`: id, name, PC_1 .. PC_<npcs>: the feature loadings (features kept);
//! - `stdev.tsv`: pc, stdev;
//! - `random_seed.txt`: R's `.Random.seed` after RunPCA, one integer per line.
//!
//! Doubles are written with 17 significant digits, so they read back exactly.

use crate::{Pca, PcaParams};
use scio::{fmt17, Json};
use std::fmt::Write as _;
use std::io;
use std::path::Path;

/// The files of the layout, in the order written.
pub const FILES: [&str; 4] = ["embeddings.tsv", "loadings.tsv", "stdev.tsv", "random_seed.txt"];

fn write(path: &Path, text: &str) -> io::Result<()> {
    std::fs::write(path, text).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
}

/// Writes the layout into `dir` (which must exist). `barcodes` are the cells; `features` the PCA
/// input's rows (id, name), of which `p.features` are the kept ones; `rng` R's `.Random.seed`
/// after RunPCA.
pub fn write_outputs(dir: &Path, barcodes: &[String], features: &[(String, String)], p: &Pca, rng: &[i32]) -> io::Result<()> {
    let nc = p.ncells;
    assert_eq!(barcodes.len(), nc);
    let pcs: Vec<String> = (1..=p.npcs).map(|k| format!("PC_{k}")).collect();
    let mut t = format!("barcode\t{}\n", pcs.join("\t"));
    for (c, b) in barcodes.iter().enumerate() {
        t.push_str(b);
        for k in 0..p.npcs {
            t.push('\t');
            t.push_str(&fmt17(p.embeddings[c + k * nc]));
        }
        t.push('\n');
    }
    write(&dir.join("embeddings.tsv"), &t)?;
    let nk = p.features.len();
    let mut t = format!("id\tname\t{}\n", pcs.join("\t"));
    for (i, &f) in p.features.iter().enumerate() {
        let _ = write!(t, "{}\t{}", features[f].0, features[f].1);
        for k in 0..p.npcs {
            t.push('\t');
            t.push_str(&fmt17(p.loadings[i + k * nk]));
        }
        t.push('\n');
    }
    write(&dir.join("loadings.tsv"), &t)?;
    let mut t = String::from("pc\tstdev\n");
    for (k, s) in p.sdev.iter().enumerate() {
        let _ = writeln!(t, "{}\t{}", pcs[k], fmt17(*s));
    }
    write(&dir.join("stdev.tsv"), &t)?;
    let seed: String = rng.iter().map(|v| format!("{v}\n")).collect();
    write(&dir.join("random_seed.txt"), &seed)
}

/// The parameters as `summary.json` records them.
pub fn params_json(params: &PcaParams) -> Json {
    Json::obj(vec![
        ("npcs", Json::Int(params.npcs as i64)),
        ("seed_use", Json::Int(params.seed_use.unwrap_or(0) as i64)),
        ("weight_by_var", Json::Bool(params.weight_by_var)),
        ("irlba_maxit", Json::Int(params.irlba.maxit as i64)),
        ("irlba_tol", Json::Num(params.irlba.tol)),
    ])
}

/// The `result` section of `summary.json` (deterministic). `n_features` is the input's rows.
pub fn result_json(p: &Pca, n_features: usize) -> Json {
    let ir = &p.irlba;
    Json::obj(vec![
        ("n_cells", Json::Int(p.ncells as i64)),
        ("n_features", Json::Int(n_features as i64)),
        ("n_pca_features", Json::Int(p.features.len() as i64)),
        ("npcs", Json::Int(p.npcs as i64)),
        ("irlba_d", Json::nums(&ir.d)),
        ("irlba_iter", Json::Int(ir.iter as i64)),
        ("irlba_mprod", Json::Int(ir.mprod as i64)),
        ("irlba_converged", Json::Bool(ir.converged)),
        ("irlba_work", Json::Int(ir.work as i64)),
        ("irlba_invariant_draws", Json::Int(ir.invariant_draws as i64)),
        ("total_variance", Json::Num(p.total_variance)),
    ])
}

/// The manifest's entries for this layout (docs/PLAN.md C3f, M13), in the directory `dir`
/// relative to the output directory (`""` or `"pca"`); their shapes are read from the files.
pub fn manifest_entries(dir: &str) -> Vec<scio::manifest::Entry> {
    use scio::manifest::{Entry, Ids};
    vec![
        Entry::new("embeddings.tsv", "pca_embeddings", "pca_embeddings", "tsv", "RunPCA's cell embeddings (u %*% diag(d)): barcode, then PC_1, PC_2, ...; Seurat's pca reduction")
            .header(true)
            .rows(Ids::named("cells", "embeddings.tsv", "barcode")),
        Entry::new("loadings.tsv", "pca_loadings", "pca_loadings", "tsv", "RunPCA's feature loadings: id, name, then PC_1, PC_2, ...; the rows are the PCA's input features (the variable features, in VariableFeatures() order)")
            .header(true)
            .rows(Ids::named("variable_features", "loadings.tsv", "name")),
        Entry::new("stdev.tsv", "pca_stdev", "pca_stdev", "tsv", "RunPCA's standard deviations (d / sqrt(cells - 1)): pc, stdev").header(true).rows(Ids::named("components", "stdev.tsv", "pc")),
        Entry::new("random_seed.txt", "rng_state", "none", "text", "R's .Random.seed after RunPCA, one integer per line"),
    ]
    .into_iter()
    .map(|e| e.under(dir))
    .collect()
}
