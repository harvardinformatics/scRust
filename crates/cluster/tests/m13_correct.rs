//! Milestone M13, K10 teacher-forced (docs/PLAN.md C3f, C5): sctransform's `correct()` from R's
//! own model parameters (`stages/vst/model_pars_fit.tsv`), cell attributes (`cell_attr.tsv`),
//! genes (`K10_genes`) and the counts, on the two base objects, against M8's portable
//! instrumented runs: the corrected counts and their log1p (by the SHA-256 of R's files' data
//! lines and R's digests of the values, `reference/k10.tsv`), `nCount_SCT` and `nFeature_SCT`
//! (R's `metadata.tsv`), correct()'s median log_umi and its residuals beyond the clip range
//! (R's `stages.json`). Its `mu` is an OpenBLAS product, so the test first checks that the kernel
//! in use is the references' (Haswell). The result does not depend on the thread count.

use scio::{Json, Npy, Table};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering")
}

fn gz(p: &Path) -> Vec<u8> {
    let path = if p.exists() { p.to_path_buf() } else { PathBuf::from(format!("{}.gz", p.display())) };
    scio::tenx::read_maybe_gz(&path).unwrap()
}

fn table(p: &Path) -> Table {
    Table::parse(&String::from_utf8(gz(p)).unwrap()).unwrap()
}

/// The `%.17g` text of the corrected counts' data lines, as R writes them, digested (the layout of
/// `reference/k10.tsv`): the banner, the size line and the entries.
fn data_lines_sha256(c: &sctransform::Corrected, log1p: bool) -> String {
    let l = sctransform::correct::Log1p::new();
    let mut h = scio::sha256::Sha256::new();
    h.update(format!("%%MatrixMarket matrix coordinate {} general\n{} {} {}\n", if log1p { "real" } else { "integer" }, c.ngenes, c.ncells, c.nnz()).as_bytes());
    for j in 0..c.ncells {
        for k in c.p[j]..c.p[j + 1] {
            let v = if log1p { scio::fmt17(l.of(c.x[k])) } else { c.x[k].to_string() };
            h.update(format!("{} {} {v}\n", c.i[k] + 1, j + 1).as_bytes());
        }
    }
    h.hex()
}

#[test]
fn correct_on_r_s_model_gives_r_s_corrected_counts() {
    for object in ["filtered_seurat_tenx_test", "filtered_seurat_emptydrops_test"] {
        let reference = fixtures().join(object).join("reference");
        rcompat::blas::require_fixture_kernel(&reference, "correct() teacher-forced (K10)");
        let input = fixtures().join(object).join("input");
        let counts = scio::sparse::read_mtx(&input).unwrap();
        let features = scio::tenx::read_features(&input).unwrap();
        let gi = Npy::parse(&gz(&reference.join("stages/K10_genes.npy"))).unwrap();
        let genes: Vec<usize> = (0..gi.nrow()).map(|i| gi.get_i32(i, 0) as usize).collect();
        let mp = table(&reference.join("stages/vst/model_pars_fit.tsv"));
        let col = |t: &Table, n: &str| t.col(n).unwrap();
        let (nc, tc, ic, lc) = (col(&mp, "name"), col(&mp, "theta"), col(&mp, "(Intercept)"), col(&mp, "log_umi"));
        let pars: Vec<[f64; 3]> = mp.rows.iter().map(|r| [scio::parse17(&r[tc]).unwrap(), scio::parse17(&r[ic]).unwrap(), scio::parse17(&r[lc]).unwrap()]).collect();
        assert!(genes.iter().zip(&mp.rows).all(|(&g, r)| features[g].name == r[nc]), "{object}: K10_genes are model_pars_fit's genes");
        let ca = table(&reference.join("stages/vst/cell_attr.tsv"));
        let log_umi: Vec<f64> = ca.rows.iter().map(|r| scio::parse17(&r[col(&ca, "log_umi")]).unwrap()).collect();
        let rows = sctransform::GeneRows::from_csc(&counts, &genes);
        let input = sctransform::ResidualInput { model_pars_fit: &pars, log_umi: &log_umi };
        let params = sctransform::SctParams::default();
        let want: BTreeMap<String, String> = table(&reference.join("k10.tsv")).rows.iter().map(|r| (r[0].clone(), r[1].clone())).collect();
        let mut first: Option<sctransform::Corrected> = None;
        for threads in [1, 3] {
            let (_, c) = sctransform::residuals_corrected(&rows, &input, &params, false, threads).unwrap();
            assert_eq!(data_lines_sha256(&c, false), want["corrected_counts_data_lines_sha256"], "{object}: the corrected counts");
            assert_eq!(data_lines_sha256(&c, true), want["log1p_data_lines_sha256"], "{object}: their log1p");
            assert_eq!(c.x_sha256(false), want["sct_counts_x_sha256"]);
            assert_eq!(c.x_sha256(true), want["sct_data_x_sha256"]);
            assert_eq!(c.median_log_umi.to_bits(), want["correct_median_log_umi"].parse::<f64>().unwrap().to_bits());
            assert_eq!(c.beyond_clip.to_string(), want["correct_residuals_beyond_clip"]);
            let meta = table(&reference.join("metadata.tsv"));
            let ncount: Vec<String> = c.n_count().iter().map(|&v| scio::fmt17(v)).collect();
            let nfeat: Vec<String> = c.n_feature().iter().map(|v| v.to_string()).collect();
            assert_eq!(ncount, meta.rows.iter().map(|r| r[col(&meta, "nCount_SCT")].clone()).collect::<Vec<_>>(), "{object}: nCount_SCT");
            assert_eq!(nfeat, meta.rows.iter().map(|r| r[col(&meta, "nFeature_SCT")].clone()).collect::<Vec<_>>(), "{object}: nFeature_SCT");
            let sc = Json::parse(&String::from_utf8(gz(&reference.join("stages/stages.json"))).unwrap()).unwrap();
            let sc = sc.get("scalars").and_then(|s| s.get("chain")).unwrap();
            assert_eq!(sc.get("vst_umi_corrected_nnz").and_then(Json::as_f64), Some(c.nnz() as f64));
            if let Some(f) = &first {
                assert_eq!(f, &c, "{object}: the thread count changed the result");
            }
            first = Some(c);
        }
    }
}
