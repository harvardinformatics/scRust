//! Milestone M15b: SCTransform's `vars.to.regress` against R 4.5.3 and Seurat 5.5.1
//! (validation/R/regress_vectors.R; fixtures tests/fixtures/rcompat/linpack_qr.tsv.gz and
//! linpack_resid.tsv.gz, made under OpenBLAS's Haswell kernel):
//!
//! - `lm(GENE ~ <vars>, qr = TRUE)$qr` (LINPACK `dqrdc2` through `Cdqrls`, `tol = 1e-07`), as
//!   Seurat's `RegressOutMatrix()` takes it once: the compact QR, `qraux`, the pivot and the rank,
//!   bit for bit, on full-rank models, a zero column, a constant column, two collinear columns, a
//!   near-constant one (rank 1), square and wide models and an outlier;
//! - `qr.resid(qr, y)` (`dqrrsd` -> `dqrsl`), which `RegressOutMatrix()` returns row by row (the
//!   R script checked that its rows equal `qr.resid()`'s), bit for bit;
//! - [`sctransform::Regression`], which the chain uses, on the same rows, at several thread counts.
//!
//! `dqrdc2` and `dqrsl` call OpenBLAS's `dnrm2`, `dscal`, `ddot` and `daxpy`, whose results depend
//! on the kernel: the tests run under the fixtures' kernel (`rcompat::blas::require_fixture_kernel`).

use rcompat::linpack::{lm_qr, QrResid};
use sctransform::Regression;
use std::collections::HashMap;
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/rcompat").join(name)
}

fn fixture(name: &str) -> Vec<Vec<String>> {
    let path = fixture_path(name);
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let text = String::from_utf8(scio::gz::gunzip(&bytes).expect("gzip")).expect("utf-8");
    text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).map(|l| l.split('\t').map(str::to_owned).collect()).collect()
}

fn f(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn fl(s: &str) -> Vec<f64> {
    if s == "-" {
        return Vec::new();
    }
    s.split(',').map(f).collect()
}

fn same(a: &[f64], b: &[f64]) -> Option<usize> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).position(|(x, y)| x.to_bits() != y.to_bits())
}

struct Case {
    n: usize,
    p: usize,
    vars: Vec<String>,
    x: Vec<f64>,
    qr: Vec<f64>,
    qraux: Vec<f64>,
    pivot: Vec<usize>,
    rank: usize,
}

fn cases() -> Vec<(String, Case)> {
    fixture("linpack_qr.tsv.gz")
        .into_iter()
        .map(|r| {
            let c = Case {
                n: r[1].parse().unwrap(),
                p: r[2].parse().unwrap(),
                vars: r[3].split(',').map(str::to_owned).collect(),
                x: fl(&r[4]),
                qr: fl(&r[5]),
                qraux: fl(&r[6]),
                pivot: r[7].split(',').map(|v| v.parse().unwrap()).collect(),
                rank: r[8].parse().unwrap(),
            };
            (r[0].clone(), c)
        })
        .collect()
}

fn kernel() {
    rcompat::blas::require_fixture_kernel(&fixture_path("linpack_qr.tsv.gz"), "linpack_qr.tsv.gz");
}

#[test]
fn lm_qr_as_r_makes_it() {
    kernel();
    let cs = cases();
    assert_eq!(cs.len(), 14);
    let mut deficient = 0;
    for (name, c) in &cs {
        let qr = lm_qr(c.x.clone(), c.n, c.p).unwrap();
        assert_eq!((qr.rank, &qr.pivot), (c.rank, &c.pivot), "{name}: rank and pivot");
        assert_eq!(same(&qr.qraux, &c.qraux), None, "{name}: qraux");
        assert_eq!(same(&qr.qr, &c.qr), None, "{name}: qr");
        if qr.rank < qr.p {
            deficient += 1;
        }
    }
    // the zero, constant, collinear, near-constant and wide models
    assert_eq!(deficient, 5);
}

#[test]
fn qr_resid_and_regress_out_matrix_as_r() {
    kernel();
    let cs: HashMap<String, Case> = cases().into_iter().collect();
    let rows = fixture("linpack_resid.tsv.gz");
    assert_eq!(rows.len(), 3 * cs.len());
    let mut by_case: HashMap<&str, Vec<(Vec<f64>, Vec<f64>)>> = HashMap::new();
    for r in &rows {
        let c = &cs[&r[0]];
        let qr = lm_qr(c.x.clone(), c.n, c.p).unwrap();
        let res = QrResid::new(&qr);
        let (y, want) = (fl(&r[2]), fl(&r[3]));
        let (mut qty, mut rsd) = (vec![0.0; c.n], vec![0.0; c.n]);
        res.resid(&y, &mut qty, &mut rsd);
        assert_eq!(same(&rsd, &want), None, "{} row {}", r[0], r[1]);
        by_case.entry(&r[0]).or_default().push((y, want));
    }
    // the chain's Regression on the same rows: the covariates are x's columns after the intercept
    for (name, rows) in by_case {
        let c = &cs[name];
        let covs: Vec<Vec<f64>> = (1..c.p).map(|j| c.x[j * c.n..(j + 1) * c.n].to_vec()).collect();
        let reg = Regression::new(c.vars.clone(), &covs, c.n).unwrap();
        let ys: Vec<f64> = rows.iter().flat_map(|(y, _)| y.clone()).collect();
        let want: Vec<f64> = rows.iter().flat_map(|(_, w)| w.clone()).collect();
        for threads in [1, 2, 3] {
            let mut got = ys.clone();
            reg.regress_rows(&mut got, c.n, threads);
            assert_eq!(same(&got, &want), None, "{name}: RegressOutMatrix, {threads} threads");
        }
    }
}
