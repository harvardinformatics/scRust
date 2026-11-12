//! sctransform's (and glmGamPoi's) functions against the installed R packages, on the synthetic
//! inputs of validation/R/sctransform_vectors.R (tests/fixtures/sctransform/).

use sctransform::glm::fit_beta_one_group;
use std::fs;
use std::path::PathBuf;

fn fixture(name: &str) -> Vec<Vec<String>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/sctransform").join(name);
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').map(str::to_owned).collect()).collect()
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

fn assert_same(got: f64, want: f64, what: &str) {
    assert!(got.to_bits() == want.to_bits(), "{what}: got {got:e} ({:016x}), R gives {want:e} ({:016x})", got.to_bits(), want.to_bits());
}

#[test]
fn fit_beta_one_group_against_glmgampoi() {
    let rows = fixture("fit_beta.tsv");
    assert_eq!(rows.len(), 24);
    let mut fallbacks = 0;
    for (k, row) in rows.iter().enumerate() {
        let y: Vec<f64> = row[0].split(',').map(|v| v.parse::<f64>().unwrap() / 2.0).collect();
        let off = fl(&row[1]);
        let (theta, start) = (f(&row[2]), f(&row[3]));
        fallbacks += (row[4] == "1") as usize;
        let what = format!("case {} (n = {}, theta {theta:e}, Newton falls back: {})", k + 1, y.len(), row[4]);
        let fit = fit_beta_one_group(&y, &off, theta, start).unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_same(fit.beta, f(&row[5]), &format!("{what}: beta"));
        assert_eq!(fit.iterations, row[6].parse::<i32>().unwrap(), "{what}: iterations");
        assert_same(fit.deviance, f(&row[7]), &format!("{what}: deviance"));
    }
    assert_eq!(fallbacks, 14);
}

#[test]
fn is_outlier_against_sctransform() {
    use sctransform::reg::is_outlier;
    let rows = fixture("is_outlier.tsv");
    assert_eq!(rows.len(), 10);
    let mut flagged = 0;
    for (k, row) in rows.iter().enumerate() {
        let (x, y) = (fl(&row[0]), fl(&row[1]));
        let want: Vec<bool> = row[2].split(',').map(|v| v == "1").collect();
        let got = is_outlier(&y, &x).unwrap_or_else(|e| panic!("case {}: {e}", k + 1));
        assert_eq!(got, want, "is_outlier case {} (n = {})", k + 1, x.len());
        flagged += want.iter().filter(|&&b| b).count();
    }
    assert!(flagged > 50);
}
