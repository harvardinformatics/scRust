//! R's own construction of DoubletFinder's SCT input and of the matrix `vst()` works on
//! (validation/R/sctransform_vectors.R, fixtures in tests/fixtures/sctransform/): on synthetic
//! counts with explicit zeros, Rust's input must store exactly R's entries (`i`, `p`, `x`), keep
//! the same genes, and give the same `median(umi@x)`, `row_gmean` and `row_var`, bit for bit.

use doubletfinder::sct::sct_input;
use scio::CscMatrix;
use sctransform::setup::{row_gmean, row_var};
use std::fs;

fn fixture(name: &str) -> Vec<Vec<String>> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures/sctransform/").to_string() + name;
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).map(|l| l.split('\t').map(str::to_string).collect()).collect()
}

fn f(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn fl(s: &str) -> Vec<f64> {
    if s == "-" { Vec::new() } else { s.split(',').map(f).collect() }
}

fn il(s: &str) -> Vec<usize> {
    if s == "-" { Vec::new() } else { s.split(',').map(|v| v.parse().unwrap()).collect() }
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn sct_input_matches_r() {
    let rows = fixture("sct_input.tsv");
    assert_eq!(rows.len(), 5);
    for (case, r) in rows.iter().enumerate() {
        let (nrow, ncol): (usize, usize) = (r[0].parse().unwrap(), r[1].parse().unwrap());
        let m = CscMatrix { nrow, ncol, p: il(&r[3]), i: il(&r[2]).into_iter().map(|v| v as u32).collect(), x: fl(&r[4]) };
        m.validate().unwrap();
        let (a, b) = (il(&r[5]), il(&r[6]));
        let pairs: Vec<(usize, usize)> = a.into_iter().zip(b).collect();
        let real: Vec<usize> = (0..ncol).collect();
        let dw = sct_input(&m, &real, &pairs);
        assert_eq!(dw.i.iter().map(|&v| v as usize).collect::<Vec<_>>(), il(&r[7]), "case {case}: SCT input i");
        assert_eq!(dw.p, il(&r[8]), "case {case}: SCT input p");
        assert_eq!(bits(&dw.x), bits(&fl(&r[9])), "case {case}: SCT input x");
        // vst's gene filter, umi[genes, ], and its statistics
        let min_cells = [3, 5, 2, 5, 1][case];
        let csr = dw.to_csr();
        let genes: Vec<usize> = (0..nrow).filter(|&g| csr.row(g).1.iter().filter(|&&v| v >= 0.01).count() >= min_cells).collect();
        assert_eq!(genes, il(&r[10]), "case {case}: kept rows");
        let umi2 = dw.select_rows(&genes);
        assert_eq!(umi2.i.iter().map(|&v| v as usize).collect::<Vec<_>>(), il(&r[11]), "case {case}: filtered i");
        assert_eq!(umi2.p, il(&r[12]), "case {case}: filtered p");
        assert_eq!(bits(&umi2.x), bits(&fl(&r[13])), "case {case}: filtered x");
        let rows2 = sctransform::GeneRows::from_csc(&dw, &genes);
        assert_eq!(rows2.stored_median().to_bits(), f(&r[14]).to_bits(), "case {case}: median(umi@x)");
        let gm: Vec<f64> = (0..genes.len()).map(|k| row_gmean(rows2.row(k).1, dw.ncol, 1.0)).collect();
        let rv: Vec<f64> = (0..genes.len()).map(|k| row_var(rows2.row(k).1, dw.ncol)).collect();
        assert_eq!(bits(&gm), bits(&fl(&r[15])), "case {case}: row_gmean");
        assert_eq!(bits(&rv), bits(&fl(&r[16])), "case {case}: row_var");
    }
}
