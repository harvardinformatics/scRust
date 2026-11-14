//! Stage K0 (docs/PLAN.md C5): `Kiss64Random` against vectors printed by RcppAnnoy 0.0.23's own
//! `kissrandom.h`, compiled with the production R env's compiler and flags (milestone M8,
//! `tests/fixtures/clustering/k0/`): `kiss()`, `index(n)` for n from 1 to 2^32 + 1, `flip()`, and
//! one stream through the mixed calls an index build makes, from the default seed every Seurat
//! build uses and from others.

use neighbors::Kiss64Random;
use std::path::PathBuf;

fn k0(name: &str) -> Vec<Vec<String>> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/k0").join(name);
    let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).map(|l| l.split('\t').map(String::from).collect()).collect()
}

fn list(s: &str) -> Vec<u64> {
    s.split(',').map(|v| v.parse().unwrap()).collect()
}

#[test]
fn kiss_values() {
    let rows = k0("kiss64_kiss.tsv");
    assert!(rows.len() >= 4);
    for r in &rows {
        let mut g = Kiss64Random::new(r[0].parse().unwrap());
        let want = list(&r[2]);
        assert_eq!(want.len(), r[1].parse::<usize>().unwrap());
        let got: Vec<u64> = (0..want.len()).map(|_| g.kiss()).collect();
        assert_eq!(got, want, "seed {}", r[0]);
    }
}

#[test]
fn index_values() {
    let rows = k0("kiss64_index.tsv");
    assert!(rows.len() >= 40);
    for r in &rows {
        let mut g = Kiss64Random::new(r[0].parse().unwrap());
        let n: u64 = r[1].parse().unwrap();
        let want = list(&r[3]);
        let got: Vec<u64> = (0..want.len()).map(|_| g.index(n as usize) as u64).collect();
        assert_eq!(got, want, "seed {}, n {n}", r[0]);
    }
}

#[test]
fn flip_values() {
    let rows = k0("kiss64_flip.tsv");
    for r in &rows {
        let mut g = Kiss64Random::new(r[0].parse().unwrap());
        let want = list(&r[2]);
        let got: Vec<u64> = (0..want.len()).map(|_| g.flip() as u64).collect();
        assert_eq!(got, want, "seed {}", r[0]);
    }
}

#[test]
fn mixed_stream() {
    let rows = k0("kiss64_mixed.tsv");
    assert_eq!(rows.len(), 2);
    for r in &rows {
        let mut g = Kiss64Random::new(r[0].parse().unwrap());
        let calls: Vec<&str> = r[1].split(',').collect();
        let want = list(&r[2]);
        assert_eq!(calls.len(), want.len());
        assert!(calls.contains(&"f"), "the stream mixes flips in");
        for (t, (c, w)) in calls.iter().zip(&want).enumerate() {
            let got = match c.strip_prefix('i') {
                Some(n) => g.index(n.parse().unwrap()) as u64,
                None => {
                    assert_eq!(*c, "f");
                    g.flip() as u64
                }
            };
            assert_eq!(got, *w, "seed {}, call {} ({c})", r[0], t + 1);
        }
    }
}
