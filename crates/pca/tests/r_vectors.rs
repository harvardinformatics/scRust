//! irlba 2.4.1 and Seurat's RowVar against the installed R packages, on the synthetic inputs of
//! validation/R/pca_vectors.R (tests/fixtures/pca/): the branches the reference data never reach.
//! irlba's products are OpenBLAS calls and the vectors were made under the portable profile's
//! kernel, so the irlba test needs `OPENBLAS_CORETYPE=Haswell` (env.sh's default).

use pca::{irlba, IrlbaParams, Products};
use rcompat::RRng;
use std::fs;
use std::path::PathBuf;

fn fixture(name: &str) -> Vec<Vec<String>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/pca").join(name);
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

fn same(got: &[f64], want: &[f64], what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: length");
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(g.to_bits() == w.to_bits(), "{what}[{k}]: got {g:e}, R gives {w:e}");
    }
}

/// The input matrix of a recipe, as pca_vectors.R's make_input() builds it.
fn make_input(recipe: &str, m: usize, n: usize, data_seed: i32) -> Vec<f64> {
    let mut rng = RRng::new(data_seed);
    let parts: Vec<&str> = recipe.split(':').collect();
    match parts[0] {
        "rnorm" => rcompat::rnorm(&mut rng, m * n, 0.0, 1.0),
        "lowrank" => {
            let r: usize = parts[1].parse().unwrap();
            let x = rcompat::rnorm(&mut rng, m * r, 0.0, 1.0);
            let y = rcompat::rnorm(&mut rng, r * n, 0.0, 1.0);
            rcompat::blas::matprod(&x, m, r, &y, r, n)
        }
        "diag" => {
            let (a, b): (usize, usize) = (parts[1].parse().unwrap(), parts[2].parse().unwrap());
            let c: f64 = parts[3].parse().unwrap();
            let mut d = vec![0.0; m * n];
            for i in 0..a {
                d[i + i * m] = 2.0;
            }
            for i in a..a + b {
                d[i + i * m] = 1.0;
            }
            let noise = rcompat::rnorm(&mut rng, m * n, 0.0, 1.0);
            if c == 0.0 {
                d
            } else {
                d.iter().zip(&noise).map(|(d, z)| d + c * z).collect()
            }
        }
        "cluster" => {
            // the first k singular values 1, 1 - gap, 1 - 2 gap, ...; then `rest` values
            // 0.5 (1 - j / rest)
            let k: usize = parts[1].parse().unwrap();
            let gap: f64 = parts[2].parse().unwrap();
            let rest: usize = parts[3].parse().unwrap();
            let mut vals: Vec<f64> = (0..k).map(|i| 1.0 - i as f64 * gap).collect();
            vals.extend((0..rest).map(|j| 0.5 * (1.0 - j as f64 / rest as f64)));
            let mut d = vec![0.0; m * n];
            for (i, v) in vals.iter().enumerate() {
                d[i + i * m] = *v;
            }
            d
        }
        other => panic!("unknown recipe {other}"),
    }
}

#[test]
fn irlba_against_r() {
    if std::env::var("OPENBLAS_CORETYPE").ok().as_deref() != Some("Haswell") {
        eprintln!("skipped: needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
        return;
    }
    let rows = fixture("irlba.tsv");
    assert_eq!(rows.len(), 19);
    let (mut tiny, mut invariant, mut expanded, mut unconverged, mut nv1) = (0, 0, 0, 0, 0);
    for row in &rows {
        let case = &row[0];
        let (m, n): (usize, usize) = (row[2].parse().unwrap(), row[3].parse().unwrap());
        let a = make_input(&row[1], m, n, row[4].parse().unwrap());
        let nv: usize = row[5].parse().unwrap();
        let params = IrlbaParams { work: (row[6] != "-").then(|| row[6].parse().unwrap()), maxit: row[7].parse().unwrap(), ..IrlbaParams::default() };
        let rng = RRng::new(row[8].parse().unwrap());
        let start_work = params.work.unwrap_or(nv + 16);
        for threads in [1, 3] {
            let p = Products::new(&a, m, n, threads);
            let r = irlba(&p, nv, &params, &mut rng.clone()).unwrap_or_else(|e| panic!("{case}: {e}"));
            let mut g = rng.clone();
            let r2 = irlba(&p, nv, &params, &mut g).unwrap();
            assert_eq!(r.iter, row[9].parse::<usize>().unwrap(), "{case}: iter");
            assert_eq!(r.mprod, row[10].parse::<usize>().unwrap(), "{case}: mprod");
            assert_eq!(r.converged, row[11] == "1", "{case}: converged");
            same(&r.d, &fl(&row[12]), &format!("{case}: d"));
            same(&r.u, &fl(&row[13]), &format!("{case}: u"));
            same(&r.v, &fl(&row[14]), &format!("{case}: v"));
            let seed: Vec<i32> = row[15].split(',').map(|v| v.parse().unwrap()).collect();
            assert_eq!(g.random_seed(), seed, "{case}: .Random.seed after irlba");
            same(&r2.d, &r.d, &format!("{case}: repeat"));
            if threads == 1 {
                eprintln!(
                    "{case}: {m} x {n}, nv {nv}: iter {}, mprod {}, converged {}, tiny {}, invariant draws {}, work {} -> {}",
                    r.iter, r.mprod, r.converged, r.tiny, r.invariant_draws, start_work, r.work
                );
                tiny += r.tiny as usize;
                invariant += (r.invariant_draws > 0) as usize;
                expanded += (r.work > start_work.min(m.min(n) - 2).max(nv + 1)) as usize;
                unconverged += (!r.converged) as usize;
                nv1 += (nv == 1) as usize;
            }
        }
    }
    // the branches the cases were chosen for: the tiny-problem svd, the near-invariant-subspace
    // draws, the clustered-spectrum expansion of the working dimension, maxit reached, nv = 1
    assert!(tiny >= 2 && invariant >= 3 && expanded >= 4 && unconverged >= 2 && nv1 >= 2, "tiny {tiny}, invariant {invariant}, expanded {expanded}, unconverged {unconverged}, nv1 {nv1}");
}

#[test]
fn rowvar_against_seurat() {
    let rows = fixture("rowvar.tsv");
    assert_eq!(rows.len(), 17);
    for row in &rows {
        let (nr, nc): (usize, usize) = (row[1].parse().unwrap(), row[2].parse().unwrap());
        let parts: Vec<&str> = row[0].split(':').collect();
        let x: Vec<f64> = match parts[0] {
            "rnorm" => {
                let mut rng = RRng::new(parts[1].parse().unwrap());
                let (scale, offset) = (f(parts[2]), f(parts[3]));
                // offset + scale * rnorm(), filled by row
                rcompat::rnorm(&mut rng, nr * nc, 0.0, 1.0).iter().map(|z| offset + scale * z).collect()
            }
            "const" => vec![f(parts[1]); nr * nc],
            other => panic!("unknown recipe {other}"),
        };
        let v: Vec<f64> = (0..nr).map(|r| pca::rowvar::row_var(&x[r * nc..(r + 1) * nc])).collect();
        same(&v, &fl(&row[3]), &format!("RowVar {} ({nr} x {nc})", row[0]));
        same(&[rcompat::sum(&v, false)], &[f(&row[4])], &format!("sum(RowVar) {}", row[0]));
    }
}
