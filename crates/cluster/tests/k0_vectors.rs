//! Stage K0 (docs/PLAN.md C5): the generators FindClusters uses, against vectors from the
//! original code (M8, `tests/fixtures/clustering/k0/`):
//!
//! - `JavaRandom.nextInt(n)` from fresh generators, as Seurat 5.5.1's own `ModularityOptimizer.cpp`
//!   prints it (seeds as the `int` FindClusters passes, negative ones included; n from 1 to
//!   2^31 − 1, powers of two included);
//! - successive `generateRandomPermutation` draws from one generator;
//! - R 4.5.3's `set.seed(1); sample(x, 1)` for n tied names, the draw GroupSingletons makes,
//!   with `.Random.seed[2]` afterwards.

use cluster::JavaRandom;
use std::path::PathBuf;

fn k0(name: &str) -> Vec<Vec<String>> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/k0").join(name);
    let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    text.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').map(String::from).collect()).collect()
}

#[test]
fn javarandom_next_int() {
    let rows = k0("javarandom_nextint.tsv");
    assert!(rows.len() >= 250);
    let (mut pow2, mut big) = (0, 0);
    for r in &rows {
        let seed: i32 = r[0].parse().unwrap();
        let n: i32 = r[1].parse().unwrap();
        let count: usize = r[2].parse().unwrap();
        let want: Vec<i32> = r[3].split(',').map(|v| v.parse().unwrap()).collect();
        assert_eq!(want.len(), count);
        let mut g = JavaRandom::from_r_seed(seed);
        let got: Vec<i32> = (0..count).map(|_| g.next_int(n)).collect();
        assert_eq!(got, want, "seed {seed}, n {n}");
        pow2 += (n & n.wrapping_neg() == n) as usize;
        big += (n > 1 << 30) as usize;
    }
    // both branches, and the rejection loop's wrap-around region
    assert!(pow2 > 50 && big > 10, "{pow2} powers of two, {big} n above 2^30");
}

#[test]
fn javarandom_permutations() {
    let rows = k0("javarandom_permutation.tsv");
    assert_eq!(rows.len(), 15);
    for r in &rows {
        let seed: i32 = r[0].parse().unwrap();
        let sizes: Vec<usize> = r[1].split(',').map(|v| v.parse().unwrap()).collect();
        let perms: Vec<Vec<u32>> = r[2].split(';').map(|p| p.split(',').map(|v| v.parse().unwrap()).collect()).collect();
        let mut g = JavaRandom::from_r_seed(seed);
        for (s, want) in sizes.iter().zip(&perms) {
            assert_eq!(&g.permutation(*s), want, "seed {seed}, size {s}");
        }
    }
}

#[test]
fn groupsingletons_sample() {
    let rows = k0("groupsingletons_sample.tsv");
    assert!(rows.len() > 200);
    for r in &rows {
        let n: usize = r[0].parse().unwrap();
        let pos: usize = r[1].parse().unwrap();
        let mti: i32 = r[3].parse().unwrap();
        let names: Vec<String> = (1..=n).map(|i| format!("c{i}")).collect();
        let mut rng = rcompat::RRng::new(1);
        let got = rcompat::sample::sample(&mut rng, &names, 1, false).unwrap();
        assert_eq!(got[0], r[2], "n {n}");
        assert_eq!(names.iter().position(|x| *x == got[0]).unwrap() + 1, pos);
        assert_eq!(rng.random_seed()[1], mti, "n {n}: .Random.seed[2]");
    }
}
