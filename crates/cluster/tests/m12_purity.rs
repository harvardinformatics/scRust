//! Milestone M12 on the committed tier-T objects (docs/PLAN.md C3e, C5: K9): bluster's
//! `neighborPurity`, its exhaustive search with knncolle's arithmetic and the certificate.
//!
//! - **Teacher-forced** (the two base objects, whose embeddings are here): R's own embeddings and
//!   labels give R's `K9_*` arrays bit for bit: each cell's 50th distance, the radius, the
//!   neighbour counts, the neighbour lists in BiocNeighbors' order, the weights, the totals, the
//!   purities and the maxima. No BLAS: any kernel, any machine.
//! - **From the counts** (the seven committed objects; the chain's other outputs are
//!   `m11_chain`'s): `neighborhood_purity` and `purity_maximum` against M8's portable runs, under
//!   the references' kernel, and the certificate: on an object with an odd number of cells the
//!   radius is the median cell's 50th distance, so that cell and its 50th neighbour lie exactly
//!   on it, and the certificate flags that pair (both ways) whenever it is within the radius;
//!   nothing else is flagged.

use cluster::chain::{self, Options};
use cluster::chain_output::Writer;
use cluster::find_clusters::Factor;
use cluster::metadata::RnaMetadata;
use cluster::purity::neighbor_purity;
use cluster::ChainParams;
use knn::purity::{Kind, Step};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
}

fn reference(object: &str) -> PathBuf {
    fixtures().join("clustering").join(object).join("reference")
}

fn gz_text(p: &Path) -> String {
    let path = if p.exists() { p.to_path_buf() } else { PathBuf::from(format!("{}.gz", p.display())) };
    String::from_utf8(scio::tenx::read_maybe_gz(&path).unwrap()).unwrap()
}

fn npy(p: &Path) -> Npy {
    let path = if p.exists() { p.to_path_buf() } else { PathBuf::from(format!("{}.gz", p.display())) };
    Npy::parse(&scio::tenx::read_maybe_gz(&path).unwrap()).unwrap()
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

/// R's labels with their levels.
fn labels(rdir: &Path) -> (Table, Factor) {
    let meta = Table::parse(&gz_text(&rdir.join("metadata.tsv"))).unwrap();
    let levels: Vec<String> = std::fs::read_to_string(rdir.join("seurat_clusters_levels.txt")).unwrap().lines().map(String::from).collect();
    let c = meta.col("seurat_clusters").unwrap();
    let codes: Vec<u32> = meta.rows.iter().map(|r| levels.iter().position(|l| *l == r[c]).unwrap() as u32).collect();
    (meta, Factor { name: "seurat_clusters".into(), levels, codes })
}

/// The median cell's pair on odd-sized objects: (the cell, its 50th neighbour), 1-based, as the
/// portable runs give them.
fn median_pair(object: &str) -> Option<(u32, u32)> {
    match object {
        "filtered_seurat_emptydrops_test" => Some((106, 136)),
        "soupx_emptydrops_test" => Some((141, 212)),
        "decontx_raw_emptydrops_test" => Some((117, 197)),
        _ => None,
    }
}

/// The certificate flags exactly the median pair, both ways, on the radius (0 ulps), labelled as
/// such.
fn check_certificate(object: &str, c: &knn::purity::Certificate) {
    match median_pair(object) {
        None => assert!(c.flags.is_empty(), "{object}: nothing flagged on an even-sized object, got {:?}", c.flags),
        Some((a, b)) => {
            let got: Vec<(Step, Kind, u32, u32, f64)> = c.flags.iter().map(|f| (f.step, f.kind, f.i + 1, f.j + 1, f.ulps())).collect();
            assert_eq!(got, vec![(Step::Radius, Kind::MedianPair, a, b, 0.0), (Step::Radius, Kind::MedianPair, b, a, 0.0)], "{object}: the median pair");
            assert_eq!(c.radius_gap_below, Some(0.0));
            assert_eq!(c.n_flagged_kind(Kind::Other), 0);
        }
    }
    assert!(c.kth_gap_ratio.is_none(), "{object}: no k-th distance has a neighbour within 2^-20 above it");
}

#[test]
fn purity_teacher_forced_on_r_s_embeddings() {
    for object in ["filtered_seurat_tenx_test", "filtered_seurat_emptydrops_test"] {
        let rdir = reference(object);
        let e = npy(&rdir.join("pca/embeddings.npy"));
        let (n, ncol) = (e.nrow(), e.ncol());
        let (meta, f) = labels(&rdir);
        for threads in [1, 3] {
            let p = neighbor_purity(&e.block_f64(0, n), n, ncol, &(1..=30).collect::<Vec<_>>(), &f, 50, true, true, threads).unwrap();
            let st = |name: &str| npy(&rdir.join("stages").join(format!("{name}.npy")));
            assert_eq!(bits(&p.dist_k), bits(&st("K9_dist50").f64s()), "{object}: 50th distances");
            let sc = Json::parse(&gz_text(&rdir.join("stages/stages.json"))).unwrap();
            let sc = sc.get("scalars").unwrap().get("chain").unwrap();
            assert_eq!(p.radius.to_bits(), sc.get("purity_radius").and_then(Json::as_f64).unwrap().to_bits(), "{object}: radius");
            assert_eq!(p.mean_neighbors().to_bits(), sc.get("purity_mean_neighbors").and_then(Json::as_f64).unwrap().to_bits(), "{object}: mean neighbours");
            let cnt = st("K9_neighbor_count");
            assert_eq!(p.neighbor_count.iter().map(|&c| c as i32).collect::<Vec<_>>(), cnt.col_i32(0, 0, n), "{object}: neighbour counts");
            let (ptr, idx) = p.neighbor_lists.as_ref().unwrap();
            let (rp, ri) = (st("K9_neighbors_ptr"), st("K9_neighbors_idx"));
            assert_eq!(ptr.iter().map(|&v| v as i32).collect::<Vec<_>>(), rp.col_i32(0, 0, n + 1), "{object}: list pointers");
            assert_eq!(idx.iter().map(|&v| v as i32 + 1).collect::<Vec<_>>(), ri.col_i32(0, 0, ri.nrow()), "{object}: lists, in BiocNeighbors' order");
            assert_eq!(bits(&p.weights), bits(&st("K9_weights").f64s()), "{object}: weights");
            assert_eq!(bits(&p.totals), bits(&st("K9_totals").f64s()), "{object}: totals");
            assert_eq!(bits(&p.purity), bits(&st("K9_purity").f64s()), "{object}: purities");
            let mx: Vec<i32> = p.maximum.iter().map(|m| m.unwrap() as i32).collect();
            assert_eq!(mx, st("K9_maximum").col_i32(0, 0, n), "{object}: maxima");
            let c = meta.col("neighborhood_purity").unwrap();
            let col: Vec<f64> = meta.rows.iter().map(|r| scio::parse17(&r[c]).unwrap()).collect();
            assert_eq!(bits(&p.purity), bits(&col), "{object}: metadata.tsv's neighborhood_purity");
            check_certificate(object, &p.certificate);
        }
    }
}

const OBJECTS: [&str; 7] = [
    "filtered_seurat_tenx_test",
    "filtered_seurat_emptydrops_test",
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn counts_dir(object: &str) -> PathBuf {
    let own = fixtures().join("clustering").join(object).join("input");
    if own.join("matrix.mtx.gz").exists() {
        own
    } else {
        fixtures().join("doubletfinder").join(object).join("input")
    }
}

#[test]
fn purity_from_the_counts_is_r_s() {
    for object in OBJECTS {
        let rdir = reference(object);
        rcompat::blas::require_fixture_kernel(&rdir, "neighborhood purity from the counts");
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m12_purity").join(object);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let input = counts_dir(object);
        let counts = scio::sparse::read_mtx(&input).unwrap();
        let mut w = Writer::new(&dir, scio::tenx::read_barcodes(&input).unwrap(), 2);
        w.features = Some(scio::tenx::read_features(&input).unwrap().iter().map(|f| (f.id.clone(), f.name.clone())).collect());
        w.rna = Some(RnaMetadata::of(&counts));
        let params = ChainParams::workflow();
        let r = chain::from_counts(counts, &params, 2, Options::default(), &mut w).unwrap();
        w.finish_labels(&r.clusters, r.silhouette.as_ref(), r.purity.as_ref()).unwrap();
        let p = r.purity.as_ref().unwrap();
        let (want, f) = labels(&rdir);
        let got = Table::read(&dir.join("metadata.tsv")).unwrap();
        let (a, b) = (want.col("neighborhood_purity").unwrap(), got.col("neighborhood_purity").unwrap());
        let differ = want.rows.iter().zip(&got.rows).filter(|(x, y)| x[a] != y[b]).count();
        assert_eq!(differ, 0, "{object}: neighborhood_purity differs in {differ} cells");
        // purity_maximum, where R kept K9_maximum
        let mx = rdir.join("stages/K9_maximum.npy");
        if mx.exists() || PathBuf::from(format!("{}.gz", mx.display())).exists() {
            let m = npy(&mx);
            let c = got.col("purity_maximum").unwrap();
            for (i, row) in got.rows.iter().enumerate() {
                assert_eq!(row[c], f.levels[m.get_i32(i, 0) as usize], "{object}: purity_maximum of cell {}", i + 1);
            }
        }
        check_certificate(object, &p.certificate);
        // summary.json's record of it
        let pj = cluster::chain_output::purity_json(p, &f.levels, &w.barcodes);
        let c = pj.get("certificate").unwrap();
        assert_eq!(c.get("flagged_pairs").and_then(Json::as_f64).unwrap() as usize, p.certificate.flags.len());
        assert_eq!(c.get("flagged_median_pair").and_then(Json::as_f64).unwrap() as usize, p.certificate.flags.len());
        if let Some(Json::List(pairs)) = c.get("pairs") {
            assert!(pairs.iter().all(|q| q.get("kind").and_then(Json::as_str) == Some("median_pair")));
        }
    }
}

/// Unweighted purities, another k, and k capped at n - 1, against a direct computation from the
/// definitions (R's arithmetic in plain loops).
#[test]
fn other_settings_by_the_definition() {
    let rdir = reference("filtered_seurat_tenx_test");
    let e = npy(&rdir.join("pca/embeddings.npy"));
    let (n, ncol) = (e.nrow(), e.ncol());
    let emb = e.block_f64(0, n);
    let (_, f) = labels(&rdir);
    let dims: Vec<usize> = (1..=10).collect();
    for (k, weighted) in [(20usize, false), (5, true), (400, true)] {
        let p = neighbor_purity(&emb, n, ncol, &dims, &f, k, weighted, false, 2).unwrap();
        let ke = k.min(n - 1);
        assert_eq!(p.k, ke);
        // the raw distances, by the definition
        let raw = |i: usize, j: usize| (0..dims.len()).fold(0.0f64, |s, d| {
            let x = emb[d * n + i] - emb[d * n + j];
            s + x * x
        });
        let mut dk = Vec::with_capacity(n);
        for i in 0..n {
            let mut v: Vec<f64> = (0..n).filter(|&j| j != i).map(|j| raw(i, j)).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            dk.push(v[ke - 1].sqrt());
        }
        assert_eq!(bits(&p.dist_k), bits(&dk), "k {k}: k-th distances");
        let radius = rcompat::quantile::median(&dk);
        assert_eq!(p.radius.to_bits(), radius.to_bits());
        let t2 = radius * radius;
        let size: Vec<usize> = (0..f.levels.len()).map(|l| f.codes.iter().filter(|&&c| c as usize == l).count()).collect();
        for i in 0..n {
            let nk = p.uclust.len();
            let mut col = vec![0.0f64; nk];
            let w = |c: u32| if weighted { 1.0 / size[c as usize] as f64 } else { 1.0 };
            col[f.codes[i] as usize] = w(f.codes[i]);
            for j in 0..n {
                if j != i && raw(i, j) <= t2 {
                    col[f.codes[j] as usize] += w(f.codes[j]);
                }
            }
            assert_eq!(bits(&p.sums[i * nk..(i + 1) * nk]), bits(&col), "k {k}: cell {i}'s sums");
            let total = cluster::purity::row_sum(&col);
            assert_eq!(p.purity[i].to_bits(), (col[f.codes[i] as usize] / total).to_bits());
        }
    }
}
