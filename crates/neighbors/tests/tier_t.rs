//! Stages K3 to K6 on tier T (docs/PLAN.md C4 and C5), teacher-forced: FindNeighbors at the
//! workflow's settings (dims 1:30, k 20, 50 trees, Euclidean, prune 1/15) on the PCA embeddings
//! of the two base objects' instrumented R runs (portable profile, milestone M8's committed
//! fixtures), against that run's own dumps:
//!
//! - K3: the items, the first n nodes of Annoy's file, equal R's float rounding of the
//!   embeddings;
//! - K4: the whole index file, byte for byte (R's `save()` of a second, identical build);
//! - K5: every cell's 20 neighbour ids and float distances;
//! - K6: `SCT_nn` and `SCT_snn` as R stored them: dimensions, entries in stored order, values
//!   bitwise.
//!
//! No BLAS is involved, so these hold on any machine, with any thread count.

use neighbors::{find_neighbors, NeighborParams};
use scio::Npy;
use std::path::{Path, PathBuf};

fn fixture(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering").join(object).join("reference")
}

fn gunzip(path: &Path) -> Vec<u8> {
    scio::gz::gunzip(&std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))).unwrap()
}

fn check_object(object: &str, threads: usize) {
    let r = fixture(object);
    let emb = Npy::read(&r.join("pca/embeddings.npy")).unwrap();
    let (n, ncol) = (emb.nrow(), emb.ncol());
    let data = emb.block_f64(0, n);
    let p = NeighborParams::workflow();
    let got = find_neighbors(&data, n, ncol, &p, threads).unwrap();
    assert_eq!(got.k, 20);

    // K3: items are the embeddings' first 30 columns rounded to float
    for i in 0..n {
        let item = got.index.get_item(i);
        for z in 0..30 {
            assert_eq!(item[z].to_bits(), (emb.get(i, z) as f32).to_bits(), "{object}: item {i}, dim {z}");
        }
    }
    // K4: the index file
    let want = gunzip(&r.join("stages/annoy.ann.gz"));
    let bytes = got.index.to_bytes();
    assert_eq!(bytes.len(), want.len(), "{object}: index size");
    if bytes != want {
        let s = got.index.node_bytes();
        let first = bytes.iter().zip(&want).position(|(a, b)| a != b).unwrap();
        panic!("{object}: annoy.ann differs from R's, first at node {} ({})", first / s, got.index.field_at(first % s));
    }
    // K5: ids (1-based in R) and distances
    let ri = Npy::parse(&gunzip(&r.join("stages/K5_nn_index.npy.gz"))).unwrap();
    let rd = Npy::parse(&gunzip(&r.join("stages/K5_nn_dist.npy.gz"))).unwrap();
    assert_eq!((ri.nrow(), ri.ncol()), (n, 20));
    for i in 0..n {
        for j in 0..20 {
            assert_eq!(got.nn_idx[i * 20 + j] + 1, ri.get_i32(i, j), "{object}: cell {i}, neighbour {j}");
            assert_eq!((got.nn_dist[i * 20 + j] as f64).to_bits(), rd.get(i, j).to_bits(), "{object}: cell {i}, distance {j}");
        }
    }
    // K6: both graphs, entries in stored order
    for (name, g) in [("nn", &got.nn), ("snn", &got.snn)] {
        let e = scio::mtx::read_entries(&r.join(format!("graphs/{name}.mtx.gz"))).unwrap();
        assert_eq!((e.nrow, e.ncol, e.vals.len()), (g.nrow, g.ncol, g.nnz()), "{object} {name}: shape");
        let mut t = 0;
        for col in 0..g.ncol {
            let (rows, vals) = g.col(col);
            for (row, v) in rows.iter().zip(vals) {
                assert_eq!((e.rows[t], e.cols[t]), (*row, col as u32), "{object} {name}: entry {t}");
                assert_eq!(e.vals[t].to_bits(), v.to_bits(), "{object} {name}: value of entry {t}");
                t += 1;
            }
        }
    }
}

#[test]
fn base_tenx_object() {
    check_object("filtered_seurat_tenx_test", 1);
}

#[test]
fn base_emptydrops_object() {
    check_object("filtered_seurat_emptydrops_test", 1);
}

#[test]
fn independent_of_threads() {
    check_object("filtered_seurat_tenx_test", 7);
    check_object("filtered_seurat_emptydrops_test", 3);
}
