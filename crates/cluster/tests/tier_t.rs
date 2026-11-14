//! Stage K7 (docs/PLAN.md C5) on the committed tier-T objects, from R's own SNN graphs
//! (teacher forcing): FindClusters at the workflow's settings, against M8's references
//! (`tests/fixtures/clustering/<object>/`):
//!
//! - every random start's iterations, modularity (bitwise) and partition, against the
//!   instrumented optimizer (`optimizer/starts.tsv`, `partitions.npy.gz`, `ids.txt`);
//! - the optimizer's result and the final ids (`reference/stages/K7_*.npy.gz`), the labels
//!   (`reference/metadata.tsv.gz`), their levels and the workflow's files (`cluster_ids.txt`,
//!   `nclusters.txt`), and the communities R's hooks recorded;
//! - the files `scrust cluster` writes, the same at 1, 3 and 7 threads.
//!
//! The two base objects carry their SNN graphs; the other tier-T objects are checked by
//! `validation/slurm/m10_cluster.sbatch` from `data/clustering/`.

use cluster::{find_clusters, ClusterParams};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

fn fixture(object: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering").join(object)
}

fn read(p: &Path) -> Vec<u8> {
    let gz = PathBuf::from(format!("{}.gz", p.display()));
    scio::tenx::read_maybe_gz(if p.exists() { p } else { &gz }).unwrap()
}

fn text(p: &Path) -> String {
    String::from_utf8(read(p)).unwrap()
}

fn npy(p: &Path) -> Npy {
    Npy::parse(&read(p)).unwrap()
}

fn check_object(object: &str) {
    let dir = fixture(object);
    let input = cluster::input::read_snn_dir(&dir.join("reference")).unwrap();
    let n = input.barcodes.len();
    let params = ClusterParams::workflow();
    let r = find_clusters(&input.snn, &params, 1, true).unwrap();
    assert_eq!(r.resolutions.len(), 1);
    let res = &r.resolutions[0];
    let o = &res.optimizer;

    // every start, against the instrumented optimizer
    let starts = Table::parse(&text(&dir.join("optimizer/starts.tsv"))).unwrap();
    assert_eq!(starts.rows.len(), o.starts.len());
    for (s, row) in starts.rows.iter().enumerate() {
        let it: usize = row[1].parse().unwrap();
        let q = scio::parse17(&row[2]).unwrap();
        assert_eq!(o.starts[s].iterations, it, "{object} start {}", s + 1);
        assert_eq!(o.starts[s].modularity.to_bits(), q.to_bits(), "{object} start {}: {} vs {}", s + 1, o.starts[s].modularity, q);
    }
    let parts = npy(&dir.join("optimizer/partitions.npy"));
    assert_eq!((parts.nrow(), parts.ncol()), (o.starts.len(), n));
    for (s, st) in o.starts.iter().enumerate() {
        let p = st.partition.as_ref().unwrap();
        for c in 0..n {
            assert_eq!(p[c] as i32, parts.get_i32(s, c), "{object} start {} cell {c}", s + 1);
        }
    }
    let ids: Vec<u32> = text(&dir.join("optimizer/ids.txt")).lines().map(|l| l.parse().unwrap()).collect();
    assert_eq!(o.ids, ids);

    // R's run: the optimizer's result, the final ids, the labels and levels
    let rr = dir.join("reference");
    let k7o = npy(&rr.join("stages/K7_optimizer_ids.npy"));
    let k7f = npy(&rr.join("stages/K7_final_ids.npy"));
    for c in 0..n {
        assert_eq!(o.ids[c] as i32, k7o.get_i32(c, 0));
        assert_eq!(res.column.labels()[c], k7f.get_i32(c, 0).to_string());
    }
    let meta = Table::parse(&text(&rr.join("metadata.tsv"))).unwrap();
    let (ccol, scol) = (meta.col("SCT_snn_res.0.8").unwrap(), meta.col("seurat_clusters").unwrap());
    for (c, row) in meta.rows.iter().enumerate() {
        assert_eq!(row[0], input.barcodes[c]);
        assert_eq!(res.column.labels()[c], row[ccol]);
        assert_eq!(r.seurat_clusters.labels()[c], row[scol]);
    }
    let lines = |name: &str| text(&rr.join(name)).lines().map(String::from).collect::<Vec<_>>();
    assert_eq!(res.column.levels, lines("SCT_snn_res.0.8_levels.txt"));
    assert_eq!(r.seurat_clusters.levels, lines("seurat_clusters_levels.txt"));
    assert_eq!(r.cluster_ids(), lines("cluster_ids.txt"));
    let stages = Json::parse(&text(&rr.join("stages/stages.json"))).unwrap();
    let sc = stages.get("scalars").and_then(|s| s.get("chain")).unwrap();
    assert_eq!(sc.get("optimizer_communities").and_then(Json::as_f64), Some(o.n_communities as f64));
    assert_eq!(sc.get("n_singletons").and_then(Json::as_f64), Some(0.0));
    assert!(res.grouped.singletons.is_empty());

    // the files, at several thread counts (one resolution here; see multiple_resolutions below)
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/test-tmp/cluster").join(object);
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let names = cluster::output::write_labels(&out, &input.barcodes, &r).unwrap();
    for name in ["cluster_ids.txt", "nclusters.txt", "SCT_snn_res.0.8_levels.txt", "seurat_clusters_levels.txt"] {
        assert_eq!(std::fs::read(out.join(name)).unwrap(), read(&rr.join(name)), "{object} {name}");
    }
    assert!(names.contains(&"metadata.tsv".to_string()));
}

#[test]
fn filtered_seurat_tenx_test() {
    check_object("filtered_seurat_tenx_test");
}

#[test]
fn filtered_seurat_emptydrops_test() {
    check_object("filtered_seurat_emptydrops_test");
}

#[test]
fn multiple_resolutions_are_independent_of_threads() {
    let input = cluster::input::read_snn_dir(&fixture("filtered_seurat_tenx_test").join("reference")).unwrap();
    let params = ClusterParams { resolutions: vec![0.2, 0.4, 0.8, 1.2, 2.0], ..ClusterParams::workflow() };
    let one = find_clusters(&input.snn, &params, 1, true).unwrap();
    // each resolution is the same as alone
    let alone = find_clusters(&input.snn, &ClusterParams::workflow(), 1, true).unwrap();
    assert_eq!(one.resolutions[2].optimizer.ids, alone.resolutions[0].optimizer.ids);
    assert_eq!(one.resolutions[2].column.levels, alone.resolutions[0].column.levels);
    assert_eq!(one.seurat_clusters.name, "seurat_clusters");
    assert_eq!(one.resolutions[4].column.name, "SCT_snn_res.2");
    let json = |r: &cluster::FindClusters| scio::json::to_file_text(&cluster::output::result_json(r));
    for t in [3, 7] {
        let many = find_clusters(&input.snn, &params, t, true).unwrap();
        assert_eq!(json(&many), json(&one), "{t} threads");
        assert_eq!(many.seurat_clusters, one.seurat_clusters);
        for (a, b) in many.resolutions.iter().zip(&one.resolutions) {
            assert_eq!(a.optimizer.starts, b.optimizer.starts);
            assert_eq!(a.column, b.column);
        }
    }
}
