//! Milestone M11 on the committed tier-T objects (docs/PLAN.md C5: K1, K2, K8; C4, "Tier T"):
//! the workflow's chain from the counts, and approxSilhouette teacher-forced on R's own
//! embeddings and labels.
//!
//! - **From the counts** (the two base objects, whose counts are here, and the five
//!   DoubletFinder inputs, whose counts are M2's fixtures): SCTransform, RunPCA, FindNeighbors,
//!   FindClusters and the silhouette widths at the workflow's settings, against M8's portable
//!   instrumented runs (`reference/`): every label, the levels and the workflow's files, every
//!   silhouette width bitwise, the counts' columns, the scalars R's hooks recorded (SCTransform's
//!   and irlba's), and for the base objects the PCA embeddings and the K3, K5, K7 and K8 arrays.
//!   SCTransform's and RunPCA's products are OpenBLAS calls and the references are portable runs,
//!   so the test first checks that the kernel in use is the references' (Haswell), and stops with
//!   both kernels named when it is not. The results do not depend on the thread count (under
//!   any kernel: that test compares Rust's runs with each other).
//! - **Teacher-forced silhouette** (the base objects): R's own embeddings and labels give R's
//!   `K8_*` arrays and `silhouette_width` column bit for bit. No BLAS: any kernel, any machine.

use cluster::chain::{self, Options};
use cluster::chain_output::Writer;
use cluster::find_clusters::Factor;
use cluster::metadata::RnaMetadata;
use cluster::{approx_silhouette, ChainParams, Observer};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

const OBJECTS: [&str; 7] = [
    "filtered_seurat_tenx_test",
    "filtered_seurat_emptydrops_test",
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
}

fn reference(object: &str) -> PathBuf {
    fixtures().join("clustering").join(object).join("reference")
}

/// The object's counts: the base objects' own, or M2's DoubletFinder input.
fn counts_dir(object: &str) -> PathBuf {
    let own = fixtures().join("clustering").join(object).join("input");
    if own.join("matrix.mtx.gz").exists() {
        own
    } else {
        fixtures().join("doubletfinder").join(object).join("input")
    }
}

fn gz_text(p: &Path) -> String {
    let path = if p.exists() { p.to_path_buf() } else { PathBuf::from(format!("{}.gz", p.display())) };
    String::from_utf8(scio::tenx::read_maybe_gz(&path).unwrap()).unwrap()
}

fn npy(p: &Path) -> Npy {
    let path = if p.exists() { p.to_path_buf() } else { PathBuf::from(format!("{}.gz", p.display())) };
    Npy::parse(&scio::tenx::read_maybe_gz(&path).unwrap()).unwrap()
}

fn out_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m11_chain").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// The whole chain on an object's counts into `dir`, with every stage dump.
fn run(object: &str, dir: &Path, threads: usize) -> cluster::Chain {
    let input = counts_dir(object);
    let barcodes = scio::tenx::read_barcodes(&input).unwrap();
    let features = scio::tenx::read_features(&input).unwrap();
    let counts = scio::sparse::read_mtx(&input).unwrap();
    let mut w = Writer::new(dir, barcodes, threads).with_dumps(scio::stages::DumpLevel::Full);
    w.features = Some(features.iter().map(|f| (f.id.clone(), f.name.clone())).collect());
    w.rna = Some(RnaMetadata::of(&counts));
    w.annoy = true;
    w.purity_lists = true;
    let params = ChainParams::workflow();
    w.neighbor_params = Some(params.neighbors.clone());
    let r = chain::from_counts(counts, &params, threads, Options { keep_partitions: false, keep_residuals: true, keep_purity_lists: true, correct: false, ..Options::default() }, &mut w).unwrap();
    w.finish_labels(&r.clusters, r.silhouette.as_ref(), r.purity.as_ref()).unwrap();
    w.finish_dumps().unwrap();
    r
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn the_chain_from_the_counts_is_r_s() {
    for object in OBJECTS {
        let rdir = reference(object);
        rcompat::blas::require_fixture_kernel(&rdir, "the clustering chain from the counts");
        let dir = out_dir(object);
        let r = run(object, &dir, 3);
        // metadata.tsv: the labels, the counts' columns, every silhouette width bitwise
        let want = Table::parse(&gz_text(&rdir.join("metadata.tsv"))).unwrap();
        let got = Table::parse(&std::fs::read_to_string(dir.join("metadata.tsv")).unwrap()).unwrap();
        assert_eq!(want.rows.len(), got.rows.len(), "{object}: cells");
        for col in ["barcode", "nCount_RNA", "nFeature_RNA", "SCT_snn_res.0.8", "seurat_clusters", "silhouette_width"] {
            let (a, b) = (want.col(col).unwrap(), got.col(col).unwrap_or_else(|| panic!("{object}: no {col}")));
            let differ = want.rows.iter().zip(&got.rows).filter(|(x, y)| x[a] != y[b]).count();
            assert_eq!(differ, 0, "{object}: {col} differs in {differ} cells");
        }
        // the levels and the workflow's files
        for f in ["seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt", "cluster_ids.txt", "nclusters.txt"] {
            assert_eq!(std::fs::read_to_string(rdir.join(f)).unwrap(), std::fs::read_to_string(dir.join(f)).unwrap(), "{object}: {f}");
        }
        // the scalars R's hooks recorded: SCTransform's, irlba's, Annoy's, the optimizer's
        let rs = Json::parse(&gz_text(&rdir.join("stages/stages.json"))).unwrap();
        let cs = Json::parse(&std::fs::read_to_string(dir.join("stages/stages.json")).unwrap()).unwrap();
        let (rs, cs) = (rs.get("scalars").unwrap().get("chain").unwrap(), cs.get("scalars").unwrap().get("chain").unwrap());
        let mut compared = 0;
        for key in [
            "vst_n_cells", "n_theta_inf", "reg_bw", "mean_cell_sum", "n_outliers", "n_all_poisson_genes", "n_smoothing_genes", "min_variance", "n_genes", "n_genes_step1",
            "n_cells_step1", "n_bins_fit", "clip_range", "pca_n_features", "irlba_iter", "irlba_mprod", "total_variance", "annoy_f", "annoy_n_nodes", "annoy_file_bytes", "knn_k",
            "knn_self_first", "knn_self_missing", "optimizer_communities", "n_singletons", "silhouette_dims",
        ] {
            let (a, b) = (rs.get(key).and_then(Json::as_f64), cs.get(key).and_then(Json::as_f64));
            assert!(a.is_some(), "{object}: R recorded no {key}");
            assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits), "{object}: scalar {key}");
            compared += 1;
        }
        assert_eq!(compared, 26);
        assert_eq!(r.clusters.n_cells, want.rows.len());
        // the base objects carry the PCA and the stage references
        if rdir.join("pca/embeddings.npy").exists() {
            let e = npy(&rdir.join("pca/embeddings.npy"));
            let t = Table::read(&dir.join("pca/embeddings.tsv")).unwrap();
            for c in 0..e.ncol() {
                let got: Vec<f64> = t.rows.iter().map(|row| scio::parse17(&row[c + 1]).unwrap()).collect();
                assert_eq!(bits(&e.col_f64(c, 0, e.nrow())), bits(&got), "{object}: embeddings PC_{}", c + 1);
            }
            for f in ["K5_nn_index.npy", "K7_optimizer_ids.npy", "K7_final_ids.npy", "K8_silhouette_other.npy"] {
                let (a, b) = (npy(&rdir.join("stages").join(f)), npy(&dir.join("stages").join(f)));
                assert_eq!((a.nrow(), a.ncol()), (b.nrow(), b.ncol()), "{object}: {f} shape");
                for c in 0..a.ncol() {
                    assert_eq!(a.col_i32(c, 0, a.nrow()), b.col_i32(c, 0, b.nrow()), "{object}: {f}");
                }
            }
            for f in ["K5_nn_dist.npy", "K8_silhouette_width.npy", "K8_centroids.npy", "K8_cluster_var.npy"] {
                let (a, b) = (npy(&rdir.join("stages").join(f)), npy(&dir.join("stages").join(f)));
                assert_eq!((a.nrow(), a.ncol()), (b.nrow(), b.ncol()), "{object}: {f} shape");
                assert_eq!(bits(&a.block_f64(0, a.nrow())), bits(&b.block_f64(0, b.nrow())), "{object}: {f}");
            }
            // K3 and K4: Annoy's items are the first nodes of R's saved index
            let ann = scio::tenx::read_maybe_gz(&rdir.join("stages/annoy.ann.gz")).unwrap();
            assert!(ann == std::fs::read(dir.join("stages/annoy.ann")).unwrap(), "{object}: annoy.ann");
            // K6: the graphs, every stored entry and value
            for g in ["nn.mtx.gz", "snn.mtx.gz"] {
                let (a, b) = (scio::mtx::read_entries(&rdir.join("graphs").join(g)).unwrap(), scio::mtx::read_entries(&dir.join("neighbors").join(g)).unwrap());
                assert_eq!((a.rows, a.cols, bits(&a.vals)), (b.rows, b.cols, bits(&b.vals)), "{object}: {g}");
            }
        }
        eprintln!("{object}: {} cells, {} clusters: identical to R's portable run", r.n_cells, r.clusters.seurat_clusters.levels.len());
    }
}

/// Rust against Rust under one kernel: any kernel.
#[test]
fn independent_of_the_thread_count() {
    let object = "soupx_emptydrops_test";
    let files = |t: usize| -> Vec<(String, Vec<u8>)> {
        let dir = out_dir(&format!("threads_{t}"));
        run(object, &dir, t);
        let mut v = Vec::new();
        for sub in ["", "pca", "neighbors", "stages"] {
            let d = dir.join(sub);
            let mut names: Vec<String> = std::fs::read_dir(&d).unwrap().filter_map(|e| e.ok()).filter(|e| e.path().is_file()).map(|e| e.file_name().to_string_lossy().into_owned()).collect();
            names.sort();
            for n in names {
                v.push((format!("{sub}/{n}"), std::fs::read(d.join(&n)).unwrap()));
            }
        }
        v
    };
    let one = files(1);
    assert!(one.len() > 40, "{} files", one.len());
    for t in [2, 5] {
        let many = files(t);
        assert_eq!(one.len(), many.len());
        for ((na, a), (nb, b)) in one.iter().zip(&many) {
            assert_eq!(na, nb);
            assert!(a == b, "{na} differs at {t} threads");
        }
    }
}

/// R's own embeddings and labels: approxSilhouette's widths, closest clusters, centroids and
/// within-cluster variances, bit for bit. BLAS-free, so on any kernel.
#[test]
fn silhouette_teacher_forced_on_r_s_embeddings() {
    for object in ["filtered_seurat_tenx_test", "filtered_seurat_emptydrops_test"] {
        let rdir = reference(object);
        let e = npy(&rdir.join("pca/embeddings.npy"));
        let (n, ncol) = (e.nrow(), e.ncol());
        let meta = Table::parse(&gz_text(&rdir.join("metadata.tsv"))).unwrap();
        let levels: Vec<String> = std::fs::read_to_string(rdir.join("seurat_clusters_levels.txt")).unwrap().lines().map(String::from).collect();
        let c = meta.col("seurat_clusters").unwrap();
        let codes: Vec<u32> = meta.rows.iter().map(|r| levels.iter().position(|l| *l == r[c]).unwrap() as u32).collect();
        let f = Factor { name: "seurat_clusters".into(), levels, codes };
        let s = approx_silhouette(&e.block_f64(0, n), n, ncol, &(1..=30).collect::<Vec<_>>(), &f).unwrap();
        assert_eq!(bits(&s.width), bits(&npy(&rdir.join("stages/K8_silhouette_width.npy")).f64s()), "{object}: widths");
        let w = meta.col("silhouette_width").unwrap();
        let col: Vec<f64> = meta.rows.iter().map(|r| scio::parse17(&r[w]).unwrap()).collect();
        assert_eq!(bits(&s.width), bits(&col), "{object}: metadata.tsv's silhouette_width");
        let o = npy(&rdir.join("stages/K8_silhouette_other.npy"));
        let other: Vec<i32> = s.other.iter().map(|x| x.unwrap() as i32).collect();
        assert_eq!(other, o.col_i32(0, 0, n), "{object}: closest other cluster");
        let cent = npy(&rdir.join("stages/K8_centroids.npy"));
        let (k, d) = (cent.nrow(), cent.ncol());
        for a in 0..k {
            for b in 0..d {
                assert_eq!(cent.get(a, b).to_bits(), s.centroids[a * d + b].to_bits(), "{object}: centroid {a}, {b}");
            }
        }
        assert_eq!(bits(&s.cluster_var), bits(&npy(&rdir.join("stages/K8_cluster_var.npy")).f64s()), "{object}: cluster_var");
    }
}

/// An observer that keeps nothing still sees every stage, in order.
#[test]
fn the_chain_reports_its_stages_in_order() {
    struct Order(Vec<String>);
    impl Observer for Order {
        fn phase(&mut self, name: &str) {
            self.0.push(name.into());
        }
    }
    // the base object's SNN graph and embeddings (R's): FindClusters and the scores only
    let rdir = reference("filtered_seurat_tenx_test");
    let e = scio::mtx::read_entries(&rdir.join("graphs/snn.mtx.gz")).unwrap();
    let snn = cluster::input::csc_from_entries(&e).unwrap();
    let emb = npy(&rdir.join("pca/embeddings.npy"));
    let mut o = Order(Vec::new());
    let r = chain::from_snn(&snn, Some((&emb.block_f64(0, emb.nrow()), emb.ncol())), &ChainParams::workflow(), 2, Options::default(), &mut o).unwrap();
    assert_eq!(o.0, vec!["cluster", "silhouette", "purity"]);
    assert!(r.silhouette.unwrap().computed);
    assert!(r.purity.unwrap().computed);
}
