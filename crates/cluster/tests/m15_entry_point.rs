//! Milestone M15 (docs/PLAN.md C3h; A3's in-memory hand-offs): the chain's library entry point,
//! [`cluster::cluster_counts`], from counts in memory to a complete typed result.
//!
//! - Counts taken or borrowed give the same result, and borrowed counts are left as they were.
//! - What [`Options`] keeps changes nothing computed: the result that keeps only the labels (as
//!   clustering on the fly does) holds the same labels, scores and stage summaries as the one
//!   that keeps every stage `scrust cluster`'s files need.
//! - The observer-only form ([`chain::from_counts`]) gives the same clustering.
//!
//! Rust against Rust, so under any OpenBLAS kernel.

use cluster::chain::{self, Clustered, NoObserver, Options};
use cluster::ChainParams;
use std::borrow::Cow;
use std::path::{Path, PathBuf};

fn input(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder").join(object).join("input")
}

/// The same clustering: its result section (deterministic: no timings) and the labels.
fn assert_same_clustering(a: &cluster::FindClusters, b: &cluster::FindClusters) {
    assert_eq!(cluster::output::result_json(a), cluster::output::result_json(b), "FindClusters' result");
    assert_eq!(a.seurat_clusters, b.seurat_clusters, "seurat_clusters");
}

/// Everything two results computed, whatever they kept.
fn assert_same_computation(a: &Clustered, b: &Clustered) {
    assert_eq!(a.cells, b.cells);
    assert_eq!(a.features, b.features);
    assert_eq!(a.rna, b.rna);
    assert_same_clustering(&a.chain.clusters, &b.chain.clusters);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    match (&a.chain.silhouette, &b.chain.silhouette) {
        (Some(x), Some(y)) => assert_eq!((bits(&x.width), &x.other), (bits(&y.width), &y.other), "silhouette"),
        (x, y) => assert_eq!(x.is_some(), y.is_some(), "silhouette computed in one only"),
    }
    match (&a.chain.purity, &b.chain.purity) {
        (Some(x), Some(y)) => assert_eq!((bits(&x.purity), &x.maximum), (bits(&y.purity), &y.maximum), "purity"),
        (x, y) => assert_eq!(x.is_some(), y.is_some(), "purity computed in one only"),
    }
    assert_eq!(a.sct.summary, b.sct.summary, "SCTransform's summary");
    assert_eq!((&a.sct.rng_exit, &a.sct.scale_data_sha256), (&b.sct.rng_exit, &b.sct.scale_data_sha256));
    assert_eq!((&a.pca.summary, &a.pca.rng_exit), (&b.pca.summary, &b.pca.rng_exit), "RunPCA's summary");
    assert_eq!((&a.neighbors.summary, a.neighbors.nnz), (&b.neighbors.summary, b.neighbors.nnz), "FindNeighbors' summary");
}

#[test]
fn the_entry_point_keeps_what_it_is_asked_to_and_computes_the_same() {
    let counts = scio::Counts::read(&input("soupx_tenx_test")).unwrap();
    let before = counts.clone();
    let params = ChainParams::workflow();
    let all = Options { keep_reductions: true, keep_sct: true, ..Options::default() };
    // borrowed: the caller keeps its counts, untouched
    let full = chain::cluster_counts(Cow::Borrowed(&counts), &params, 2, all, &mut NoObserver).unwrap();
    assert_eq!(counts, before, "borrowed counts are left as they were");
    let (rows, p) = full.pca.kept.as_ref().expect("RunPCA kept");
    assert_eq!((p.ncells, p.npcs), (300, 50));
    assert!(!rows.is_empty() && p.features.len() <= rows.len());
    let nb = full.neighbors.kept.as_ref().expect("FindNeighbors kept");
    assert_eq!((nb.nn.nnz(), nb.snn.nnz()), full.neighbors.nnz);
    let sct = full.sct.kept.as_ref().expect("SCTransform's model kept");
    assert!(sct.rows.x.is_empty() && sct.residuals.y.is_none(), "the kept genes' counts and the residual matrix are freed");
    assert_eq!(sct.residuals.scale_data.len(), sct.residuals.top_features.len() * 300);
    assert!(sct.corrected.is_none() && full.sct.sct_columns.is_none(), "no corrected counts unless asked");
    assert_eq!(full.labels().len(), 300);
    assert_eq!(full.labels(), &full.chain.clusters.seurat_clusters);
    // taken, keeping only the labels (clustering on the fly), on another thread count
    let lean = chain::cluster_counts(Cow::Owned(counts.clone()), &params, 3, Options::default(), &mut NoObserver).unwrap();
    assert!(lean.pca.kept.is_none() && lean.neighbors.kept.is_none() && lean.sct.kept.is_none());
    assert_same_computation(&full, &lean);
    // the corrected counts, kept with the model
    let exported = chain::cluster_counts(Cow::Borrowed(&counts), &params, 1, Options { correct: true, ..all }, &mut NoObserver).unwrap();
    let c = exported.sct.kept.as_ref().unwrap().corrected.as_ref().expect("the corrected counts");
    assert_eq!(exported.sct.sct_columns, Some((c.n_count(), c.n_feature())));
    assert_same_clustering(&exported.chain.clusters, &full.chain.clusters);
    // the observer-only form
    let legacy = chain::from_counts(counts.matrix.clone(), &params, 2, Options::default(), &mut NoObserver).unwrap();
    assert_same_clustering(&legacy.clusters, &full.chain.clusters);
    // the scores only when asked
    let mut none = params.clone();
    none.stats.silhouette = false;
    none.stats.purity = false;
    let bare = chain::cluster_counts(Cow::Borrowed(&counts), &none, 2, Options::default(), &mut NoObserver).unwrap();
    assert!(bare.chain.silhouette.is_none() && bare.chain.purity.is_none());
    assert_same_clustering(&bare.chain.clusters, &full.chain.clusters);
}

/// The workflow's PCA guard, with doubletfinder.R's context in its message.
#[test]
fn the_guard_names_its_context() {
    let e = chain::pca_guard_in(40, 50, Some("doubletfinder")).unwrap_err().to_string();
    assert!(e.starts_with("doubletfinder: [LOW_QUALITY_SAMPLE] only 40 cell(s) remain - too few to compute 50 principal components"), "{e}");
    assert!(chain::pca_guard_in(51, 50, Some("doubletfinder")).is_ok());
    // on the counts: 40 cells, the preset's guard
    let counts = scio::Counts::read(&input("soupx_tenx_test")).unwrap();
    let keep: Vec<bool> = (0..counts.n_cells()).map(|j| j < 40).collect();
    let few = counts.keep_cells(&keep);
    let e = chain::cluster_counts(Cow::Owned(few), &ChainParams::workflow(), 1, Options { guard_context: Some("doubletfinder"), ..Options::default() }, &mut NoObserver).unwrap_err().to_string();
    assert!(e.starts_with("doubletfinder: [LOW_QUALITY_SAMPLE] only 40 cell(s) remain"), "{e}");
}
