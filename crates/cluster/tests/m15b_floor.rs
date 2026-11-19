//! Milestone M15b (docs/PLAN.md A5b, C5's K14): the floor of runs that regressed covariates out of
//! `scale.data` (SCTransform's `vars.to.regress`), `[clustering.floor.kernel_regressed]`.
//!
//! - **When it applies:** only when both runs record the same non-empty `vars.to.regress` (R's
//!   reference runs in `result.regression.vars`, scrust's in `run.parameters.sctransform`);
//!   every other comparison keeps `[clustering.floor.kernel]`, as before M15b.
//! - **What it holds:** R's own differences across kernels with the regression, the maxima
//!   validation/python/m15b_floor.py measured on M15b's 134 cases (R's production runs against
//!   its portable runs, both directions), and no more; they exceed `[clustering.floor.kernel]`
//!   where the regression makes `scale.data` kernel-dependent.

use cluster::validate::{floor, floor_for, floor_regressed, regressed};
use scio::grade::{self, Class, Outcome, Stats};
use scio::Json;

fn r_run(vars: &[&str]) -> Json {
    let v: Vec<String> = vars.iter().map(|v| format!("{v:?}")).collect();
    Json::parse(&format!("{{\"result\": {{\"regression\": {{\"vars\": [{}], \"rank\": 2}}}}}}", v.join(", "))).unwrap()
}

fn scrust_run(vars: &[&str]) -> Json {
    let v: Vec<String> = vars.iter().map(|v| format!("{v:?}")).collect();
    let sct = if vars.is_empty() { "{\"seed_use\": 12345}".to_string() } else { format!("{{\"seed_use\": 12345, \"vars_to_regress\": [{}]}}", v.join(", ")) };
    Json::parse(&format!("{{\"run\": {{\"parameters\": {{\"sctransform\": {sct}}}}}, \"result\": {{}}}}")).unwrap()
}

#[test]
fn the_regressed_floor_only_when_both_runs_regressed_the_same() {
    let name = |r: Option<&Json>, c: &Json| floor_for(r, c).0.source;
    assert_eq!(regressed(Some(&r_run(&["percent.mt", "nCount_RNA"]))), vec!["percent.mt", "nCount_RNA"]);
    assert_eq!(regressed(Some(&scrust_run(&["percent.mt"]))), vec!["percent.mt"]);
    assert!(regressed(Some(&scrust_run(&[]))).is_empty() && regressed(None).is_empty());
    let (f, v) = floor_for(Some(&r_run(&["percent.mt", "nCount_RNA"])), &scrust_run(&["percent.mt", "nCount_RNA"]));
    assert_eq!((f.source.as_str(), v), ("clustering.floor.kernel_regressed", Some(vec!["percent.mt".to_string(), "nCount_RNA".to_string()])));
    // no regression, another list, another order, one side only, or no reference: as before M15b
    for (r, c) in [(r_run(&[]), scrust_run(&[])), (r_run(&["percent.mt"]), scrust_run(&["percent.mt", "nCount_RNA"])),
                   (r_run(&["nCount_RNA", "percent.mt"]), scrust_run(&["percent.mt", "nCount_RNA"])), (r_run(&["percent.mt"]), scrust_run(&[])),
                   (r_run(&[]), scrust_run(&["percent.mt"]))] {
        assert_eq!(name(Some(&r), &c), "clustering.floor.kernel");
    }
    assert_eq!(name(None, &scrust_run(&["percent.mt"])), "clustering.floor.kernel");
    // a run of scrust as the reference (both scrust's records)
    assert_eq!(name(Some(&scrust_run(&["score"])), &scrust_run(&["score"])), "clustering.floor.kernel_regressed");
}

/// One measured array: the maxima over M15b's 134 cases of R's production runs against its
/// portable runs (validation/python/m15b_floor.py, data/clustering/m15b/floor/r_vs_r.tsv), each
/// statistic at its own maximum, so the worst any case showed.
fn measured(n: usize, n_differ: usize, rows_differ: usize, max_abs: f64, max_scaled: f64) -> Stats {
    Stats { n, n_differ, rows_differ, max_abs, max_scaled, max_rel: 0.0, flipped: 0, incomparable: false }
}

#[test]
fn the_regressed_floor_holds_r_against_r_and_no_more() {
    let tolerance = grade::mode(Some("Haswell"), Some("SapphireRapids"));
    let (regressed_floor, unregressed_floor) = (floor_regressed(), floor());
    let graded = |f: &scio::grade::Floor, key: &str, s: &Stats| grade::grade(s, Class::BlasDependent, &tolerance, f, key);
    let cases: Vec<(&str, Stats)> = vec![
        // scale.data: max_abs on soupx_emptydrops_4500871lane ids_ncount; the fraction (0.7145) on
        // posthocfilt_threshold_scdblfinder_decontx_raw_tenx_test ids_ncount (43,639 of 61,074)
        ("S2d_scale_data", measured(61_074, 43_639, 702, 5.577760475716786e-13, 2.5636277230434938e-14)),
        ("S3_embeddings_dims_1_30", measured(2_382_270, 2_373_288, 79_409, 6.181126721571673e-10, 2.541942065390956e-12)),
        ("S3_embeddings_dims_31_50", measured(3_260, 3_260, 163, 5.462276971124336e-06, 6.558590028828964e-07)),
        ("S3_loadings", measured(39_650, 39_536, 793, 1.6660416942912426e-07, 5.947046484191046e-07)),
        ("S3_sdev", measured(50, 48, 48, 2.617905892066119e-13, 4.959296748297619e-14)),
        ("total_variance", measured(1, 1, 1, 1.1368683772161603e-13, 2.1158087588688924e-16)),
        ("K8_silhouette_width", measured(79_409, 79_366, 79_366, 1.945443806050662e-12, 2.0765184218349937e-12)),
        // purity: 34 cells over 18 cases (at most 2 in one), the largest difference 0.0258
        ("K9_purity", measured(247, 34, 34, 0.0258032275131006, 0.0258032275131006)),
    ];
    for (key, s) in &cases {
        assert_eq!(graded(&regressed_floor, key, s), Outcome::WithinFloor, "{key}: R's own difference across kernels, with the regression");
        // beyond [clustering.floor.kernel], which the regression's runs therefore do not use
        assert!(matches!(graded(&unregressed_floor, key, s), Outcome::Differ(_)), "{key}: within the unregressed floor");
    }
    // and no more: past any bound, a difference
    let past = |s: &Stats, f: &dyn Fn(&mut Stats)| {
        let mut t = s.clone();
        f(&mut t);
        t
    };
    let sd = &cases[0].1;
    assert!(matches!(graded(&regressed_floor, "S2d_scale_data", &past(sd, &|t| t.max_abs = 5.7e-13)), Outcome::Differ(m) if m.contains("max_abs")));
    assert!(matches!(graded(&regressed_floor, "S2d_scale_data", &past(sd, &|t| t.n_differ = t.n * 73 / 100)), Outcome::Differ(m) if m.contains("max_fraction_differing")));
    let pur = &cases[7].1;
    assert!(matches!(graded(&regressed_floor, "K9_purity", &past(pur, &|t| t.rows_differ = 35)), Outcome::Differ(m) if m.contains("cells_differing")));
    assert!(matches!(graded(&regressed_floor, "K9_purity", &past(pur, &|t| t.max_abs = 0.027)), Outcome::Differ(m) if m.contains("max_abs")));
    // an array the regressed floor does not list stays exact (it stands alone)
    assert!(matches!(graded(&regressed_floor, "K5_nn_dist", &cases[1].1), Outcome::Differ(m) if m.contains("not in the floor")));
    // the same kernel: exact, whatever the floor
    let exact = grade::mode(Some("SapphireRapids"), Some("SapphireRapids"));
    assert!(matches!(grade::grade(&cases[0].1, Class::BlasDependent, &exact, &regressed_floor, "S2d_scale_data"), Outcome::Differ(_)));
}
