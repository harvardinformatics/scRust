//! Milestone M6 on the tier-T fixtures (docs/PLAN.md B5 and B6): the whole of DoubletFinder in
//! Rust, from the counts, the cluster labels and the seed alone, on the five test objects.
//!
//! - Against R's direct reference run (`reference/`): `cells.tsv`, `sweep_stats.tsv` and
//!   `bcmvn.tsv` byte for byte, and every field of `summary.json`'s `result` section, including
//!   the digest of R's final `.Random.seed`.
//! - Against R's instrumented run (`stages/`): every stage array the fixtures hold, S1 to S7
//!   (one check per array and context), every `.Random.seed` checkpoint and their sequence, and
//!   every scalar; S6 and S7 recomputed from R's S5. The random-number state flows through the
//!   seven contexts from Rust's own draws: nothing is read from R's run.
//! - Against the saved workflow output (`saved/`): the doublet set, the kept cells' pANN, the
//!   column names.
//! - `singlets/` (M6b): the input's counts without the doublets, every stored entry as it is,
//!   the saved output's kept cells in order, the input's feature lines. (The saved objects'
//!   own counts are compared on `data/`, through `validation/slurm/m6_doubletfinder.sbatch`.)
//! - The results do not depend on the number of threads, `singlets/` included.
//!
//! SCTransform's and RunPCA's products are OpenBLAS calls and the fixtures are the portable
//! runs, so the tests need `OPENBLAS_CORETYPE=Haswell` (env.sh's default).

use doubletfinder::chain::{rng_digest, RustSource};
use doubletfinder::input::{read_counts, read_input};
use doubletfinder::pann::NeighbourDumps;
use doubletfinder::record::{DumpLevel, Recorder};
use doubletfinder::run::{run, Output, Params};
use doubletfinder::singlets::write_singlets;
use doubletfinder::validate::{compare, compare_saved, Check, Status};
use scio::{tenx, Json, Table};
use std::path::{Path, PathBuf};

/// `singlets/` is the input's counts restricted to the cells the run keeps, in input order, each
/// column exactly as stored; the cells are the saved workflow output's kept cells; the feature
/// lines are the input's.
fn check_singlets(object: &str, out: &Output, dir: &Path) {
    let input = fixture(object).join("input");
    let w = write_singlets(&input, out, dir, "scrust test", 2).unwrap();
    let written = dir.join("singlets");
    let m = scio::sparse::read_mtx(&input).unwrap();
    let got = scio::sparse::read_mtx(&written).unwrap();
    let keep: Vec<usize> = (0..m.ncol).filter(|&j| !out.doublet[j]).collect();
    assert_eq!((got.nrow, got.ncol, got.nnz()), (m.nrow, keep.len(), w.nnz), "{object}");
    let mut at = 0;
    for (k, &j) in keep.iter().enumerate() {
        let (rows, x) = m.col(j);
        let (grows, gx) = got.col(k);
        assert_eq!(rows, grows, "{object}: column {j}");
        assert!(x.iter().zip(gx).all(|(a, b)| a.to_bits() == b.to_bits()), "{object}: column {j}");
        at += x.len();
    }
    assert_eq!(at, got.nnz());
    let kept = Table::read(&fixture(object).join("saved/cells.tsv")).unwrap();
    let saved: Vec<String> = kept.rows.iter().map(|r| r[0].clone()).collect();
    assert_eq!(tenx::read_barcodes(&written).unwrap(), saved, "{object}: kept cells");
    let lines = |d: &Path| tenx::read_lines(&tenx::find(d, "features.tsv").unwrap()).unwrap();
    assert_eq!(lines(&written), lines(&input), "{object}: features");
    assert_eq!(w.field, tenx::read_mtx_header(&input).unwrap().field, "{object}");
}

const OBJECTS: [&str; 5] = [
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn fixture(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder").join(object)
}

fn haswell() -> bool {
    let ok = std::env::var("OPENBLAS_CORETYPE").ok().as_deref() == Some("Haswell");
    if !ok {
        eprintln!("skipped: the tier-T fixtures need OPENBLAS_CORETYPE=Haswell (source env.sh)");
    }
    ok
}

fn out_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m6_tier_t").join(name);
    if d.exists() {
        std::fs::remove_dir_all(&d).unwrap();
    }
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// The whole run on one fixture's input; with `dumps`, the stage arrays go to `<dir>/stages`.
/// Returns the output and the generator's final state.
fn run_object(object: &str, threads: usize, dumps: Option<&Path>) -> (Output, Vec<i32>) {
    let input_dir = fixture(object).join("input");
    let (input, _) = read_input(&input_dir, &input_dir.join("seurat_clusters.tsv"), None).unwrap();
    let counts = read_counts(&input_dir, &input).unwrap();
    let params = Params {
        threads,
        keep_stages: dumps.is_some(),
        neighbour_dumps: NeighbourDumps { digests: dumps.is_some(), positions: false },
        ..Params::default()
    };
    let mut source = RustSource::new(counts, params.seed, params.pcs, threads);
    if let Some(dir) = dumps {
        source = source.with_recorder(Recorder::new(&dir.join("stages"), DumpLevel::Light));
    }
    let out = run(&input, &params, &mut source, &mut |_| {}).unwrap();
    let state = source.rng_state();
    if let Some(r) = source.take_recorder() {
        r.finish(&out).unwrap();
    }
    (out, state)
}

fn failures(object: &str, checks: &[Check]) -> usize {
    let bad: Vec<&Check> = checks.iter().filter(|c| c.status == Status::Differ).collect();
    for c in &bad {
        eprintln!("{object}: {} {}: {} of {} differ; {}", c.stage, c.name, c.n_differ, c.n, c.detail);
    }
    bad.len()
}

#[test]
fn every_stage_from_the_counts_is_r_s() {
    if !haswell() {
        return;
    }
    let mut bad = 0;
    for object in OBJECTS {
        let dir = out_dir(object);
        let (out, state) = run_object(object, 3, Some(&dir));
        out.write_tables(&dir).unwrap();
        // the run's OpenBLAS kernel, as `scrust doubletfinder` records it in its fingerprint
        let run_info = Json::obj(vec![("fingerprint", Json::obj(vec![("openblas_corename", Json::Str(rcompat::blas::corename()))]))]);
        let summary = Json::obj(vec![("result", out.result_json(Some(rng_digest(&state)))), ("run", run_info)]);
        std::fs::write(dir.join("summary.json"), scio::json::to_file_text(&summary)).unwrap();

        let f = fixture(object);
        let report = compare(&f.join("reference"), &f.join("stages"), &dir, 2).unwrap();
        bad += failures(object, &report.checks);
        let saved = compare_saved(&f.join("saved"), &dir).unwrap();
        bad += failures(object, &saved.checks);
        assert_eq!(saved.checks.iter().filter(|c| c.status == Status::Identical).count(), 6, "{object}: saved-output checks");
        // the fixtures hold no exported singlets: that check is not applicable here
        assert!(saved.checks.iter().any(|c| c.name == "singlets" && c.status == Status::NotApplicable), "{object}");
        check_singlets(object, &out, &dir);

        let identical = |stage: &str| report.checks.iter().filter(|c| c.stage == stage && c.status == Status::Identical).count();
        // every result field, the final RNG state's digest included
        assert_eq!(identical("result"), 20, "{object}: result fields");
        assert_eq!(identical("environment"), 1, "{object}: the OpenBLAS kernel");
        // every RNG checkpoint of R's run, in R's order
        let rng_rows = Table::read(&f.join("stages/rng.tsv")).unwrap().rows.len();
        let rng_same = report.checks.iter().filter(|c| c.name.contains(": RNG ") && c.status == Status::Identical).count();
        assert_eq!(rng_same, rng_rows, "{object}: RNG checkpoints");
        assert!(report.checks.iter().any(|c| c.name == "checkpoint sequence" && c.status == Status::Identical), "{object}");
        // every array the fixture holds, in every context: S1 to S7
        for want in [
            "final: S1_doublets", "sweep_pN0.05: S1_doublets", "final: S2a_cells", "sweep_pN0.3: S2a_gene_stats", "final: S2a_genes_step1",
            "final: S2b_glm", "sweep_pN0.2: S2b_model_pars", "final: S2c_model_pars_fit", "final: S2c_smoothing_genes",
            "final: S2d_residual_stats", "final: S2d_variable_features", "final: S2d_pca_features", "final: S3_embeddings",
            "sweep_pN0.15: S3_loadings", "final: S3_sdev", "final: S3_irlba_d", "final: S4_digest", "sweep_pN0.25: S5_pann_sweep",
            "final: S5_pann_final", "S6_gkde", "final: S7_doublet", "final: irlba_iter", "final: total_variance", "sweep_pN0.1: reg_bw",
            "cells.tsv", "sweep_stats.tsv", "bcmvn.tsv",
        ] {
            assert!(report.checks.iter().any(|c| c.name == want && c.status == Status::Identical), "{object}: {want} not compared");
        }
        // what the fixtures leave out (per-context matrices, S6_bkde) is the only thing missing
        for c in report.checks.iter().filter(|c| c.status == Status::Missing) {
            assert!(c.detail.contains("reference's directory"), "{object}: {} {} missing: {}", c.stage, c.name, c.detail);
        }
    }
    assert_eq!(bad, 0, "checks differ (see stderr)");
}

#[test]
fn independent_of_the_thread_count() {
    if !haswell() {
        return;
    }
    let object = "soupx_tenx_test";
    let (one, s1) = run_object(object, 1, None);
    for threads in [2, 7] {
        let (many, s) = run_object(object, threads, None);
        assert_eq!(one.cells_table(), many.cells_table(), "{threads} threads");
        assert_eq!(one.sweep_stats_table(), many.sweep_stats_table(), "{threads} threads");
        assert_eq!(one.bcmvn_table(), many.bcmvn_table(), "{threads} threads");
        assert_eq!(s1, s, "{threads} threads: final RNG state");
    }
    // the reference run's own files
    let f = fixture(object).join("reference");
    let dir = out_dir("threads");
    one.write_tables(&dir).unwrap();
    for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
        assert!(std::fs::read(f.join(file)).unwrap() == std::fs::read(dir.join(file)).unwrap(), "{file}");
    }
    // singlets/, byte for byte, whatever the number of threads writing it
    let input = fixture(object).join("input");
    let files = |threads: usize| -> Vec<Vec<u8>> {
        let d = out_dir(&format!("threads_singlets_{threads}"));
        write_singlets(&input, &one, &d, "scrust test", threads).unwrap();
        ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"].iter().map(|f| std::fs::read(d.join("singlets").join(f)).unwrap()).collect()
    };
    let first = files(1);
    for threads in [2, 7] {
        assert!(files(threads) == first, "{threads} threads: singlets/ differs");
    }
}
