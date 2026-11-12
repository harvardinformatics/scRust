//! Milestone M5 on the tier-T fixtures (docs/PLAN.md B5 and B6): Rust's RunPCA (irlba) against the
//! instrumented R runs' S3 dumps in all 7 contexts of each of the five test objects, after Rust's
//! own SCTransform (the chained run: the counts, R's S1 draws and R's RNG state at SCTransform's
//! entry); and, with the tier-P data in `data/`, teacher-forced on R's own PCA input.
//! irlba's products are OpenBLAS calls and the fixtures are the portable runs, so the tests need
//! `OPENBLAS_CORETYPE=Haswell` (env.sh's default).

use doubletfinder::pca::{check_chained, check_teacher};
use doubletfinder::sct::{context, sct_contexts, RngCheckpoints, Scalars};
use doubletfinder::stages::StageDir;
use doubletfinder::validate::{Check, Status};
use scio::sparse::read_mtx;
use std::path::{Path, PathBuf};

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
        eprintln!("skipped: S3 needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
    }
    ok
}

fn report(object: &str, checks: &[Check]) -> usize {
    let bad: Vec<&Check> = checks.iter().filter(|c| c.status != Status::Identical).collect();
    for c in &bad {
        eprintln!("{object}: {} {} {:?} ({} of {} differ) {}", c.stage, c.name, c.status, c.n_differ, c.n, c.detail);
    }
    bad.len()
}

#[test]
fn s3_after_rust_sctransform_is_exact() {
    if !haswell() {
        return;
    }
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in OBJECTS {
        let dir = fixture(object);
        let counts = read_mtx(&dir.join("input")).unwrap();
        let stages = StageDir::open(&dir.join("stages")).unwrap();
        let rng = RngCheckpoints::open(&stages).unwrap();
        let scalars = Scalars::open(&stages).unwrap();
        let ctxs = sct_contexts(&stages);
        assert_eq!(ctxs.len(), 7, "{object}: contexts");
        let mut checks = Vec::new();
        for ctx in &ctxs {
            let c = context(&stages, &rng, &counts, ctx).unwrap();
            // thread counts that split nothing and everything give the same bits
            let threads = if ctx == "final" { 4 } else { 1 };
            checks.extend(check_chained(&stages, &rng, &scalars, &c, &params, threads).unwrap());
        }
        assert_eq!(checks.len(), 13 * 7, "{object}: checks");
        failures += report(object, &checks);
    }
    assert_eq!(failures, 0, "S3 checks differ (see stderr)");
}

/// Teacher-forced on R's own PCA input (`S2d_scale_data.<ctx>.npy`, in `data/` only), every
/// portable instrumented run of the tier-T objects, when `SCRUST_DATA_TESTS` is set.
#[test]
fn s3_teacher_forced_on_data() {
    if std::env::var_os("SCRUST_DATA_TESTS").is_none() || !haswell() {
        return;
    }
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in OBJECTS {
        let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/doubletfinder").join(object);
        for run in ["instrumented-portable-sapphirerapids", "instrumented-portable-cascadelake"] {
            let sdir = base.join("runs").join(run).join("stages");
            if !sdir.join("index.tsv").exists() {
                continue;
            }
            let counts = read_mtx(&base.join("input")).unwrap();
            let stages = StageDir::open(&sdir).unwrap();
            let rng = RngCheckpoints::open(&stages).unwrap();
            let scalars = Scalars::open(&stages).unwrap();
            let mut checks = Vec::new();
            for ctx in sct_contexts(&stages) {
                let c = context(&stages, &rng, &counts, &ctx).unwrap();
                let (ch, source) = check_teacher(&stages, &rng, &scalars, &c, &params, 2).unwrap();
                assert_eq!(source, "R", "{object} {run} {ctx}: R's PCA input");
                checks.extend(ch);
            }
            failures += report(&format!("{object} {run}"), &checks);
        }
    }
    assert_eq!(failures, 0, "teacher-forced S3 checks differ (see stderr)");
}

/// In R's own runs, irlba draws its start vector and nothing more: R's state at irlba's exit is
/// its state at irlba's entry advanced by `rnorm(nfeatures)`, in every context of every
/// instrumented run (tier T from the fixtures; tier P from `data/` when `SCRUST_DATA_TESTS` is
/// set). So the near-invariant-subspace branches, whose draws the port reproduces on R vectors
/// (crates/pca/tests/r_vectors.rs), are never reached by the reference data.
#[test]
fn r_irlba_draws_only_its_start_vector() {
    let mut dirs: Vec<PathBuf> = OBJECTS.iter().map(|o| fixture(o).join("stages")).collect();
    if std::env::var_os("SCRUST_DATA_TESTS").is_some() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/doubletfinder");
        for object in std::fs::read_dir(&data).into_iter().flatten().flatten() {
            for run in std::fs::read_dir(object.path().join("runs")).into_iter().flatten().flatten() {
                let s = run.path().join("stages");
                if run.file_name().to_string_lossy().starts_with("instrumented-") && s.join("index.tsv").exists() {
                    dirs.push(s);
                }
            }
        }
    }
    let mut contexts = 0;
    for dir in &dirs {
        let stages = StageDir::open(dir).unwrap();
        let rng = RngCheckpoints::open(&stages).unwrap();
        let scalars = Scalars::open(&stages).unwrap();
        for ctx in sct_contexts(&stages) {
            let nf = scalars.get(&ctx, "pca_n_features").unwrap() as usize;
            let mut g = rng.rng(&ctx, "irlba_entry").unwrap();
            let _ = rcompat::rnorm(&mut g, nf, 0.0, 1.0);
            assert_eq!(Some(g.random_seed()), rng.state(&ctx, "irlba_exit"), "{} {ctx}: irlba drew more than rnorm({nf})", dir.display());
            contexts += 1;
        }
    }
    eprintln!("{contexts} contexts in {} runs: irlba drew rnorm(nfeatures) only", dirs.len());
    assert!(contexts >= 35);
}
