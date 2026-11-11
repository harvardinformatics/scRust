//! Milestone M4 on the tier-T fixtures (docs/PLAN.md B5 and B6): Rust's SCTransform against the
//! instrumented R run's stage dumps, in all 7 contexts of each of the five test objects, each
//! sub-stage teacher-forced on R's own inputs: S2a from the counts, R's S1 draws and R's RNG state
//! at SCTransform's entry.

use doubletfinder::sct::{check_s2a, context, sct_contexts, RngCheckpoints};
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

fn report(object: &str, checks: &[Check]) -> usize {
    let bad: Vec<&Check> = checks.iter().filter(|c| c.status != Status::Identical).collect();
    for c in &bad {
        eprintln!("{object}: {} {} {:?} ({} of {} differ) {}", c.stage, c.name, c.status, c.n_differ, c.n, c.detail);
    }
    bad.len()
}

#[test]
fn s2a_setup_is_exact() {
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in OBJECTS {
        let dir = fixture(object);
        let counts = read_mtx(&dir.join("input")).unwrap();
        let stages = StageDir::open(&dir.join("stages")).unwrap();
        let rng = RngCheckpoints::open(&stages).unwrap();
        let ctxs = sct_contexts(&stages);
        assert_eq!(ctxs.len(), 7, "{object}: contexts");
        let mut checks = Vec::new();
        for ctx in &ctxs {
            let c = context(&stages, &rng, &counts, ctx).unwrap();
            let (_, _, ch) = check_s2a(&stages, &rng, &c, &params).unwrap();
            checks.extend(ch);
        }
        failures += report(object, &checks);
    }
    assert_eq!(failures, 0, "S2a checks differ (see stderr)");
}

/// Tier-P data (gitignored `data/`), when present: every instrumented run of `object`.
fn data_runs(object: &str) -> Vec<(PathBuf, PathBuf)> {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/doubletfinder").join(object);
    let mut out = Vec::new();
    if let Ok(rd) = std::fs::read_dir(base.join("runs")) {
        let mut runs: Vec<PathBuf> = rd.filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| {
            p.file_name().unwrap().to_string_lossy().starts_with("instrumented-") && p.join("stages/index.tsv").exists()
        }).collect();
        runs.sort();
        for r in runs {
            out.push((base.join("input"), r.join("stages")));
        }
    }
    out
}

/// S2a on the small production object (gene sampling, no cell sampling), when `data/` is present
/// and `SCRUST_DATA_TESTS` is set (a few seconds on one core).
#[test]
fn s2a_setup_small_production_object() {
    if std::env::var_os("SCRUST_DATA_TESTS").is_none() {
        return;
    }
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in ["soupx_tenx_SRR18256555", "soupx_tenx_10xmousebrain"] {
        let runs = data_runs(object);
        if runs.is_empty() {
            continue;
        }
        let counts = read_mtx(&runs[0].0).unwrap();
        for (_, st) in runs.iter().take(1) {
            let stages = StageDir::open(st).unwrap();
            let rng = RngCheckpoints::open(&stages).unwrap();
            let mut checks = Vec::new();
            for ctx in sct_contexts(&stages) {
                let c = context(&stages, &rng, &counts, &ctx).unwrap();
                let (_, _, ch) = check_s2a(&stages, &rng, &c, &params).unwrap();
                checks.extend(ch);
            }
            let n_ok = checks.iter().filter(|c| c.status == Status::Identical).count();
            eprintln!("{object} {}: {n_ok} of {} S2a checks identical", st.display(), checks.len());
            failures += report(object, &checks);
        }
    }
    assert_eq!(failures, 0, "S2a checks differ (see stderr)");
}

#[test]
fn s2b_fits_are_exact() {
    use doubletfinder::sct::{check_s2b, RS2a};
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in OBJECTS {
        let dir = fixture(object);
        let counts = read_mtx(&dir.join("input")).unwrap();
        let stages = StageDir::open(&dir.join("stages")).unwrap();
        let rng = RngCheckpoints::open(&stages).unwrap();
        let mut checks = Vec::new();
        for ctx in sct_contexts(&stages) {
            let c = context(&stages, &rng, &counts, &ctx).unwrap();
            let r = RS2a::read(&stages, &ctx).unwrap();
            let (_, ch) = check_s2b(&stages, &rng, &c, &r, &params, 4).unwrap();
            checks.extend(ch);
        }
        failures += report(object, &checks);
    }
    assert_eq!(failures, 0, "S2b checks differ (see stderr)");
}

#[test]
fn s2c_regularisation_is_exact() {
    use doubletfinder::sct::{check_s2c, Scalars, RS2a};
    let params = sctransform::SctParams::default();
    let mut failures = 0;
    for object in OBJECTS {
        let dir = fixture(object);
        let counts = read_mtx(&dir.join("input")).unwrap();
        let stages = StageDir::open(&dir.join("stages")).unwrap();
        let rng = RngCheckpoints::open(&stages).unwrap();
        let scalars = Scalars::open(&stages).unwrap();
        let mut checks = Vec::new();
        for ctx in sct_contexts(&stages) {
            let c = context(&stages, &rng, &counts, &ctx).unwrap();
            let r = RS2a::read(&stages, &ctx).unwrap();
            let (_, ch) = check_s2c(&stages, &scalars, &c, &r, &params).unwrap();
            checks.extend(ch);
        }
        failures += report(object, &checks);
    }
    assert_eq!(failures, 0, "S2c checks differ (see stderr)");
}
