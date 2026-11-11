//! Milestone M3 on the tier-T fixtures (docs/PLAN.md B5 and B6): the DoubletFinder core on R's
//! own sweep subsample and PCA embeddings (teacher forcing) must reproduce the reference run
//! exactly: `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv` byte for byte, the `result` section
//! of `summary.json` field by field, and every S4-S7 stage block the fixtures hold, bit for bit;
//! S6 and S7 also teacher-forced on R's own S5. Results must not depend on the thread count.

use doubletfinder::pann::NeighbourDumps;
use doubletfinder::run::{run, Input, Output, Params};
use doubletfinder::stages::{stage_dumps, StageSource};
use doubletfinder::validate::{compare, Status};
use scio::{tenx, Json, Table};
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

fn input(object: &str) -> Input {
    let dir = fixture(object).join("input");
    let barcodes = tenx::read_barcodes(&dir).unwrap();
    let n_genes = tenx::read_features(&dir).unwrap().len();
    let t = Table::read(&dir.join("seurat_clusters.tsv")).unwrap();
    assert_eq!(t.header, ["barcode", "seurat_clusters"]);
    assert!(t.rows.iter().zip(&barcodes).all(|(r, b)| &r[0] == b));
    let clusters = t.rows.iter().map(|r| r[1].clone()).collect();
    let levels = std::fs::read_to_string(dir.join("seurat_clusters_levels.txt")).unwrap();
    Input { barcodes, clusters, cluster_levels: levels.lines().map(String::from).collect(), n_genes }
}

fn run_object(object: &str, threads: usize, dumps: bool) -> Output {
    let params = Params {
        threads,
        keep_stages: dumps,
        neighbour_dumps: NeighbourDumps { digests: dumps, positions: false },
        ..Params::default()
    };
    let mut source = StageSource::new(&fixture(object).join("stages")).unwrap();
    run(&input(object), &params, &mut source, &mut |_| {}).unwrap()
}

fn out_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m3_tier_t").join(name);
    if d.exists() {
        std::fs::remove_dir_all(&d).unwrap();
    }
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn reproduces_the_reference_runs() {
    for object in OBJECTS {
        let out = run_object(object, 3, true);
        let dir = out_dir(object);
        out.write_tables(&dir).unwrap();
        stage_dumps(&out).write(&dir.join("stages")).unwrap();
        let summary = Json::obj(vec![("result", out.result_json(None))]);
        std::fs::write(dir.join("summary.json"), scio::json::to_file_text(&summary)).unwrap();

        let reference = fixture(object).join("reference");
        for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
            assert!(
                std::fs::read(reference.join(file)).unwrap() == std::fs::read(dir.join(file)).unwrap(),
                "{object}: {file} differs from R's"
            );
        }
        let report = compare(&reference, &fixture(object).join("stages"), &dir, 2).unwrap();
        for c in &report.checks {
            assert_ne!(c.status, Status::Differ, "{object}: {} {}: {}", c.stage, c.name, c.detail);
        }
        // everything the fixtures hold was compared: files, every result field but the RNG
        // digest, the S4 digests, S5, S6_gkde and S7, and the teacher-forced S6 and S7
        let identical: Vec<&str> = report.checks.iter().filter(|c| c.status == Status::Identical).map(|c| c.name.as_str()).collect();
        for want in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv", "pK", "nExp", "k", "pANN_column", "S4_digest.npy",
                     "S5_pann_sweep.npy", "S5_pann_final.npy", "S6_gkde.npy", "S7_doublet.npy",
                     "S6_gkde (teacher-forced on R's S5)", "BCreal (teacher-forced on R's S5)", "S7_doublet (teacher-forced on R's S5)"] {
            assert!(identical.contains(&want), "{object}: {want} not compared");
        }
        assert_eq!(report.checks.iter().filter(|c| c.stage == "result" && c.status == Status::Identical).count(), 19, "{object}");
        // the saved workflow output's columns
        let saved = std::fs::read_to_string(fixture(object).join("saved/saved.json")).unwrap();
        assert!(saved.contains(&format!("\"pANN_column\": \"{}\"", out.pann_column)), "{object}");
    }
}

#[test]
fn independent_of_the_thread_count() {
    let object = "cellbender_fromraw_test";
    let one = run_object(object, 1, true);
    for threads in [2, 5, 16] {
        let many = run_object(object, threads, true);
        assert_eq!(one.cells_table(), many.cells_table(), "{threads} threads");
        assert_eq!(one.sweep_stats_table(), many.sweep_stats_table(), "{threads} threads");
        assert_eq!(one.bcmvn_table(), many.bcmvn_table(), "{threads} threads");
        let (a, b) = (one.sweep.as_ref().unwrap(), many.sweep.as_ref().unwrap());
        for (x, y) in a.iterations.iter().zip(&b.iterations) {
            assert_eq!(x.digests, y.digests);
        }
        assert_eq!(one.final_run.digests, many.final_run.digests);
    }
    // without dumps, the final run skips the sort of the neighbours: same pANN
    let fast = run_object(object, 4, false);
    assert_eq!(one.cells_table(), fast.cells_table());
}
