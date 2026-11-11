//! Stage S4 (docs/PLAN.md B5), teacher-forced on R's own PCA embeddings (stage S3): for every
//! sweep iteration and the final run of the tier-T objects, the first `ind` (sweep) or k + 1
//! (final) order positions of every real cell must equal R's. Checked through R's per-cell
//! digest (stages/S4_digest.npy, committed) and, when data/ is present, against the full
//! positions R dumped (S4_positions.<ctx>.npy).

use knn::{neighbours, Digest, Order, Points};
use scio::{Json, Npy, Table};
use std::path::{Path, PathBuf};

const OBJECTS: [&str; 5] = [
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn block(index: &Table, file: &str, ctx: &str) -> (usize, usize) {
    let r = index.rows.iter().find(|r| r[0] == file && r[1] == ctx).unwrap_or_else(|| panic!("{file} {ctx}"));
    (r[2].parse().unwrap(), r[3].parse().unwrap())
}

#[test]
fn s4_teacher_forced_on_r_embeddings() {
    for object in OBJECTS {
        let stages = root().join("tests/fixtures/doubletfinder").join(object).join("stages");
        let data_stages = root().join("data/doubletfinder").join(object).join("runs/instrumented-portable-sapphirerapids/stages");
        let index = Table::read(&stages.join("index.tsv")).unwrap();
        let json = Json::parse(&std::fs::read_to_string(stages.join("stages.json")).unwrap()).unwrap();
        let emb = Npy::read(&stages.join("S3_embeddings.npy")).unwrap();
        let digests = Npy::read(&stages.join("S4_digest.npy")).unwrap();
        let contexts: Vec<String> = match json.get("contexts").unwrap() {
            Json::Vector(v) => v.iter().map(|c| c.as_str().unwrap().to_string()).collect(),
            _ => panic!(),
        };
        for ctx in &contexts {
            let scalars = json.get("scalars").unwrap().get(ctx).unwrap();
            let (estart, n) = block(&index, "S3_embeddings.npy", ctx);
            let (dstart, n_real) = block(&index, "S4_digest.npy", ctx);
            let l = if ctx == "final" {
                scalars.get("k").unwrap().as_f64().unwrap() as usize + 1
            } else {
                scalars.get("ind").unwrap().as_f64().unwrap() as usize
            };
            let data = emb.block_f64(estart, n);
            let points = Points::new(&data, n, emb.ncol());
            let dg = Digest::new(l);
            let got = neighbours(&points, n_real, l, Order::Sorted, 4, |_, pos| (dg.of(pos), pos.to_vec()));
            for (i, ((h1, h2), _)) in got.iter().enumerate() {
                let want = (digests.get(dstart + i, 0), digests.get(dstart + i, 1));
                assert_eq!((*h1 as f64, *h2 as f64), want, "{object} {ctx} cell {i}: S4 digest");
            }
            if let Ok(positions) = Npy::read(&data_stages.join(format!("S4_positions.{ctx}.npy"))) {
                assert_eq!((positions.nrow(), positions.ncol()), (l.min(n), n_real));
                for (i, (_, pos)) in got.iter().enumerate() {
                    let want: Vec<i32> = positions.col_i32(i, 0, positions.nrow());
                    let have: Vec<i32> = pos.iter().map(|&p| p as i32 + 1).collect();
                    assert_eq!(have, want, "{object} {ctx} cell {i}: S4 positions");
                }
            }
        }
    }
}
