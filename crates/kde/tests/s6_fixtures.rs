//! Stage S6 (docs/PLAN.md B5), teacher-forced on R's own pANN: for every (pN, pK) of the
//! tier-T sweeps, R's S5 pANN vector goes through `bkde`, `approxfun` on
//! `seq(min, max, length.out = n)` and the bimodality coefficient, and each result must equal
//! R's bit for bit: the density at the grid (S6_bkde, kept in data/ when present), the values
//! passed to `bimodality_coefficient` (S6_gkde, committed) and BCreal (reference/sweep_stats.tsv).

use kde::{bimodality_coefficient, bkde, density_function};
use rcompat::arith::seq_length_out;
use scio::{Npy, Table};
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

/// (file, ctx) -> (row_start, nrow) from a stage index.
fn blocks(index: &Table, file: &str) -> Vec<(String, usize, usize)> {
    index
        .rows
        .iter()
        .filter(|r| r[0] == file)
        .map(|r| (r[1].clone(), r[2].parse().unwrap(), r[3].parse().unwrap()))
        .collect()
}

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

#[test]
fn s6_teacher_forced_on_r_pann() {
    for object in OBJECTS {
        let fixture = root().join("tests/fixtures/doubletfinder").join(object);
        let stages = fixture.join("stages");
        let index = Table::read(&stages.join("index.tsv")).unwrap();
        let s5 = Npy::read(&stages.join("S5_pann_sweep.npy")).unwrap();
        let gkde = Npy::read(&stages.join("S6_gkde.npy")).unwrap();
        // the full bkde dump of the same run, when the data directory is present
        let data_bkde = root()
            .join("data/doubletfinder")
            .join(object)
            .join("runs/instrumented-portable-sapphirerapids/stages/S6_bkde.npy");
        let bkde_r = Npy::read(&data_bkde).ok();
        let sweep = Table::read(&fixture.join("reference/sweep_stats.tsv")).unwrap();
        let s5_blocks = blocks(&index, "S5_pann_sweep.npy");
        let g_blocks = blocks(&index, "S6_gkde.npy");
        let b_blocks = blocks(&index, "S6_bkde.npy");
        let npk = s5.ncol();
        assert_eq!(g_blocks.len(), s5_blocks.len() * npk, "{object}");
        assert_eq!(sweep.rows.len(), g_blocks.len(), "{object}");
        for (item, (ctx, start, len)) in g_blocks.iter().enumerate() {
            let (pn_i, pk_i) = (item / npk, item % npk);
            let (s5_ctx, s5_start, s5_len) = &s5_blocks[pn_i];
            assert!(ctx.starts_with(&format!("summary:pN_{}_pK_", &s5_ctx["sweep_pN".len()..])), "{object} {ctx}");
            let pann = s5.col_f64(pk_i, *s5_start, *s5_len);
            let est = bkde(&pann).unwrap();
            if let Some(b) = &bkde_r {
                let (bctx, bstart, blen) = &b_blocks[item];
                assert_eq!(bctx, ctx);
                assert!(same_bits(&est.x, &b.col_f64(0, *bstart, *blen)), "{object} {ctx}: bkde grid");
                assert!(same_bits(&est.y, &b.col_f64(1, *bstart, *blen)), "{object} {ctx}: bkde density");
            }
            let (lo, hi) = pann.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
            let xs = seq_length_out(lo, hi, pann.len() as f64).unwrap();
            let g = density_function(&est).unwrap().eval_all(&xs);
            assert!(same_bits(&g, &gkde.col_f64(0, *start, *len)), "{object} {ctx}: gkde");
            let bc = bimodality_coefficient(&g);
            let want: f64 = sweep.rows[item][2].parse().unwrap();
            assert_eq!(bc.to_bits(), want.to_bits(), "{object} {ctx}: BCreal {bc:e} vs {want:e}");
        }
    }
}
