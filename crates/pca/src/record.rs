//! Stage dumps of one `RunPCA()` call (stage S3, and the PCA's input under S2d), in the layout
//! of the instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`, whose RunPCA hooks the
//! clustering runs reuse). Shared by DoubletFinder's runs (M6) and the clustering chain (M11).

use crate::Pca;
use scio::stages::{DumpLevel, StageRecorder};
use scio::Json;
use std::io;

/// Records one RunPCA call in context `ctx`, as R's hooks in `RunPCA.default` and `irlba` do: the
/// RNG checkpoints `pca_entry`, `irlba_entry`, `irlba_exit` and `pca_exit`; `S2d_pca_features`
/// (the input rows PrepDR keeps, as rows of the counts: `input_rows[k]` is input row `k`'s);
/// at level standard and above `S2d_scale_data` (the matrix RunPCA.default receives: the kept
/// rows of `input`, which is the row-major `input_rows.len() x ncells` matrix given to
/// [`crate::run_pca`]); `S3_embeddings`, `S3_loadings`, `S3_sdev`, `S3_irlba_d`; and the scalars
/// `pca_n_features`, `total_variance`, `irlba_iter`, `irlba_mprod`. `entry` and `exit` are the
/// generator's states at RunPCA's entry and exit; `pca_entry` is left out when the entry state
/// is not known (a chain started from SCTransform's output: RunPCA reseeds, so it does not
/// matter).
pub fn record(rec: &mut StageRecorder, ctx: &str, input_rows: &[usize], input: &[f64], p: &Pca, entry: Option<&[i32]>, exit: &[i32]) -> io::Result<()> {
    if let Some(e) = entry {
        rec.rng(ctx, "S3", "pca_entry", e.to_vec());
    }
    rec.rng(ctx, "S3", "irlba_entry", p.rng_irlba_entry.clone());
    rec.rng(ctx, "S3", "irlba_exit", exit.to_vec());
    rec.rng(ctx, "S3", "pca_exit", exit.to_vec());
    let standard = rec.level >= DumpLevel::Standard;
    let dir = rec.dir.clone();
    let w = &mut rec.writer;
    let feats: Vec<i32> = p.features.iter().map(|&k| input_rows[k] as i32).collect();
    w.add_i32("S2d_pca_features", ctx, feats.clone(), feats.len(), 1, "S2d");
    if standard {
        // the matrix RunPCA.default receives: scale.data's rows that PrepDR keeps
        let n = p.ncells;
        if p.features.len() == input_rows.len() {
            w.write_big_f64_rows(&dir, "S2d_scale_data", ctx, input, feats.len(), n, "S2d")?;
        } else {
            let kept: Vec<f64> = p.features.iter().flat_map(|&k| input[k * n..(k + 1) * n].iter().copied()).collect();
            w.write_big_f64_rows(&dir, "S2d_scale_data", ctx, &kept, feats.len(), n, "S2d")?;
        }
    }
    w.add_f64("S3_embeddings", ctx, p.embeddings.clone(), p.ncells, p.npcs, "S3");
    w.add_f64("S3_loadings", ctx, p.loadings.clone(), p.features.len(), p.npcs, "S3");
    w.add_f64("S3_sdev", ctx, p.sdev.clone(), 1, p.npcs, "S3");
    w.add_f64("S3_irlba_d", ctx, p.irlba.d.clone(), 1, p.npcs, "S3");
    for (k, v) in [
        ("pca_n_features", Json::Int(feats.len() as i64)),
        ("total_variance", Json::Num(p.total_variance)),
        ("irlba_iter", Json::Int(p.irlba.iter as i64)),
        ("irlba_mprod", Json::Int(p.irlba.mprod as i64)),
    ] {
        rec.scalar(ctx, k, v);
    }
    Ok(())
}
