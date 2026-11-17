//! Stage dumps in the layout of the instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`):
//! reading R's (teacher forcing, plan A5b) and writing Rust's, so that the two compare block by
//! block.
//!
//! A stage directory holds `.npy` arrays and `index.tsv`, one row per (array, context): `file`,
//! `ctx`, `row_start`, `nrow`, `ncol`, `dtype`, `stage`. Arrays of one quantity from all
//! contexts are stacked along rows; large matrices have one file per context. Contexts:
//! `sweep_pN<pN>` (pN as `paste()` writes it), `summary:<sweep result name>`, `final`.
//! Cells are 0-based rows of the embeddings (real cells, then doublets); order positions are R's
//! 1-based ones.

use crate::pann::Embeddings;
use crate::run::{EmbeddingSource, Output};
use crate::DfError;
use rcompat::as_character;
use scio::Npy;
use std::path::Path;

pub use scio::stages::{Block, StageDir, StageWriter};

/// Teacher forcing: R's own embeddings (stage S3) and sweep subsample (stage S1) from an
/// instrumented run's stage directory.
pub struct StageSource {
    stages: StageDir,
    embeddings: Option<Npy>,
}

impl StageSource {
    pub fn new(dir: &Path) -> Result<StageSource, DfError> {
        Ok(StageSource { stages: StageDir::open(dir)?, embeddings: None })
    }

    fn embeddings(&mut self, ctx: &str) -> Result<Embeddings, DfError> {
        if self.embeddings.is_none() {
            self.embeddings = Some(self.stages.read("S3_embeddings.npy")?);
        }
        let b = self
            .stages
            .block("S3_embeddings.npy", ctx)
            .ok_or_else(|| DfError::Input(format!("{}: no S3 embeddings for context {ctx}", self.stages.dir.display())))?
            .clone();
        let a = self.embeddings.as_ref().unwrap();
        Embeddings::new(a.block_f64(b.row_start, b.nrow), b.nrow, b.ncol)
    }
}

impl EmbeddingSource for StageSource {
    fn sweep_cells(&mut self, n: usize) -> Result<Vec<usize>, DfError> {
        match self.stages.block("S1_subsample.npy", "sweep") {
            None => Ok((0..n).collect()),
            Some(b) => {
                let b = b.clone();
                let a = self.stages.read("S1_subsample.npy")?;
                let cells: Vec<usize> = a.col_i32(0, b.row_start, b.nrow).into_iter().map(|i| i as usize).collect();
                if cells.iter().any(|&i| i >= n) {
                    return Err(DfError::Input("S1_subsample: cell index out of range".into()));
                }
                Ok(cells)
            }
        }
    }

    fn sweep_embeddings(&mut self, _pn_index: usize, pn: f64, _real_cells: &[usize], _n_doublets: usize) -> Result<Embeddings, DfError> {
        self.embeddings(&format!("sweep_pN{}", as_character(pn)))
    }

    fn final_embeddings(&mut self, _pn: f64, _n_doublets: usize) -> Result<Embeddings, DfError> {
        // the sweep's embeddings are not needed any more
        self.embeddings("final").inspect(|_| self.embeddings = None)
    }
}

/// Positions per cell (0-based), as R's L x n matrix of 1-based positions, column-major.
fn positions_matrix(per_cell: &[Vec<u32>]) -> (Vec<i32>, usize) {
    let l = per_cell.first().map_or(0, |p| p.len());
    let mut m = Vec::with_capacity(l * per_cell.len());
    for p in per_cell {
        debug_assert_eq!(p.len(), l);
        m.extend(p.iter().map(|&x| x as i32 + 1));
    }
    (m, l)
}

fn digest_columns(d: &[(u64, u64)]) -> Vec<f64> {
    let mut v: Vec<f64> = d.iter().map(|x| x.0 as f64).collect();
    v.extend(d.iter().map(|x| x.1 as f64));
    v
}

/// Rust's stage arrays S4-S7 for a run, in R's layout and with R's names: `S4_digest` (and
/// `S4_positions.<ctx>` when kept), `S5_pann_sweep`, `S5_pann_final`, `S6_bkde`, `S6_gkde`,
/// `S7_doublet`. S1-S3 (doublets, SCTransform, PCA) are not Rust's in teacher-forced runs and
/// are not written; a full Rust run records them with [`crate::record::Recorder`].
pub fn stage_dumps(out: &Output) -> StageWriter {
    let mut w = StageWriter::default();
    add_outputs(&mut w, out);
    w
}

/// Adds the stage arrays S4-S7 of `out` to `w` (see [`stage_dumps`]).
pub fn add_outputs(w: &mut StageWriter, out: &Output) {
    if let Some(sweep) = &out.sweep {
        for (it, pn) in sweep.iterations.iter().zip(&sweep.pn) {
            let ctx = format!("sweep_pN{}", pn.label);
            if let Some(d) = &it.digests {
                w.add_f64("S4_digest", &ctx, digest_columns(d), d.len(), 2, "S4");
            }
            if let Some(p) = &it.positions {
                let (m, l) = positions_matrix(p);
                w.add_big_i32("S4_positions", &ctx, m, l, p.len(), "S4");
            }
            let npk = it.pann.len();
            let mut v = Vec::with_capacity(it.n_real * npk);
            for col in &it.pann {
                v.extend_from_slice(col);
            }
            w.add_f64("S5_pann_sweep", &ctx, v, it.n_real, npk, "S5");
        }
        for (s, name) in sweep.stats.iter().zip(out.sweep_names()) {
            let ctx = format!("summary:{name}");
            if let Some(b) = &s.bkde {
                let mut v = b.x.clone();
                v.extend_from_slice(&b.y);
                w.add_f64("S6_bkde", &ctx, v, b.x.len(), 2, "S6");
            }
            if let Some(g) = &s.gkde {
                w.add_f64("S6_gkde", &ctx, g.clone(), g.len(), 1, "S6");
            }
        }
    }
    let f = &out.final_run;
    if let Some(d) = &f.digests {
        w.add_f64("S4_digest", "final", digest_columns(d), d.len(), 2, "S4");
    }
    if let Some(p) = &f.positions {
        let (m, l) = positions_matrix(p);
        w.add_big_i32("S4_positions", "final", m, l, p.len(), "S4");
    }
    w.add_f64("S5_pann_final", "final", f.pann.clone(), f.pann.len(), 1, "S5");
    w.add_i32("S7_doublet", "final", out.doublet.iter().map(|&d| d as i32).collect(), out.doublet.len(), 1, "S7");
}
